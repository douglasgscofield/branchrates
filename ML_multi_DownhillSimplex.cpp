#include "ML_multi_DownhillSimplex.h"


void ML_multi_DownhillSimplex::print_data_structures() const {
    print_data_structures(simplex, funk_vals);
}

void ML_multi_DownhillSimplex::print_data_structures(const mat_ratep_type& p,
                                                     const v_ratep_type& y) const
{
    std::cout << "y:funk_vals[0.." << y.size()-1 << "] p:simplex[0.."
              << p.size()-1 << "][0.." << p[0].size()-1 << "]" << std::endl;
    for (int i = 0; i < p.size(); ++i) {
        std::cout << " v:funk_vals[" << i << "]=" << y[i] << " p:simplex=[";
        for (int j = 0; j < p[i].size(); ++j) {
            std::cout << p[i][j] << " ";
        }
        std::cout << "]" << std::endl;
    }
    std::cout << std::endl;
}


void ML_multi_DownhillSimplex::start_amoeba(int& nfunk)
{
    nfunk = 0;

    amoeba(simplex, 
           funk_vals, 
           get_function_tolerance(),
           &ML_multi::bounded_eval_neg_log_likelihood_at,
           // &ML_multi::eval_neg_log_likelihood_at,
           nfunk);

    // Now set the ratepvector and Likelihood parameters to reflect
    // the minimum found at vertex 0 of the simplex
    v_ratep_type parameters(simplex[0]);
    assert(parameters.size() == branch_rate_manager.get_ratepvector().size());
    branch_rate_manager.get_ratepvector().assign(parameters);
    compute();

    if (DEBUG_START_AMOEBA) {
        std::cout << "end of start_amoeba(): ";
        std::cout << " ftol=" << get_gradient_tolerance();
        std::cout << " nfunk=" << nfunk;
        std::cout << " minimum at simplex[0]=" << get_neg_log_likelihood();
        std::cout << std::endl << std::endl << "parm\tval" << std::endl;
        for (int i = 0; i < branch_rate_manager.get_ratepvector().size(); ++i) {
            std::cout << branch_rate_manager.get_ratepvector()[i].get_name()
                      << "\t"
                      << branch_rate_manager.get_ratepvector()[i].get_ratep()
                      << std::endl;
        }
        std::cout << std::endl;
        if (DEBUG_PRINT_DATA_STRUCTURES)
            print_data_structures();
    }
}


void ML_multi_DownhillSimplex::initialize_simplex(const v_ratep_type& p0,
                                                  const double s_len,
                                                  ptr_eval_func funk)
{
    const int n_parameters = p0.size();
    const int num_ratep = branch_rate_manager.get_ratepvector().size();
    assert(n_parameters == num_ratep);
    simplex.resize(n_parameters + 1);
    funk_vals.resize(n_parameters + 1);
    simplex[0] = p0;
    funk_vals[0] = (this->*funk)(simplex[0]);
    for (int i = 1; i < (n_parameters + 1); ++i) {
        simplex[i].resize(n_parameters);
        simplex[i] = p0;
        // now generate P_i by simulating offsets by unit vectors
        simplex[i][i-1] += (1.0 * s_len);
        funk_vals[i] = (this->*funk)(simplex[i]);
    }
    if (DEBUG_START_AMOEBA) {
        std::cout << "end of initialize_simplex(): " << std::endl;
        std::cout << "i\tratep[i].name\tratep[i].ratep\tp0[i]" << std::endl;
        for (int i = 0; i < branch_rate_manager.get_ratepvector().size(); ++i) {
            std::cout << i
                      << "\t"
                      << branch_rate_manager.get_ratepvector()[i].get_name()
                      << "\t"
                      << branch_rate_manager.get_ratepvector()[i].get_ratep()
                      << "\t"
                      << p0[i]
                      << std::endl;
        }
        std::cout << std::endl;
    }
}


void ML_multi_DownhillSimplex::amoeba(mat_ratep_type& p, 
                                      v_ratep_type& y,
                                      const double ftol,
                                      ptr_eval_func funk,
                                      int& nfunk
                                     )
{
    const double TINY=1.0e-10;
    int i, ihi, ilo, inhi, j;
    double rtol, ysave, ytry;
    int mpts = p.size();
    int ndim = p[0].size();
    v_ratep_type psum(ndim);
    get_psum(p, psum);
    for (;;) {
        ilo = 0;
        ihi = y[0] > y[1] ? (inhi = 1, 0) : (inhi = 0, 1);
        for (i = 0; i < mpts; ++i) {
            if (y[i] <= y[ilo]) ilo = i;
            if (y[i] > y[ihi]) {
                inhi = ihi;
                ihi = i;
            } else if (y[i] > y[inhi] && i != ihi) inhi = i;
        }
        rtol = 2.0 * std::abs(y[ihi] - y[ilo]) / 
               (std::abs(y[ihi]) + std::abs(y[ilo]) + TINY);
        if (rtol < ftol) {
            SWAP(y[0], y[ilo]);
            for (i = 0; i < ndim; ++i) SWAP(p[0][i], p[ilo][i]);
            break;
        }
        if (nfunk >= get_NMAX()) {
            if (CONFIG_DIE_ON_NMAX_EXCEEDED) {
                std::cerr << "amoeba: NMAX " << get_NMAX() << " exceeded "
                          << nfunk << std::endl;
                assert(false);
            }
            // otherwise, put the lowest at vertex 0 and return to try again
            SWAP(y[0], y[ilo]);
            for (i = 0; i < ndim; ++i) SWAP(p[0][i], p[ilo][i]);
            break;
        }
        nfunk += 2;
        ytry = amotry(p, y, psum, funk, ihi, -1.0);
        if (ytry <= y[ilo])
            ytry = amotry(p, y, psum, funk, ihi, 2.0);
        else if (ytry >= y[inhi]) {
            ysave = y[ihi];
            ytry = amotry(p, y, psum, funk, ihi, 0.5);
            if (ytry >= ysave) {
                for (i = 0; i < mpts; ++i) {
                    if (i != ilo) {
                        for (j = 0; j < ndim; ++j)
                            p[i][j] = psum[j] = 0.5*(p[i][j] + p[ilo][j]);
                        if (DEBUG_BOUNDS_TRACE) bounds_trace(p[i], "amoeba");
                        y[i] = (this->*funk)(psum);
                    }
                }
                nfunk += ndim;
                get_psum(p, psum);
            }
        } else --nfunk;
        if (DEBUG_AMOEBA) {
            if (! DEBUG_MONITOR_X10 || (nfunk % 10) == 0) {
                std::cout << "amoeba(): nfunk=" << nfunk;
                std::cout << " y[ilo]=" << y[ilo];
                std::cout << " y[ihi]=" << y[ihi];
                std::cout << std::endl;
            }
        }
    }
}


double ML_multi_DownhillSimplex::amotry(mat_ratep_type& p, 
                                        v_ratep_type& y,
                                        v_ratep_type& psum,
                                        ptr_eval_func funk,
                                        const int ihi,
                                        const double fac
                                       )
{
    int j, n_adjusted;
    double fac1, fac2, ytry;
    int ndim = p[0].size();
    v_ratep_type ptry(ndim);
    fac1 = (1.0 - fac) / ndim;
    fac2 = fac1 - fac;
    for (j = 0; j < ndim; ++j)
        ptry[j] = psum[j]*fac1 - p[ihi][j]*fac2;
    if (CONFIG_BOUNDS_ADJUST) {
        bounds_adjust(ptry, n_adjusted);
    }
    if (DEBUG_BOUNDS_TRACE) bounds_trace(ptry, "amotry");
    ytry = (this->*funk)(ptry);
    if (ytry < y[ihi]) {
        y[ihi] = ytry;
        for (j = 0; j < ndim; ++j) {
            psum[j] += ptry[j] - p[ihi][j];
            p[ihi][j] = ptry[j];
        }
        if (CONFIG_BOUNDS_ADJUST) {
            bounds_adjust(ptry, n_adjusted);
        }
        if (DEBUG_BOUNDS_TRACE) bounds_trace(p[ihi], "amotry");
    }
    if (DEBUG_AMOTRY) {
        std::cout << "amotry(): ytry=" << ytry << std::endl;
        if (DEBUG_PRINT_DATA_STRUCTURES)
            print_data_structures(p, y);
    }
    return(ytry);
}


void ML_multi_DownhillSimplex::get_psum(const mat_ratep_type& p,
                                        v_ratep_type& psum)
{
    int i, j;
    double sum;
    const int mpts = p.size();
    const int ndim = p[0].size();
    assert(ndim == psum.size());
    for (j = 0; j < ndim; ++j) {
        for (sum = 0.0, i = 0; i < mpts; ++i) {
            sum += p[i][j];
        }
        psum[j] = sum;
    }
}

