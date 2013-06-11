#include "ML_multi.h"

const double ML_multi::chisq_095_1df        = 3.84145882069412;
const double ML_multi::lrt_095_interval     = ML_multi::chisq_095_1df / 2.0;
const double ML_multi::lrt_profile_interval = ML_multi::lrt_095_interval * 2.0;
const double ML_multi::interval_parm_multiple_bounds[2] = {0.001, 100};
const double ML_multi::interval_parm_steps = 100;


double ML_multi::bounded_eval_neg_log_likelihood_at(const v_ratep_type& p)
{
    const int n_parameters = p.size();
    v_ratep_type newp(p);
    int n_adjusted;

    bounds_adjust(newp, n_adjusted);
    if (DEBUG_BOUNDS && n_adjusted == n_parameters) {
        std::cerr << "bounded_eval: all parameters violated bounds" << std::endl;
    }
    return(eval_neg_log_likelihood_at(newp));
}


double ML_multi::eval_neg_log_likelihood_at(const v_ratep_type& p)
{
    branch_rate_manager.get_ratepvector().assign(p);
    compute();
    // return(get_neg_log_likelihood());
    return(get_neg_log_likelihood());
}


void ML_multi::eval_neg_log_likelihood_gradient_at(const v_ratep_type& in,
                                                   v_ratep_type& out)
{
    assert(in.size() == out.size());
    RatePVector& rpv = branch_rate_manager.get_ratepvector();
    rpv.assign(in);
    compute();
    //double neg_logL = get_neg_log_likelihood();
    for (size_t i = 0; i < in.size(); ++i) {
        double old_parm = rpv[i].get_ratep();
        // An "old-style" numerical derivative
        // rpv[i].update_ratep(bounds[i].delta);
        // compute();
        // double neg_logL_grad = get_neg_log_likelihood();
        // double grad = (neg_logL_grad - neg_logL) / bounds[i].delta;
        // out[i] = grad;

        // An alternative, using what we've learned for calculating
        // numerical derivatives.
        ////////double delta = bounds[i].delta;
        double delta = Likelihood::delta;
        rpv[i].set_ratep(old_parm - delta);
        compute();
        double neg_logL_m_delta = get_neg_log_likelihood();
        rpv[i].set_ratep(old_parm + delta);
        compute();
        double neg_logL_p_delta = get_neg_log_likelihood();
        double grad = (neg_logL_p_delta - neg_logL_m_delta) /
                      (2.0 * delta);
        out[i] = grad;
        rpv[i].set_ratep(old_parm);
    }
}


#define   DEBUG_LNSRCH   0

void ML_multi::lnsrch(const v_ratep_type& xold,
                      const double fold,
                      const v_ratep_type& g,
                      v_ratep_type& p,
                      v_ratep_type& x,
                      double& f,
                      const double stpmax,
                      bool& check,
                      ptr_eval_func func
                      //double func(const v_ratep_type& p)
                     )
{
    // from Press et al 2002, pp. 389-390
    const double ALF = 1.0e-4, TOLX = std::numeric_limits<double>::epsilon();
    int i;
    double a, alam, alam2 = 0.0, alamin, b, disc, f2 = 0.0;
    double rhs1, rhs2, slope, sum, temp, test, tmplam;

    int n = xold.size();
    check = false;
    sum = 0.0;
    for (i = 0; i < n; ++i) sum += p[i]*p[i];
    sum = std::sqrt(sum);
    if (sum > stpmax)
        for (i = 0; i < n; ++i) p[i] *= stpmax/sum;
    slope = 0.0;
    for (i = 0; i < n; ++i)
        slope += g[i]*p[i];
    if (slope >= 0.0) {
        std::cerr << "lnsrch(): Roundoff problem" << std::endl;
    }
    if (DEBUG_LNSRCH) {
        std::cerr << "slope=" << slope << std::endl;
    }
    assert(slope < 0.0);
    test = 0.0;
    for (i = 0; i < n; ++i) {
        temp = std::abs(p[i]) / MAX(std::abs(xold[i]), 1.0);
        if (temp > test) test = temp;
    }
    alamin = TOLX / test;
    alam = 1.0;
    for (;;) {
        for (i = 0; i < n; ++i) x[i] = xold[i] + alam*p[i];
        if (DEBUG_LNSRCH) {
            std::cout << "alam=" << alam;
            for (int i = 0; i < n; ++i) {
                std::cout << "      xold[" << i << "]=" << xold[i];
                std::cout << " x[" << i << "]=" << x[i];
                std::cout << " p[" << i << "]=" << p[i];
                std::cout << std::endl;
            }
        }
        for (i = 0; i < n; ++i) {
            // Go through all the parameters, check for out of range.  If a
            // parameter exceeds an upper or lower bound, then peg the
            // parameter at the bound.
            if (x[i] < bounds[i].lowerbound) x[i] = bounds[i].lowerbound;
            if (x[i] > bounds[i].upperbound) x[i] = bounds[i].upperbound;
        }
        if (DEBUG_LNSRCH) {
            std::cout << "alam=" << alam;
            for (int i = 0; i < n; ++i) {
                std::cout << "      xold[" << i << "]=" << xold[i];
                std::cout << " x[" << i << "]=" << x[i];
                std::cout << " p[" << i << "]=" << p[i];
                std::cout << std::endl;
            }
        }
        f = (this->*func)(x);   // f = func(x);
        if (alam < alamin) {
            for (i = 0; i < n; ++i) x[i] = xold[i];
            check = true;
            return;
        } else if (f <= fold + ALF*alam*slope) {
            return;
        } else {
            if (alam == 1.0) {
                tmplam = -slope / (2.0 * (f - fold - slope));
            } else {
                rhs1 = f - fold - alam*slope;
                rhs2 = f2 - fold - alam2*slope;
                a = (rhs1/(alam*alam) - rhs2/(alam2*alam2)) / (alam - alam2);
                b = (-alam2*rhs1/(alam*alam) + alam*rhs2/(alam2*alam2)) /
                    (alam - alam2);
                if (a == 0.0) {
                    tmplam = -slope/(2.0*b);
                } else {
                    disc = b*b - 3.0*a*slope;
                    if (disc < 0.0) tmplam = 0.5*alam;
                    else if (b <= 0.0) tmplam = (-b + std::sqrt(disc))/(3.0*a);
                    else tmplam = -slope / (b + std::sqrt(disc));
                }
                if (tmplam > 0.5*alam)
                    tmplam = 0.5*alam;
            }
        }
        alam2 = alam;
        f2 = f;
        alam = MAX(tmplam, 0.1*alam);
    }
}


///////////////////////////
// bounds-checking of parameters
///////////////////////////
void ML_multi::bounds_adjust(v_ratep_type& p, int& n_adjusted) const
{
    n_adjusted = 0;
    for (size_t i = 0; i < p.size(); ++i) {
        if (p[i] < bounds[i].lowerbound) {
            if (DEBUG_BOUNDS) {
                std::cerr << "bounds_adjust: p[" << i << "]=" << p[i]
                          << " is less than lowerbound[" << i << "]="
                          << bounds[i].lowerbound << ", adjusted." << std::endl;
            }
            p[i] = bounds[i].lowerbound;
            ++n_adjusted;
        } else if (p[i] > bounds[i].upperbound) {
            if (DEBUG_BOUNDS) {
                std::cerr << "bounds_adjust: p[" << i << "]=" << p[i]
                          << " is greater than upperbound[" << i << "]="
                          << bounds[i].upperbound << ", adjusted." << std::endl;
            }
            p[i] = bounds[i].upperbound;
            ++n_adjusted;
        }
    }
}


void ML_multi::bounds_trace(const v_ratep_type& p, const std::string& loc) const
{
    int n_adjusted = 0;
    for (size_t i = 0; i < p.size(); ++i) {
        if (p[i] < bounds[i].lowerbound) {
            std::cerr << loc << ": p[" << i << "]=" << p[i]
                      << " is less than lowerbound[" << i << "]="
                      << bounds[i].lowerbound << std::endl;
            ++n_adjusted;
        } else if (p[i] > bounds[i].upperbound) {
            std::cerr << loc << ": p[" << i << "]=" << p[i]
                      << " is greater than upperbound[" << i << "]="
                      << bounds[i].upperbound << std::endl;
            ++n_adjusted;
        }
    }
    if (n_adjusted != 0) {
        std::cerr << loc << ": a total of " << n_adjusted
                  << " bounds violations." << std::endl;
        assert(n_adjusted == 0);
    }
}


///////////////////////////
// Surface evaluation routines
///////////////////////////
void ML_multi::surface_range_print(const std::string& ratep1,
                                   const std::string& ratep2,
                                   const double lower,
                                   const double upper,
                                   const double stepsize,
                                   std::ostream& os)
{
    RatePVector& rpv = branch_rate_manager.get_ratepvector();
    int rid1 = rpv.get_ratep_id(ratep1);
    int rid2 = rpv.get_ratep_id(ratep2);
    surface_range_print(rid1, rid2, lower, upper, stepsize, os);
}

void ML_multi::surface_range_print(const int rid1,
                                   const int rid2,
                                   const double lower,
                                   const double upper,
                                   const double stepsize,
                                   std::ostream& os)
{
    os << "ratep1.name" << "\t" << "ratep1.val";
    os << "\t" << "ratep2.name" << "\t" << "ratep2.val";
    os << "\t" << "neg.log.L" << std::endl;
    RatePVector& rpv = branch_rate_manager.get_ratepvector();
    for (double d1 = lower; d1 <= upper; d1 += stepsize) {
        for (double d2 = lower; d2 <= upper; d2 += stepsize) {
            rpv[rid1].set_ratep(d1);
            rpv[rid2].set_ratep(d2);
            compute();
            os << rpv[rid1].get_name() << "\t" << d1;
            os << "\t" << rpv[rid2].get_name() << "\t" << d2;
            os << "\t" << get_neg_log_likelihood();
            os << std::endl;
        }
    }
}


void ML_multi::surface_print(const std::string& ratep1,
                             const std::string& ratep2,
                             std::ostream& os,
                             const double lowerdelta,
                             const double upperdelta,
                             const double stepsize)
{
    RatePVector& rpv = branch_rate_manager.get_ratepvector();
    int rid1 = rpv.get_ratep_id(ratep1);
    int rid2 = rpv.get_ratep_id(ratep2);
    surface_print(rid1, rid2, os, lowerdelta, upperdelta, stepsize);
}


void ML_multi::surface_print(const int rid1,
                             const int rid2,
                             std::ostream& os,
                             const double lowerdelta,
                             const double upperdelta,
                             const double stepsize)
{
    os << "ratep1.name" << "\t" << "ratep1.val";
    os << "\t" << "ratep2.name" << "\t" << "ratep2.val";
    os << "\t" << "neg.log.L" << std::endl;
    RatePVector& rpv = branch_rate_manager.get_ratepvector();
    double old_parm1 = rpv[rid1].get_ratep();
    double old_parm2 = rpv[rid2].get_ratep();
    for (double d1 = lowerdelta; d1 <= upperdelta; d1 += stepsize) {
        for (double d2 = lowerdelta; d2 <= upperdelta; d2 += stepsize) {
            double new_parm1 = old_parm1 + d1;
            rpv[rid1].set_ratep(new_parm1);
            double new_parm2 = old_parm2 + d2;
            rpv[rid2].set_ratep(new_parm2);
            compute();
            os << rpv[rid1].get_name() << "\t" << new_parm1;
            os << "\t" << rpv[rid2].get_name() << "\t" << new_parm2;
            os << "\t" << get_neg_log_likelihood();
            os << std::endl;
        }
    }
}


//////////////////////////
// Print table of partial likelihoods across each parameter
//////////////////////////

void ML_multi::partials_print(std::ostream& os,
                              const double lowerdelta,
                              const double upperdelta,
                              const double stepsize)
{
    os << "ratep.name" << "\t" << "val" << "\t" << "neg.log.L" << std::endl;
    RatePVector& rpv = branch_rate_manager.get_ratepvector();
    for (size_t i = 0; i < rpv.size(); ++i) {
        double old_parm = rpv[i].get_ratep();
        for (double d = lowerdelta; d < upperdelta; d+= stepsize) {
            double new_parm = old_parm + d;
            rpv[i].set_ratep(new_parm);
            compute();
            os << rpv[i].get_name();
            os << "\t" << new_parm;
            os << "\t" << get_neg_log_likelihood();
            os << std::endl;
        }
        rpv[i].set_ratep(old_parm);
    }
}


//////////////////////////
// Generate intervals for parameters in which the likelihood >=
// max_log_likelihood - diff_log_likelihood.
//////////////////////////

const ML_multi::mat_ratep_type&
    ML_multi::get_log_likelihood_intervals(const double diff_log_likelihood)
{
    compute();
    const double max_log_likelihood = get_log_likelihood();
    double       parm_log_likelihood;
    
    double       old_parm;         // MLE value of current parameter
    const double parm_fraction = 0.001;
    double       parm_adjustment;  // parm_fraction*old_parm
    int          n_adjustments;
    const int    max_n_adjustments = 10000;
    double       new_parm;         // adjusted parameter value

    RatePVector& rpv = branch_rate_manager.get_ratepvector();
    
    static mat_ratep_type intervals(rpv.size());

    for (size_t i = 0; i < rpv.size(); ++i) {
        intervals[i].resize(2);  // [0] is lower, [1] is upper

        old_parm = rpv[i].get_ratep();
        // now calculate the interval in which we'll do the profile
        parm_adjustment = old_parm * parm_fraction;

        if (DEBUG_INTERVALS1) {
            std::cout << "ML_multi::get_log_likelihood_intervals for parm "
                      << i << " val=" << old_parm << " log.L="
                      << max_log_likelihood << std::endl;
        }

        parm_log_likelihood = max_log_likelihood;
        n_adjustments = 0;
        // compute likelihoods, moving down in the interval
        for (new_parm = old_parm;
             parm_log_likelihood >= max_log_likelihood - diff_log_likelihood
             && n_adjustments <= max_n_adjustments
             && new_parm >= old_parm * interval_parm_multiple_bounds[0];
             new_parm -= parm_adjustment, ++n_adjustments) {
            rpv[i].set_ratep(new_parm);
            compute();
            parm_log_likelihood = get_log_likelihood();
            if (DEBUG_INTERVALS2) {
                std::cout << "parm " << i << " low[" << n_adjustments
                          << "]: parm=" << new_parm << " logL="
                          << parm_log_likelihood << std::endl;
            }
        }
        intervals[i][0] = new_parm;
        if (DEBUG_INTERVALS1) {
            std::cout << "parm " << i << " interval: initial=" << old_parm
                      << " max.logL=" << max_log_likelihood
                      << " low[" << n_adjustments << "]=" << new_parm
                      << " low.logL=" << parm_log_likelihood << std::endl;
        }

        // compute likelihoods, moving up in the interval
        parm_log_likelihood = max_log_likelihood;
        n_adjustments = 0;
        for (new_parm = old_parm;
             parm_log_likelihood >= max_log_likelihood - diff_log_likelihood
             && n_adjustments <= max_n_adjustments
             && new_parm <= old_parm * interval_parm_multiple_bounds[0];
             new_parm += parm_adjustment, ++n_adjustments) {
            // compute likelihoods, moving down in the interval
            rpv[i].set_ratep(new_parm);
            compute();
            parm_log_likelihood = get_log_likelihood();
            if (DEBUG_INTERVALS2) {
                std::cout << "parm " << i << " high[" << n_adjustments
                          << "]: parm=" << new_parm << " logL="
                          << parm_log_likelihood << std::endl;
            }
        }
        intervals[i][1] = new_parm;
        if (DEBUG_INTERVALS1) {
            std::cout << "parm " << i << " interval: initial=" << old_parm
                      << " max.logL=" << max_log_likelihood
                      << " high[" << n_adjustments << "]=" << new_parm
                      << " high.logL=" << parm_log_likelihood << std::endl;
        }

        rpv[i].set_ratep(old_parm);
    }
    return(intervals);
}


//////////////////////////
// Print table of profile likelihoods for each parameter, across the range
// of values from each parameter in which the likelihood >=
// max_log_likelihood - diff_log_likelihood.
//////////////////////////

void ML_multi::profile_likelihoods_print(std::ostream& os,
                                         const double  diff_log_likelihood)
{
    compute();
    //const double max_log_likelihood = get_log_likelihood();
    
    const double parm_interval_stepsize = 0.01;
    const mat_ratep_type intervals =
                     get_log_likelihood_intervals(diff_log_likelihood);

    os << "ratep.name" << "\t" << "val" << "\t" << "log.L" << std::endl;
    os.precision(15);
    RatePVector& rpv = branch_rate_manager.get_ratepvector();
    for (size_t i = 0; i < rpv.size(); ++i) {
        if (DEBUG_PROFILE) {
        }
        double old_parm = rpv[i].get_ratep();
        double parm_interval = (intervals[i][1] - intervals[i][0]) *
                               parm_interval_stepsize;
        for (double new_parm = intervals[i][0];
             new_parm < intervals[i][1] + parm_interval;
             new_parm += parm_interval) {
            rpv[i].set_ratep(new_parm);
            compute();
            os << rpv[i].get_name();
            os << "\t" << new_parm;
            os << "\t" << get_log_likelihood();
            os << std::endl;
        }
        rpv[i].set_ratep(old_parm);
    }
}


