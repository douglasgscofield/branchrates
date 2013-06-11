#ifndef ML_MULTI_DOWNHILLSIMPLEX_H
#define ML_MULTI_DOWNHILLSIMPLEX_H

// This implementation of the Downhill Simplex method is due to
// Numerical Recipes in C++, for which I have purchased a single-copy
// license and am therefore allowed to distribute in compiled form in a
// noncommercial program.  The Gnu Scientific Library includes an
// implementation of the Downhill Simplex method, and I intend to switch
// to that so that I can distribute source.

#include "ML_multi.h"
#include "RatePVector.h"
#include "HistoryMultiParameter.h"

#include <cmath>
#include <limits>
#include <iostream>
#include <cassert>

#define  CONFIG_BOUNDS_ADJUST         1
#define  CONFIG_DIE_ON_NMAX_EXCEEDED  0

#define  DEBUG_FIND_MAXIMUM           1
#define  DEBUG_START_AMOEBA           1
#define  DEBUG_AMOEBA                 1
#define  DEBUG_MONITOR_X10            1
#define  DEBUG_AMOTRY                 0
#define  DEBUG_BOUNDS_TRACE           0
#define  DEBUG_PRINT_DATA_STRUCTURES  0

class ML_multi_DownhillSimplex : public ML_multi {
    public:
        virtual void   print_data_structures() const;
        virtual void   print_data_structures(const mat_ratep_type& p,
                                             const v_ratep_type& y) const;

    protected:
        mat_ratep_type simplex;
        v_ratep_type   funk_vals;

        long           NMAX;
        int            num_restart;
        const int      max_restart;
        double         scale_length;

        double         get_scale_length() const;
        void           set_scale_length(double arg);
        long           get_NMAX() const;
        void           set_NMAX(long arg);
        void           find_maximum();
        void           start_amoeba(int& nfunk);
        void           initialize_simplex(const v_ratep_type& p0,
                                          const double s_len,
                                          ptr_eval_func funk);
        void           get_psum(const mat_ratep_type& p,
                                v_ratep_type& psum);
        void           amoeba(mat_ratep_type& p, 
                              v_ratep_type& y,
                              const double ftol,
                              ptr_eval_func funk,
                              int& nfunk
                             );
        double         amotry(mat_ratep_type& p, 
                              v_ratep_type& y,
                              v_ratep_type& psum,
                              ptr_eval_func funk,
                              const int ihi,
                              const double fac
                             );
        /*****
        void           amoeba(ptr_eval_func funk,
                              int& nfunk
                             );
        double         amotry(v_ratep_type& psum,
                              ptr_eval_func funk,
                              const int ihi,
                              const double fac
                             );
        *****/

    public:
        ML_multi_DownhillSimplex(const TraitMatrix& tm, PhyloTree& pt,
                                 BranchRateManager& brm);
        void           maximize(const v_ratep_type& lb,
                                const v_ratep_type& ub,
                                const double acc,
                                const double f_tol = 0.0001,
                                const double s_len = 0.01,
                                const int n_restart = 0);
        void           maximize(const int n_parameters,
                                const double lb,
                                const double ub,
                                const double acc,
                                const double f_tol = 0.0001,
                                const double s_len = 0.01,
                                const int n_restart = 0);
};

///////////////////////////////
/////////////////////////////// Likelihood methods
///////////////////////////////

inline ML_multi_DownhillSimplex::ML_multi_DownhillSimplex(const TraitMatrix& tm,
                                                          PhyloTree& pt, 
                                                          BranchRateManager& brm)
                                                          // const double grad_tol)
    : ML_multi(tm, pt, brm), NMAX(10000), max_restart(50)
{
    if (branch_rate_manager.get_ratepvector().size() < 2) {
        std::cerr << "ML_multi_DownhillSimplex: at least 2 parameters required"
                  << std::endl;
        assert(branch_rate_manager.get_ratepvector().size() >= 2);
    }
}


inline void ML_multi_DownhillSimplex::maximize(const v_ratep_type& lb,
                                               const v_ratep_type& ub,
                                               const double acc,
                                               const double f_tol,
                                               const double s_len,
                                               const int n_restart)
{
    assert(lb.size() == ub.size());
    int num_ratep = branch_rate_manager.get_ratepvector().size();
    double delta = determine_machine_delta();
    assert(num_ratep == (int)lb.size());
    assert(num_ratep > 1);
    set_function_tolerance(f_tol);
    set_scale_length(s_len);
    bounds.resize(num_ratep);
    for (int i = 0; i < num_ratep; ++i) {
        bounds[i].lowerbound = lb[i];
        bounds[i].upperbound = ub[i];
        bounds[i].accuracy = acc;
        bounds[i].delta = delta;
    }
    num_restart = n_restart ? n_restart : max_restart;
    find_maximum();
}


inline void ML_multi_DownhillSimplex::maximize(const int n_parameters,
                                               const double lb, 
                                               const double ub,
                                               const double acc,
                                               const double f_tol,
                                               const double s_len,
                                               const int n_restart)
{
    int num_ratep = branch_rate_manager.get_ratepvector().size();
    assert(n_parameters == num_ratep);
    assert(num_ratep > 0);
    // assert(num_ratep > 1);
    double delta = determine_machine_delta();
    set_function_tolerance(f_tol);
    set_scale_length(s_len);
    bounds.resize(num_ratep);
    set_bounds(num_ratep, lb, ub, delta, acc);
    num_restart = n_restart ? n_restart : max_restart;
    find_maximum();
}


inline void ML_multi_DownhillSimplex::find_maximum() {
    int nfunk = 0;  // number of function evaluations, set within start_amoeba()
    int restart = 0;
    if (DEBUG_FIND_MAXIMUM)
        std::cout << "find_maximum(): " << num_restart << " restarts" << std::endl;
    if (branch_rate_manager.get_ratepvector().get_initialized() == false) {
        assign_midpoint();
    }
    compute();
    set_scale_firstval(Likelihood::get_neg_log_likelihood());
    v_ratep_type p0(branch_rate_manager.get_ratepvector().get_v_ratep());
    while(++restart <= num_restart && nfunk != 1) {
        if (DEBUG_FIND_MAXIMUM) std::cout << " initializing simplex...";
        initialize_simplex(p0,
                           get_scale_length(),
                           &ML_multi::bounded_eval_neg_log_likelihood_at);
        if (DEBUG_FIND_MAXIMUM) std::cout << " start_amoeba()..." << std::endl;
        start_amoeba(nfunk);
        if (nfunk >= get_NMAX()) {
            if (DEBUG_FIND_MAXIMUM)
                std::cout << "start_amoeba() terminated due to exceeding NMAX"
                          << std::endl;
        }
        if (DEBUG_FIND_MAXIMUM)
            std::cout << "DONE start_amoeba().  restarts remaining="
                      << num_restart - restart << std::endl << std::endl;
        p0 = branch_rate_manager.get_ratepvector().get_v_ratep();
    }
    if (DEBUG_FIND_MAXIMUM)
        std::cout << "number of restarts performed: "
                  << restart << std::endl;
}


inline double ML_multi_DownhillSimplex::get_scale_length() const {
    return(scale_length);
}


inline void ML_multi_DownhillSimplex::set_scale_length(double arg) {
    scale_length = arg;
}


inline long ML_multi_DownhillSimplex::get_NMAX() const {
    return(NMAX);
}


inline void ML_multi_DownhillSimplex::set_NMAX(long arg) {
    NMAX = arg;
}


#endif  // ML_MULTI_DOWNHILLSIMPLEX_H

