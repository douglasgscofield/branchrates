#ifndef ML_MULTI_H
#define ML_MULTI_H

#include "Likelihood.h"
#include "RatePVector.h"
#include "HistoryMultiParameter.h"

#include <cmath>
#include <iostream>
#include <cassert>

#define  DEBUG_BOUNDS        0
#define  DEBUG_INTERVALS1    1
#define  DEBUG_INTERVALS2    0
#define  DEBUG_PROFILE       1

class ML_multi : public Likelihood {
    public:
        static const double    chisq_095_1df;
        static const double    lrt_095_interval;
        static const double    lrt_profile_interval;
        static const double    interval_parm_multiple_bounds[2];
        static const double    interval_parm_steps;

            // function for evaluating likelihood, and a pointer to same
        double          eval_neg_log_likelihood_at(const v_ratep_type& p);
        double          bounded_eval_neg_log_likelihood_at(const v_ratep_type& p);
        typedef double (ML_multi::* ptr_eval_func) (const v_ratep_type&);

        // function for evaluating gradient, and a pointer to same
        void            eval_neg_log_likelihood_gradient_at(
                                 const v_ratep_type& in,
                                 v_ratep_type& out);
        typedef void (ML_multi::* ptr_eval_gradient_func) (const v_ratep_type&,
                                                           v_ratep_type&);

    protected:
        virtual void    find_maximum() = 0;

        double          gradient_tolerance;
        double          function_tolerance;
        double          scale_firstval;
        virtual void    assign_midpoint();
        void            bounds_adjust(v_ratep_type& p, int& n_adjusted) const;
        void            bounds_trace(const v_ratep_type& p,
                                     const std::string& loc = "bounds_trace") const;
        double          get_scale_firstval() const;
        void            set_scale_firstval(double arg);
        double          scale_by_firstval(double arg) const;
        // HistoryMultiParameter   history;


    protected:
        void            lnsrch(const v_ratep_type& xold,
                               const double fold,
                               const v_ratep_type& g,
                               v_ratep_type& p,
                               v_ratep_type& x,
                               double& f,
                               const double stpmax,
                               bool& check,
                               ptr_eval_func func
                              );

    public:
        ML_multi(const TraitMatrix& tm, PhyloTree& pt,
                 BranchRateManager& brm,
                 const double prior_prob_state0 = 0.5,
                 double grad_tol = 0.0000001);

        virtual void    maximize(const v_ratep_type& lb,
                                 const v_ratep_type& ub,
                                 const double acc = 0.0000001);
        virtual void    maximize(const int n_parameters,
                                 const double lb,
                                 const double ub,
                                 const double acc = 0.0000001);

        virtual void    partials_print(std::ostream& os = std::cout,
                                       const double lowerdelta = -0.01,
                                       const double upperdelta = +0.01,
                                       const double stepsize = 0.001);

        const mat_ratep_type&
                        get_log_likelihood_intervals(
                                 const double diff_log_likelihood =
                                                             lrt_095_interval);

        void            profile_likelihoods_print(
                                 std::ostream& os = std::cout,
                                 const double diff_log_likelihood =
                                                         lrt_profile_interval);

        virtual void    surface_range_print(const std::string& ratep1,
                                            const std::string& ratep2,
                                            const double lower,
                                            const double upper,
                                            const double stepsize,
                                            std::ostream& os = std::cout);

        virtual void    surface_range_print(const int rid1,
                                            const int rid2,
                                            const double lower,
                                            const double upper,
                                            const double stepsize,
                                            std::ostream& os = std::cout);

        virtual void    surface_print(const std::string& ratep1,
                                      const std::string& ratep2,
                                      std::ostream& os = std::cout,
                                      const double lowerdelta = -0.01,
                                      const double upperdelta = +0.01,
                                      const double stepsize = 0.001);

        virtual void    surface_print(const int rid1,
                                      const int rid2,
                                      std::ostream& os = std::cout,
                                      const double lowerdelta = -0.01,
                                      const double upperdelta = +0.01,
                                      const double stepsize = 0.001);

        virtual void    set_bounds(const int n_parameters,
                                   const double lb,
                                   const double ub,
                                   const double del,
                                   const double acc = 0.0000001);
        double          get_gradient_tolerance() const;
        void            set_gradient_tolerance(double arg);
        double          get_function_tolerance() const;
        void            set_function_tolerance(double arg);
};

///////////////////////////////
/////////////////////////////// ML_multi methods
///////////////////////////////

inline ML_multi::ML_multi(const TraitMatrix& tm,
                            PhyloTree& pt, 
                            BranchRateManager& brm,
                            const double prior_prob_state0,
                            const double grad_tol)
    : Likelihood(tm, pt, brm, prior_prob_state0), gradient_tolerance(grad_tol),
      scale_firstval(-9999999.9999999)
{
    /* empty */
}

inline void ML_multi::set_bounds(const int n_parameters,
                                 const double lb,
                                 const double ub,
                                 const double del,
                                 const double acc)
{
    for (int i = 0; i < n_parameters; ++i) {
        bounds[i].lowerbound = lb;
        bounds[i].upperbound = ub;
        bounds[i].accuracy = acc;
        bounds[i].delta = del;
    }
}

inline double ML_multi::scale_by_firstval(double arg) const {
    return(arg / get_scale_firstval());
}

inline void ML_multi::maximize(const v_ratep_type& lb,
                               const v_ratep_type& ub,
                               const double acc)
{
    assert(lb.size() == ub.size());
    int num_ratep = branch_rate_manager.get_ratepvector().size();
    assert(num_ratep == lb.size());
    assert(num_ratep > 1);
    double delta = determine_machine_delta();
    bounds.resize(num_ratep);
    for (int i = 0; i < num_ratep; ++i) {
        bounds[i].lowerbound = lb[i];
        bounds[i].upperbound = ub[i];
        bounds[i].accuracy = acc;
        bounds[i].delta = delta;
    }
    find_maximum();
}


inline void ML_multi::maximize(const int n_parameters,
                               const double lb, 
                               const double ub,
                               const double acc) {
    int num_ratep = branch_rate_manager.get_ratepvector().size();
    assert(n_parameters == num_ratep);
    assert(num_ratep > 0);
    // assert(num_ratep > 1);
    double delta = determine_machine_delta();
    bounds.resize(num_ratep);
    set_bounds(num_ratep, lb, ub, delta, acc);
    find_maximum();
}

inline void ML_multi::assign_midpoint() {
    // establish starting points for all parameters
    assert(branch_rate_manager.get_ratepvector().get_initialized() == false);

    int n_parm = bounds.size();
    v_ratep_type initials(n_parm);
    for (int i = 0; i < n_parm; ++i) {
        initials[i] = (bounds[i].upperbound + bounds[i].lowerbound) / 2.0;
    }
    branch_rate_manager.get_ratepvector().assign(initials);
}

inline double ML_multi::get_gradient_tolerance() const {
    return(gradient_tolerance);
}


inline void ML_multi::set_gradient_tolerance(double arg) {
    gradient_tolerance = arg;
}


inline double ML_multi::get_function_tolerance() const {
    return(function_tolerance);
}


inline void ML_multi::set_function_tolerance(double arg) {
    function_tolerance = arg;
}


inline double ML_multi::get_scale_firstval() const {
    return(scale_firstval);
}


inline void ML_multi::set_scale_firstval(double arg) {
    scale_firstval = arg;
}


#endif  // ML_MULTI_H


