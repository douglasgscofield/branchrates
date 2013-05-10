#ifndef ML_MULTI_QUASINEWTON_H
#define ML_MULTI_QUASINEWTON_H

#include "ML_multi.h"
#include "RatePVector.h"
#include "HistoryMultiParameter.h"

#include <cmath>
#include <limits>
#include <iostream>
#include <cassert>

#define  DEBUG_START_DFPMIN  0

class ML_multi_QuasiNewton : public ML_multi {
    protected:
        void          find_maximum();
        void          start_dfpmin();
        void          dfpmin(v_ratep_type& p,
                             const double gtol,
                             int& iter,
                             double& fret,
                             ptr_eval_func func,
                             ptr_eval_gradient_func dfunc
                             // double func(const v_ratep_type& p),
                             // void dfunc(const v_ratep_type& in,
                             //            v_ratep_type& out)
                            );

    public:
        ML_multi_QuasiNewton(const TraitMatrix& tm, PhyloTree& pt,
                             BranchRateManager& brm,
                             const double grad_tol = 0.0000001);
        double        get_neg_log_likelihood() const;
};

///////////////////////////////
/////////////////////////////// Likelihood methods
///////////////////////////////

inline ML_multi_QuasiNewton::ML_multi_QuasiNewton(const TraitMatrix& tm,
                                                  PhyloTree& pt, 
                                                  BranchRateManager& brm, 
                                                  const double grad_tol)
    : ML_multi(tm, pt, brm, grad_tol)
{
    /* empty */
}

inline double ML_multi_QuasiNewton::get_neg_log_likelihood() const {
    // return(scale_log(Likelihood::get_neg_log_likelihood()));
    return(scale_by_firstval(Likelihood::get_neg_log_likelihood()));
}


inline void ML_multi_QuasiNewton::find_maximum() {
    assign_midpoint();
    start_dfpmin();
}


inline void ML_multi_QuasiNewton::start_dfpmin()
{
    int iter;
    double fret;
    v_ratep_type parameters(branch_rate_manager.get_ratepvector().get_v_ratep());

    dfpmin(parameters, 
           get_gradient_tolerance(), 
           iter, 
           fret, 
           &ML_multi::eval_neg_log_likelihood_at,
           &ML_multi::eval_neg_log_likelihood_gradient_at
           );

    // Now set the ratepvector and Likelihood parameters to reflect
    // the minimum-finding result
    branch_rate_manager.get_ratepvector().assign(parameters);
    compute();

    if (DEBUG_START_DFPMIN) {
        std::cout << "end of start_dfpmin():";
        std::cout << " gradient_tolerance=" << get_gradient_tolerance(),
        std::cout << " iterations=" << iter;
        std::cout << " minimum=" << fret;
        std::cout << std::endl;
        std::cout << "parm\tval" << std::endl;
        for (int i = 0; i < branch_rate_manager.get_ratepvector().size(); ++i) {
            std::cout << branch_rate_manager.get_ratepvector()[i].get_name()
                      << "\t"
                      << branch_rate_manager.get_ratepvector()[i].get_ratep()
                      << std::endl;
        }
    }
}


#endif  // ML_MULTI_QUASINEWTON_H

