#ifndef ML_MULTI_POWELL_H
#define ML_MULTI_POWELL_H

#include "ML_multi.h"
#include "RatePVector.h"
#include "HistoryMultiParameter.h"

#include <cmath>
#include <limits>
#include <iostream>
#include <cassert>

#define  DEBUG_START_POWELL  0

class ML_multi_Powell : public ML_multi {
    public:
        typedef double (ML_multi_Powell::* ptr_eval_1_func) (const double);

    private:
        int             ncom;
        v_ratep_type    pcom;
        v_ratep_type    xicom;
        ptr_eval_func   nrfunc;
        double          f1dim(const double x);
        void            mnbrak(double& ax,
                               double& bx,
                               double& cx,
                               double& fa,
                               double& fb,
                               double& fc,
                               ptr_eval_1_func f
                              );
        double          brent(const double ax,
                              const double bx,
                              const double cx,
                              ptr_eval_1_func f,
                              const double tol,
                              double& xmin
                             );
        void            shft3(double& a, double& b, double& c, const double d) const;

    protected:
        void            find_maximum();
        void            start_powell();
        void            initialize_directions(const v_ratep_type& p,
                                              mat_ratep_type& dir) const;
        void            powell(v_ratep_type& p,
                               mat_ratep_type& xi,
                               const double ftol,
                               int& iter,
                               double& fret,
                               ptr_eval_func func
                              );
        void            linmin(v_ratep_type& p,
                               v_ratep_type& xi,
                               double& fret,
                               ptr_eval_func func
                              );

    public:
        ML_multi_Powell(const TraitMatrix& tm, 
                        PhyloTree& pt,
                        BranchRateManager& brm);
};

///////////////////////////////
/////////////////////////////// Likelihood methods
///////////////////////////////

inline ML_multi_Powell::ML_multi_Powell(const TraitMatrix& tm,
                                        PhyloTree& pt, 
                                        BranchRateManager& brm) 
    : ML_multi(tm, pt, brm)
{
    /* empty */
}


inline void ML_multi_Powell::find_maximum() {
    ///////////////////////////////////////////////
    ///////////////////////////////////////////////
    ///////////////////////////////////////////////
    ///////////////////////////////////////////////
    ///////////////////////////////////////////////
    ///////////////////////////////////////////////
    ///////////////////////////////////////////////
    ///////////////////////////////////////////////
    ///////////////////////////////////////////////
    ///////////////////////////////////////////////
    assign_midpoint();
    start_powell();
}

inline void ML_multi_Powell::initialize_directions(const v_ratep_type& p,
                                                   mat_ratep_type& dir) const
{
    ///////////////////////////////////////////////
    ///////////////////////////////////////////////
    ///////////////////////////////////////////////
    ///////////////////////////////////////////////
    ///////////////////////////////////////////////
    ///////////////////////////////////////////////
    ///////////////////////////////////////////////
    ///////////////////////////////////////////////
    ///////////////////////////////////////////////
    ///////////////////////////////////////////////
}


inline void ML_multi_Powell::start_powell()
{
    int iter;
    double fret;
    v_ratep_type parameters(branch_rate_manager.get_ratepvector().get_v_ratep());
    mat_ratep_type directions;

    ///////////////////////////////////////////////
    ///////////////////////////////////////////////
    ///////////////////////////////////////////////
    ///////////////////////////////////////////////
    ///////////////////////////////////////////////
    ///////////////////////////////////////////////
    ///////////////////////////////////////////////
    ///////////////////////////////////////////////
    ///////////////////////////////////////////////
    ///////////////////////////////////////////////
    initialize_directions(parameters,
                          directions);

    powell(parameters, 
           directions,
           get_function_tolerance(),
           iter, 
           fret, 
           &ML_multi::eval_neg_log_likelihood_at
           );

    // Now set the ratepvector and Likelihood parameters to reflect
    // the minimum-finding result
    branch_rate_manager.get_ratepvector().assign(parameters);
    compute();

    if (DEBUG_START_POWELL) {
        std::cout << "end of start_powell():";
        std::cout << " function_tolerance=" << get_function_tolerance(),
        std::cout << " iterations=" << iter;
        std::cout << " minimum=" << fret;
        std::cout << std::endl;
        std::cout << "parm\tval" << std::endl;
        for (size_t i = 0; i < branch_rate_manager.get_ratepvector().size(); ++i) {
            std::cout << branch_rate_manager.get_ratepvector()[i].get_name()
                      << "\t"
                      << branch_rate_manager.get_ratepvector()[i].get_ratep()
                      << std::endl;
        }
    }
}


inline void ML_multi_Powell::shft3(double& a,
                              double& b, 
                              double& c, 
                              const double d) const {
    a = b;
    b = c;
    c = d;
}


#endif  // ML_MULTI_POWELL_H

