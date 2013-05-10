#ifndef ML_SINGLE_H
#define ML_SINGLE_H

#include "Likelihood.h"
#include "RatePVector.h"
#include "HistorySingleParameter.h"

#include <cmath>
#include <iostream>
#include <cassert>

class ML_single : public Likelihood {
    protected:
        virtual void             find_maximum() = 0;
        virtual double           adjust() = 0;
        virtual bool             converged();

        HistorySingleParameter   history;

    public:
        ML_single(const TraitMatrix& tm, PhyloTree& pt,
                  BranchRateManager& brm);

        virtual void             maximize(const double lower,
                                          const double upper,
                                          const double accuracy);
};

///////////////////////////////
/////////////////////////////// Likelihood methods
///////////////////////////////

inline ML_single::ML_single(const TraitMatrix& tm,
                            PhyloTree& pt, 
                            BranchRateManager& brm)
    : Likelihood(tm, pt, brm)
{
    /* empty */
}

inline void ML_single::maximize(const double lb, const double ub,
                                const double acc) {
    int num_ratep = branch_rate_manager.get_ratepvector().size();
    assert(num_ratep == 1);
    bounds.resize(num_ratep);
    bounds[0].lowerbound = lb;
    bounds[0].upperbound = ub;
    bounds[0].accuracy = acc;
    bounds[0].delta = determine_machine_delta();
    find_maximum();
}


#endif  // ML_SINGLE_H


