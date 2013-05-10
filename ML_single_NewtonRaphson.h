#ifndef ML_SINGLE_NEWTONRAPHSON_H
#define ML_SINGLE_NEWTONRAPHSON_H

#include "ML_single.h"

#include <cmath>
#include <iostream>
#include <cassert>

///////////////////////
// Debug flags 
///////////////////////
#define MAXDEBUG                         0
#define MAXDEBUG_ADJUSTSINGLE            0
#define MAXDEBUG_HISTORY_AT_CONVERGENCE  1


class ML_single_NewtonRaphson : public ML_single {
    protected:
        void                     find_maximum();
        double                   adjust();

    public:
        ML_single_NewtonRaphson(const TraitMatrix& tm, PhyloTree& pt,
                                BranchRateManager& brm);
};

///////////////////////////////
/////////////////////////////// Likelihood methods
///////////////////////////////

inline ML_single_NewtonRaphson::ML_single_NewtonRaphson(const TraitMatrix& tm,
                                                    PhyloTree& pt, 
                                                    BranchRateManager& brm)
    : ML_single(tm, pt, brm)
{
    /* empty */
}

#endif  // ML_NEWTONRAPHSON_SINGLE_H


