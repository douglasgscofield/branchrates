#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#include "PhyloTree.h"
#include "TraitMatrix.h"
#include "BranchRateManager.h"
#include "RatePVector.h"

#include <cmath>
#include <iostream>
#include <cassert>
#include <limits>

///////////////////////
// Templates in lieu of macros
///////////////////////
template<class T>
inline const T SQR(const T a) {
    return(a*a); 
}

template<class T>
inline const T MAX(const T& a, const T& b) {
    return(b > a ? (b) : (a)); 
}

template<class T>
inline void SWAP(T& a, T& b) {
    T dum = a; a = b; b = dum;
}

template<class T>
inline const T SIGN(const T& a, const T& b) {
    return(b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));
}

///////////////////////
// Debug flags common to all likelihood routines
///////////////////////
#define DEBUG_LIKELIHOOD_TRAIT        0
#define DEBUG_LIKELIHOOD_CALCS        0
#define DEBUG_LIKELIHOOD_COMPUTE      0


class Likelihood {
    public:
        typedef RatePVector::v_ratep_type 
                                 v_ratep_type;
        typedef std::vector<v_ratep_type>
                                 mat_ratep_type;

        double delta;
        
        class RatePBounds {
            public:
                double lowerbound;
                double upperbound;
                double accuracy;
                double delta;
        };
        typedef std::vector<RatePBounds>
                                 ratepbounds_vector_type;

        ratepbounds_vector_type  bounds;

    protected:
        const TraitMatrix&       trait_matrix;
        PhyloTree&               tree;
        BranchRateManager&       branch_rate_manager;
        std::vector<double>      prior_prob_root_state;
        double                   log_likelihood;
        void                     set_log_likelihood(double log_L);

    public:
        Likelihood(const TraitMatrix& tm, 
                   PhyloTree& pt,
                   BranchRateManager& brm,
                   const double prior_prob_state0 = 0.5);

        virtual void             set_prior_prob_root_state(const double prior);
        double                   determine_machine_delta() const;
        void                     compute();
        double                   compute_for_single_parameter(double parm);
        double                   L_trait(int trait) const ;
        double                   L_trait_root_state(int trait,
                                           TraitMatrix::trait_type state) const;
        double                   L_trait_node_state(int trait,
                                           const PhyloTreeNode* node,
                                           TraitMatrix::trait_type state) const;
        double                   get_log_likelihood() const;
        virtual double           get_neg_log_likelihood() const;
        double                   get_likelihood() const;
        double                   get_single_parameter() const;

        double                   scale_none(double arg) const;
        double                   scale_log(double arg) const;

        void                     print(std::ostream& os = std::cout) const;
        friend std::ostream&     operator<<(std::ostream& os,
                                            const Likelihood& lhood);
};

///////////////////////////////
/////////////////////////////// Likelihood methods
///////////////////////////////

inline double Likelihood::get_single_parameter() const {
    RatePVector& rpv = branch_rate_manager.get_ratepvector();
    assert(rpv.size() == 1);
    return(rpv[0].get_ratep());
}

inline double Likelihood::compute_for_single_parameter(double parm) {
    branch_rate_manager.set_unidir_rate(parm);
    compute();
    return(get_log_likelihood());
}

inline double Likelihood::scale_log(double arg) const {
    return(std::log(arg));
}

inline double Likelihood::scale_none(double arg) const {
    return(arg);
}

inline double Likelihood::get_neg_log_likelihood() const {
    return(-get_log_likelihood());
}

inline double Likelihood::get_log_likelihood() const {
    return(log_likelihood);
}

inline void Likelihood::set_log_likelihood(double log_L) {
    log_likelihood = log_L;
}

inline double Likelihood::get_likelihood() const {
    return(std::exp(get_log_likelihood()));
}

inline void Likelihood::set_prior_prob_root_state(const double prior) {
    prior_prob_root_state.resize(TraitMatrix::max_num_states);
    prior_prob_root_state[0] = prior;
    double step = (1.0 - prior) / double(prior_prob_root_state.size() - 1);
    for (int i = 1; i < prior_prob_root_state.size(); ++i) {
        prior_prob_root_state[i] = step;
    }
}

inline Likelihood::Likelihood(const TraitMatrix& tm, PhyloTree& pt,
                       BranchRateManager& brm,
                       const double prior_prob_state0) :
                       delta(determine_machine_delta()),
                       trait_matrix(tm), tree(pt), branch_rate_manager(brm),
                       log_likelihood(-99999.99999)
{

    // establish prior probabilities of trait states at the root node
    
    set_prior_prob_root_state(prior_prob_state0);
    // prior_prob_root_state[0] = 0.99872;
    // prior_prob_root_state[1] = 0.00128;
    // for (int i = 0; i < prior_prob_root_state.size(); ++i) {
    //     // this scheme assigned equal probabilities to each state
    //     double total_prob = 1.0;
    //     prior_prob_root_state[i] = total_prob /
    //                                    double(prior_prob_root_state.size());
    // }
}


#endif  // LIKELIHOOD_H

