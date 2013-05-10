#include "Likelihood.h"

///////////////////////////
// Compute the likelihood of a single state of a particular trait across all
// taxa, for a subtree starting at a particular node.  Felsenstein (2004) calls
// this a conditional likelihood for the subtree, in that it's conditioned on
// the state of interest.
//
// If we're at a leaf node, the likelihood is simply the probability that the
// trait state at the leaf node is the same as the state of interest.  This
// will be 1.0 if true or 0.0 if false.
//
// If we're an internal node, then we proceed to cycle through all possible
// states that each of our descendant nodes can carry, compute the conditional
// likelihood for each of those state-descendant combinations via recursive
// calls, weight them by the probability of a transition from this node's state
// of interest to each of the possible descendant states, and then return that
// likelihood.
///////////////////////////
double Likelihood::L_trait_node_state(int trait, const PhyloTreeNode* node,
                                      TraitMatrix::trait_type state) const {
    if (DEBUG_LIKELIHOOD_CALCS) {
        std::clog << "inside Likelihood:L_trait_node_state("
                  << "trait=" << trait << ", node_id=" << node->get_node_id()
                  << " (" << ((void*)node) << ")" << ", state=" << state << ")"
                  << std::endl;
    }

    double thistrait_node_state_L;
    
    // Leaf node

    if (node->get_is_leaf() == true) {
        // int col = trait_matrix.map_taxon_id_to_col(node->get_taxon_id());
        // Must assure that .vals[node->get_leaf_id()] is the same column that
        // is found by
        // .vals[trait_matrix.map_taxon_id_to_col(node->get_taxon_id())]
        if (trait_matrix.matrix[trait].vals[node->get_leaf_id()] == state) {
            thistrait_node_state_L = 1.0;
        } else {
            thistrait_node_state_L = 0.0;
        }
        if (DEBUG_LIKELIHOOD_CALCS) {
            std::clog << "   Likelihood to be returned="
                      << thistrait_node_state_L << std::endl;
        }
        return(thistrait_node_state_L);
    }
    
    // Internal node

    thistrait_node_state_L = 1.0;
    for (int d = 0; d < node->descendants.size(); ++d) {
        double partial_sum = 0.0;
        for (int s = 0; s < TraitMatrix::max_num_states; ++s) {
            double p_trans = node->descendants[d]->get_p_transition(state, s);
            double descendant_L = L_trait_node_state(trait, 
                                                     node->descendants[d], s);
            partial_sum += (p_trans * descendant_L);
        }
        thistrait_node_state_L *= partial_sum;
    }

    if (DEBUG_LIKELIHOOD_CALCS) {
        std::clog << "   Likelihood to be returned="
                  << thistrait_node_state_L << std::endl;
    }
    return(thistrait_node_state_L);
}


///////////////////////////
// Compute the conditional likelihood of a single state of a single trait
// across all taxa, for the entire tree starting at the root.
///////////////////////////
double Likelihood::L_trait_root_state(int trait,
                                      TraitMatrix::trait_type state) const {
    if (DEBUG_LIKELIHOOD_CALCS) {
        std::clog << "inside Likelihood:L_trait_root_state(trait="
                  << trait << ", state=" << state << ")" << std::endl;
    }
    
    double thistrait_root_state_L;
    
    thistrait_root_state_L = L_trait_node_state(trait, &tree, state);

    if (DEBUG_LIKELIHOOD_CALCS) {
        std::clog << "   Likelihood to be returned=" << thistrait_root_state_L
                  << std::endl;
    }
    return(thistrait_root_state_L);
}


///////////////////////////
// Compute the likelihood of a particular trait across all taxa,
// for the entire tree.  Each possible trait state is weighted at
// the root by the values in Likelihood.prior_prob_root_state[state].
// This is the method called to compute a likelihood value for
// a single trait.
///////////////////////////
double Likelihood::L_trait(int trait) const {
    if (DEBUG_LIKELIHOOD_TRAIT || DEBUG_LIKELIHOOD_CALCS) {
        std::clog << "inside Likelihood::L_trait(trait=" << trait << ")"
                  << std::endl;
    }

    double thistrait_L = 0.0;
    
    for (int state = 0; state < TraitMatrix::max_num_states; ++state) {
        double thistrait_state_root_L = prior_prob_root_state[state] *
                                        L_trait_root_state(trait, state);
        thistrait_L += thistrait_state_root_L;
    }

    if (DEBUG_LIKELIHOOD_CALCS) {
        std::clog << "   Likelihood to be returned=" << thistrait_L << std::endl;
    }
    return(thistrait_L);
}


///////////////////////////
// Compute the likelihood of the trait_matrix given the tree and transition
// probabilities.  This is the method called to compute a likelihood value.
///////////////////////////
void Likelihood::compute() {
    if (DEBUG_LIKELIHOOD_COMPUTE || DEBUG_LIKELIHOOD_CALCS) {
        std::clog << "inside Likelihood::compute()" << std::endl;
    }
    if (trait_matrix.matrix.size() == 0) {
        std::cerr << "Likelihood::compute():  No characters" << std::endl;
        exit(1);
    }
    
    double log_L = 0.0;
    double thistrait_root_log_L;

    branch_rate_manager.calculate_p_transition();

    for (int trait = 0; trait < trait_matrix.matrix.size(); ++trait) {

        if (trait_matrix.matrix[trait].get_freq() > 0) {

            thistrait_root_log_L = trait_matrix.matrix[trait].get_freq() *
                                   std::log(L_trait(trait));
                            
            log_L += thistrait_root_log_L;
        }
    }
    
    if (DEBUG_LIKELIHOOD_COMPUTE || DEBUG_LIKELIHOOD_CALCS) {
        std::clog << "   Likelihood computed=" << log_L << std::endl;
    }
    set_log_likelihood(log_L);
}


///////////////////////////
// Utility functions including Likelihood::print() and ostream& operator<<.
///////////////////////////
double Likelihood::determine_machine_delta() const {
    double delta = -999999.999999;
    if (std::numeric_limits<double>::is_specialized == true) {
        delta = std::sqrt(std::numeric_limits<double>::epsilon());
        // delta = std::pow(std::numeric_limits<double>::epsilon(), 0.5);
        // delta = std::pow(std::numeric_limits<double>::epsilon(), 0.66666666667);
    } else {
        delta = 1.49012e-08;
    }
    return(delta);
}

void Likelihood::print(std::ostream& os) const {
    os << "log_likelihood=" << log_likelihood;
    os << " trait_matrix=" << ((void*)&trait_matrix);
    os << " tree=" << ((void*)&tree);
    os << " branch_rate_manager=" << ((void*)&branch_rate_manager);
    os << " prior_prob_root_state[0.." << (prior_prob_root_state.size()-1)
       << "]={";
    for (int i = 0; i < prior_prob_root_state.size(); ++i) {
        os << " " << prior_prob_root_state[i];
    }
    os << "}";
    os << std::endl;
    RatePVector& rpv = branch_rate_manager.get_ratepvector();
    os << "number of parameters=" << rpv.size();
    os << std::endl;
    rpv.print(os);
}



std::ostream& operator<<(std::ostream& os, const Likelihood& lhood) {
    os << "Likelihood:";
    lhood.print(os);
    os << std::endl;
}


