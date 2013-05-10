#ifndef SIMULATE_H
#define SIMULATE_H

#include "TraitMatrix.h"
#include "TaxonMatrix.h"
#include "PhyloTree.h"
#include "GSL.h"

#include <cassert>
#include <vector>
#include <string>
#include <iostream>


class Simulate {

    ////////////////////////////////////// beginning sim_characters()
    //
    // Data structures and methods for simulation of character state vectors.
    // Instantiate a Simulate object, then use sim_characters().  It clears the
    // trait matrix in TraitMatrix.matrix via resize(0) prior to filling it, so
    // the same TraitMatrix can be reused.  It does not clear the taxon
    // information or anything else, only the trait values.

    private:
        class sim_characters_config {
            public:
                GSL         gsl;
                double      prior_prob_state0;
                int         num_character_states;
                int         cur_node;
                std::vector<TraitMatrix::trait_type>
                            nodestates;
                std::vector<TraitMatrix::trait_type>
                            leafstates;
        };
        void                sim_characters_subtree(
                                    sim_characters_config& c,
                                    const PhyloTreeNode* node);
    public:
        void                run_sim_characters(TraitMatrix& trait_matrix,
                                               const PhyloTree& tree,
                                               const int sim_internal_tree,
                                               const int sim_ntaxa,
                                               const std::vector<int>&
                                                       sim_nchars,
                                               const double prior_state0,
                                               const double forw_base,
                                               const double back_relmag,
                                               const char sim_description[]);
        void                sim_characters(
                                    TraitMatrix& trait_matrix, 
                                    const PhyloTree& tree, 
                                    int num,
                                    const double prior_state0);
        static const std::string& 
                            form_filename_sim_characters(
                                    const int sim_internal_tree,
                                    const int sim_ntaxa,
                                    const int sim_nchars,
                                    const double prior_state0,
                                    const double forw_base,
                                    const double back_relmag,
                                    const char sim_description[]);

    //
    ////////////////////////////////////// end sim_characters()
};

#endif // SIMULATE_H


