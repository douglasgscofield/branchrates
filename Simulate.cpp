#include "Simulate.h"

const std::string& Simulate::form_filename_sim_characters(
                                  const int sim_internal_tree,
                                  const int sim_ntaxa,
                                  const int sim_nchars,
                                  const double prior_state0,
                                  const double forw_base,
                                  const double back_relmag,
                                  const char sim_description[])
{
    static std::string filename;
    std::ostringstream buf;
    buf << "sim";
    buf << "_tree=" << sim_internal_tree;
    buf << "_tx=" << sim_ntaxa;
    buf << "_nch=" << sim_nchars;
    buf << "_forwbase=" << forw_base;
    buf << "_backrelmag=" << back_relmag;
    buf << "_" << sim_description;
    buf << "_prior=" << prior_state0;
    buf << ".traits";
    filename = buf.str();
    return(filename);
}

void Simulate::run_sim_characters(TraitMatrix& trait_matrix,
                                  const PhyloTree& tree,
                                  const int sim_internal_tree,
                                  const int sim_ntaxa,
                                  const std::vector<int>& sim_nchars,
                                  const double prior_state0,
                                  const double forw_base,
                                  const double back_relmag,
                                  const char sim_description[])
{
    for (int i = 0; i < sim_nchars.size(); ++i) {
        std::cout << "sim " << sim_nchars[i]
                  << " chars on internal tree " << sim_internal_tree;
        std::cout << " for " << sim_ntaxa << " taxa,";
        std::cout << " forw_base=" << forw_base;
        std::cout << " back_relmag=" << back_relmag;
        std::cout << " ...";
        sim_characters(trait_matrix, tree, sim_nchars[i], prior_state0);
        const std::string& filename_chars = 
                form_filename_sim_characters(sim_internal_tree,
                        sim_ntaxa, sim_nchars[i],
                        prior_state0, forw_base, back_relmag, sim_description);
        std::cout << " writing " << filename_chars << "...";
        trait_matrix.write_file(filename_chars);
        std::cout << " done." << std::endl;
    }
}

void Simulate::sim_characters(TraitMatrix& trait_matrix, 
                              const PhyloTree& tree, 
                              int num,
                              const double prior_state0)
{
    assert(trait_matrix.get_taxa() != NULL);
    assert(tree.get_num_nodes() > 0);
    assert(tree.get_num_leaves() <= (sizeof(int)*8) - 1);

    int                     num_char_config = 0x1 << tree.get_num_leaves();
    sim_characters_config   c;
    
    c.prior_prob_state0 = prior_state0;  // probability of state 0 in root
    c.num_character_states = 2;

    // set up trait_matrix with columns in ascending leaf_id order
    TaxonMatrix::taxon_id_vector_type v_tid;
    v_tid.assign(tree.get_num_leaves(), TaxonMatrix::notset);
    for (int l = 0; l < tree.get_num_leaves(); ++l) {
        TaxonMatrix::taxon_id_type tid = 
                          trait_matrix.get_taxa()->taxon_id_for_leaf_id(l);
        assert(tid != TaxonMatrix::notset);
        v_tid[l] = tid;
    }
    trait_matrix.set_taxon_ids(v_tid);

    trait_matrix.matrix.resize(0);  // clear it, if it's been used before
    trait_matrix.matrix.resize(num_char_config);

    for (int s = 0; s < num_char_config; ++s) {
        trait_matrix.matrix[s].set_id(s);
        trait_matrix.matrix[s].set_freq(0);
        trait_matrix.matrix[s].vals.assign(tree.get_num_leaves(), -1);
        for (int l = 0; l < tree.get_num_leaves(); ++l) {
            int bitpos = tree.get_num_leaves() - l - 1;
            trait_matrix.matrix[s].vals[l] = (s & (0x1 << bitpos)) >> bitpos;
        }
    }

    c.nodestates.assign(tree.get_num_nodes(), -1);
    c.leafstates.assign(tree.get_num_leaves(), -1);
    
    for (int i = 0; i < num; ++i) {
        c.nodestates[tree.get_node_id()] =
                (c.gsl.ran3() < c.prior_prob_state0) ? 0 : 1;
        sim_characters_subtree(c, &tree);
        
        // figure out which state to assign to
        int states_index = 0;
        for (int l = 0; l < tree.get_num_leaves(); ++l) {
            // this assert() is only good if valid states are 0 or 1
            assert(c.leafstates[l] == 0 || c.leafstates[l] == 1);
            int bitpos = tree.get_num_leaves() - l - 1;
            states_index |= (c.leafstates[l] << bitpos);
        }
        trait_matrix.matrix[states_index].increment_freq();
    }
}


void Simulate::sim_characters_subtree(sim_characters_config& c, 
                                      const PhyloTreeNode* node)
{
    int thisstate;
    if (node->get_is_root() == false) {
        // we're not the root, so compute a simulated state
        // if (c.gsl.ran3() < c.probtrans[node->get_node_id()]
        //                    [c.nodestates[node->get_ancestor()->get_node_id()]]
        //                    [0]) {
        if (c.gsl.ran3() < node->get_p_transition(
                     c.nodestates[node->get_ancestor()->get_node_id()], 0)) {
            thisstate = 0;
        } else {
            thisstate = 1;
        }
        c.nodestates[node->get_node_id()] = thisstate;
    }
    if (node->get_is_leaf() == true) {
        c.leafstates[node->get_leaf_id()] = thisstate;
        return;
    }
    for (int i = 0; i < node->descendants.size(); ++i) {
        sim_characters_subtree(c, node->descendants[i]);
    }
}


