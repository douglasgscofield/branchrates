#include "RatePVector.h"
#include "RatePMap.h"
#include "PhyloTree.h"
#include "TaxonMatrix.h"
#include "TraitMatrix.h"
#include "BranchRateManager.h"
#include "ML_single_NewtonRaphson.h"
#include "ML_multi_QuasiNewton.h"
#include "ML_multi_DownhillSimplex.h"
#include "Simulate.h"

#include <string>
#include <iostream>
#include <fstream>

extern void sim_chars(TraitMatrix& trait_matrix, const PhyloTree& tree, int num);

int main() {
    TaxonMatrix         taxon_matrix;
    PhyloTree           tree;
    BranchRateManager   branch_rate_manager;
    RatePMap            ratep_map;
    int                 restarts = 50;

    // establish taxa
    // taxon_matrix.read_internal_data();
    taxon_matrix.read_file("koonin.taxa");
    // exit(0);
    
    // establish tree
    // since I haven't implemented a tree reader yet (I just haven't gotten
    // around to it yet, it's a trivial recursive problem), load the tree
    // from its internal representation.  internal tree 1 is the tree
    // found in the file koonin.tree, but with branchlengths multiplied
    // by 100
    tree.read_internal_tree(1);  // internal tree 1 is the koonin.tree
    tree.set_taxon_matrix(taxon_matrix);
    tree.set_branchratemanager(branch_rate_manager);
    tree.allocate_branchrate_ids();
    // std::cout << tree;
    // std::cout << taxon_matrix;

    TraitMatrix trait_matrix;
    trait_matrix.set_taxa(taxon_matrix);
    trait_matrix.read_file("koonin.traits");
    // std::cout << trait_matrix;
    // exit(0);

    // map rate parameters to branches
    // ratep_map.read_file("koonin.map");
    // ratep_map.read_file("koonin.all1rate.map");
    ratep_map.read_file("koonin.hybrid.map");
    // ratep_map.read_file("koonin.2rate.map");
    // ratep_map.read_file("koonin.1rate.map");
    std::cout << ratep_map;
    branch_rate_manager.allocate_ratep_from_map(ratep_map);
    // branch_rate_manager.read_ratep_init_file("koonin.hybrid.map.init");
    branch_rate_manager.read_ratep_init_file("testing.init");
    std::cout << tree;
    // exit(0);


    // now compute likelihood of tree

    //////////////////////
    // QuasiNewton
    // ML_multi_QuasiNewton ML(trait_matrix, tree, branch_rate_manager);
    // ML.maximize(branch_rate_manager.get_ratepvector().size(),
    //             0.0000000000000001,
    //             0.1);
    // std::cout << ML;
    
    /////////////////////
    // DownhillSimplex
    ML_multi_DownhillSimplex ML(trait_matrix, tree, branch_rate_manager);
    ML.set_prior_prob_root_state(0.99);
    restarts = 15;
    ML.maximize(branch_rate_manager.get_ratepvector().size(),  // num params
                0.000000000000000000001,   // lower bound of param range
                0.1,                       // upper bound of param range
                0.0000001,                 // accuracy
                0.000001,                  // function tolerance (see code)
                0.00001,                   // scale length (see code)
                restarts);                 // number of restarts
    std::cout << ML;
    ML.profile_likelihoods_print(std::cout);


    //ML.surface_range_print("ratef1", "rateb1", 0.001, 0.1, 0.002);


    /////////////////////
    // single_NewtonRaphson
//    if (branch_rate_manager.get_ratepvector().size() == 1) {
//        ML_single_NewtonRaphson ML_s(trait_matrix, tree, branch_rate_manager);
//        ML_s.maximize(0.0000000000000001, 0.001, 0.00000000001);
//        std::cout << "ML_s.maximize() stops at " << ML_s.get_single_parameter()
//                  << " log-L=" << ML_s.get_log_likelihood() << std::endl;
//    }
//
    // ML.partials_print();
    // ML.surface_range_print("rate6", "rate1", 0.01, 0.2, 0.01);

    // ML.maximize(0.001, 0.2, 0.00000001);
    // std::cout << "ML.maximize() stops at " << ML.get_single_parameter()
    //           << std::endl;

}



    // simulate characters

    // branch_rate_manager.set_bidir_rates(fr, br);
//    RatePVector& rpv = branch_rate_manager.get_ratepvector();
//    double fr0 = 0.05;
//    double br0 = 0.025;
//    double fr1 = 0.04;
//    double br1 = 0.015;
//    double fr2 = 0.08;
//    double br2 = 0.015;
//    int rid = rpv.get_ratep_id("rate1");
//    rpv[rid].set_ratep(fr0);
//    rid = rpv.get_ratep_id("rate2");
//    rpv[rid].set_ratep(br0);
//    rid = rpv.get_ratep_id("rate3");
//    rpv[rid].set_ratep(fr1);
//    rid = rpv.get_ratep_id("rate4");
//    rpv[rid].set_ratep(br1);
//    rid = rpv.get_ratep_id("rate5");
//    rpv[rid].set_ratep(fr2);
//    rid = rpv.get_ratep_id("rate6");
//    rpv[rid].set_ratep(br2);
//    branch_rate_manager.calculate_p_transition();
//    int num_sim_reps = 200000;
//    sim_chars(trait_matrix, tree, num_sim_reps);
    // std::cout << trait_matrix;
    // std::cout << tree;
    // std::exit(0);

//
//    Likelihood::v_ratep_type parms(2);
//    Likelihood::v_ratep_type grad(2);
//    std::cout << "br\tlogL\tgrad0\tgrad1" << std::endl;
//    for (double br = 0.015; br < 0.035; br += 0.0005) {
//        ////branch_rate_manager.set_unidir_rate(br);
//        //branch_rate_manager.set_bidir_rates(0.75, br);
//        //ML.compute();
//        //std::cout << br << "\t" << ML.get_log_likelihood() << std::endl;
//        parms[0] = 0.05;
//        parms[1] = br;
//        double neglogL = ML.eval_neg_log_likelihood_at(parms);
//        ML.eval_neg_log_likelihood_gradient_at(parms, grad);
//        std::cout << br << "\t" << neglogL << "\t" << grad[0] << "\t" << grad[1];
//        std::cout << std::endl;
//    }


