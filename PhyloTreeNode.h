#ifndef PHYLOTREENODE_H
#define PHYLOTREENODE_H

#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <cassert>

#include "TaxonMatrix.h"
#include "BranchRateManager.h"

//////////////////////
// Below these comments are data structures which define an internal tree
// that is used to establish a tree for testing purposes
//////////////////////
#define   MAX_init_tree_descendants   2
#define   indent_block                "|   "

struct s_init_tree {
    int          nodenum;
    int          ancestor;
    int          descendant[MAX_init_tree_descendants];
    const char*  name;
    double       branchlength;  /* length of branch to ancestor */
    int          is_leaf;
};

static std::vector<s_init_tree*> internal_trees;

/************************************
static s_init_tree rogozin_tree[] = {
    {0,  -1, { 1,  2}, "Root",          0.0,      0},
    {1,   0, {-1, -1}, "Arth",          0.1,      1},
    {2,   0, { 3,  4}, "IntNonPlant",   0.293312, 0},
    {3,   2, { 5,  6}, "IntFungiPlas",  0.099932, 0},
    {4,   2, { 7,  8}, "IntMetazoa",    0.042289, 0},
    {5,   3, {-1, -1}, "Scpo",          0.108095, 1},
    {6,   3, { 9, 10}, "IntSacPlas",    0.133098, 0},
    {7,   4, {-1, -1}, "Cael",          0.158298, 1},
    {8,   4, {11, 12}, "IntArthHuman",  0.075582, 0},
    {9,   6, {-1, -1}, "Plfa",          0.610389, 1},
    {10,  6, {-1, -1}, "Sace",          0.097062, 1},
    {11,  8, {-1, -1}, "Hosa",          0.187505, 1},
    {12,  8, {13, 14}, "IntArth",       0.060409, 0},
    {13, 12, {-1, -1}, "Anga",          0.096329, 1},
    {14, 12, {-1, -1}, "Drme",          0.088562, 1},
    {-1, -1, {-1, -1}, NULL,           -1.0,      0},
};
************************************/

/************************************/
static s_init_tree internaltree0_small[] = {
    {  0, -1, { 1,  2}, "root",  0.0, 0},
    {  1,  0, {-1, -1}, "a",     4.0, 1},
    {  2,  0, { 3,  4}, "int1",  1.0, 0},
    {  3,  2, {-1, -1}, "b",     2.0, 1},
    {  4,  2, {-1, -1}, "c",     3.0, 1},
    { -1, -1, {-1, -1}, NULL,   -1.0, 0},
};
/************************************/

/************************************/
static s_init_tree internaltree1_rogozin[] = {
    {0,  -1, { 1,  2}, "Root",          0.0,      0},
    {1,   0, {-1, -1}, "Plasmodium",    61.0, 1},
    {2,   0, { 3,  4}, "NonPlasmodium",  1.2,  0},
    {3,   2, {-1, -1}, "Arabidopsis",   29.3, 1},
    {4,   2, { 5,  6}, "FungiPlusMetazoa", 4.8, 0},
    {5,   4, { 7,  8}, "Fungi",         10.0, 0},
    {6,   4, { 9, 10}, "Metazoa",        4.2, 0},
    {7,   5, {-1, -1}, "Saccharomyces", 23.0, 1},
    {8,   5, {-1, -1}, "Schizosaccharomyces", 10.8, 1},
    {9,   6, {-1, -1}, "Celegans",      15.8, 1},
    {10,  6, {11, 12}, "NonCelegansMetazoans", 7.5, 0},
    {11, 10, {-1, -1}, "Human",         18.8, 1},
    {12, 10, {13, 14}, "Arthropods",     6.0, 0},
    {13, 12, {-1, -1}, "Anopheles",      9.6, 1},
    {14, 12, {-1, -1}, "Drosophila",     8.9, 1},
    {-1, -1, {-1, -1}, NULL,           -1.0,      0},
};
/************************************/

/************************************/
static s_init_tree internaltree2_threetaxa[] = {
    {  0, -1, { 1,  2}, "n0",  0.0, 0},
    {  1,  0, {-1, -1}, "v1", 30.0, 1},
    {  2,  0, { 3,  4}, "n4", 15.0, 0},
    {  3,  2, {-1, -1}, "v2", 15.0, 1},
    {  4,  2, {-1, -1}, "v3", 15.0, 1},
    { -1, -1, {-1, -1}, NULL, -1.0, 0},
};
/************************************/

/************************************/
static s_init_tree internaltree3_threetaxa[] = {
    {  0, -1, { 1,  2}, "n0",  0.0, 0},
    {  1,  0, {-1, -1}, "v1", 30.0, 1},
    {  2,  0, { 3,  4}, "n4", 30.0, 0},
    {  3,  2, {-1, -1}, "v2", 30.0, 1},
    {  4,  2, {-1, -1}, "v3", 15.0, 1},
    { -1, -1, {-1, -1}, NULL, -1.0, 0},
};
/************************************/

/************************************/
static s_init_tree internaltree4_threetaxa[] = {
    {  0, -1, { 1,  2}, "n0",   0.0, 0},
    {  1,  0, {-1, -1}, "v1", 100.0, 1},
    {  2,  0, { 3,  4}, "n4",  50.0, 0},
    {  3,  2, {-1, -1}, "v2",  50.0, 1},
    {  4,  2, {-1, -1}, "v3",  50.0, 1},
    { -1, -1, {-1, -1}, NULL,  -1.0, 0},
};
/************************************/

/************************************/
static s_init_tree internaltree5_sixtaxa[] = {
    {  0, -1, { 5, 10}, "n0",    0.0, 0},
    {  1,  3, {-1, -1}, "v1",   50.0, 1},
    {  2,  3, {-1, -1}, "v2",   50.0, 1},
    {  3,  5, { 1,  2}, "n7",   50.0, 0},
    {  4,  5, {-1, -1}, "v3",  100.0, 1},
    {  5,  0, { 3,  4}, "n8",   50.0, 0},
    {  6,  8, {-1, -1}, "v4",   50.0, 1},
    {  7,  8, {-1, -1}, "v5",   50.0, 1},
    {  8, 10, { 6,  7}, "n9",   50.0, 0},
    {  9, 10, {-1, -1}, "v6",  100.0, 1},
    { 10,  0, { 8,  9}, "n10",  50.0, 0},
    { -1, -1, {-1, -1}, NULL,   -1.0, 0},
};
/************************************/

//////////////////////
// Above these comments are data structures which define an internal tree
// that is used to establish a tree for testing purposes
//////////////////////

class PhyloTreeNode {
    private:  // PhyloTree does not need access here, PhyloTreeNode's methods
              // will handle all the tree-related operations
        std::string          name;
        std::string          description;
        std::string          comment;
        TaxonMatrix::taxon_id_type
                             taxon_id;
        int                  node_id;
        int                  leaf_id;
        bool                 is_leaf;
        bool                 is_root;
        double               branchlength;
        static const int     notset = -1;

        PhyloTreeNode*       ancestor;

        BranchRateManager*   branchratemanager;
        int                  branchrate_id;
        const BranchRateManager::p_trans_type*
                             ptr_p_trans;

        void                 set_branchrate_id(int brid);
        
    protected:
        void                 read_internal_tree(const s_init_tree* table, int n);
        int                  count_nodes() const;
        int                  count_leaves() const;
        int                  count_internal_nodes() const;
        int                  count_branches() const;


        void                 set_branchlength(double arg);
        void                 set_name(const std::string& arg);
        void                 set_description(const std::string& arg);
        void                 append_comment(const std::string& arg);
        void                 set_taxon_id(TaxonMatrix::taxon_id_type tid);
        void                 set_is_root(bool arg);
        void                 set_is_leaf(bool arg);
        void                 set_ancestor(PhyloTreeNode* arg);

    public:
        PhyloTreeNode();
        ~PhyloTreeNode();

        typedef std::vector<PhyloTreeNode*>
                             descendants_type;
        descendants_type     descendants;

        const std::string&   get_name() const;
        const std::string&   get_description() const;
        const std::string&   get_comment() const;
        int                  get_node_id() const;
        int                  get_leaf_id() const;
        TaxonMatrix::taxon_id_type
                             get_taxon_id() const;
        bool                 get_is_leaf() const;
        bool                 get_is_root() const;
        double               get_branchlength() const;

        PhyloTreeNode*       get_ancestor() const;
        PhyloTreeNode*       find_root();
        void                 add_descendant(PhyloTreeNode* node);
        int                  num_descendants() const;
        bool                 has_descendants() const;

        void                 find_taxa_for_leaves(TaxonMatrix& tm);

        BranchRateManager*   get_branchratemanager() const;
        void                 set_branchratemanager(BranchRateManager& b);
        int                  get_branchrate_id() const;
        void                 allocate_branchrate_ids();
        
        double               get_p_transition(int from, int to) const;
        double               get_rate_forward() const;
        double               get_rate_backward() const;

        static const std::string&
                             form_newick_label(const std::string& name);
        void                 basic_header_print(std::ostream& os = std::cout) const;
        void                 print(std::ostream& os = std::cout) const;
        void                 full_print(std::ostream& os = std::cout,
                                        const int indent = 0) const;
        void                 brief_print(std::ostream& os = std::cout,
                                         const int indent = 0) const;
        void                 newick_print(std::ostream& os = std::cout) const;

        friend std::ostream& operator<<(std::ostream& os,
                                        const PhyloTreeNode& node);
};

///////////////////////////
/////////////////////////// PhyloTreeNode methods
///////////////////////////

inline PhyloTreeNode::PhyloTreeNode() : name(""), description(""), comment(""),
    taxon_id(notset), node_id(notset), leaf_id(notset), is_leaf(false),
    is_root(false), branchlength(0.0), ancestor(NULL), branchratemanager(NULL),
    branchrate_id(notset) {
    // initialization of internal_trees[] vector
    internal_trees.push_back(internaltree0_small);
    internal_trees.push_back(internaltree1_rogozin);
    internal_trees.push_back(internaltree2_threetaxa);
    internal_trees.push_back(internaltree3_threetaxa);
    internal_trees.push_back(internaltree4_threetaxa);
    internal_trees.push_back(internaltree5_sixtaxa);
    //internal_trees.push_back(sixtaxa_tree1);
    //internal_trees.push_back(sixtaxa_tree2);
}

inline PhyloTreeNode::~PhyloTreeNode() {
    if (get_is_root()) {
        // delete all my descendants
    }
}

inline int PhyloTreeNode::get_node_id() const {
    return(node_id);
}

inline int PhyloTreeNode::get_leaf_id() const {
    return(leaf_id);
}

inline double PhyloTreeNode::get_p_transition(int from, int to) const {
    //return(ptr_p_trans[from][to]);
    //return(((*ptr_p_trans)[from])[to]);
    //return(ptr_p_trans->at(from).at(to));
    return(branchratemanager->get_p_transition(branchrate_id, from, to));
}

inline double PhyloTreeNode::get_rate_forward() const {
    return(branchratemanager->get_rate_forward(branchrate_id));
}

inline double PhyloTreeNode::get_rate_backward() const {
    return(branchratemanager->get_rate_backward(branchrate_id));
}

inline void PhyloTreeNode::set_branchratemanager(BranchRateManager& b) {
    branchratemanager = &b;
}

inline BranchRateManager* PhyloTreeNode::get_branchratemanager() const {
    return(branchratemanager);
}

inline void PhyloTreeNode::set_branchrate_id(int brid) {
    branchrate_id = brid;
}

inline int PhyloTreeNode::get_branchrate_id() const {
    return(branchrate_id);
}

inline PhyloTreeNode* PhyloTreeNode::find_root() {
    PhyloTreeNode* a = this;
    while (a != NULL && !a->get_is_root()) a = a->get_ancestor();
    return(a);
}

inline double PhyloTreeNode::get_branchlength() const {
    assert(branchratemanager == NULL ||
           branchrate_id < 0 ||
           (branchratemanager->get_branchlength(branchrate_id) == branchlength));
    return(branchlength);
}

inline void   PhyloTreeNode::set_branchlength(double arg) {
    branchlength = arg;
}

inline const std::string& PhyloTreeNode::get_name() const {
    return(name);
}

inline void   PhyloTreeNode::set_name(const std::string& arg) {
    name = arg;
}

inline const std::string& PhyloTreeNode::get_description() const {
    return(description);
}

inline void   PhyloTreeNode::set_description(const std::string& arg) {
    description = arg;
}

inline const std::string& PhyloTreeNode::get_comment() const {
    return(comment);
}

inline void PhyloTreeNode::append_comment(const std::string& arg) {
    comment += arg;
}

inline TaxonMatrix::taxon_id_type PhyloTreeNode::get_taxon_id() const
{
    return(taxon_id);
}

inline void PhyloTreeNode::set_taxon_id(TaxonMatrix::taxon_id_type tid)
{
    taxon_id = tid;
}

inline bool   PhyloTreeNode::get_is_leaf() const {
    return(is_leaf);
}

inline void   PhyloTreeNode::set_is_leaf(bool arg) {
    is_leaf = arg;
}

inline bool   PhyloTreeNode::get_is_root() const {
    return(is_root);
}

inline void   PhyloTreeNode::set_is_root(bool arg) {
    is_root = arg;
}

inline PhyloTreeNode* PhyloTreeNode::get_ancestor() const {
    return(ancestor);
}

inline void PhyloTreeNode::set_ancestor(PhyloTreeNode* arg) {
    ancestor = arg;
}

inline void PhyloTreeNode::add_descendant(PhyloTreeNode* node) {
    node->set_ancestor(this);
    descendants.push_back(node);
}

inline int PhyloTreeNode::num_descendants() const {
    return((int)descendants.size());
}

inline bool PhyloTreeNode::has_descendants() const {
    return(!descendants.empty());
}

inline void PhyloTreeNode::print(std::ostream& os) const {
    full_print(os);
    // newick_print(os);
}

#endif // PHYLOTREENODE_H


