#ifndef PHYLOTREE_H
#define PHYLOTREE_H

#include "PhyloTreeNode.h"
#include "TaxonMatrix.h"

#include <iostream>
#include <cassert>


class PhyloTree : public PhyloTreeNode {
    private:
        TaxonMatrix*  taxon_matrix;
        int           num_nodes;
        int           num_internal_nodes;
        int           num_leaves;
        int           num_branches;
        
        // This is the root of the tree, so use all the members of
        // PhyloTreeNode to manage as if it were root, and as a tree node
        // in general.

        // information on prior probabilities of trait states in the root
        // extensions of the algorithm could include this in each node, but
        // we're not concerned with that now

    public:
        PhyloTree();
        PhyloTree(const int treenum);
        ~PhyloTree();

        void                 count();
        int                  get_num_nodes() const;
        int                  get_num_internal_nodes() const;
        int                  get_num_leaves() const;
        int                  get_num_branches() const;

        void                 read_internal_tree(const int treenum = 0);
        inline TaxonMatrix*  get_taxon_matrix();
        inline void          set_taxon_matrix(TaxonMatrix& tm);
        void                 find_taxa_for_leaves();

        inline void          print(std::ostream& os = std::cout) const;

        friend std::ostream& operator<<(std::ostream& os, const PhyloTree& node);
};

///////////////////////////////
/////////////////////////////// PhyloTree methods
///////////////////////////////

inline PhyloTree::PhyloTree() :
        PhyloTreeNode(), taxon_matrix(NULL) {
    /* empty */
}

inline PhyloTree::PhyloTree(const int treenum) :
        PhyloTreeNode(), taxon_matrix(NULL) {
    read_internal_tree(treenum);
}

inline PhyloTree::~PhyloTree() {
}

inline int PhyloTree::get_num_nodes() const {
    return(num_nodes);
}

inline int PhyloTree::get_num_internal_nodes() const {
    return(num_internal_nodes);
}

inline int PhyloTree::get_num_leaves() const {
    return(num_leaves);
}

inline int PhyloTree::get_num_branches() const {
    return(num_branches);
}

inline void PhyloTree::find_taxa_for_leaves() {
    TaxonMatrix* ptm = get_taxon_matrix();
    assert(NULL != ptm);
    PhyloTreeNode::find_taxa_for_leaves(*ptm);
}

inline void PhyloTree::read_internal_tree(const int treenum) {
    // read one of the internal_trees[]
    PhyloTreeNode::read_internal_tree(internal_trees[treenum], 0);
    count();
}

inline void PhyloTree::count() {
    num_nodes = count_nodes();
    num_leaves = count_leaves();
    num_internal_nodes = count_internal_nodes();
    num_branches = count_branches();
}

inline TaxonMatrix* PhyloTree::get_taxon_matrix() {
    return(taxon_matrix);
}

inline void PhyloTree::set_taxon_matrix(TaxonMatrix& tm) {
    taxon_matrix = &tm;
    find_taxa_for_leaves();
}

inline void PhyloTree::print(std::ostream& os) const {
    os << "PhyloTree.taxon_matrix=" << ((void*)taxon_matrix);
    os << " num_nodes=" << num_nodes;
    os << " num_leaves=" << num_leaves;
    os << " num_internal_nodes=" << num_internal_nodes;
    os << " num_branches=" << num_branches;
    os << std::endl;
    PhyloTreeNode::print(os);
}

#endif // PHYLOTREE_H


