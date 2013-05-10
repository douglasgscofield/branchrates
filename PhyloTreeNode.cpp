#include "PhyloTreeNode.h"

void PhyloTreeNode::find_taxa_for_leaves(TaxonMatrix& tm) {
    assert(tm.taxa.empty() == false);
    if (get_is_leaf() == true) {
        TaxonMatrix::taxon_id_type tid = tm.get_taxon_id(name);
        if (tid < 0) {
            std::cerr << "find_taxa_for_leaves: no taxon_id for taxon "
                      << name << std::endl;
            std::cerr << tm;
        }
        assert(tid >= 0);
        assert(get_leaf_id() >= 0);
        tm.taxa[tid].set_leaf_id(get_leaf_id());
        set_taxon_id(tid);
    } else {
        for (int i = 0; i < descendants.size(); ++i) {
            descendants[i]->find_taxa_for_leaves(tm);
        }
    }
}

void PhyloTreeNode::read_internal_tree(const s_init_tree* table, int n) {
    static int leafnum = -1;
    set_name(table[n].name);
    set_branchlength(table[n].branchlength);
    node_id = table[n].nodenum;
    if (0 == n)
        set_is_root(true);
    set_is_leaf(table[n].is_leaf);
    if (get_is_leaf())
        leaf_id = ++leafnum;
    set_ancestor(NULL);
    for (int i = 0; i < MAX_init_tree_descendants; ++i) {
        int d;
        if ((d = table[n].descendant[i]) > 0) {
            assert(! get_is_leaf());
            PhyloTreeNode* pdesc = new PhyloTreeNode;
            pdesc->read_internal_tree(table, d);
            add_descendant(pdesc);
        }
    }
}

int PhyloTreeNode::count_nodes() const {
    int n = 1;  // one for me
    for (int i = 0; i < descendants.size(); ++i) {
        // and ones for my homies
        n += descendants[i]->count_nodes();
    }
    return(n);
}

int PhyloTreeNode::count_leaves() const {
    int n = (get_is_leaf() == true) ? 1 : 0;
    for (int i = 0; i < descendants.size(); ++i) {
        n += descendants[i]->count_leaves();
    }
    return(n);
}

int PhyloTreeNode::count_internal_nodes() const {
    int n = (get_is_leaf() == false) ? 1 : 0;
    for (int i = 0; i < descendants.size(); ++i) {
        n += descendants[i]->count_internal_nodes();
    }
    return(n);
}

int PhyloTreeNode::count_branches() const {
    int n = 0;
    for (int i = 0; i < descendants.size(); ++i) {
        ++n;
        n += descendants[i]->count_branches();
    }
    return(n);
}

void PhyloTreeNode::allocate_branchrate_ids() {
    assert(branchratemanager != NULL);
    assert(get_branchrate_id() < 0);
    if (get_is_root() == false) {
        int brid = branchratemanager->allocate_new_branchrate_id(get_name(),
                                                          get_branchlength());
        set_branchrate_id(brid);
        
        //ptr_p_trans = branchratemanager->get_addr_p_trans(brid);
        ptr_p_trans = NULL;
    }
    for (int i = 0; i < descendants.size(); ++i) {
        descendants[i]->set_branchratemanager(*branchratemanager);
        descendants[i]->allocate_branchrate_ids();
    }
}


std::ostream& operator<<(std::ostream& os, const PhyloTreeNode& node) {
    os << "PhyloTreeNode: ";
    node.print(os);
    os << std::endl;
    return(os);
}

void PhyloTreeNode::basic_header_print(std::ostream& os) const {
    os << "<PhyloTreeNode";
    os << " name=:" << name << ":";
    os << " taxon_id=" << taxon_id;
    os << " nid=" << node_id;
    os << " lid=" << leaf_id;
    if (get_is_root())
        os << " ROOT";
    if (get_is_leaf())
        os << " LEAF";
    os << " branchrate_id=" << branchrate_id;
    os << " branchlength=" << branchlength;
    os << ">";
}

void PhyloTreeNode::brief_print(std::ostream& os, const int indent) const{
    for (int i = 0; i < indent; ++i) os << indent_block;
    basic_header_print(os);
    if (0 <= branchrate_id) {
        assert(get_is_root() == false);
        os << " ptr_p_trans=" << (void*)ptr_p_trans;
        //os << " ptr_p_trans->size()=" << ptr_p_trans->size();
        //os << " (*ptr_p_trans).size()=" << (*ptr_p_trans).size();
        //os << " (*ptr_p_trans)[1][1]=" << (*ptr_p_trans).at(1).at(1);
        //os << " ptr_p_trans->at(0).at(0)=" << ptr_p_trans->at(0).at(0);
        //os << " ptr_p_trans->at(0)[1]=" << ptr_p_trans->at(0)[1];
        //os << " ptr_p_trans->at(0)[2]=" << ptr_p_trans->at(0)[2];
        //os << " (*ptr_p_trans)[1][1]=" << (*ptr_p_trans)[1][1];
        //os << " ptr_p_trans->size()=" << ptr_p_trans->size();
        os << " rate_forward=" << get_rate_forward();
        os << " rate_backward=" << get_rate_backward();
        os << " p_trans[0][0,1]=[" << get_p_transition(0, 0) << ", "
           << get_p_transition(0, 1) << "]";
        os << " p_trans[1][0,1]=[" << get_p_transition(1, 0) << ", "
           << get_p_transition(1, 1) << "]";
    }
    os << std::endl;
    for (int i = 0; i < descendants.size(); ++i) {
        descendants[i]->brief_print(os, indent + 1);
    }
    os << std::endl;
}

void PhyloTreeNode::full_print(std::ostream& os, const int indent) const {
    static const int annotate_description = 1;
    static const int annotate_comment = 1;
    static const int annotate_branchrate_entries = 1;
    
    // first line: general node information
    for (int i = 0; i < indent; ++i) os << indent_block;
    basic_header_print(os);
    os << std::endl;
    
    // second line: potentially long strings
    if (annotate_description || annotate_comment) {
        for (int i = 0; i < indent; ++i)
            os << indent_block;
        os << "<---";
        if (annotate_description)
            os << "description=:" << description << ":";
        if (annotate_comment)
            os << " comment=:" << comment << ":";
        os << ">" << std::endl;
    }
    
    // third line: rate information
    for (int i = 0; i < indent; ++i)
        os << indent_block;
    os << "<---";
    os << "branchratemanager=" << ((void*)branchratemanager);
    os << " branchrate_id=" << branchrate_id << " ";
    if (annotate_branchrate_entries)
        branchratemanager->print(branchrate_id, os);
    os << ">" << std::endl;
    
    for (int i = 0; i < descendants.size(); ++i) {
        descendants[i]->full_print(os, indent + 1);
    }
}


void PhyloTreeNode::newick_print(std::ostream& os) const {
    static const int annotate_id = 1;
    static const int annotate_root = 1;
    static const int annotate_leaf = 1;
    static const int annotate_description = 1;
    static const int annotate_comment = 0;
    // need to check against Newick grammar
    if (!get_is_leaf()) os << "(";
    for (int i = 0; i < descendants.size(); ++i) {
        descendants[i]->newick_print(os);
    }
    if (!get_is_leaf()) {
        os << ")";
    }
    os << form_newick_label(name) << ":" << branchlength;
    if (annotate_root && get_is_root()) os << "[ROOT]";
    if (annotate_leaf && get_is_leaf()) os << "[LEAF]";
    if (annotate_id) {
        os << "[taxon_id=" << taxon_id << "]";
        os << "[node_id=" << node_id << "]";
        if (get_is_leaf()) os << "[leaf_id=" << leaf_id << "]";
    }
    if (annotate_description && !description.empty()) {
        os << "[DESCRIPTION=" << description << "]";
    }
    if (annotate_comment) {
        os << "[COMMENT BEGIN]" << comment << "[COMMENT END]";
    }
}

const std::string& PhyloTreeNode::form_newick_label(const std::string& name) {
    static std::string ans;
    ans = "";
    bool quote_needed = false, blank_seen = false;
    for (int i = 0; i < name.length(); ++i) {
        switch(name[i]) {
            case ' ':
                blank_seen = true;
                break;
            case '\'': case '(': case ')': case ';':
            case '[':  case ']': case ':': case ',':
                quote_needed = true;
                i = name.length();  // terminate loop
                break;
        }
    }
    if (!quote_needed) {  //
        if (!blank_seen)
            return(name);
        else {
            for (int i = 0; i < name.length(); ++i)
                ans += (name[i] == ' ') ? '_' : name[i];
            return(ans);
        }
    }
    ans = '\'';
    for (int i = 0; i < name.length(); ans += name[i++]) {
        if (name[i] == '\'')
            ans += '\'';
    }
    ans += '\'';
    return(ans);
}


