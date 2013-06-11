#ifndef BRANCHRATEMANAGER_H
#define BRANCHRATEMANAGER_H

#include "RatePVector.h"
#include "TraitMatrix.h"
#include "RatePMap.h"

#include <string>
#include <cassert>
#include <sstream>
#include <valarray>

#define   indent_block     "|   "

#define   DEBUG_READ_RATEP_INIT_FILE  1

class PhyloTreeNode;

class BranchRateManager {
    public:
        typedef std::vector< std::vector<double> > p_trans_type;

    private:
        const std::string            name_prefix;
        const std::string            name_forward;
        const std::string            name_backward;
        static const int             notset = -1;

        class BranchRate {
            private:
                int                  this_branchrate_id;
                std::string          name;

                PhyloTreeNode*       tip_node;
                double               branchlength;

                int                  forward_id;
                int                  backward_id;
                p_trans_type         p_trans;
                bool                 ratep_allocated;

            public:
                BranchRate();
                ~BranchRate();

                void                 set_this_branchrate_id(int brid);
                int                  get_this_branchrate_id() const;
                void                 set_name(const std::string& n);
                const std::string&   get_name() const;
                bool                 get_ratep_allocated() const;
                void                 set_ratep_allocated(bool arg);
                void                 set_forward_id(int fid);
                int                  get_forward_id() const;
                void                 set_backward_id(int bid);
                int                  get_backward_id() const;
                void                 set_branchlength(double arg);
                double               get_branchlength() const;
                double               get_p_trans(int from, int to) const;
                void                 set_p_trans(int from, int to, double p);
                const p_trans_type&  get_ref_p_trans() const;
                void                 set_tip_node(PhyloTreeNode* node);
                PhyloTreeNode*       get_tip_node() const;
                void                 print(std::ostream& os = std::cout) const;

        };

        std::vector<BranchRate>      branchrates;  // manage entries in a RatePVector

        RatePVector                  ratepv;  // actual rate parameters in a vector

        // private copy-constructor, to help detect pass-by-value errors
        BranchRateManager(const BranchRateManager& brm);

        void                         set_forward_id(int brid, int fid);
        void                         set_backward_id(int brid, int bid);

    public:
        BranchRateManager();
        ~BranchRateManager();

        RatePVector&         get_ratepvector();
        void                 allocate_ratep_from_map(RatePMap& map);
        void                 read_ratep_init_file(const std::string& filename);

        int                  allocate_new_branchrate_id(const std::string& n,
                                                        double bl);
        bool                 get_ratep_allocated(int brid) const;
        int                  get_forward_id(int brid) const;
        int                  get_backward_id(int brid) const;
        double               get_rate_forward(int brid) const;
        double               get_rate_backward(int brid) const;
        double               get_branchlength(int brid) const;
        void                 calculate_p_transition();
        void                 set_unidir_rate(double r);
        void                 set_unidir_rate(int brid, double r);
        void                 set_bidir_rates(double r_forw, double r_back);
        void                 set_bidir_rates(int brid, double r_forw, double r_back);
        double               get_p_transition(int brid, int from, int to) const;
        void                 set_p_transition(int brid, int from, int to, double p);
        const p_trans_type*  get_addr_p_trans(int brid) const;
        PhyloTreeNode*       get_tip_node(int brid) const;
        void                 print(int brid, std::ostream& os = std::cout) const;
        void                 print(std::ostream& os = std::cout) const;
        friend std::ostream& operator<<(std::ostream& os,
                                        const BranchRateManager& brm);
};

//////////////////////////////
////////////////////////////// BranchRateManager::BranchRate methods
//////////////////////////////

inline const BranchRateManager::p_trans_type& 
    BranchRateManager::BranchRate::get_ref_p_trans() const {

    //const BranchRateManager::p_trans_type* p = &p_trans;
    //std::cout << "get_ref_p_trans(): p_trans_type* p=" << ((void*)p) << std::endl;
    
    return(p_trans);
}



inline void BranchRateManager::BranchRate::print(std::ostream& os) const {
    os << "[" << this_branchrate_id << "]";
    os << " name=:" << get_name() << ":";
    os << " ratep_allocated=" << (get_ratep_allocated() == true ? "true" : "false");
    os << " forward_id=" << get_forward_id();
    os << " backward_id=" << get_backward_id();
    os << " branchlength=" << get_branchlength();
    os << " &p_trans=" << ((void*)&p_trans);
    os << " p_trans.size()=" << p_trans.size();
}

inline double BranchRateManager::BranchRate::get_p_trans(int from, int to) const {
    return(p_trans[from][to]);
}

inline void BranchRateManager::BranchRate::set_p_trans(int from, int to,
                                                       double p) {
    p_trans[from][to] = p;
}

inline bool BranchRateManager::BranchRate::get_ratep_allocated() const {
    return(ratep_allocated);
}

inline void BranchRateManager::BranchRate::set_ratep_allocated(bool arg) {
    ratep_allocated = arg;
}

inline void BranchRateManager::BranchRate::set_tip_node(PhyloTreeNode* node) {
    tip_node = node;
}

inline PhyloTreeNode* BranchRateManager::BranchRate::get_tip_node() const
{
    return(tip_node);
}

inline void BranchRateManager::BranchRate::set_name(const std::string& n)
{
    name = n;
}

inline const std::string& BranchRateManager::BranchRate::get_name() const
{
    return(name);
}

inline void BranchRateManager::BranchRate::set_branchlength(double bl)
{
    branchlength = bl;
}

inline double BranchRateManager::BranchRate::get_branchlength() const
{
    return(branchlength);
}

inline void BranchRateManager::BranchRate::set_this_branchrate_id(int brid)
{
    this_branchrate_id = brid;
}

inline int BranchRateManager::BranchRate::get_this_branchrate_id() const
{
    return(this_branchrate_id);
}

inline int BranchRateManager::BranchRate::get_forward_id() const
{
    return(forward_id);
}

inline void BranchRateManager::BranchRate::set_forward_id(int id)
{
    forward_id = id;
}

inline int BranchRateManager::BranchRate::get_backward_id() const
{
    return(backward_id);
}

inline void BranchRateManager::BranchRate::set_backward_id(int id)
{
    backward_id = id;
}

inline BranchRateManager::BranchRate::BranchRate() : this_branchrate_id(notset),
    name(""), forward_id(notset), backward_id(notset), ratep_allocated(false) {
    p_trans.resize(TraitMatrix::max_num_states);
    for (size_t i = 0; i < p_trans.size(); ++i) {
        p_trans[i].resize(TraitMatrix::max_num_states);
        for (size_t j = 0; j < p_trans[i].size(); ++j) {
            p_trans[i][j] = 100.0 + double(i)*10.0 + double(j);
        }
    }
}

inline BranchRateManager::BranchRate::~BranchRate() {
    /* EMPTY */
}

////////////////////////////////////////
//////////////////////////////////////// BranchRateManager methods
////////////////////////////////////////

inline const BranchRateManager::p_trans_type*
    BranchRateManager::get_addr_p_trans(int brid) const {
    const BranchRateManager::p_trans_type* p = &(branchrates[brid].get_ref_p_trans());
    std::cout << std::endl << "get_addr_p_trans(): p_trans_type* p=" << ((void*)p) << std::endl;
    return(&(branchrates[brid].get_ref_p_trans()));
}

inline RatePVector& BranchRateManager::get_ratepvector() {
    return(ratepv);
}

inline BranchRateManager::BranchRateManager() : name_prefix("branch_"),
    name_forward("forward_"), name_backward("backward_")
{
    /* EMPTY */
}

inline BranchRateManager::~BranchRateManager() {
    /* EMPTY */
}

inline int BranchRateManager::allocate_new_branchrate_id(const std::string& n,
                                                         double bl) {
    int end = branchrates.size();
    branchrates.resize(end + 1);
    if (n == "") {
        std::stringstream ss;
        ss << name_prefix << end;
        branchrates[end].set_name(ss.str());
    } else {
        branchrates[end].set_name(n);
    }
    branchrates[end].set_branchlength(bl);
    branchrates[end].set_this_branchrate_id(end);
    return(end);
}

inline void BranchRateManager::set_forward_id(int brid, int fid) {
    branchrates[brid].set_forward_id(fid);
}

inline int BranchRateManager::get_forward_id(int brid) const {
    return(branchrates[brid].get_forward_id());
}

inline void BranchRateManager::set_backward_id(int brid, int bid) {
    branchrates[brid].set_backward_id(bid);
}

inline int BranchRateManager::get_backward_id(int brid) const {
    return(branchrates[brid].get_backward_id());
}

inline bool BranchRateManager::get_ratep_allocated(int brid) const
{
    return(branchrates[brid].get_ratep_allocated());
}

inline double BranchRateManager::get_branchlength(int brid) const {
    return(branchrates[brid].get_branchlength());
}

inline double BranchRateManager::get_rate_forward(int brid) const {
    return( ((BranchRateManager*)this)->ratepv[branchrates[brid].get_forward_id()].get_ratep() );
//    return(    ratepv[branchrates[brid].get_forward_id()].get_ratep()    );
}

inline double BranchRateManager::get_rate_backward(int brid) const {
    return( ((BranchRateManager*)this)->ratepv[branchrates[brid].get_backward_id()].get_ratep() );
//    return(    ratepv[branchrates[brid].get_backward_id()].get_ratep()    );
}

inline double BranchRateManager::get_p_transition(int brid, int from,
                                                        int to) const {
    return(branchrates[brid].get_p_trans(from, to));
}

inline void BranchRateManager::set_p_transition(int brid, int from, int to,
                                                double p) {
    branchrates[brid].set_p_trans(from, to, p);
}

inline PhyloTreeNode* BranchRateManager::get_tip_node(int brid) const {
    return(branchrates[brid].get_tip_node());
}

#endif  // BRANCHRATEMANAGER_H


