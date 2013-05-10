#ifndef RATEPMAP_H
#define RATEPMAP_H

#include "RatePVector.h"
//#include "PhyloTreeNode.h"

#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cassert>

class RatePMap {
    private:
        class MapEntry {
            public:
                std::string             branchname;
                //PhyloTreeNode*          node_for_branch;
                bool                    is_forward;
                bool                    is_backward;
                std::string             ratename;
                int                     ratep_id;
                MapEntry(const std::string& b, const std::string& r);
        };
        std::vector<MapEntry*> map;
        //RatePVector*           ratep_vector;
        //bool                   tree_node_decorated;
        //bool                   ratep_allocated;
        
    public:
        RatePMap();
        RatePMap(const std::string& filename);
        ~RatePMap();

        void                   read_file(const std::string& filename);

        int                    size() const;
        const std::string&     get_branchname(int i) const;
        bool                   get_is_forward(int i) const;
        bool                   get_is_backward(int i) const;
        const std::string&     get_ratename(int i) const;
        int                    get_ratep_id(int i) const;
        void                   set_ratep_id(int i, int rid);

        void                   map_forward_rate(const std::string& b,
                                                const std::string& r);
        void                   map_backward_rate(const std::string& b,
                                                 const std::string& r);
        //void                   decorate_with_tree_nodes(PhyloTreeNode* node);
        //void                   allocate_ratep();
        void                   is_ok() const;

        void                   print(std::ostream& os = std::cout) const;
        friend std::ostream&   operator<<(std::ostream& os, const RatePMap& rpm);
};

///////////////////////////////
///////////////////////////////  RatePMap methods
///////////////////////////////

inline int RatePMap::size() const {
    return(map.size());
}

inline RatePMap::RatePMap() {
    /* EMPTY */
}

inline RatePMap::RatePMap(const std::string& filename) {
    read_file(filename);
}

inline RatePMap::MapEntry::MapEntry(const std::string& b, const std::string& r)
    : branchname(b), ratename(r), is_forward(false), is_backward(false),
      ratep_id(-1)
{
    /* EMPTY */
}

inline void RatePMap::map_forward_rate(const std::string& b,
                                       const std::string& r)
{
    MapEntry* p = new MapEntry(b, r);
    p->is_forward = true;
    map.push_back(p);
}

inline void RatePMap::map_backward_rate(const std::string& b,
                                        const std::string& r)
{
    MapEntry* p = new MapEntry(b, r);
    p->is_backward = true;
    map.push_back(p);
}

inline const std::string& RatePMap::get_branchname(int i) const {
    return(map[i]->branchname);
}

inline bool RatePMap::get_is_forward(int i) const {
    return(map[i]->is_forward);
}

inline bool RatePMap::get_is_backward(int i) const {
    return(map[i]->is_backward);
}

inline const std::string& RatePMap::get_ratename(int i) const {
    return(map[i]->ratename);
}

inline int RatePMap::get_ratep_id(int i) const {
    return(map[i]->ratep_id);
}

inline void RatePMap::set_ratep_id(int i, int rid) {
    map[i]->ratep_id = rid;
}

#endif  // RATEPMAP_H

