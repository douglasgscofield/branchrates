#include "RatePMap.h"

RatePMap::~RatePMap()
{
    for (size_t i = 0; i < map.size(); ++i) {
        delete map[i];
    }
}


void RatePMap::read_file(const std::string& filename)
{
    std::ifstream file(filename.c_str());
    assert(file);
    const int bufsz = 4096;
    char buf[bufsz+1];
    int numline = 0;
    const int maxfields = 3;
    while (file.eof() == false) {
        ++numline;
        file.getline(buf, bufsz, '\n');
        std::istringstream inp(buf);
        std::vector<std::string> field(maxfields);
        int nf = 0;
        std::string w;
        while (inp >> w) {
            field[nf++] = w;
            if (nf >= maxfields)
                break;
        }
        if (nf == 0)
            continue;
        if (nf != maxfields) {
            std::cerr << filename.c_str() << ":" << numline 
                        << " need three fields" << std::endl;
            assert(false);
        }
        if (field[0] == "forward") {
            map_forward_rate(field[1], field[2]);
        } else if (field[0] == "backward") {
            map_backward_rate(field[1], field[2]);
        } else {
            std::cerr << filename.c_str() << ":" << numline 
                        << " " << field[0] << " is unrecognized first field"
                        << std::endl;
            assert(false);
        }
    }
}


std::ostream& operator<<(std::ostream& os, const RatePMap& rpm)
{
    os << "RatePMap: " << std::endl;
    rpm.print(os);
    os << std::endl;
    return(os);
}


void RatePMap::print(std::ostream& os) const
{
    //os << "  ratep_vector=" << ratep_vector;
    //os << "  tree_node_decorated=" << (tree_node_decorated ? "true" : "false");
    //os << "  ratep_allocated=" << (ratep_allocated ? "true" : "false");
    //os << std::endl;
    for (size_t i = 0; i < map.size(); ++i) {
        os << "  branch=" << map[i]->branchname;
        if (map[i]->is_forward)
            os << " forward";
        if (map[i]->is_backward)
            os << " backward";
        // os << "  is_forward=" << (map[i]->is_forward ? "true" : "false");
        // os << "  is_backward=" << (map[i]->is_backward ? "true" : "false");
        os << "  ratename=" << map[i]->ratename;
        os << "  ratep_id=" << map[i]->ratep_id;
        os << std::endl;
    }
}


void RatePMap::is_ok() const
{
    // check for duplicate entries for a branch rate
    for (size_t i = 0; i < map.size()-1; ++i) {
        for (size_t j = i+1; j < map.size(); ++j) {
            if (map[i]->branchname == map[j]->branchname) {
                std::cerr << "duplicate RatePMap entries for " << map[i]->branchname
                     << std::endl;
            }
        }
    }
}

//void RatePMap::decorate_with_tree_nodes(PhyloTreeNode* node) {
//    assert(node != NULL);
//    int i;
//    for (i = 0; i < map.size(); ++i) {
//        if (node->get_name() == map[i]->branchname) {
//            map[i]->node_for_branch = node;
//            break;
//        }
//    }
//    if (i == map.size()) {
//        // we didn't find the name of the node
//        std::cerr << "RatePMap::decorate: couldn't find a name for " << node
//            << std::endl;
//    }
//    for (int i = 0; i < node->descendants.size(); ++i) {
//        RatePMap::decorate_with_tree_nodes(node->descendants[i]);
//    }
//    tree_node_decorated = true;
//}

