#include "PhyloTree.h"

std::ostream& operator<< (std::ostream& os, const PhyloTree& node) {
    os << "PhyloTree: ";
    node.print(os);
    os << std::endl;
    return(os);
}

