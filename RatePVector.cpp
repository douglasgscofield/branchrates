#include "RatePVector.h"

void RatePVector::print(std::ostream& os) const {
    os << "name_prefix=:" << name_prefix << ":";
    os << " prev_in_use=" << (prev_in_use == true ? "true" : "false")
       << std::endl;
    for (int i = 0; i < ratepvec.size(); ++i) {
        os << "    [" << i << "] ";
        ratepvec[i].print(os);
        os << std::endl;
    }
}

std::ostream& operator<< (std::ostream& os, const RatePVector& rpv) {
    os << "RatePVector:" << std::endl;
    rpv.print(os);
    return(os);
}


