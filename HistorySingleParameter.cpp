#include "HistorySingleParameter.h"

std::ostream& operator<<(std::ostream& os, const HistorySingleParameter& h) {
    os << "HistorySingleParameter[" << (h.size()-1) << "..0]" << std::endl;
    h.print(os);
    return(os);
}

/////////////////////////
///////////////////////// HistorySingleParameter methods
/////////////////////////

void HistorySingleParameter::print(std::ostream& os) const {
    for (int i = (size()-1); i >= 0; --i) {
        os << "   [" << i << "] ";
        (((HistorySingleParameter*)this)->history)[i].print(os);
    }
}

void HistorySingleParameter::simple_print(std::ostream& os) const {
    os << "parameter" << "\t" << "logL" << std::endl;
    for (int i = (size()-1); i >= 0; --i) {
        os << i << "\t";
        (((HistorySingleParameter*)this)->history)[i].simple_print(os);
        os << std::endl;
    }
}

/////////////////////////
///////////////////////// HistorySingleParameter::HistoryEntry methods
/////////////////////////

inline void HistorySingleParameter::HistoryEntry::print(std::ostream& os) const {
    os << " delta=" << delta;
    os << " parameter[-delta,-0.5delta,0,+0.5delta,+delta]=["
       << parameter_m_delta << ", "
       // << parameter_m_05delta << ", "
       << parameter << ", "
       // << parameter_p_05delta << ", "
       << parameter_p_delta << "]";
    os << " epsilon=" << epsilon;
    os << " logL[-delta,-0.5delta,0,+0.5delta,+delta]=["
       << logL_m_delta << ", "
       // << logL_m_05delta << ", "
       << logL << ", "
       // << logL_p_05delta << ", "
       << logL_p_delta << "]";
    os << " f_lo=" << f_lo;
    os << " f_hi=" << f_hi;
    os << " f_parameter=" << f_parameter;
    os << " fp_parameter=" << fp_parameter;
    os << std::endl;
}

inline void HistorySingleParameter::HistoryEntry::simple_print(
                                                   std::ostream& os) const {
    os << parameter << "\t" <<  logL;
}


