#ifndef HISTORY_MULTI_PARAMETER_H
#define HISTORY_MULTI_PARAMETER_H

#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>

class HistoryMultiParameter {
    public:
        class HistoryEntry {
            public:
                double parameter;         double delta;
                double lowerbound;        double upperbound;
                double parameter_m_delta; // double parameter_m_05delta;
                double parameter_p_delta; // double parameter_p_05delta;
                double epsilon;       
                double logL_m_delta;      // double logL_m_05delta;
                double logL;
                double logL_p_delta;      // double logL_p_05delta;
                double f_lo;              double f_hi;
                double f_parameter;       double fp_parameter;

                HistoryEntry(double pa, double ep, double de, double lo, 
                             double up, double lL_m_de,
                             // double lL_m_05de,
                             double lL,
                             // double lL_p_05de,
                             double lL_p_de);

                void           print(std::ostream& os = std::cout) const;
                void           simple_print(std::ostream& os = std::cout) const;
        };
        
        std::vector<HistoryEntry>
                               history;

        void                   add(double pa, double ep, double de, double lo,
                                   double up, double lL_m_de,
                                   // double lL_m_05de,
                                   double lL,
                                   // double lL_p_05de,
                                   double lL_p_de);
        void                   clear();
        int                    size() const;
        void                   print(std::ostream& os = std::cout) const;
        void                   simple_print(std::ostream& os = std::cout) const;
        HistoryEntry&          operator[](int i);
        friend std::ostream&   operator<<(std::ostream& os,
                                          const HistoryMultiParameter& h);
};

/////////////////////////
///////////////////////// HistoryMultiParameter::HistoryEntry methods
/////////////////////////

inline HistoryMultiParameter::HistoryEntry::HistoryEntry(double pa, double ep,
        double de, double lo, double up, double lL_m_de,
        // double lL_m_05de,
        double lL,
        // double lL_p_05de,
        double lL_p_de) : 
        parameter(pa), delta(de), lowerbound(lo), upperbound(up), 
        parameter_m_delta(parameter-delta), 
        // parameter_m_05delta(parameter-0.5*delta),
        parameter_p_delta(parameter+delta), 
        // parameter_p_05delta(parameter+0.5*delta),
        epsilon(ep), logL_m_delta(lL_m_de),
        // logL_m_05delta(lL_m_05de),
        logL(lL),
        // logL_p_05delta(lL_p_05de),
        logL_p_delta(lL_p_de), 
        f_lo((logL - logL_m_delta)/delta), f_hi((logL_p_delta - logL)/delta), 
        // f_parameter((logL_p_05delta - logL_m_05delta)/delta),
        f_parameter((logL_p_delta - logL_m_delta)/(2.0*delta)),
        fp_parameter((f_hi - f_lo)/delta)
{
    /* empty */
}

inline HistoryMultiParameter::HistoryEntry&
                       HistoryMultiParameter::operator[](int i) {
    return(history[i]);
}

/////////////////////////
///////////////////////// HistoryMultiParameter methods
/////////////////////////

inline void HistoryMultiParameter::add(double pa, double ep, double de,
                                        double lo, double up, double lL_m_de,
                                        // double lL_m_05de,
                                        double lL,
                                        // double lL_p_05de,
                                        double lL_p_de)
{
    HistoryEntry* he = new HistoryEntry(pa, ep, de, lo, up, lL_m_de,
                                        // lL_m_05de,
                                        lL,
                                        // lL_p_05de,
                                        lL_p_de);
    history.push_back(*he);
}

inline void HistoryMultiParameter::clear() {
    history.clear();
}

inline int HistoryMultiParameter::size() const {
    return(history.size());
}

#endif  // HISTORY_MULTI_PARAMETER_H


