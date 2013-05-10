#include "ML_single_NewtonRaphson.h"

//////////////////////
// Maximize the likelihood when dealing with a single parameter
//////////////////////
void ML_single_NewtonRaphson::find_maximum() {

    history.clear();

    RatePVector& rpv = branch_rate_manager.get_ratepvector();
    assert(rpv.size() == 1);
    int thisparm = 0;
    Likelihood::RatePBounds& b = bounds[thisparm];
    
    b.delta = determine_machine_delta();
    // double b.delta = 0.0000001;
    double parameter = (b.lowerbound + 
                        b.upperbound) / 2.0;
    double logL_m_delta = compute_for_single_parameter(parameter - 
                                                       b.delta);
    // double logL_m_05delta = compute_for_single_parameter(parameter -
    //                                                 (0.5*b.delta));
    double logL = compute_for_single_parameter(parameter);
    // double logL_p_05delta = compute_for_single_parameter(parameter +
    //                                                 (0.5*b.delta));
    double logL_p_delta = compute_for_single_parameter(parameter + 
                                                  b.delta);

    history.add(parameter, b.accuracy, b.delta, 
                b.lowerbound, b.upperbound, 
                logL_m_delta,
                // logL_m_05delta,
                logL,
                // logL_p_05delta,
                logL_p_delta);

    do {
        parameter = adjust();
        if (parameter < b.lowerbound || parameter > b.upperbound) {
            std::cout << "warning: exceeded bound" << std::endl;
        }
        logL_m_delta = compute_for_single_parameter(parameter - b.delta);
        // logL_m_05delta = compute_for_single_parameter(parameter - (0.5*b.delta));
        logL = compute_for_single_parameter(parameter);
        // logL_p_05delta = compute_for_single_parameter(parameter + (0.5*b.delta));
        logL_p_delta = compute_for_single_parameter(parameter + b.delta);
        history.add(parameter, b.accuracy, b.delta, b.lowerbound, 
                    b.upperbound, logL_m_delta,
                    // logL_m_05delta,
                    logL,
                    // logL_p_05delta,
                    logL_p_delta);
        if (MAXDEBUG) {
            std::cout << "History complete:" << std::endl;
            // history.print();
            history.simple_print();
        }
    } while (! converged());
    
    if (MAXDEBUG_HISTORY_AT_CONVERGENCE) {
        history.simple_print();
    }

    int hsz = history.size();
    if (compute_for_single_parameter(b.lowerbound) > history[hsz - 1].logL) {
        std::cerr << "log-Likelihood(b.lowerbound) > log-Likelihood(converged)"
                  << std::endl;
    }
    if (compute_for_single_parameter(b.upperbound) > history[hsz - 1].logL) {
        std::cerr << "log-Likelihood(b.upperbound) > log-Likelihood(converged)"
                  << std::endl;
    }
    rpv[thisparm].set_ratep(history[hsz - 1].parameter);

    set_log_likelihood(history[hsz - 1].logL);
}


double ML_single_NewtonRaphson::adjust() {
    int hsz = history.size();
    assert(hsz > 0);
    HistorySingleParameter::HistoryEntry& hist = history[hsz - 1];

    double newparameter, dparameter, sign_dparameter;

    sign_dparameter = (hist.f_parameter < 0.0) ? -1.0 : +1.0;

    double dparameter_numerator = ((hist.logL_p_delta - hist.logL) * hist.delta);
    double dparameter_denominator = (hist.logL_p_delta - (2.0 * hist.logL) + 
                                     hist.logL_m_delta);

    if (dparameter_numerator == 0.0) {
        // The numeric derivative to the right of parameter is 0
        if (hist.logL - hist.logL_m_delta == 0.0) {
            // The numeric derivative to the left of parameter is 0, so we're in a
            // zero-slope neighborhood, and we don't know where to go.  Issue a
            // warning, return the original parameter, and let converged_single()
            // stop us.
            std::cerr << "adjust_single() in zero-slope neighborhood with parameter="
                      << hist.parameter << " and delta=" << hist.delta 
                      << ", are the bounds" << "[" << hist.lowerbound << ", "
                      << hist.upperbound << "] wrong?" << std::endl;
        }
        newparameter = hist.parameter;
    } else if (dparameter_denominator == 0.0) {
        // We won't converge correctly in this neighborhood following N-R, as
        // the 2nd derivative of the likelihood surface is 0; so, adjust parameter
        // by following the slope upwards.  We may find a suitable neighborhood
        // in the future.
        newparameter = (sign_dparameter > 0) ? hist.parameter_p_delta 
                                             : hist.parameter_m_delta;
    } else {
        dparameter = std::abs(dparameter_numerator / dparameter_denominator);
        newparameter = hist.parameter + (sign_dparameter * dparameter);
    }
    
    if (MAXDEBUG_ADJUSTSINGLE) {
        std::cout << "adjust():" << std::endl;
        std::cout << "parameter=" << hist.parameter;
        std::cout << " f(parameter)=" << hist.f_parameter;
        std::cout << " f'(parameter)=" << hist.fp_parameter;
        std::cout << " f(parameter)/f'(parameter)=" << 
                     hist.f_parameter/hist.fp_parameter;
        std::cout << " sign_dparameter=" << sign_dparameter;
        std::cout << " dparameter=" << dparameter;
        std::cout << " newparameter=" << newparameter;
        std::cout << std::endl;
    }
    if (newparameter < hist.lowerbound) {
        std::cerr << "adjust_rate() wants to return a parameter="
                  << newparameter << " less than lowerbound=" 
                  << hist.lowerbound << std::endl;
        newparameter = hist.lowerbound;
    } else if (newparameter > hist.upperbound) {
        std::cerr << "adjust_rate() wants to return a parameter="
                  << newparameter << " greater than upperbound=" 
                  << hist.upperbound << std::endl;
        newparameter = hist.upperbound;
    }
    return (newparameter);
}


