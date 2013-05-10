#include "ML_single.h"


bool ML_single::converged() {
    int hsz = history.size();
    if (hsz <= 2) {
        return (false);
    }

    HistorySingleParameter::HistoryEntry& hist_2prev = history[hsz - 3];
    HistorySingleParameter::HistoryEntry& hist_prev = history[hsz - 2];
    HistorySingleParameter::HistoryEntry& hist = history[hsz - 1];
    // if we met the convergence criteria for the last two steps
    double diff_prev = (hist_2prev.parameter - hist_prev.parameter);
    double diff = (hist_prev.parameter - hist.parameter);
    if (std::abs(diff_prev) < hist_prev.epsilon &&
        std::abs(diff) < hist.epsilon) {
        return (true);
    } else {
        return (false);
    }
}

