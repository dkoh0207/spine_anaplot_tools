/**
 * @file category_nue.h
 * @brief Header file for definitions of analysis variables which can be used
 * to categorize interactions in the nue analysis.
 */

#ifndef CATEGORY_NUE_H
#define CATEGORY_NUE_H

#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRInteractionDLP.h"
#include "sbnanaobj/StandardRecord/SRInteractionTruthDLP.h"

#include "include/utilities.h"

namespace category::nue
{
    double category(const caf::SRInteractionTruthDLPProxy & obj)
    {
        double cat(8);
        if(cuts::nue::signal_1eNp(obj)) cat = 0;
        else if(cuts::nue::is_nue_cc(obj)) cat = 1;
        else if(cuts::nue::is_nue_nc(obj)) cat = 2;
        else if(cuts::nue::is_numu_cc(obj)) cat = 3;
        else if(cuts::nue::is_numu_nc(obj)) cat = 4;
        else if(cuts::cosmic(obj)) cat = 5;
        return cat;
    }
}

#endif // CATEGORY_NUE_H