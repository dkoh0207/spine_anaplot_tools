/**
 * @file cuts_nue.h
 * @brief Header file for definitions of analysis cuts specific to the nue
 * analysis.
 * @details This file contains definitions of analysis cuts which can be used
 * to select interactions specific to the nue analysis. The cuts are
 * intended to be used in conjunction with the generic cuts defined in cuts.h.
 * @author mueller@fnal.gov
*/
#ifndef CUTS_NUE_H
#define CUTS_NUE_H

#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>

#include "include/utilities.h"
#include "include/bnb_nue/variables_nue.h"

/**
 * @namespace cuts::nue
 * @brief Namespace for organizing cuts specific to the nue analysis.
 * @details This namespace is intended to be used for organizing cuts which act
 * on interactions specific to the nue analysis. Each cut is implemented as
 * a function which takes an interaction object as an argument and returns a
 * boolean. The function should be templated on the type of interaction object if
 * the cut is intended to be used on both true and reconstructed interactions.
 * @note The namespace is intended to be used in conjunction with the cuts
 * namespace, which is used for organizing generic cuts which act on interactions.
 */
namespace cuts::nue
{
    /**
     * @brief Apply a 1eNp topological (final state) cut.
     * @details The interaction must have a topology matching 1muNp as defined by
     * the conditions in the @ref utilities::count_primaries() function.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has a 1eNp topology.
     * @note This cut is intended to be used for the nue analysis.
     */
    template<class T>
        bool topological_1eNp_cut(const T & obj)
        {
            std::vector<uint32_t> c(utilities::count_primaries_nue(obj));
            return c[0] == 0 && c[1] == 1 && c[2] == 0 && c[3] == 0 && c[4] >= 1;
        }

    /**
     * @brief Apply a fiducial volume, containment, flash time (BNB), and 1eNp
     * topological cut (logical "and" of each).
     * @details This function applies a fiducial volume, containment, flash time
     * (BNB), and 1eNp topological cut on the interaction using the logical "and"
     * of each previously defined cut.  
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * flash time, and 1eNp topological cut.
     * @note This cut is intended to be used for the nue analysis.
     */
    template<class T>
        bool all_1eNp_cut(const T & obj) { return fiducial_cut<T>(obj) && containment_cut<T>(obj) && flash_cut<T>(obj) && topological_1eNp_cut<T>(obj); }

    /**
     * @brief Apply a cut to select muon neutrinos
     * @details This function applies a cut to select
     * interactions with a muon neutrino
     * in the final state as defined by the
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has a muon neutrino.
     * @note This cut is intended to be used for the nue analysis.
     */
    bool is_numu(const caf::SRInteractionTruthDLPProxy & obj) { return neutrino(obj) && (obj.pdg_code == 14 || obj.pdg_code == -14); }
    
    /**
     * @brief Apply a cut to select electron neutrinos
     * @details This function applies a cut to select
     * interactions with an electron neutrino
     * in the final state as defined by the
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has an electron neutrino.
     * @note This cut is intended to be used for the nue analysis.
     */
    bool is_nue(const caf::SRInteractionTruthDLPProxy & obj) { return neutrino(obj) && (obj.pdg_code == 12 || obj.pdg_code == -12); }

    /**
     * @brief Apply a cut to select nue CC interactions
     * @details This function applies a cut to select
     * interactions with a nue CC interaction
     * in the final state as defined by the
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has a nue CC interaction.
     * @note This cut is intended to be used for the nue analysis.
     */
    bool is_nue_cc(const caf::SRInteractionTruthDLPProxy & obj) { return is_nue(obj) && (obj.current_type == 0); }

    /**
     * @brief Apply a cut to select nue NC interactions
     * @details This function applies a cut to select
     * interactions with a nue NC interaction
     * in the final state as defined by the
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has a nue NC interaction.
     * @note This cut is intended to be used for the nue analysis.
     */
    bool is_nue_nc(const caf::SRInteractionTruthDLPProxy & obj) { return is_nue(obj) && (obj.current_type == 1); }

    /**
     * @brief Apply a cut to select numu CC interactions
     * @details This function applies a cut to select
     * interactions with a numu CC interaction
     * in the final state as defined by the
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has a numu CC interaction.
     * @note This cut is intended to be used for the nue analysis.
     */
    bool is_numu_cc(const caf::SRInteractionTruthDLPProxy & obj) { return is_numu(obj) && (obj.current_type == 0); }

    /**
     * @brief Apply a cut to select numu NC interactions
     * @details This function applies a cut to select
     * interactions with a numu NC interaction
     * in the final state as defined by the
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has a numu NC interaction.
     * @note This cut is intended to be used for the nue analysis.
     */
    bool is_numu_nc(const caf::SRInteractionTruthDLPProxy & obj) { return is_numu(obj) && (obj.current_type == 1); }

    /**
     * @brief Apply a cut to select the 1eNp signal.
     * @details This function applies a cut on the final state, fiducial volume,
     * and containment of the interaction. This is the "true" 1eNp signal.
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * and 1eNp topological cut.
     * @note This cut is intended to be used for the nue analysis for
     * defining the signal.
     */
    bool signal_1eNp(const caf::SRInteractionTruthDLPProxy & obj) { return is_nue(obj) && fiducial_cut(obj) && containment_cut(obj) && topological_1eNp_cut(obj); }

    /**
     * @brief Apply a cut to select the 1eNp non-signal.
     * @details This function applies a cut on the final state, fiducial volume,
     * and containment of the interaction. This is the "true" 1eNp non-signal
     * (1eNp topology, but not signal).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * and 1eNp topological cut.
     * @note This cut is intended to be used for the nue analysis for
     * defining the non-signal.
     */
    bool non_signal_1eNp(const caf::SRInteractionTruthDLPProxy & obj) { return is_nue(obj) && !(fiducial_cut(obj) && containment_cut(obj) && topological_1eNp_cut(obj)); }

    template <class T>
        bool shower_dedx_cut(const T & obj)
        {
            double dedx = vars::nue::leading_shower_start_dedx(obj);
            if (dedx > SHOWER_DEDX_CUT_VAL)
            {
                return false;
            }
            else
            {
                return true;
            }
        }
    
    template <class T>
        bool shower_vertex_distance_cut(const T & obj)
        {
            double distance = vars::nue::leading_shower_vertex_distance(obj);
            if (distance > SHOWER_VERTEX_DISTANCE_CUT_VAL)
            {
                return false;
            }
            else
            {
                return true;
            }
        }

    template <class T>
        bool shower_directional_spread_cut(const T & obj)
        {
            double spread = vars::nue::leading_shower_directional_spread(obj);
            if (spread > SHOWER_DIRECTIONAL_SPREAD_CUT_VAL)
            {
                return false;
            }
            else
            {
                return true;
            }
        }

    template <class T>
        bool shower_axial_spread_cut(const T & obj)
        {
            double spread = vars::nue::leading_shower_axial_spread(obj);
            if (spread < SHOWER_AXIAL_SPREAD_CUT_VAL)
            {
                return false;
            }
            else
            {
                return true;
            }
        }

    template <class T>
        bool shower_start_straightness_cut(const T & obj)
        {
            double straightness = vars::nue::leading_shower_start_straightness(obj);
            if (straightness < SHOWER_START_STRAIGHTNESS_CUT_VAL)
            {
                return false;
            }
            else
            {
                return true;
            }
        }

    template <class T>
        bool all_shower_cut(const T & obj)
        {
            return shower_dedx_cut(obj) && shower_vertex_distance_cut(obj) && shower_directional_spread_cut(obj) && shower_axial_spread_cut(obj) && shower_start_straightness_cut(obj);
        }

    template <class T>
        bool complete_1eNp_cut(const T & obj)
        {
            return all_1eNp_cut(obj) && all_shower_cut(obj);
        }
}
#endif // CUTS_NUE_H