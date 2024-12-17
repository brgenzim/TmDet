// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <unordered_map>
#include <string>

/**
 * @brief namespace of tmdet types
 * @namespace Tmdet
 * @namespace Types
 */
namespace Tmdet::Types {

    /**
     * @brief description of a chain type
     */
    struct Chain {
        /**
         * @brief name of chain type
         */
        std::string name;

        /**
         * @brief description of chain type
         */
        std::string description;

        /**
         * @brief check if two chain types are equal
         * 
         * @param other 
         * @return true 
         * @return false 
         */
        bool operator == (const Chain& other) const {
            return name == other.name;
        }

        /**
         * @brief check is chain type is alpha
         * 
         * @return true 
         * @return false 
         */
        bool isAlpha() const {
            return name == "alpha";
        }

        /**
         * @brief check is chain type is beta
         * 
         * @return true 
         * @return false 
         */
        bool isBeta() const {
            return name == "beta";
        }
    };

    namespace ChainType {
        const Chain ALPHA = {
            "alpha", 
            "Chain containing alpha helical transmembrane segment(s)"
        };
        const Chain BETA = {
            "beta", 
            "Chain containing beta barrel transmembrane domain"
        };
        const Chain NON_TM = {
            "non_tm",
            "Chain without any transmembrane region"
        };
        const Chain LOW_RES = {
            "low_res",
            "Chain only with backbone atoms"
        };
        const Chain UNK = {
            "unknown",
            "Chain type is not determined"
        };
        const Chain NOT_SELECTED = {
            "not_selected",
            "Chain is unselected"
        };
    }

    const std::unordered_map<std::string, const Chain> Chains = {
        { "alpha", ChainType::ALPHA },
        { "beta", ChainType::BETA },
        { "non_tm", ChainType::NON_TM },
        { "low_res", ChainType::LOW_RES},
        { "unknown", ChainType::UNK},
        { "not_selected", ChainType::NOT_SELECTED}
    };

}
