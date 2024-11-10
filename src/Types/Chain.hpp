#ifndef __TMDET_TYPES_CHAIN__
#define __TMDET_TYPES_CHAIN__

#include <unordered_map>
#include <string>

namespace Tmdet::Types {

    struct Chain {
        std::string name;
        std::string description;

        bool operator == (const Chain& other) const {
            return name == other.name;
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

#endif
