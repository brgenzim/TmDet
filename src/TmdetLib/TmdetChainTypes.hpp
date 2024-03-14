#ifndef __UNITMP_TMDETLIB_TMDET_CHAIN_TYPES__
#define __UNITMP_TMDETLIB_TMDET_CHAIN_TYPES__

#include <unordered_map>

namespace UniTmp::TmdetLib {

    struct TmdetChainType {
        char *name;
        char *description;
    };

    namespace TmdetChainTypes {
        const TmdetChainType ALPHA {
            "alpha", 
            "Chain containing alpha helical transmembrane segment(s)"
        };
        const TmdetChainType BETA {
            "beta", 
            "Chain containing beta barrel transmembrane domain"
        };
        const TmdetChainType NONTM {
            "non_tm",
            "Chain without any transmembrane region"
        };

        const std::unordered_map<const char*, TmdetChainType> all {
            {"alpha", ALPHA},
            {"beta", BETA},
            {"non_tm", NONTM}
        };
    };

}

#endif
