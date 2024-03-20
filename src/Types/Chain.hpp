#ifndef __TMDET_TYPES_CHAIN__
#define __TMDET_TYPES_CHAIN__

#include <unordered_map>
#include <string>

using namespace std;

namespace Tmdet::Types {

    struct Chain {
        string name;
        string description;
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
        const Chain NONTM = {
            "non_tm",
            "Chain without any transmembrane region"
        };
    }

    const unordered_map<string, const Chain> Chains = {
        { "alpha", ChainType::ALPHA },
        { "beta", ChainType::BETA },
        { "non_tm", ChainType::NONTM }
    };

}

#endif
