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

    const unordered_map<string, Chain> Chains = {
        { "alpha", {
                "alpha", 
                "Chain containing alpha helical transmembrane segment(s)"
            }
        },
        { "beta", {
                "beta", 
                "Chain containing beta barrel transmembrane domain"
            }
        },
        { "non_tm", {
                "non_tm",
                "Chain without any transmembrane region"
            }
        }
    };

}

#endif
