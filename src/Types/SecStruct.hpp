#ifndef __TMDET_TYPES_SECSTRUCT__
#define __TMDET_TYPES_SECSTRUCT__

#include <unordered_map>
#include <string>

using namespace std;

namespace Tmdet::Types {

    struct SecStruct {
        string name;
        char code;
    };

    const unordered_map<char, SecStruct> SecStructs = {
        { 'H', {
                "Helix", 
                'H',
            }
        },
        { 'G', {
                "GHelix", 
                'G',
            }
        },
        { 'I', {
                "IHelix", 
                'I',
            }
        },
        { 'T', {
                "Turn",
                'T',
            }
        },
        { 'B', {
                "Bend",
                'B',
            }
        },
        { 'E', {
                "Extended",
                'E',
            }
        },
        { 'S', {
                "SBend",
                'S',
            }
        },
        { '-', {
                "Unknown",
                '-',
            }
        },
    };
}

#endif
