#ifndef __TMDET_TYPES_REGION__
#define __TMDET_TYPES_REGION__

#include <unordered_map>
#include <string>

using namespace std;

namespace Tmdet::Types {

    struct Region {
        string name;
        char code;
        string description;
    };

    const unordered_map<char, Region> Regions = {
        { 'M', {
                "transmembrane", 
                'M',
                "Region crossing the membrane"
            }
        },
        { 'H', {
                "Alpha helical transmembrane", 
                'H',
                "Alpha helical region crossing the membrane"
            }
        },
        { 'B', {
                "Beta barrel strain", 
                'B',
                "Beta strain crossing the membrane"
            }
        },
        { '1', {
                "side1",
                '1',
                "Side one of the membrane"
            }
        },
        { '2', {
                "side2",
                '2',
                "Other side of the membrane"
            }
        },
        { 'L', {
                "loop",
                'L',
                "Re-entrant membrane loop"
            }
        },
        { 'F', {
                "interfacial helix",
                'F',
                "Interfacial helix on membrane surface"
            }
        },
        { 'I', {
                "inside",
                'I',
                "Beta barrel inside element"
            }
        },
        { 'N', {
                "chain begin",
                'N',
                "N-terminal chain segment"
            }
        },
        { 'C', {
                "chain end",
                'C',
                "C-terminal chain segment"
            }
        },
        { 'U', {
                "unkown",
                'U',
                "Unknow localization"
            }
        }
    };
}

#endif
