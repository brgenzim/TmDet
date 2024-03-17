#ifndef __TMDET_TYPES_MEMBRANE__
#define __TMDET_TYPES_MEMBRANE__

#include <unordered_map>
#include <string>

using namespace std;

namespace Tmdet::Types {

    struct Membrane {
        string name;
        string description;
    };

    const unordered_map<string, Membrane> Membranes = {
        { "Plain", {
                "Plain", 
                "Simple plain membrane used most of membrane proteins"
            }
        },
        { "Curved", {
                "Curved",
                "Curved membrane represented by a radius"
            }
        },
        { "Double", {
                "Double",
                "Double membrane for bacterial proteins"
            }
        }
    };

}

#endif
