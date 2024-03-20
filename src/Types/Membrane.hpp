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

    namespace MembraneType {
        const Membrane PLAIN = {
            "Plain", 
            "Simple plain membrane used most of membrane proteins"
        };
        const Membrane CURVED = {
            "Curved",
            "Curved membrane represented by a radius"
        };
    }

    const unordered_map<string, const Membrane> Membranes = {
        { "Plain", MembraneType::PLAIN },
        { "Curved", MembraneType::CURVED }
    };

}

#endif
