#pragma once

#include <unordered_map>
#include <string>

namespace Tmdet::Types {

    struct Membrane {
        std::string name;
        std::string description;
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

    const std::unordered_map<std::string, const Membrane> Membranes = {
        { "Plain", MembraneType::PLAIN },
        { "Curved", MembraneType::CURVED }
    };

}
