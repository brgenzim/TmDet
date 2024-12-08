#pragma once

#include <unordered_map>
#include <string>

/**
 * @brief namespace of types
 * @namespace Tmdet
 * @namespace Types
 */
namespace Tmdet::Types {

    struct Membrane {
        std::string name;
        std::string description;

        bool operator==(Membrane other) {
            return name == other.name;
        }

        bool isPlane() const {
            return name == "Plane";
        }
    };

    namespace MembraneType {
        const Membrane PLANE = {
            "Plane", 
            "Simple plain membrane used most of membrane proteins"
        };
        const Membrane BLENDED = {
            "Blended",
            "Blended membrane represented by one radius"
        };
    }

    const std::unordered_map<std::string, const Membrane> Membranes = {
        { "Plane", MembraneType::PLANE },
        { "Blended", MembraneType::BLENDED }
    };

}
