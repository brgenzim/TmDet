// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <unordered_map>
#include <string>

/**
 * @brief namespace of tmdet types
 * @namespace Tmdet
 * @namespace Types
 */
namespace Tmdet::Types {

    /**
     * @brief definition of membrane type
     */
    struct Membrane {
        /**
         * @brief name of the membrane type
         */
        std::string name;

        /**
         * @brief description of the membrane type
         */
        std::string description;

        /**
         * @brief check if two membrane types are equal
         * 
         * @param other 
         * @return true 
         * @return false 
         */
        bool operator==(Membrane other) {
            return name == other.name;
        }

        /**
         * @brief check if the membrane type is plane
         * 
         * @return true 
         * @return false 
         */
        bool isPlane() const {
            return name == "Plane";
        }
    };

    namespace MembraneType {
        const Membrane PLANE = {
            "Plane", 
            "Simple plain membrane used most of membrane proteins"
        };
        const Membrane CURVED = {
            "Curved",
            "Curved membrane represented by one radius"
        };
    }

    const std::unordered_map<std::string, const Membrane> Membranes = {
        { "Plane", MembraneType::PLANE },
        { "Curved", MembraneType::CURVED }
    };

}
