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
     * @brief definition of oligomer type
     */
    struct Oligomer {
        /**
         * @brief name of the oligomer type
         */
        std::string name;

        /**
         * @brief description of oligomer type
         */
        std::string description;
    };

    namespace OligomerType {
        const Oligomer MONOMER = {
            "Monomer", 
            "There is only one chain in the protein"
        };
        const Oligomer HOMO_OLIGOMER = {
            "HomoOligomer",
            "There is only one chain type but several times"
        };
        const Oligomer HETERO_OLIGOMER = {
            "HeteroOligomer",
            "There several different chains in the protein"
        };
        const Oligomer HOMO_HETERO_OLIGOMER = {
            "HomoHeteroOligomer",
            "There several different chains in the protein but multiple times"
        };
        const Oligomer HETERO_WITH_HOMO_OLIGOMER = {
            "HeteroWithHomoOligomer",
            "There are several different chains, but one of them appears more than two times"
        };
    }

    const std::unordered_map<std::string, const Oligomer> Oligomers = {
        { "Monomer", OligomerType::MONOMER },
        { "HomoOligomer", OligomerType::HOMO_OLIGOMER },
        { "HeteroOligomer", OligomerType::HETERO_OLIGOMER },
        { "HomoHeteroOligomer", OligomerType::HOMO_HETERO_OLIGOMER }
    };

}
