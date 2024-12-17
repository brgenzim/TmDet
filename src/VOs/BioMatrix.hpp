// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <string>
#include <vector>
#include <VOs/Matrix.hpp>

/**
 * @brief namespace for value objects
 */
namespace Tmdet::VOs {

    /**
     * @brief the assembly in the membrane
     */
    struct BioMatrix {

        /**
         * @brief transformation and chain ids
         */
        std::vector<Matrix> matrices;

        /**
         * @brief deleted chains (e.g. antibodies, nanobodies, biologically 
         *        not relevant assemblies, etc)
         */
        std::vector<std::string> deletedChainIds;
    };
}
