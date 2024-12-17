// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <string>
#include <VOs/Chain.hpp>

/**
 * @brief namespace for tmdet data transfer objects
 *
 * @namespace Tmdet
 * @namespace DTOs
 */
namespace Tmdet::DTOs {

    struct Dssp {

        /**
         * @brief Get secondary structure of a chain as string
         * 
         * @param chain 
         * @return std::string 
         */
        static std::string getSecondaryStructure(const Tmdet::VOs::Chain& chain);
    };
}
