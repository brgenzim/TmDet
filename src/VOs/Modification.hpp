// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <string>

/**
 * @brief namespace for value objects
 */
namespace Tmdet::VOs {

    /**
     * @brief description of modifications
     */
    struct Modification {

        /**
         * @brief date of the modification
         */
        std::string date;

        /**
         * @brief short description of changes
         */
        std::string descr;
    };
}
