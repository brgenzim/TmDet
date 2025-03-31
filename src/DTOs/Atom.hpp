// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <string>
#include <gemmi/model.hpp>
#include <VOs/Atom.hpp>

/**
 * @brief namespace for tmdet data transfer objects
 *
 * @namespace Tmdet
 * @namespace DTOs
 */
namespace Tmdet::DTOs {

    /**
     * @brief writing atom data
     */
    struct Atom {

        /**
         * @brief string representation of the atom
         * 
         * @param atom
         */
        static std::string toString(const Tmdet::VOs::Atom& atom);
    };
}
