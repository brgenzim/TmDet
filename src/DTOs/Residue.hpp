// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <string>
#include <gemmi/model.hpp>
#include <VOs/Residue.hpp>

/**
 * @brief namespace for tmdet data transfer objects
 *
 * @namespace Tmdet
 * @namespace DTOs
 */
namespace Tmdet::DTOs {

    /**
     * @brief reading and writing residue data and parsing
     *        it to Residue Value Object
     */
    struct Residue {

        /**
         * @brief get data from gemmi residue
         * 
         * @param residue 
         * @param chainIdx 
         * @param residueIdx 
         * @return Tmdet::VOs::Residue 
         */
        static Tmdet::VOs::Residue get(gemmi::Residue& residue, int chainIdx, int residueIdx);

        /**
         * @brief string representation of the residue
         * 
         * @param residueVO 
         */
        static std::string toString(const Tmdet::VOs::Residue& residueVO);
    };
}
