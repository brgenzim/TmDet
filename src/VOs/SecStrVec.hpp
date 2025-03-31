// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <string>
#include <gemmi/math.hpp>
#include <Types/SecStruct.hpp>

/**
 * @brief namespace for value objects
 */
namespace Tmdet::VOs {

    /**
     * @brief description of a secondary structure vector
     */
    struct SecStrVec {
        /**
         * @brief type of the secondary structure element
         */
        Tmdet::Types::SecStruct type;

        /**
         * @brief three dimensional coordinate of the begin of the vector
         */
        gemmi::Vec3 begin;

        /**
         * @brief three dimensional coordinate of the end of the vector
         */
        gemmi::Vec3 end;

        /**
         * @brief index of the chain in the protein value object
         */
        int chainIdx;

        /**
         * @brief index of the first residues in the chain value object
         */
        int begResIdx;

        /**
         * @brief index of the last residues in the chain value object
         */
        int endResIdx;

        /**
         * @brief beta sheet index
         */
        int sheetIdx = -1;

        /**
         * @brief beta barrel index
         */
        int barrelIdx = -1;
    };
}
