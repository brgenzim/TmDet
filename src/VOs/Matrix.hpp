// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <string>
#include <VOs/TMatrix.hpp>

/**
 * @brief namespace for value objects
 */
namespace Tmdet::VOs {

    /**
     * @brief transformation matrix and source chain identifiers
     *        for generating Bio Matrix
     */
    struct Matrix {

        /**
         * @brief identifier of the matrix
         */
        int id;

        /**
         * @brief chain identifier of the chain copied
         */
        std::string sourceChainId;

        /**
         * @brief chain identifier of the new chain
         */
        std::string newChainId;

        /**
         * @brief transformation matrix for generating coordinates of
         *        new chain from the source chain
         */
        TMatrix tmatrix;
    };
}
