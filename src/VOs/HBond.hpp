// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

namespace Tmdet::VOs {

    /**
     * @brief simple representation of a hydrogen bond
     */
    struct HBond {

        /**
         * @brief calculated energy of the hydrogen bond
         */
        double energy = 1e30;

        /**
         * @brief chain index of the akceptor residue
         */
        int toChainIdx = -1;

        /**
         * @brief residue index of the akceptor
         */
        int toResIdx = -1;

        int fromResIdx = -1;
    };
}
