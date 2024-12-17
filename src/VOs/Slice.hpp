// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

namespace Tmdet::VOs {

    /**
     * @brief description of a slice
     */
    struct Slice {
        
        /**
         * @brief apolar surface
         */
        double apol = 0.0;

        /**
         * @brief outside water accessible surface of the atoms in the slice
         */
        double surf = 0.0;

        /**
         * @brief ratio of straight element in the slice
         */
        double straight = 0.0;

        /**
         * @brief number of C alpha atoms in the slice
         */
        int numCa = 0;
        
        /**
         * @brief number of sec structure ends
         */
        int ssEnd = 0;

        /**
         * @brief raw, not smoothed q value
         */
        double rawQ = 0.0;

        /**
         * @brief calculated qValue for the slice
         */
        double qValue = 0.0;
    };
}
