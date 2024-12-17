// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <Types/Region.hpp>

/**
 * @brief namespace for value objects
 */
namespace Tmdet::VOs {

    struct SeqPosition {
        int authId;
        char authIcode;
        int labelId;
        int idx;
    };

    /**
     * @brief description of a region in sequence
     */
    struct Region {

        /**
         * @brief start of the region in Tmdet ValueObject Chain
         */
        SeqPosition beg;

        /**
         * @brief end of the region in Tmdet ValueObject Chain
         */
        SeqPosition end;
        
        /**
         * @brief region type
         */
        Tmdet::Types::Region type = Tmdet::Types::RegionType::UNK;

    };
}
