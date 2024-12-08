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
