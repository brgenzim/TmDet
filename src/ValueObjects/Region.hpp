#pragma once

#include <string>
#include <Types/Region.hpp>

/**
 * @brief namespace for value objects
 */
namespace Tmdet::ValueObjects {

    /**
     * @brief description of a region in sequence
     */
    struct Region {

        /**
         * @brief start of the region in Tmdet ValueObject Chain
         */
        unsigned int beg = 0;

        /**
         * @brief end of the region in Tmdet ValueObject Chain
         */
        unsigned int end = 0;
        
        /**
         * @brief region type
         */
        Tmdet::Types::Region type = Tmdet::Types::RegionType::UNK;
    };
}
