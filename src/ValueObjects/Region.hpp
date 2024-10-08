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
         * @brief start of the region in sequence position
         */
        int beg;

        /**
         * @brief start of the region in PDB numbering
         */
        int rbeg;

        /**
         * @brief insertion code in PDB structure of 
         *        region start position
         */
        std::string begi;

        /**
         * @brief end of the region in sequence position
         */
        int end;
        
        /**
         * @brief end of the region in PDB numbering
         */
        int rend;
        
        /**
         * @brief insertion code in PDB structure of 
         *        region end position
         */
        std::string endi;
        
        /**
         * @brief region type
         */
        Tmdet::Types::Region type;
    };
}
