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
        int beg = 0;

        /**
         * @brief start of the region in PDB by auth
         */
        int beg_auth_seq_id = 0;

        /**
         * @brief insertion code of start residue
         */
        char beg_auth_seq_icode = ' ';

        /**
         * @brief start of the region in PDB by label
         */
        int beg_label_seq_id = 0;

        /**
         * @brief end of the region in sequence position
         */
        int end = 0;
        
        /**
         * @brief end of the region in PDB by auth
         */
        int end_auth_seq_id = 0;
        
        /**
         * @brief insertion code of start residue
         */
        char end_auth_seq_icode = ' ';

        /**
         * @brief end of the region in PDB by label
         */
        int end_label_seq_id = 0;
        
        /**
         * @brief region type
         */
        Tmdet::Types::Region type = Tmdet::Types::RegionType::UNK;
    };
}
