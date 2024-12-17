// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <string>
#include <vector>
#include <gemmi/metadata.hpp>
#include <gemmi/model.hpp>
#include <VOs/Chain.hpp>

/**
 * @brief namespace for tmdet data transfer objects
 *
 * @namespace Tmdet
 * @namespace DTOs
 */
namespace Tmdet::DTOs {

    /**
     * @brief reading and parsing gemmi chain data
     *        to Chain Value Object
     */
    struct Chain {

        /**
         * @brief get chain value object from gemmi chain data
         * 
         * @param protein 
         * @param chain 
         * @param entityId 
         * @param chainIdx 
         * @return Tmdet::VOs::Chain 
         */
        static Tmdet::VOs::Chain get(gemmi::Structure& protein, gemmi::Chain& chain, int chainIdx);

        /**
         * @brief string representation of the chain
         * 
         * @param chain 
         * @return std::string 
         */
        static std::string toString(const Tmdet::VOs::Chain& chain);

        /**
         * @brief Get the entity index of a chain
         * 
         * @param entities 
         * @param chain 
         * @return int 
         */
        static int getEntityIdx(const std::vector<gemmi::Entity> entities, const std::string entityId);
    };
}
