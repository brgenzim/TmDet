#pragma once

#include <string>
#include <vector>
#include <gemmi/metadata.hpp>
#include <gemmi/model.hpp>
#include <ValueObjects/Chain.hpp>

/**
 * @brief namespace for data transfer objects
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
         * @param entityIdx
         * @return Tmdet::ValueObjects::Chain 
         */
        static Tmdet::ValueObjects::Chain get(gemmi::Structure& protein, gemmi::Chain& chain, int chainIdx, int entityIdx);

        /**
         * @brief string representation of the chain
         * 
         * @param chain 
         * @return std::string 
         */
        static std::string toString(const Tmdet::ValueObjects::Chain& chain);

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
