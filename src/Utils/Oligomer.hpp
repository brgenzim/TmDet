// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <vector>
#include <string>
#include <gemmi/model.hpp>
#include <Types/Oligomer.hpp>

/**
 * @brief namespace for Tmdet
 * @namespace Tmdet
 */
namespace Tmdet {
    
    /**
    * @brief namespace for various utilities
    * @namespace Utils
    */
    namespace Utils {

        struct Oligomer {

            static std::vector<gemmi::Entity> getPolymerEntities(const gemmi::Structure& structure) {
                std::vector<gemmi::Entity> ret;
                for (const auto& entity : structure.entities) {
                    if (entity.entity_type == gemmi::EntityType::Polymer) {
                        // WARNING: entity.subchains contains label asym ids and not auth asym id (that is used in chain name)
                        ret.emplace_back(entity);
                    }
                }
                return ret;
            }

            static std::vector<gemmi::Entity> getHomoOligomerEntities(const gemmi::Structure& structure) {
                auto entities = getPolymerEntities(structure);
                std::vector<gemmi::Entity> result;
                for (auto& entity : entities) {
                    if (entity.subchains.size() > 1) {
                        result.emplace_back(entity);
                    }
                }
                return result;
            }

            static bool isEntityOligomerized(const std::string& entityId, const gemmi::Structure& structure) {
                for (const auto& entity: getHomoOligomerEntities(structure)) {
                    if (entity.name == entityId) {
                        return true;
                    }
                }
                return false;
            }

            static bool isMonomer(const std::vector<gemmi::Entity>& entities) {
                return (entities.size() == 1 && entities[0].subchains.size() == 1);
            }

            static bool isHomoOligomer(const std::vector<gemmi::Entity>& entities) {
                return (entities.size() == 1 && entities[0].subchains.size() > 1);
            }

            static bool isHomoHeteroOligomer(const std::vector<gemmi::Entity>& entities) {
                if (entities.size() > 1) {
                    for( const auto& entity: entities ) {
                        if (entities[0].subchains.size() != entity.subchains.size()) {
                            return false;
                        }
                    }
                    return (entities[0].subchains.size() > 1);
                }
                return false;
            }

            static bool isHeteroWithHomoOligomer(const std::vector<gemmi::Entity>& entities) {
                if (entities.size() > 1) {
                    auto max = entities[0].subchains.size();
                    for( const auto& entity: entities ) {
                        if ( max < entity.subchains.size()) {
                            return true;
                        }
                    }
                }
                return false;
            }

            static Tmdet::Types::Oligomer getOligomerType(gemmi::Structure& structure) {
                Tmdet::Types::Oligomer ret = Tmdet::Types::OligomerType::HETERO_OLIGOMER;
                const auto& entities = getPolymerEntities(structure);
                if (isMonomer(entities)) {
                    ret = Tmdet::Types::OligomerType::MONOMER;
                }
                else if (isHomoOligomer(entities)) {
                    ret = Tmdet::Types::OligomerType::HOMO_OLIGOMER;
                }
                else if (isHomoHeteroOligomer(entities)) {
                    ret = Tmdet::Types::OligomerType::HOMO_HETERO_OLIGOMER;
                }
                else if (isHeteroWithHomoOligomer(entities)) {
                    ret = Tmdet::Types::OligomerType::HETERO_WITH_HOMO_OLIGOMER;
                }
                return ret;
            }
        };
    };
}
