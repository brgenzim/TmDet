#ifndef __TMDET_UTILS_OLIGOMER__
#define __TMDET_UTILS_OLIGOMER__

#include <array>
#include <string>
#include <any>
#include <gemmi/model.hpp>
#include <Types/Oligomer.hpp>

using namespace std;

namespace Tmdet::Utils {

    struct _chains {
        string id; // entity id from _entity category
        int n;
        vector<string> chids;
    };

    struct Oligomer {

        static vector<_chains> getNumberOfChains(gemmi::Structure structure) {
            vector<_chains> ret;
            for (const auto& entity : structure.entities) {
                if (entity.entity_type == gemmi::EntityType::Polymer) {
                    // WARNING: entity.subchains contains label asym ids!
                    ret.emplace_back(_chains(entity.name, entity.subchains.size(), entity.subchains));
                }
            }

            return ret;
        }

        // TODO: maybe it is unneccessary
        // static vector<vector<string>> getHomoOligomerChains(gemmi::Structure structure) {
        //     auto entities = getNumberOfChains(structure);
        //     vector<vector<string>> result;
        //     for (auto& chains : entities) {
        //         if (chains.chids.size() > 1) {
        //             result.emplace_back(chains.chids);
        //         }
        //     }
        //     return result;
        // }

        static bool isEntityOligomerized(string entityId, vector<_chains>& chains) {
            bool result = false;
            auto ch = chains.begin();
            for (; ch != chains.end(); ch++) {
                if (ch->id == entityId) {
                    result = ch->n > 1;
                    break;
                }
            }
            return result;
        }

        static bool isMonomer(vector<_chains>& chains) {
            return (chains.size() == 1 && chains[0].n == 1);
        }

        static bool isHomoOligomer(vector<_chains>& chains) {
            return (chains.size() == 1 && chains[0].n > 1);
        }

        static bool isHomoHeteroOligomer(vector<_chains>& chains) {
            if (chains.size() > 1) {
                for( const auto& chain: chains ) {
                    if (chains[0].n != chain.n) {
                        return false;
                    }
                }
                return (chains[0].n > 1);
            }
            return false;
        }

        static bool isHeteroWithHomoOligomer(vector<_chains>& chains) {
            int max = 0;
            if (chains.size() > 1) {
                for( const auto& chain: chains ) {
                    if (chain.n > max) {
                        max = chain.n;
                    }
                }
            }
            return (max > 1);
        }

        static Tmdet::Types::Oligomer getOligomerType(gemmi::Structure structure) {
            Tmdet::Types::Oligomer ret = Tmdet::Types::OligomerType::HETERO_OLIGOMER;
            vector<_chains> chains = getNumberOfChains(structure);
            if (isMonomer(chains)) {
                ret = Tmdet::Types::OligomerType::MONOMER;
            }
            else if (isHomoOligomer(chains)) {
                ret = Tmdet::Types::OligomerType::HOMO_OLIGOMER;
            }
            else if (isHomoHeteroOligomer(chains)) {
                ret = Tmdet::Types::OligomerType::HOMO_HETERO_OLIGOMER;
            }
            else if (isHeteroWithHomoOligomer(chains)) {
                ret = Tmdet::Types::OligomerType::HETERO_WITH_HOMO_OLIGOMER;
            }
            return ret;
        }

    };
}
#endif