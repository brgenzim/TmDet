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
        int id;
        int n;
        vector<string> chids;
    };

    struct Oligomer {
        
        static vector<_chains> getNumberOfChains(gemmi::cif::Document doc) {
            vector<_chains> ret;
            for (auto row : doc.sole_block().find("_entity.",{"id", "type", "pdbx_number_of_molecules"})) {
                if (row[1] == "polymer") {
                    ret.emplace_back(_chains({stoi(row[0]),stoi(row[2]),{}}));
                }
            }
            return ret;
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

        static Tmdet::Types::Oligomer getOligomerType(gemmi::cif::Document doc) {
            Tmdet::Types::Oligomer ret = Tmdet::Types::OligomerType::HETERO_OLIGOMER;
            vector<_chains> chains = getNumberOfChains(doc);
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