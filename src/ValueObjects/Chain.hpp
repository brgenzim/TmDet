#ifndef __TMDET_VALUE_OBJECTS_CHAIN__
#define __TMDET_VALUE_OBJECTS_CHAIN__

#include <string>
#include <vector>
#include <Types/Chain.hpp>
#include <ValueObjects/Region.hpp>
#include <ValueObjects/Residue.hpp>
#include <gemmi/model.hpp>

namespace Tmdet::ValueObjects {

    struct Chain {
        std::string id;
        std::string entityId;
        bool selected = true;
        int numtm = 0;
        std::string seq;
        std::vector<Tmdet::ValueObjects::Region> regions;
        Tmdet::Types::Chain type = Tmdet::Types::ChainType::UNK;
        gemmi::Chain& gemmi;
        std::vector<Residue> residues;
        int idx = 0;
        int length = 0;

        explicit Chain(gemmi::Chain& _gemmi) : 
            id(_gemmi.name), 
            gemmi(_gemmi) {
                for (const auto &residue : gemmi.residues) {
                    if (!residue.atoms.empty()) {
                        this->entityId = residue.entity_id;
                        break;
                    }
                }
        }

        ~Chain()=default;
    };
}

#endif
