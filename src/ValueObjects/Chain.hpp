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
        bool selected;
        int numtm;
        std::string seq;
        std::vector<Tmdet::ValueObjects::Region> regions;
        Tmdet::Types::Chain type;
        gemmi::Chain& gemmi;
        std::vector<Residue> residues;
        int idx;
        int length;

        explicit Chain(gemmi::Chain& _gemmi) : gemmi(_gemmi) {
            this->id = gemmi.name;
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
