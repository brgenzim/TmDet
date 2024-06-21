#ifndef __TMDET_VALUE_OBJECTS_CHAIN__
#define __TMDET_VALUE_OBJECTS_CHAIN__

#include <string>
#include <vector>
#include <Types/Chain.hpp>
#include <ValueObjects/Region.hpp>
#include <ValueObjects/Residue.hpp>
#include <gemmi/model.hpp>

using namespace std;

namespace Tmdet::ValueObjects {

    struct Chain {
        string id;
        string entityId;
        bool selected;
        int numtm;
        string seq;
        vector<Tmdet::ValueObjects::Region> regions;
        Tmdet::Types::Chain type;
        gemmi::Chain& gemmi;
        vector<Residue> residues;
        int idx;
        int length;
        Chain(gemmi::Chain& _gemmi) : gemmi(_gemmi) {
            this->id = gemmi.name;
            for (auto residue : gemmi.residues) {
                if (residue.atoms.size() > 0) {
                    this->entityId = residue.entity_id;
                    break;
                }
            }
        }
        ~Chain() {}
    };
}

#endif