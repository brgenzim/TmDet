#ifndef __TMDET_VALUE_OBJECTS_RESIDUE__
#define __TMDET_VALUE_OBJECTS_RESIDUE__

#include <string>
#include <vector>
#include <any>
#include <unordered_map>
#include <Types/SecStruct.hpp>
#include <ValueObjects/Atom.hpp>
#include <ValueObjects/HBond.hpp>
#include <gemmi/model.hpp>

namespace Tmdet::ValueObjects {

    struct Residue {
        gemmi::Residue& gemmi;
        std::vector<Atom> atoms;
        double surface;
        Tmdet::Types::SecStruct ss = Tmdet::Types::SecStructs.at('-');
        HBond hbond1;
        HBond hbond2;
        // Represents order of residues after alignment.
        // Use this to measure residues distance in the sequence.
        int idx;
        int chainIdx;
        std::unordered_map<std::string,std::any> temp;

        explicit Residue(gemmi::Residue& _gemmi) : gemmi(_gemmi) {}

        ~Residue()=default;

        int resn() const {
            return gemmi.seqid.num.value;
        }
    };
}

#endif
