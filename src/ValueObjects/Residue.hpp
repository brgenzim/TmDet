#ifndef __TMDET_VALUE_OBJECTS_RESIDUE__
#define __TMDET_VALUE_OBJECTS_RESIDUE__

#include <string>
#include <vector>
#include <any>
#include <iostream>
#include <unordered_map>
#include <Types/Residue.hpp>
#include <Types/SecStruct.hpp>
#include <ValueObjects/Atom.hpp>
#include <ValueObjects/Chain.hpp>
#include <ValueObjects/HBond.hpp>
#include <gemmi/model.hpp>

namespace Tmdet::ValueObjects {

    struct Chain;

    struct Residue {


        gemmi::Residue& gemmi;
        std::vector<Atom> atoms;
        double surface;
        Tmdet::Types::Residue type;
        Tmdet::Types::SecStruct ss = Tmdet::Types::SecStructs.at('-');
        HBond hbond1;
        HBond hbond2;
        // Represents order of residues after alignment.
        // Use this to measure residues distance in the sequence.
        int idx;
        int chainIdx;
        Tmdet::ValueObjects::Chain& parentChain;
        std::unordered_map<std::string,std::any> temp;

        explicit Residue(gemmi::Residue& _gemmi, Chain& chainVO) : 
            gemmi(_gemmi),
            parentChain(chainVO) {
                if (Tmdet::Types::Residues.contains(_gemmi.name)) {
                    type = Tmdet::Types::Residues.at(_gemmi.name);
                }
        }
        int resn() const;
        void setType();

        bool hasAllSideChainAtoms() const;
        bool hasAllAtoms() const;
    };
}

#endif
