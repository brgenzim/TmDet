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

using namespace std;

namespace Tmdet::ValueObjects {

    struct Residue {
        gemmi::Residue& gemmi;
        vector<Atom> atoms;
        double surface;
        Tmdet::Types::SecStruct ss = Tmdet::Types::SecStructs.at('-');
        HBond hbond1;
        HBond hbond2;
        int idx;
        int chainIdx;
        unordered_map<string,any> temp;

        Residue(gemmi::Residue& _gemmi) : gemmi(_gemmi) {}

        // copy constructor
        Residue(const Residue& other) :
            gemmi(other.gemmi),
            atoms(other.atoms),
            surface(other.surface),
            ss(other.ss),
            hbond1(other.hbond1),
            hbond2(other.hbond2),
            idx(other.idx),
            chainIdx(other.chainIdx),
            temp(other.temp) {}

        ~Residue() {}
        int resn() {
            return gemmi.seqid.num.value;
        }
    };
}

#endif
