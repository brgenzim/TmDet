#include <string>

#include <ValueObjects/Residue.hpp>

namespace Tmdet::ValueObjects {

    

    int Residue::resn() const {
        return gemmi.seqid.num.value;
    }

    bool Residue::hasAllSideChainAtoms() const {
        if (gemmi.name == "GLY") {
            return true;
        }
        unsigned int numSideChainAtoms = 0;
        for(const auto& [key,data]: type.atoms) {
            numSideChainAtoms += (data.bb?0:1);
        }
        return numSideChainAtoms >= (type.atoms.size() - 5);
    }

    bool Residue::hasAllAtoms() const {
        return atoms.size() >= (type.atoms.size() - 1);
    }
}