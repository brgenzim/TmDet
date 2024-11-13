#include <string>

#include <ValueObjects/Residue.hpp>

namespace Tmdet::ValueObjects {

    int Residue::resn() const {
        return gemmi.seqid.num.value;
    }

    void Residue::setNumberOfAtoms() {
        for(const auto& atom: atoms) {
            if (type.atoms.contains(atom.gemmi.name)) {
                if (type.atoms.at(atom.gemmi.name).bb) {
                    nba++;
                }
                else {
                    nsa++;
                }
            }
        }
    }

    bool Residue::hasAllSideChainAtoms() const {
        return (nsa == type.nsa);
    }

    bool Residue::hasOnlyBackBoneAtoms() const {
        return (nsa == 0 && nba > 0);
    }

    bool Residue::hasAllAtoms() const {
        return atoms.size() >= (type.atoms.size() - 1);
    }

    gemmi::Vec3 Residue::centre() const {
        gemmi::Vec3 ret(0,0,0);
        for(const auto& atom: atoms) {
            ret += atom.gemmi.pos;
        }
        ret /= (double)atoms.size();
        return ret;
    }
}