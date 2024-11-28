#include <string>

#include <ValueObjects/Residue.hpp>
#include <gemmi/elem.hpp>

namespace Tmdet::ValueObjects {

    void Residue::setProperties(std::string name) {
        if (Tmdet::Types::Residues.contains(name)) {
            type = Tmdet::Types::Residues.at(name);
        }
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
        return (nsa < type.nsa && nba > 0);
    }

    bool Residue::hasAllAtoms() const {
        return atoms.size() >= (type.atoms.size() - 1);
    }

    const gemmi::Atom* Residue::getCa() const {
        const gemmi::Atom* ca = gemmi.get_ca();
        if (ca == nullptr) {
            ca = gemmi.find_atom("CB", '*', gemmi::El::C);
        }
        if (ca == nullptr) {
            ca = gemmi.get_c();
        }
        if (ca == nullptr) {
            ca = gemmi.get_n();
        }
        return ca;
    }

    void Residue::transform(Tmdet::ValueObjects::TMatrix& tmatrix) {
        for(auto& atom: atoms) {
            atom.transform(tmatrix);
        }
    }

}