#include <vector>
#include <any>
#include <array>
#include <iostream>
#include <cmath>
#include <limits>
#include <gemmi/model.hpp>
#include <Utils/Fragment.hpp>

namespace Tmdet::Utils {

    void Fragment::run() {
        for (auto& chainVO: tmdetVO.chains) {
            for (auto& residueVO: chainVO.residues) {
                residueVO.temp.insert({"fragment",std::any_cast<int>(-1)});
                residueVO.temp.insert({"neighbors",std::any_cast<std::vector<_cr>>(getNeighbors(residueVO))});
            }
        }
        int fr = 0;
        for (auto& chainVO: tmdetVO.chains) {
            for (auto& residueVO: chainVO.residues) {
                if (std::any_cast<int>(residueVO.temp.at("fragment")) == -1 &&
                    std::any_cast<std::vector<_cr>>(residueVO.temp.at("neighbors")).size() > 0 ) {
                    addToFragment(residueVO,fr++);
                }
            }
        }
    }

    void Fragment::addToFragment(Tmdet::ValueObjects::Residue& residueVO, int fr) {
        residueVO.temp.at("fragment") = std::any_cast<int>(fr);
        std::cout << "addToFragment" << fr << ": " << residueVO.chainIdx << "," << residueVO.idx << std::endl;
        for (auto n: std::any_cast<std::vector<_cr>>(residueVO.temp.at("neighbors"))) {
            if (std::any_cast<int>(tmdetVO.chains[n.chain_idx].residues[n.residue_idx].temp.at("fragment")) == -1) {
                addToFragment(tmdetVO.chains[n.chain_idx].residues[n.residue_idx], fr);
            }
        }
    }

    std::vector<_cr> Fragment::getNeighbors(Tmdet::ValueObjects::Residue& residueVO) {
        std::vector<_cr> ret;
        std::vector<_cr> empty;
        const gemmi::Atom* ca_atom = residueVO.gemmi.get_ca();
        if (ca_atom) {
            for (auto mark: tmdetVO.neighbors.find_neighbors(*ca_atom, 3, 8)) {
                if (tmdetVO.chains[mark->chain_idx].residues[mark->residue_idx].atoms[mark->atom_idx].gemmi.name == "CA") {
                    _cr cr;
                    cr.chain_idx = mark->chain_idx;
                    cr.residue_idx = mark->residue_idx;
                    ret.emplace_back(cr);
                }
            }
        }
        return (ret.size()>4?ret:empty);
    }

    bool Fragment::sameChain(gemmi::Chain* chain, int chainIdx) {
        return (chain->name == tmdetVO.chains[chainIdx].id );
    }

}
