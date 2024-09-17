#include <vector>
#include <any>
#include <array>
#include <iostream>
#include <cmath>
#include <limits>
#include <gemmi/model.hpp>
#include <Utils/Fragment.hpp>
#include <Utils/Graph.hpp>

namespace Tmdet::Utils {

#define NODE_ID(c,r) any_cast<int>(tmdetVO.chains[c].residues[r].temp.at("node_index"))

    void Fragment::run() {
        auto crs = getCAlphaNetwork();
        auto clusters = createFragments(crs.size());
        writeBackFragmentInfoToStructure(clusters, crs);
        freeTempValues();
    }

    std::vector<_cr> Fragment::getNeighbors(Tmdet::ValueObjects::Residue& residueVO) {
        std::vector<_cr> ret;
        std::vector<_cr> empty;
        const gemmi::Atom* ca_atom = residueVO.gemmi.get_ca();
        if (ca_atom) {
            for (auto mark: tmdetVO.neighbors.find_neighbors(*ca_atom, 3, 9)) {
                if (tmdetVO.chains[mark->chain_idx].residues[mark->residue_idx].atoms[mark->atom_idx].gemmi.name == "CA") {
                    _cr cr;
                    cr.chain_idx = mark->chain_idx;
                    cr.residue_idx = mark->residue_idx;
                        ret.emplace_back(cr);
                }
            }
        }
        return (ret.size()>2?ret:empty);
    }

    bool Fragment::sameChain(gemmi::Chain* chain, int chainIdx) {
        return (chain->name == tmdetVO.chains[chainIdx].id );
    }

    
}
