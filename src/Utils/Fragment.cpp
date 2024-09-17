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

    std::vector<_cr> Fragment::getCAlphaNetwork() {
        int node_idx = 0;
        std::vector<_cr> crs;
        for (auto& chainVO: tmdetVO.chains) {
            for (auto& residueVO: chainVO.residues) {
                auto neighbors = getNeighbors(residueVO);
                
                residueVO.temp.insert({"neighbors",std::any_cast<std::vector<_cr>>(neighbors)});
                if (neighbors.size() > 0) {
                    residueVO.temp.insert({"node_index",std::any_cast<int>(node_idx++)});
                    _cr cr;
                    cr.chain_idx = residueVO.chainIdx;
                    cr.residue_idx = residueVO.idx;
                    crs.emplace_back(cr);
                }
            }
        }
        return crs;
    }

    std::vector<std::vector<int>> Fragment::createFragments(int size) {
        Graph G(size);
        for (auto& chainVO: tmdetVO.chains) {
            for (auto& residueVO: chainVO.residues) {
                auto neighbors = std::any_cast<std::vector<_cr>>(residueVO.temp.at("neighbors"));
                for (const auto& cr: neighbors) {
                    int u = any_cast<int>(residueVO.temp.at("node_index"));
                    if (tmdetVO.chains[cr.chain_idx].residues[cr.residue_idx].temp.contains("node_index")) {
                        int v = NODE_ID(cr.chain_idx,cr.residue_idx);
                        G.addEdge(u,v);
                    }
                }
            }
        }
        std::vector<int> res = G.optim();
        std::vector<std::vector<int>> ret(G.getNumClusters());
        for( int i=0; i<size; i++) {
            ret[res[i]].push_back(i);
        }
        return ret;
    }

    void Fragment::writeBackFragmentInfoToStructure(std::vector<std::vector<int>> clusters, std::vector<_cr> crs) {
        int cl_idx = 0;
        for (auto& cluster: clusters) {
            std::cout << "select cl" << cl_idx << ", (";
            bool first = true;
            for(auto& nodeIdx: cluster) {
                tmdetVO.chains[crs[nodeIdx].chain_idx].residues[crs[nodeIdx].residue_idx].temp.insert({"fragment",std::any_cast<int>(cl_idx)});
                if (!first) {
                    std::cout << " | ";
                }
                std::cout << "(c. " << tmdetVO.chains[crs[nodeIdx].chain_idx].id << " & ";
                std::cout << "i. " << tmdetVO.chains[crs[nodeIdx].chain_idx].residues[crs[nodeIdx].residue_idx].resn() << ")";
                first = false;
            }
            std::cout << ")" << std::endl;
            cl_idx++;
        }
    }

    void Fragment::freeTempValues() {
        for (auto& chainVO: tmdetVO.chains) {
            for (auto& residueVO: chainVO.residues) {
                residueVO.temp.erase("neighbors");
                if (residueVO.temp.contains("node_index")) {
                    residueVO.temp.erase("node_index");
                }
            }
        }
    }

}
