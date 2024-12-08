#include <vector>
#include <any>
#include <array>
#include <iostream>
#include <cmath>
#include <limits>
#include <gemmi/model.hpp>
#include <VOs/Protein.hpp>
#include <VOs/Residue.hpp>
#include <Utils/Fragment.hpp>
#include <Utils/Graph.hpp>

namespace StructVO = Tmdet::VOs;

namespace Tmdet::Utils {

#define NODE_ID(c,r) any_cast<int>(proteinVO.chains[c].residues[r].temp.at("node_index"))

    /**
     * @brief main entry point for fragment generation
     *        the resulted fragment ids are inserted
     *        to the proteinVO redidues temp map.
     * @return void
     */
    void Fragment::run() {
        auto crs = getCAlphaNetwork();
        auto clusters = createFragments(crs.size());
        writeBackFragmentInfoToStructure(clusters, crs);
        freeTempValues();
    }

    std::vector<_cr> Fragment::getNeighbors(const StructVO::Residue& residueVO) {
        std::vector<_cr> ret;
        std::vector<_cr> empty;
        const gemmi::Atom* ca_atom = residueVO.gemmi.get_ca();
        if (ca_atom) {
            for (auto mark: proteinVO.neighbors.find_neighbors(*ca_atom, 3, 9)) {
                if (proteinVO.chains[mark->chain_idx].residues[mark->residue_idx].atoms[mark->atom_idx].gemmi.name == "CA") {
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
        for (auto& chainVO: proteinVO.chains) {
            for (auto& residueVO: chainVO.residues) {
                auto neighbors = getNeighbors(residueVO);
                
                residueVO.temp.try_emplace("neighbors",std::any_cast<std::vector<_cr>>(neighbors));
                if ( !neighbors.empty() ) {
                    residueVO.temp.try_emplace("node_index",std::any_cast<int>(node_idx++));
                    _cr cr;
                    cr.chain_idx = residueVO.chainIdx;
                    cr.residue_idx = residueVO.idx;
                    crs.emplace_back(cr);
                }
            }
        }
        return crs;
    }

    std::vector<std::vector<int>> Fragment::createFragments(const unsigned long size) {
        Graph G(size);
        for (auto& chainVO: proteinVO.chains) {
            for (auto& residueVO: chainVO.residues) {
                auto neighbors = std::any_cast<std::vector<_cr>>(residueVO.temp.at("neighbors"));
                for (const auto& cr: neighbors) {
                    int u = any_cast<int>(residueVO.temp.at("node_index"));
                    if (proteinVO.chains[cr.chain_idx].residues[cr.residue_idx].temp.contains("node_index")) {
                        int v = NODE_ID(cr.chain_idx,cr.residue_idx);
                        G.addEdge(u,v);
                    }
                }
            }
        }
        std::vector<unsigned int> res = G.optim();
        std::vector<std::vector<int>> ret(G.getNumClusters());
        for(unsigned long i=0; i<size; i++) {
            ret[res[i]].push_back(i);
        }
        return ret;
    }

    void Fragment::writeBackFragmentInfoToStructure(std::vector<std::vector<int>> clusters, std::vector<_cr> crs) {
        int cl_idx = 0;
        for (auto& cluster: clusters) {
            //std::cout << "select cl" << cl_idx << ", (";
            //bool first = true;
            for(auto& nodeIdx: cluster) {
                proteinVO.chains[crs[nodeIdx].chain_idx].residues[crs[nodeIdx].residue_idx].temp.insert({"fragment",std::any_cast<int>(cl_idx)});
               /* if (!first) {
                    std::cout << " | ";
                }
                std::cout << "(c. " << proteinVO.chains[crs[nodeIdx].chain_idx].id << " & ";
                std::cout << "i. " << proteinVO.chains[crs[nodeIdx].chain_idx].residues[crs[nodeIdx].residue_idx].resn() << ")";
                first = false;*/
            }
            //std::cout << ")" << std::endl;
            cl_idx++;
        }
    }

    void Fragment::freeTempValues() {
        for (auto& chainVO: proteinVO.chains) {
            for (auto& residueVO: chainVO.residues) {
                residueVO.temp.erase("neighbors");
                if (residueVO.temp.contains("node_index")) {
                    residueVO.temp.erase("node_index");
                }
            }
        }
    }

}
