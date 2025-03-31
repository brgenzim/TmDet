// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#include <vector>
#include <any>
#include <array>
#include <iostream>
#include <cmath>
#include <limits>
#include <gemmi/model.hpp>
#include <Config.hpp>
#include <System/Logger.hpp>
#include <VOs/Protein.hpp>
#include <VOs/Residue.hpp>
#include <Utils/Fragment.hpp>

namespace Tmdet::Utils {



    /**
     * @brief main entry point for fragment generation
     *        the resulted fragment ids are inserted
     *        to the protein redidues temp map.
     * @return void
     */
    int Fragment::run() {
        numFragments=0;
        nr = protein.numberOfSelectedResidues();
        for(int i=0; i<nr; i++) {
            contactMap.push_back(std::vector<bool>(nr,false));
        }
        setContactMap();
        createFragments();
        freeTempValues();
        return numFragments;
    }

    void Fragment::setContactMap() {
        nr = protein.numberOfSelectedResidues();
        int cm_index=0;

        protein.eachSelectedResidue(
            [&](Tmdet::VOs::Residue& residue) -> void {
                residue.temp.try_emplace("fragment",std::any(-1));
                residue.temp.try_emplace("cm_index",std::any(cm_index));
                residue.temp.try_emplace("neighbors",any_cast<std::vector<_cr>>(getNeighbors(residue)));
                cm_index++;
            }
        );

        protein.eachSelectedResidue(
            [&](Tmdet::VOs::Residue& residue) -> void {
                for(auto& cr: any_cast<std::vector<_cr>>(residue.temp["neighbors"])) {
                    if (protein.chains[cr.chainIdx].selected
                        && protein.chains[cr.chainIdx].residues[cr.residueIdx].selected) {
                        auto& neighbor = protein.chains[cr.chainIdx].residues[cr.residueIdx];
                        contactMap[any_cast<int>(residue.temp["cm_index"])][any_cast<int>(neighbor.temp["cm_index"])] = true;
                    }
                }
            }
        );
    }

    std::vector<_cr> Fragment::getNeighbors(const Tmdet::VOs::Residue& residue) {
        std::vector<_cr> ret;
        std::vector<_cr> empty;
        const gemmi::Atom* ca_atom = residue.gemmi.get_ca();
        bool check=false;
        if (ca_atom) {
            for (auto mark: protein.neighbors.find_neighbors(*ca_atom, 3, 9)) {
                if (protein.chains[mark->chain_idx].selected
                    && protein.chains[mark->chain_idx].residues[mark->residue_idx].selected
                    && protein.chains[mark->chain_idx].residues[mark->residue_idx].atoms[mark->atom_idx].gemmi.name == "CA") {
                    ret.emplace_back(mark->chain_idx,mark->residue_idx);
                    check |= (std::abs(residue.idx - protein.chains[mark->chain_idx].residues[mark->residue_idx].idx) > 3);
                }
            }
        }
        return (check?ret:empty);
    }

    void Fragment::createFragments() {
        protein.eachSelectedResidue(
            [&](Tmdet::VOs::Residue& residue) -> void {
                if (any_cast<int>(residue.temp["fragment"]) == -1 
                    && any_cast<std::vector<_cr>>(residue.temp["neighbors"]).size() > 0) {
                    setFragment(residue);
                    numFragments++;
                }
            }
        );
    }

    void Fragment::setFragment(Tmdet::VOs::Residue& residue) {
        residue.temp["fragment"] = std::any(numFragments);
        for(auto& cr: any_cast<std::vector<_cr>>(residue.temp["neighbors"])) {
            auto& neighbor = protein.chains[cr.chainIdx].residues[cr.residueIdx];
            if (enableMove(residue,neighbor)) {
                setFragment(neighbor);
            }
        }
    }

    bool Fragment::enableMove(Tmdet::VOs::Residue& from, Tmdet::VOs::Residue& to) {
        if (!to.selected || !to.temp.contains("fragment") 
            || !to.temp.contains("cm_index") || !from.temp.contains("cm_index")) {
            return false;
        }
        if (any_cast<int>(to.temp["fragment"]) != -1) {
            return false;
        }
        if (from.chainIdx == to.chainIdx
            && std::abs(from.authId - to.authId) <5) {
                return true;
        }
        int numContacts = checkRegionForContact(from,to,-4,-1,-4,-1)
            + checkRegionForContact(from,to,-4,-1,1,4)
            + checkRegionForContact(from,to,1,4,-4,-1)
            + checkRegionForContact(from,to,1,4,1,4);
        return (numContacts>0);
    }

    int Fragment::checkRegionForContact(Tmdet::VOs::Residue& from, Tmdet::VOs::Residue& to, int fb, int fe, int tb, int te) {
        for (int i=fb; i<=fe; i++) {
            int k = any_cast<int>(from.temp["cm_index"]) + i;
            if (k>=0 && k<nr) {
                for (int j=tb; j<=te; j++) {
                    int l = any_cast<int>(to.temp["cm_index"]) + j;
                    if (l>=0 && l<nr) {
                        if (contactMap[k][l]) {
                            return 1;
                        }
                    }
                }
            }
        }
        return 0;
    }

    void Fragment::freeTempValues() {
        protein.eachSelectedChain(
            [&](Tmdet::VOs::Chain& chain) -> void {
                for (auto& residue: chain.residues) {
                    residue.temp.erase("neighbors");
                }
            }
        );
    }

}
