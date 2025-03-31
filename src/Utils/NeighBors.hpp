// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <vector>
#include <VOs/CR.hpp>
#include <VOs/Protein.hpp>
#include <VOs/Residue.hpp>

/**
 * @brief namespace for utils
 */
namespace Tmdet::Utils {

    /**
     * @brief storing Ca neighbors
     */
    struct NeighBors {
        static void store(Tmdet::VOs::Protein& protein) {
            protein.eachSelectedResidue(
                [&](Tmdet::VOs::Residue& residue) -> void {
                    std::vector<Tmdet::VOs::CR> neighbors;
                    auto ca =residue.getCa();
                    if (ca != nullptr) {
                        for(auto m : protein.neighbors.find_neighbors(*ca, 2, 6.5)) {
                            if (protein.chains[m->chain_idx].selected
                                && protein.chains[m->chain_idx].residues[m->residue_idx].selected
                                && protein.chains[m->chain_idx].residues[m->residue_idx].atoms[m->atom_idx].gemmi.name == "CA"
                                && (residue.chainIdx != m->chain_idx || 
                                    (residue.chainIdx == m->chain_idx 
                                        && std::abs(residue.labelId - protein.chains[m->chain_idx].residues[m->residue_idx].labelId) > 2))) {
                                neighbors.emplace_back(m->chain_idx,m->residue_idx);
                            }
                        }
                    }
                    residue.temp.try_emplace("ca_neighbors",neighbors);
                }
            );
        }

        static std::vector<Tmdet::VOs::CR> get(Tmdet::VOs::Residue& residue) {
            std::vector<Tmdet::VOs::CR> empty;
            return (residue.temp.contains("ca_neighbors")?any_cast<std::vector<Tmdet::VOs::CR>>(residue.temp["ca_neighbors"]):empty);
        }
    };
}
