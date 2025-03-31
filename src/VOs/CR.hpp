// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

/**
 * @brief namespace for value objects
 */
namespace Tmdet::VOs {

#define pcr(p,c,r) p.chains[c].residues[r]

    /**
     * @brief simple struct to store chain and residue index
     */
    struct CR {

        /**
         * @brief chain index
         */
        int chainIdx;

        /**
         * @brief residue index
         */
        int residueIdx;

    };
}
