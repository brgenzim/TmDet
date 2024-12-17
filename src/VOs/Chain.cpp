// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#include <VOs/Chain.hpp>
#include <VOs/Residue.hpp>
#include <VOs/TMatrix.hpp>

namespace Tmdet::VOs {

    void Chain::transform(Tmdet::VOs::TMatrix& tmatrix) {
        eachResidue(
            [&](Tmdet::VOs::Residue& residue) -> void {
                residue.transform(tmatrix);
            }
        );
    }

}
