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
