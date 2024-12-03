#include <ValueObjects/Chain.hpp>
#include <ValueObjects/Residue.hpp>
#include <ValueObjects/TMatrix.hpp>

namespace Tmdet::ValueObjects {

    void Chain::transform(Tmdet::ValueObjects::TMatrix& tmatrix) {
        eachResidue(
            [&](Tmdet::ValueObjects::Residue& residue) -> void {
                residue.transform(tmatrix);
            }
        );
    }

}
