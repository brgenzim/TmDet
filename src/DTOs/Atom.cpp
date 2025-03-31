// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#include <string>
#include <format>
#include <DTOs/Atom.hpp>
#include <VOs/Atom.hpp>

namespace Tmdet::DTOs {

    std::string Atom::toString(const Tmdet::VOs::Atom& atom) {
        return std::format(R"(
        ATOM idx:{: >6d} name:{} x:{:8.3f} y:{:8.3f} z:{:8.3f} surf:{:8.3f} out:{:8.3f})",
            atom.idx, atom.gemmi.name, atom.gemmi.pos.x, atom.gemmi.pos.y, atom.gemmi.pos.z,
            atom.surface, atom.outSurface);
    }
}
