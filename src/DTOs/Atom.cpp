#include <string>
#include <format>
#include <DTOs/Atom.hpp>
#include <ValueObjects/Atom.hpp>

namespace Tmdet::DTOs {

    std::string Atom::toString(const Tmdet::ValueObjects::Atom& atom) {
        return std::format(R"(
        ATOM idx:{} name:{} x:{} y:{} z:{} surf:{} out:{})",
            atom.idx, atom.gemmi.name, atom.gemmi.pos.x, atom.gemmi.pos.y, atom.gemmi.pos.z,
            atom.surface, (atom.temp.contains("outside")?std::to_string(any_cast<double>(atom.temp.at("outside"))):""));
    }
}
