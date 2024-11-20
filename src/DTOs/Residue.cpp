#include <string>
#include <format>
#include <gemmi/model.hpp>
#include <DTOs/Atom.hpp>
#include <DTOs/Residue.hpp>
#include <ValueObjects/Atom.hpp>
#include <ValueObjects/Residue.hpp>

namespace Tmdet::DTOs {

    Tmdet::ValueObjects::Residue Residue::get(gemmi::Residue& residue, int chainIdx, int residueIdx) {
        auto residueVO = Tmdet::ValueObjects::Residue(residue);
        residueVO.chainIdx = chainIdx;
        residueVO.idx = residueIdx;
        residueVO.authId = 0; //todo
        residueVO.labelId = 0; //todo
        residueVO.authIcode = ' '; //todo
        int atomIdx = 0;
        for( auto& atom: residue.atoms) {
            auto atomVO = Tmdet::ValueObjects::Atom(atom);
            atomVO.chainIdx = chainIdx;
            atomVO.residueIdx = residueIdx;
            atomVO.idx = atomIdx++;
            residueVO.atoms.emplace_back(atomVO);
        }
        residueVO.setProperties(residue.name);
        return residueVO;
    }

    Tmdet::ValueObjects::Residue Residue::unk(int chainIdx, int residueIdx, std::string name) {
        gemmi::Residue residue;
        auto residueVO = Tmdet::ValueObjects::Residue(residue);
        residueVO.chainIdx = chainIdx;
        residueVO.idx = residueIdx;
        residueVO.authId = -1;
        residueVO.labelId = -1;
        residueVO.authIcode = ' ';
        residueVO.setProperties(name);
        return residueVO;
    }

    std::string Residue::toString(const Tmdet::ValueObjects::Residue& residue) {
        std::string atoms = "";
        for(const auto& atom: residue.atoms) {
            atoms += Tmdet::DTOs::Atom::toString(atom);
        }
        return std::format(R"(
    RESIDUE authId: {} labelId{} a3: {} a1: {} ss: {} surface: {}
{})", 
            residue.authId, residue.labelId, residue.type.name, residue.type.a1,
            residue.ss.code, residue.surface, atoms);
    }
}
