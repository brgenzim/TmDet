#include <string>
#include <format>
#include <gemmi/model.hpp>
#include <DTOs/Atom.hpp>
#include <DTOs/Residue.hpp>
#include <VOs/Atom.hpp>
#include <VOs/Residue.hpp>

namespace Tmdet::DTOs {

    Tmdet::VOs::Residue Residue::get(gemmi::Residue& residue, int chainIdx, int residueIdx) {
        auto residueVO = Tmdet::VOs::Residue(residue);
        residueVO.chainIdx = chainIdx;
        residueVO.idx = residueIdx;
        residueVO.authId = (int)residue.seqid.num;
        residueVO.labelId = (int)residue.label_seq;
        residueVO.authIcode = residue.seqid.icode;
        int atomIdx = 0;
        for(auto& atom: residue.atoms) {
            auto atomVO = Tmdet::VOs::Atom(atom);
            atomVO.chainIdx = chainIdx;
            atomVO.residueIdx = residueIdx;
            atomVO.idx = atomIdx++;
            residueVO.atoms.emplace_back(atomVO);
        }
        residueVO.setProperties(residue.name);
        return residueVO;
    }

    std::string Residue::toString(const Tmdet::VOs::Residue& residue) {
        std::string atoms = "";
        for(const auto& atom: residue.atoms) {
            atoms += Tmdet::DTOs::Atom::toString(atom);
        }
        std::string temp = "";
        if (residue.temp.contains("fragment")) {
            temp += std::format(R"(TEMP fragmentId: {})",
                any_cast<int>(residue.temp.at("fragment")));
        }
        return std::format(R"(
    RESIDUE idx:{: >6d} authId:{: >6d} labelId:{: >6d} a3:{} a1:{} ss:{} surface:{:8.3f} outSurface:{:8.3f}{} temp:{})", 
            residue.idx,residue.authId, residue.labelId, residue.type.name, 
            residue.type.a1, residue.ss.code, residue.surface, residue.outSurface, temp, atoms);
    }
}
