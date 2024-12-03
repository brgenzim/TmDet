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
        residueVO.authId = (int)residue.seqid.num;
        residueVO.labelId = (int)residue.label_seq;
        residueVO.authIcode = residue.seqid.icode;
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
        /*std::string temp = "";
        for(const auto &[key,value]: residue.temp) {
            temp += std::format(R"(
        TEMP key: {} value: {}
            )",key,any_cast<double>(value));
        }*/
        return std::format(R"(
    RESIDUE idx:{: >6d} authId:{: >6d} labelId:{: >6d} a3:{} a1:{} ss:{} surface:{:8.3f} outSurface:{:8.3f}{})", 
            residue.idx,residue.authId, residue.labelId, residue.type.name, 
            residue.type.a1, residue.ss.code, residue.surface, residue.outSurface, atoms);
    }
}
