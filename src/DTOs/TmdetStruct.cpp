#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <gemmi/cifdoc.hpp>
#include <gemmi/model.hpp>
#include <gemmi/modify.hpp>
#include <gemmi/polyheur.hpp>
#include <gemmi/align.hpp>
#include <gemmi/seqalign.hpp>
#include <gemmi/seqtools.hpp>
#include <gemmi/to_cif.hpp>
#include <gemmi/to_mmcif.hpp>
#include <DTOs/TmdetStruct.hpp>
#include <Utils/Alignment.hpp>
#include <DTOs/Xml.hpp>
#include <Types/Protein.hpp>
#include <Types/Residue.hpp>
#include <Types/SecStruct.hpp>
#include <ValueObjects/TmdetStruct.hpp>
#include <ValueObjects/Chain.hpp>
#include <ValueObjects/Residue.hpp>

namespace Tmdet::DTOs {

    void TmdetStruct::readXml(Tmdet::ValueObjects::TmdetStruct& tmdetVO, const std::string& path) {
        Xml xml;
        xml.read(path);
        tmdetVO.tmp = xml.getTmp();
        tmdetVO.code = xml.getCode();
        tmdetVO.date = xml.getCreateDate();
        tmdetVO.modifications = xml.getModifications();
        tmdetVO.qValue = xml.getQvalue();
        tmdetVO.type = Tmdet::Types::Proteins.at(xml.getTmtype());
        tmdetVO.spres = xml.getSpres();
        tmdetVO.pdbkwres = xml.getPdbkwres();
        tmdetVO.bioMatrix = xml.getBioMatrix();
        tmdetVO.membranes = xml.getMembranes();
        xml.getChains(tmdetVO.chains);
    }

    void TmdetStruct::writeXml(Tmdet::ValueObjects::TmdetStruct& tmdetVO, const std::string& path) {
        Xml xml;
        xml.create();
        xml.setTmp(tmdetVO.tmp);
        xml.setCode(tmdetVO.code);
        xml.setCreateDate(tmdetVO.date);
        xml.setModifications(tmdetVO.modifications);
        xml.setQvalue(tmdetVO.qValue);
        xml.setTmtype(tmdetVO.type.name);
        xml.setSpres(tmdetVO.spres);
        xml.setPdbkwres(tmdetVO.pdbkwres);
        xml.setBioMatrix(tmdetVO.bioMatrix);
        xml.setMembranes(tmdetVO.membranes);
        xml.setChains(tmdetVO.chains);
        xml.write(path);
    }

    void TmdetStruct::writeCif(const Tmdet::ValueObjects::TmdetStruct& tmdetVO, const std::string& path) {
        std::ofstream outCif(path);
        gemmi::cif::WriteOptions options(gemmi::cif::Style::Pdbx);
        gemmi::cif::Document document = make_mmcif_document(tmdetVO.gemmi);
        gemmi::cif::write_cif_to_stream(outCif, document, options);
    }

    void TmdetStruct::parse(Tmdet::ValueObjects::TmdetStruct& tmdetVO) {
        remove_hydrogens(tmdetVO.gemmi.models[0]);
        remove_ligands_and_waters(tmdetVO.gemmi.models[0]);
        remove_alternative_conformations(tmdetVO.gemmi.models[0]);
        // keep only the first model
        tmdetVO.gemmi.models.resize(1);

        // Fill residue gaps in chains where there is usable entity sequence data
        Tmdet::Utils::Alignment::alignResidues(tmdetVO);
        int chainIdx = 0;
        for(auto& chain: tmdetVO.gemmi.models[0].chains) {
            auto chainVO = Tmdet::ValueObjects::Chain(chain);
            chainVO.idx = chainIdx++;
            int residueIdx = 0;
            for( auto& residue: chain.residues) {
                auto residueVO = Tmdet::ValueObjects::Residue(residue,chainVO);
                residueVO.chainIdx = chainVO.idx;
                residueVO.idx = residueIdx++;
                residueVO.surface = 0.0;
                int atomIdx = 0;
                for( auto& atom: residue.atoms) {
                    auto atomVO = Tmdet::ValueObjects::Atom(atom,residueVO,chainVO);
                    atomVO.chainIdx = chainVO.idx;
                    atomVO.residueIdx = residueVO.idx;
                    atomVO.idx = atomIdx++;
                    atomVO.surface = 0.0;
                    residueVO.atoms.emplace_back(atomVO);
                }
                chainVO.residues.emplace_back(residueVO);
            }
            chainVO.length = residueIdx;
            tmdetVO.chains.emplace_back(chainVO);
        }

        tmdetVO.neighbors = gemmi::NeighborSearch(tmdetVO.gemmi.models[0], tmdetVO.gemmi.cell, 9);
        tmdetVO.neighbors.populate();
    }

    void TmdetStruct::print(std::ostream& os, const Tmdet::ValueObjects::TmdetStruct& tmdetVO) {
        for(const auto& chainVO: tmdetVO.chains) {
            printChain(os, chainVO);
        }
    }

    void TmdetStruct::printChain(std::ostream& os, const Tmdet::ValueObjects::Chain& chainVO) {
        for(const auto& residueVO: chainVO.residues) {
            os << "CHAIN " << chainVO.gemmi.name << " " << chainVO.length << std::endl;
            printResidue(os, residueVO);
        }
    }

    void TmdetStruct::printResidue(std::ostream& os, const Tmdet::ValueObjects::Residue& residueVO) {
        os << "\tRESIDUE " << residueVO.idx << ":" << residueVO.resn() << "(" << residueVO.gemmi.name << ") ";
        os << residueVO.surface << " " << residueVO.ss.code << std::endl;
        if ( residueVO.temp.contains("fragment") ) {
            os << "\t\tTEMP: fragment: " << any_cast<int>(residueVO.temp.at("fragment")) << std::endl;
        }
        for(const auto& atomVO: residueVO.atoms) {
            printAtom(os, atomVO);
        }
    }

    void TmdetStruct::printAtom(std::ostream& os, const Tmdet::ValueObjects::Atom& atomVO) {
        os << "\t\tATOM " << atomVO.idx << ": " << atomVO.gemmi.name << " ";
        os << atomVO.gemmi.pos.x << " " << atomVO.gemmi.pos.y << " " << atomVO.gemmi.pos.z << " ";
        Tmdet::Types::Residue residueType = Tmdet::Types::ResidueType::getResidue(atomVO.parentResidue.gemmi.name);
        double vdw = 0.0;
        if (residueType.atoms.contains(atomVO.gemmi.name)) {
            vdw = residueType.atoms.at(atomVO.gemmi.name).atom.vdw;
        } else {
            vdw = Types::AtomType::DEFAULT_VDW;
        }
        os << atomVO.surface << " " << vdw;
        if (atomVO.temp.contains("outside")) {
            os << " out: " << any_cast<double>(atomVO.temp.at("outside"));
        }
        os << std::endl;
    }

    std::vector<std::string> TmdetStruct::getChainSequence(
        const Tmdet::ValueObjects::TmdetStruct& tmdetVO, const gemmi::Chain& chain) {

        std::vector<std::string> sequence;
        auto entityId = chain.residues[0].entity_id;
        for (const auto& entity : tmdetVO.gemmi.entities) {
            if (entity.entity_type == gemmi::EntityType::Polymer
                && entity.name == entityId) {

                sequence = entity.full_sequence;
                break;
            }
        }
        return sequence;
    }

}
