#include <string>
#include <vector>
#include <iostream>
#include <gemmi/model.hpp>
#include <gemmi/modify.hpp>
#include <gemmi/polyheur.hpp>
#include <DTOs/TmdetStruct.hpp>
#include <Utils/Xml.hpp>
#include <Types/Protein.hpp>
#include <Types/Residue.hpp>
#include <ValueObjects/TmdetStruct.hpp>

using namespace std;

namespace Tmdet::DTOS {

    void TmdetStruct::readXml(Tmdet::ValueObjects::TmdetStruct& tmdetVO, string path) {
        Tmdet::Utils::Xml xml;
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

    void TmdetStruct::writeXml(Tmdet::ValueObjects::TmdetStruct& tmdetVO, string path) {
        Tmdet::Utils::Xml xml;
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

    void TmdetStruct::parse(Tmdet::ValueObjects::TmdetStruct& tmdetVO) {
        remove_hydrogens(tmdetVO.gemmi.models[0]);
        remove_ligands_and_waters(tmdetVO.gemmi.models[0]);
        remove_alternative_conformations(tmdetVO.gemmi.models[0]);
        tmdetVO.neighbors = NeighborSearch(tmdetVO.gemmi.models[0], tmdetVO.gemmi.cell, 9);
        tmdetVO.neighbors.populate();
        int chainIdx = 0;
        for(auto& chain: tmdetVO.gemmi.models[0].chains) {
            Tmdet::ValueObjects::Chain chainVO = Tmdet::ValueObjects::Chain(chain);
            chainVO.idx = chainIdx++;
            int residueIdx = 0;
            for( auto& residue: chain.residues) {
                Tmdet::ValueObjects::Residue residueVO = Tmdet::ValueObjects::Residue(residue);
                residueVO.chainIdx = chainVO.idx;
                residueVO.idx = residueIdx++;
                residueVO.surface = 0.0;
                int atomIdx = 0;
                for( auto& atom: residue.atoms) {
                    Tmdet::ValueObjects::Atom atomVO = Tmdet::ValueObjects::Atom(atom);
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
    }

    void TmdetStruct::out(Tmdet::ValueObjects::TmdetStruct& tmdetVO) {
        for(auto& chain: tmdetVO.chains) {
            cout << "CHAIN " << chain.gemmi.name << " " << chain.length << endl;
            for( auto& residue: chain.residues) {
                cout << "\tRESIDUE " << residue.idx << ":" << residue.resn() << "(" << residue.gemmi.name << ") ";
                cout << residue.surface << " " << residue.ss.code << endl;
                for( auto& atom: residue.atoms) {
                    cout << "\t\tATOM " << atom.idx << ": " << atom.gemmi.name << " ";
                    cout << atom.gemmi.pos.x << " " << atom.gemmi.pos.y << " " << atom.gemmi.pos.z << " ";
                    cout << atom.surface << " " << atom.gemmi.element.vdw_r() << endl;
                }
            }
        }
    }
}
