#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
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

        // Fill residue gaps in chains where there is usable entity sequence data
        for(auto& chain: tmdetVO.gemmi.models[0].chains) {
            auto sequence = getChainSequence(tmdetVO, chain);
            if (sequence.empty() || chain.residues.empty()) {
                // no supporting information to do gap fix
                // or there is no residues for the iteration
                continue;
            }

            int labelSeq = 1;
            std::vector<gemmi::Residue> newChainResidues;
            for( auto& residue: chain.residues) {
                // insert residues, when there is difference between entity_poly and atom sites seq
                auto newResidues = simpleResidueGapFill(chain, residue, labelSeq, sequence);
                if (newResidues.size() > 0) {
                    labelSeq += newResidues.size();
                    newChainResidues.insert(newChainResidues.end(), newResidues.begin(), newResidues.end());
                }
                labelSeq++;
            }
            chain.append_residues(newChainResidues);
            std::sort(chain.residues.begin(), chain.residues.end(), compareResidues);

            // If chain sequence longer than atom sequence
            size_t residueCount = chain.residues.size();
            if (sequence.size() > residueCount) {
                newChainResidues.resize(0);
                int labelSeqNum = residueCount;
                const int seqLength = sequence.size();
                int nextSeqNum = chain.residues[residueCount - 1].seqid.num.value + 1;
                while (labelSeqNum < seqLength) {
                    string residueName = sequence[labelSeqNum];
                    auto newResidue = createResidue(nextSeqNum, labelSeqNum, residueName, chain.name);
                    newChainResidues.emplace_back(*newResidue);
                    labelSeqNum++;
                    nextSeqNum++;
                }
                chain.append_residues(newChainResidues);
            }
        }

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

        tmdetVO.neighbors = NeighborSearch(tmdetVO.gemmi.models[0], tmdetVO.gemmi.cell, 9);
        tmdetVO.neighbors.populate();

    }

    bool TmdetStruct::compareResidues(const gemmi::Residue& res1, const gemmi::Residue& res2) {
        return res1.label_seq.value < res2.label_seq.value;
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
                    cout << atom.surface << " " << Tmdet::Types::Residues.at(residue.gemmi.name).atoms.at(atom.gemmi.name).atom.vdw;
                    if (atom.temp.find("outside") != atom.temp.end()) {
                        cout << " out: " << any_cast<double>(atom.temp.at("outside"));
                    }
                    cout << endl;
                }
            }
        }
    }

    std::vector<string> TmdetStruct::getChainSequence(
        const Tmdet::ValueObjects::TmdetStruct& tmdetVO, const gemmi::Chain& chainVO) {

        std::vector<string> sequence;
        for (const auto& entity : tmdetVO.gemmi.entities) {
            auto begin = entity.subchains.begin();
            auto end = entity.subchains.end();
            if (entity.entity_type == gemmi::EntityType::Polymer
                && std::find(begin, end, chainVO.name) != end) {

                sequence = entity.full_sequence;
                break;
            }
        }
        return sequence;
    }

    std::vector<gemmi::Residue> TmdetStruct::simpleResidueGapFill(gemmi::Chain& chain, gemmi::Residue& residue,
        int labelSeq, const std::vector<string> sequence) {

        const int expectedIndex = residue.label_seq.value;
        int gapLength = expectedIndex - labelSeq;
        std::vector<gemmi::Residue> newResidues;
        while (gapLength > 0) {
            int newSeqId = residue.seqid.num.value - gapLength;
            auto residueForInsert = createResidue(newSeqId, labelSeq, sequence[labelSeq - 1], chain.name);
            newResidues.emplace_back(*residueForInsert);
            labelSeq++;
            gapLength--;
        }
        return newResidues;
    }

    gemmi::Residue* TmdetStruct::createResidue(int seqNum, int labelSeqNum, string name, string chainName) {
        auto residue = new gemmi::Residue();
        residue->seqid.num.value = seqNum;
        residue->label_seq.value = labelSeqNum;
        residue->name = name;
        residue->subchain = chainName;
        return residue;
    }
}

