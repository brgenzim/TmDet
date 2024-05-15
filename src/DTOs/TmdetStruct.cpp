#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <gemmi/model.hpp>
#include <gemmi/modify.hpp>
#include <gemmi/polyheur.hpp>
#include <gemmi/align.hpp>
#include <gemmi/seqalign.hpp>
#include <gemmi/seqtools.hpp>
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
        alignResidues(tmdetVO);

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
        const Tmdet::ValueObjects::TmdetStruct& tmdetVO, const gemmi::Chain& chain) {

        std::vector<string> sequence;
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
        // if (newResidues.empty()) {
        //     cout << residue.label_seq.value << ":" << residue.name;
        // } else {
        //     for (auto& res : newResidues) {
        //         cout << res.label_seq.value << ":" << res.name << ", ";
        //     }
        // }
        // cout << endl;
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

    #include <iostream>
    #include <vector>
    #include <utility>
    #include <cctype>
    #include <cstdlib>

    void TmdetStruct::alignResidues(const Tmdet::ValueObjects::TmdetStruct& tmdetVO) {

        for(auto& chain: tmdetVO.gemmi.models[0].chains) {

            auto sequence = getChainSequence(tmdetVO, chain);

            if (sequence.empty() || chain.residues.empty()) {
                // no supporting information to do gap fix
                // or there is no residues for the iteration
                continue;
            }

            std::vector<string> chainResidues;
            for( auto& residue: chain.residues) {
                chainResidues.emplace_back(residue.name);
            }


            const gemmi::AlignmentScoring* scoring = gemmi::AlignmentScoring::simple();
            // cout << "one letterSEQ: " << gemmi::one_letter_code(sequence) << endl;

            const std::vector<int> gapOpening;
            auto alignmentResult2 = gemmi::align_string_sequences(sequence, chainResidues, gapOpening, scoring);

            auto withGaps2 = alignmentResult2.add_gaps(gemmi::one_letter_code(chainResidues) , 2);
            string cigar = alignmentResult2.cigar_str();
            // cout << "CIGAR['" << chain.name << "']: " << cigar
            //     << " identity: " << alignmentResult2.calculate_identity() << endl;
            // cout << "with gaps2:    " << withGaps2 << endl;
            // exit(1);

            std::vector<gemmi::Residue>  newChainResidues;
            int labelSeqNum = 1;
            for (auto& oneLetterSeq : withGaps2) {

                if (oneLetterSeq != '-') {
                    labelSeqNum++;
                    continue;
                }

                gemmi::Residue* newResidue = createResidue(0, labelSeqNum, sequence[labelSeqNum - 1], chain.name);
                newChainResidues.emplace_back(*newResidue);
                labelSeqNum++;
            }
            chain.append_residues(newChainResidues);
            std::sort(chain.residues.begin(), chain.residues.end(), compareResidues);
        }
    }

    std::vector<std::pair<char, int>> TmdetStruct::parseCIGAR(const std::string& cigar) {
        std::vector<std::pair<char, int>> result;

        int length = 0;
        for (char c : cigar) {
            if (std::isdigit(c)) {
                length = length * 10 + (c - '0');
            } else {
                result.emplace_back(c, length);
                length = 0;
            }
        }

        return result;
    }

}

