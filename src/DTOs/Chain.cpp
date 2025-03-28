// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#include <string>
#include <format>
#include <vector>
#include <gemmi/align.hpp>
#include <gemmi/metadata.hpp>
#include <gemmi/model.hpp>
#include <gemmi/seqtools.hpp>
#include <Config.hpp>
#include <DTOs/Chain.hpp>
#include <DTOs/Dssp.hpp>
#include <DTOs/Residue.hpp>
#include <Helpers/String.hpp>
#include <System/Logger.hpp>
#include <System/Environment.hpp>
#include <Types/Residue.hpp>
#include <VOs/Protein.hpp>
#include <VOs/Chain.hpp>
#include <VOs/Residue.hpp>

namespace Tmdet::DTOs {

    /**
     * @brief It calculates one-letter sequence from chain residues.
     */
    std::string getOneLetterCodeForEmptySequence(gemmi::Chain& chain) {

        std::string oneLetterSequence;
        for (const auto& residue : chain.residues) {
            auto residueType = Tmdet::Types::ResidueType::getResidue(residue.name);
            oneLetterSequence += residueType.a1;
        }
        return oneLetterSequence;
    }

    Tmdet::VOs::Chain Chain::get(gemmi::Structure& protein, gemmi::Chain& chain, int chainIdx) {
        Tmdet::VOs::Chain chainVO;
        chainVO.id = chain.name;
        chainVO.selected = false;
        auto poly = chain.get_polymer();
        if (poly) {
            chainVO.entityId = chain.residues[0].entity_id;
            chainVO.labelId = poly.subchain_id();
            chainVO.entityIdx = getEntityIdx(protein.entities,chainVO.labelId);
            if (chainVO.entityIdx != -1 
                && protein.entities[chainVO.entityIdx].polymer_type == gemmi::PolymerType::PeptideL) {
                chainVO.selected = true;
                auto sequence = protein.entities[chainVO.entityIdx].full_sequence;
                if (sequence.size() > 0) {
                    chainVO.seq = gemmi::one_letter_code(sequence);
                } else {
                    chainVO.seq = getOneLetterCodeForEmptySequence(chain);
                }
                chainVO.idx = chainIdx;
                int residueIdx = 0;
                for(auto& residue: chain.residues) {
                    chainVO.residues.emplace_back(Tmdet::DTOs::Residue::get(residue,chainIdx,residueIdx++));
                }
                chainVO.length = residueIdx;
                if (chainVO.length > 0) {
                    detectResolution(chainVO);
                }
                else {
                    chainVO.selected = false;
                }
            }
        }
        return chainVO;
    }

    std::string Chain::toString(const Tmdet::VOs::Chain& chain) {
        std::string residues = "";
        for(const auto& residue: chain.residues) {
            residues += Tmdet::DTOs::Residue::toString(residue);
        }
        return std::format(R"(
CHAIN idx:{} authId:{} labelId:{} entityId:{} entityIdx:{} length:{} selected:{}
    Sequence: 
{}
    Secondary structure:
{}
{})",
            chain.idx, chain.id, chain.labelId, chain.entityId,
            chain.entityIdx, chain.length,(chain.selected?"True":"False"),
            Tmdet::Helpers::String::formatSequence(chain.seq, 50, 10,"\t\t"),
            Tmdet::Helpers::String::formatSequence(Tmdet::DTOs::Dssp::getSecondaryStructure(chain), 50, 10,"\t\t"),
            residues
        );
    }

    int Chain::getEntityIdx(const std::vector<gemmi::Entity> entities, const std::string entityId) {

        int idx=0;
        for (const auto& entity : entities) {
            if (entity.entity_type == gemmi::EntityType::Polymer) {
                for(const auto& subchain: entity.subchains) {
                    if (subchain == entityId) {
                        return idx;
                    }
                }
            }
            idx++;
        }
        WARN_LOG("Could not find entity: >>{}<<",entityId);
        return -1;
    }

    void Chain::detectResolution(Tmdet::VOs::Chain& chainVO) {
        int nr = 0;
        int nb = 0;
        for (const auto& residue : chainVO.residues) {
            nr += (residue.hasAllSideChainAtoms()?1:0);
            nb += (residue.hasOnlyBackBoneAtoms()?1:0);
        }
        if (nb > nr) {
            chainVO.type = Tmdet::Types::ChainType::LOW_RES;
        }
        if (nr+nb < std::stoi(environment.get("TMDET_MIN_NUMBER_OF_RESIDUES_IN_CHAIN",DEFAULT_TMDET_MIN_NUMBER_OF_RESIDUES_IN_CHAIN))) {
            chainVO.selected = false;
        }
    }

}