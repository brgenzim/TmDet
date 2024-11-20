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
#include <ValueObjects/Protein.hpp>
#include <ValueObjects/Chain.hpp>
#include <ValueObjects/Residue.hpp>

namespace Tmdet::DTOs {
    
    Tmdet::ValueObjects::Chain Chain::get(gemmi::Structure& protein, gemmi::Chain& chain, std::string entityId, int chainIdx) {
        Tmdet::ValueObjects::Chain chainVO;
        chainVO.id = chain.name;
        chainVO.entityId = entityId;
        chainVO.entityIdx = getEntityIdx(protein.entities,entityId);
        auto sequence = protein.entities[chainVO.entityIdx].full_sequence;
        chainVO.seq = gemmi::one_letter_code(sequence);
        chainVO.idx = chainIdx;
        auto poly = chain.get_polymer();
        auto alignRes = gemmi::align_sequence_to_polymer(
                    sequence,
                    poly,
                    protein.entities[chainVO.entityIdx].polymer_type);
        chainVO.length = alignRes.match_string.length();
        chainVO.labelId = poly.subchain_id();
        int residueIdx = 0;
        for (auto i=0; i<chainVO.length; i++) {
            if (alignRes.match_string[i] == '|') {
                chainVO.gemmiResidueIndeces.push_back(residueIdx);
                chainVO.residues.emplace_back(Tmdet::DTOs::Residue::get(chain.residues[residueIdx++],chainIdx,i));
            }
            else {
                chainVO.residues.emplace_back(Tmdet::DTOs::Residue::unk(chainIdx,i,sequence[i]));
            }
        }
        return chainVO;
    }

    std::string Chain::toString(const Tmdet::ValueObjects::Chain& chain) {
        std::string residues = "";
        for(const auto& residue: chain.residues) {
            residues += Tmdet::DTOs::Residue::toString(residue);
        }
        return std::format(R"(
CHAIN id: {} entityId: {} entityIdx: {} length: {} selected: {}
    Sequence: {}
    Secondary structure: {}
{})",
            chain.id, chain.entityId, chain.entityIdx, chain.length,
            (chain.selected?"True":"False"),
            Tmdet::Helpers::String::formatSequence(chain.seq, 50, 10),
            Tmdet::Helpers::String::formatSequence(Tmdet::DTOs::Dssp::getSecondaryStructure(chain), 50, 10),
            residues
        );
    }

    int Chain::getEntityIdx(const std::vector<gemmi::Entity> entities, const std::string entityId) {

        int idx=0;
        for (const auto& entity : entities) {
            if (entity.entity_type == gemmi::EntityType::Polymer
                && entity.name == entityId) {
                    return idx;
            }
            idx++;
        }
        WARN_LOG("Could not find entity: {}",entityId);
        return -1;
    }

}