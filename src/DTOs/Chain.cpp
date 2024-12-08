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
#include <VOs/Protein.hpp>
#include <VOs/Chain.hpp>
#include <VOs/Residue.hpp>

namespace Tmdet::DTOs {
    
    Tmdet::VOs::Chain Chain::get(gemmi::Structure& protein, gemmi::Chain& chain, int chainIdx) {
        Tmdet::VOs::Chain chainVO;
        chainVO.id = chain.name;
        chainVO.selected = false;
        auto poly = chain.get_polymer();
        if (poly) {
            chainVO.entityId = poly.subchain_id(); //chain.residues[0].entity_id;
            chainVO.entityIdx = getEntityIdx(protein.entities,chainVO.entityId);
            if (chainVO.entityIdx != -1) {
                chainVO.selected = true;
                auto sequence = protein.entities[chainVO.entityIdx].full_sequence;
                chainVO.seq = gemmi::one_letter_code(sequence);
                chainVO.idx = chainIdx;
        
        
                chainVO.labelId = poly.subchain_id();
                int residueIdx = 0;
                for(auto& residue: chain.residues) {
                    chainVO.residues.emplace_back(Tmdet::DTOs::Residue::get(residue,chainIdx,residueIdx++));
                }
                chainVO.length = residueIdx;
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
            DEBUG_LOG("Entity: name:{}<< type:{}<<",entity.name,(entity.entity_type==gemmi::EntityType::Polymer?"yes":"no"));
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

}