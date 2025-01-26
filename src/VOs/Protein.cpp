// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#include <string>
#include <map>
#include <gemmi/cifdoc.hpp>
#include <gemmi/gz.hpp>
#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/to_mmcif.hpp>
#include <gemmi/to_cif.hpp>
#include <gemmi/mmread.hpp>
#include <gemmi/model.hpp>
#include <gemmi/cifdoc.hpp>
#include <gemmi/polyheur.hpp>
#include <gemmi/align.hpp>
#include <Version.hpp>
#include <Config.hpp>
#include <System/FilePaths.hpp>
#include <System/Date.hpp>
#include <System/Logger.hpp>
#include <Types/Protein.hpp>
#include <VOs/Protein.hpp>
#include <Utils/Md5.hpp>

namespace Tmdet::VOs {

    void Protein::getStructure(const std::string& inputPath) {

        logger.debug("Processing protein.getStructure()");
        Tmdet::System::FilePaths::isCif(inputPath)?getCifStructure(inputPath):getEntStructure(inputPath);
        logger.debug(" Processed protein.getStructure()");
    }

    void Protein::getEntStructure(const std::string& inputPath) {
        gemmi = gemmi::read_pdb(gemmi::MaybeGzipped(inputPath));

        gemmi::setup_entities(gemmi);
        gemmi::assign_label_seq_id(gemmi,true);
        document = gemmi::make_mmcif_document(gemmi);
        setupPolymerNames();
    }

    void Protein::getCifStructure(const std::string& inputPath) {
        version = Tmdet::version();
        date = Tmdet::System::Date::get();
        document = gemmi::cif::read(gemmi::MaybeGzipped(inputPath));
        gemmi = gemmi::make_structure(std::move(document));
        setupPolymerNames();
    }

    void Protein::setupPolymerNames() {
        for (auto& entity : gemmi.entities) {
            if (entity.entity_type == gemmi::EntityType::Polymer) {
                polymerNames[entity.name] = "DESCRIPTION: N/A";
            }
        }
        for (gemmi::cif::Block& block : document.blocks) {
            for (auto entity : block.find("_entity.", {"id", "pdbx_description"})) {
                polymerNames[entity[0]] = entity[1];
            }
        }
        for(auto& [key, val]: polymerNames) {
            DEBUG_LOG("polymerNames[{}]: {}",key,val);
        }
    }

    void Protein::notTransmembrane() {
        version = (version==""?Tmdet::version():version);
        modifications.emplace_back(
            Tmdet::System::Date::get(),
            (std::string)"Not transmembrane protein"
        );
        clear();
    }

    void Protein::clear() {
        tmp = false;
        qValue = 0.0;
        type = Tmdet::Types::ProteinType::SOLUBLE;
        membranes.clear();
        eachChain(
            [&](Tmdet::VOs::Chain& chain) -> void {
                chain.regions.clear();
            }
        );
        DEBUG_LOG("Protein data cleared: {}",tmp);
    }

    std::string Protein::hash() const {
        std::string raw = code;
        for(const auto& chain: chains) {
            if (chain.selected) {
                for (const auto& residue: chain.residues) {
                    if (residue.selected) {
                        raw += residue.type.a1;
                    }
                }
            }
        }
        return Tmdet::Utils::Md5::getHash(raw);
    }

    void Protein::unSelectAll() {
        for (auto& c: chains) {
            c.selected = false;
        }
    }

    int Protein::searchChainById(const std::string& id) const {
        for (auto& c: chains) {
            if (c.selected && c.id == id) {
                return c.idx;
            }
        }
        WARN_LOG("Chain {} not found, or not selected",id);
        return -1;
    }

    int Protein::searchChainByLabId(const std::string& id) const {
        for (auto& c: chains) {
            if (c.selected && c.labelId == id) {
                return c.idx;
            }
        }
        WARN_LOG("Chain {} not found, or not selected",id);
        return -1;
    }

    gemmi::Vec3 Protein::centre() {
        gemmi::Vec3 ret(0,0,0);
        int n = 0;
        for(auto& chain: chains) {
            if (chain.selected) {
                for (auto& residue: chain.residues) {
                    if (residue.selected) {
                        for (auto& atom: residue.atoms) {
                            ret+=atom.gemmi.pos;
                            n++;
                        }
                    }
                }
            }
        }
        if (n>0) {
            ret /= n;
        }
        DEBUG_LOG("Mass centre: {}:{}:{}",ret.x,ret.y,ret.z);
        return ret;
    }

    void Protein::transform() {
        eachChain(
            [&](Tmdet::VOs::Chain& chain) -> void {
                chain.transform(tmatrix);
            }
        );
        for(auto& ssVec: secStrVecs) {
            tmatrix.transform(ssVec.begin);
            tmatrix.transform(ssVec.end);
        }
    
    }

}
