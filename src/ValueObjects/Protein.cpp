#include <string>
#include <map>
#include <gemmi/cifdoc.hpp>
#include <gemmi/gz.hpp>
#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/model.hpp>
#include <gemmi/cifdoc.hpp>
#include <Version.hpp>
#include <Config.hpp>
#include <System/Date.hpp>
#include <Types/Protein.hpp>
#include <ValueObjects/Protein.hpp>
#include <Utils/Md5.hpp>

namespace Tmdet::ValueObjects {

    void Protein::getStructure(const std::string& inputPath) {
        document = gemmi::cif::read(gemmi::MaybeGzipped(inputPath));
        gemmi = gemmi::make_structure(std::move(document));

        const auto& block = document.blocks[0];

        //
        // Fill polymerNames map
        //

        for (auto& entity : gemmi.entities) {
            if (entity.entity_type == gemmi::EntityType::Polymer) {
                polymerNames[entity.name] = "DESCRIPTION: N/A";
            }
        }

        auto entityLoop = block.find_loop_item("_entity.id")->loop;
        int loopLength = entityLoop.length();
        int idCol = entityLoop.find_tag("_entity.id");
        int descriptionCol = entityLoop.find_tag("_entity.pdbx_description");
        for (int row = 0; row < loopLength; row++) {
            const std::string& id = entityLoop.val(row, idCol);
            if (polymerNames.contains(id)) {
                std::string name(entityLoop.val(row, descriptionCol));
                if (name.starts_with("'") && name.ends_with("'")) {
                    name = name.substr(1, name.length() - 2);
                }
                polymerNames[id] = name;
            }
        }
    }

    void Protein::notTransmembrane() {
        tmp = false;
        version = (version==""?Tmdet::version():version);
        modifications.emplace_back(
            Tmdet::System::Date::get(),
            (std::string)"Not transmembrane protein"
        );
        type = Tmdet::Types::ProteinType::SOLUBLE;
        bioMatrix.matrices.clear();
        bioMatrix.deletedChainIds.clear();
        membranes.clear();
        chains.clear();
    }

    void Protein::clear() {
        tmp = false;
        version = Tmdet::version();
        date = Tmdet::System::Date::get();
        modifications.clear();
        qValue = 0.0;
        type = Tmdet::Types::ProteinType::SOLUBLE;
        bioMatrix.matrices.clear();
        bioMatrix.deletedChainIds.clear();
        membranes.clear();
        chains.clear();
    }

    std::string Protein::hash() const {
        std::string raw = code;
        for(const auto& c: chains) {
            if (c.selected) {
                raw += c.seq;
            }
        }
        void* sig = hashing::md5::hash(raw);
        return hashing::md5::sig2hex(sig);
    }

    void Protein::unSelectAll() {
        for (auto& c: chains) {
            c.selected = false;
        }
    }

    int Protein::searchChainById(const std::string& id) const {
        for (auto& c: chains) {
            if (c.id == id) {
                return c.idx;
            }
        }
        logger.error("Chain {} not found",id);
        return -1;
    }

    int Protein::searchChainByLabId(const std::string& id) const {
        for (auto& c: chains) {
            if (c.labId == id) {
                return c.idx;
            }
        }
        logger.error("Chain {} not found",id);
        return -1;
    }

    gemmi::Vec3 Protein::centre() {
        gemmi::Vec3 ret(0,0,0);
        int n = 0;
        for(auto& chain: chains) {
            if (chain.selected) {
                ret += chain.centre();
                n++;
            }
        }
        if (n>0) {
            ret /= n;
        }
        return ret;
    }
}
