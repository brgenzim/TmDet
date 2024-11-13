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
            if (c.selected && c.id == id) {
                return c.idx;
            }
        }
        logger.error("Chain {} not found, or not selected",id);
        return -1;
    }

    int Protein::searchChainByLabId(const std::string& id) const {
        for (auto& c: chains) {
            if (c.selected && c.labId == id) {
                return c.idx;
            }
        }
        logger.error("Chain {} not found, or not selected",id);
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

    template<typename T>
    void Protein::eachChain(T* obj, void (T::*func)(Tmdet::ValueObjects::Chain& chain)) {
        for(auto& chain: chains) {
            (obj->*func)(chain);
        }
    }

    template<typename T>
    void Protein::eachSelectedChain(T* obj, void (T::*func)(Tmdet::ValueObjects::Chain& chain)) {
        for(auto& chain: chains) {
            if (chain.selected) {
                (obj->*func)(chain);
            }
        }
    }

}
