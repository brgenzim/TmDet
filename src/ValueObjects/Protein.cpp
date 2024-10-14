#include <string>
#include <gemmi/cifdoc.hpp>
#include <gemmi/gz.hpp>
#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <Version.hpp>
#include <System/Date.hpp>
#include <Types/Protein.hpp>
#include <ValueObjects/Protein.hpp>
#include <Utils/Md5.hpp>

namespace Tmdet::ValueObjects {

    void Protein::getStructure(const std::string& inputPath) {
        document = gemmi::cif::read(gemmi::MaybeGzipped(inputPath));
        gemmi = gemmi::make_structure(std::move(document));
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

    int Protein::searchChainById(std::string& id) const {
        for (auto& c: chains) {
            if (c.id == id) {
                return c.idx;
            }
        }
        return -1;
    }
}
