#include <string>
#include <gemmi/cifdoc.hpp>
#include <gemmi/gz.hpp>
#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <DTOs/TmdetStruct.hpp>

namespace Tmdet::ValueObjects {

    TmdetStruct get(const std::string &inputPath) {
        gemmi::Structure pdb;
        gemmi::cif::Document document = gemmi::cif::read(gemmi::MaybeGzipped(inputPath));
        pdb = gemmi::make_structure(std::move(document));
        auto tmdetVO = Tmdet::ValueObjects::TmdetStruct(pdb, document);
        Tmdet::DTOS::TmdetStruct::parse(tmdetVO);
        return tmdetVO;
    }
}
