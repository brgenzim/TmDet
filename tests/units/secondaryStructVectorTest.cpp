#include <iostream>
#include <sstream>
#include <string>
#include <gemmi/cif.hpp>
#include <gemmi/gz.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/model.hpp>
#include <Services/ConfigurationService.hpp>
#include <ValueObjects/TmdetStruct.hpp>
#include <DTOs/TmdetStruct.hpp>
#include <Utils/Dssp.hpp>
#include <Utils/SecStrVec.hpp>

int main() {
    gemmi::Structure pdb;
    Tmdet::Services::ConfigurationService::init();
    auto basePath = Tmdet::Services::ConfigurationService::getValue(Tmdet::Services::ConfigurationService::Keys::PDB_DIRECTORY);
    auto inputPath(basePath);
    inputPath += "/af/1afo_updated.cif.gz";

    //input is mandatory
    gemmi::cif::Document document = gemmi::cif::read(gemmi::MaybeGzipped(inputPath));
    pdb = gemmi::make_structure(std::move(document));
    Tmdet::ValueObjects::TmdetStruct tmdetVO = Tmdet::ValueObjects::TmdetStruct(pdb, document);
    tmdetVO.inputPath = inputPath;
    Tmdet::DTOS::TmdetStruct::parse(tmdetVO);

    Tmdet::Utils::Dssp dssp = Tmdet::Utils::Dssp(tmdetVO);
    dssp.calcDsspOnStructure();
    dssp.writeDsspOnStructure();

    Tmdet::Utils::SecStrVec secStructVectors;
    secStructVectors.define(tmdetVO);

    return 0;
}

std::string vec2string(gemmi::Vec3 vector) {
    std::stringstream result;
    result << "[ " << vector.x << ", "
        << vector.y << ", "
        << vector.z << " ]";
    return result.str();
}
