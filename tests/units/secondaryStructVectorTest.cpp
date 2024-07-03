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
    //
    // Test intersections of plane and vectors
    int numBoth = 0;
    int numDown = 0;
    int numUp = 0;

    Tmdet::ValueObjects::Membrane membrane;
    membrane.origo = gemmi::Vec3(0, 0, 0);
    membrane.normal = gemmi::Vec3(1, 0, 0);
    membrane.h = 7;
    membrane.type = Tmdet::Types::MembraneType::PLAIN;
    secStructVectors.numCross(membrane, numBoth, numUp, numDown);
    std::cout << "numBoth: " << numBoth << " numUp: " << numUp << " numDown: " << std::endl;

    numBoth = numDown = numUp = 0;

    return 0;
}

std::string vec2string(gemmi::Vec3 vector) {
    std::stringstream result;
    result << "[ " << vector.x << ", "
        << vector.y << ", "
        << vector.z << " ]";
    return result.str();
}
