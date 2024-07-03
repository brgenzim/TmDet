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

void assertTrue(std::string caseName, bool condition);

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

    // Test case 1
    {
        membrane.origo = gemmi::Vec3(0, 0, 0);
        membrane.normal = gemmi::Vec3(1, 0, 0);
        membrane.h = 7;
        membrane.type = Tmdet::Types::MembraneType::PLAIN;
        secStructVectors.numCross(membrane, numBoth, numUp, numDown);
        assertTrue("Case1: numBoth == 2", numBoth == 2 && numUp == 0 && numDown == 0);
    }

    numBoth = numDown = numUp = 0;

    // Test case 2
    {
        membrane.origo = gemmi::Vec3(0, 0, 0);
        membrane.normal = gemmi::Vec3(1, 0, 0);
        membrane.h = 20; // between two planes
        membrane.type = Tmdet::Types::MembraneType::PLAIN;
        secStructVectors.numCross(membrane, numBoth, numUp, numDown);
        assertTrue("Case2: numBoth == 0", numBoth == 0 && numUp == 0 && numDown == 0);
    }

    // Test case 3
    {
        // membrane planes are under helicies
        membrane.origo = gemmi::Vec3(-21, 0, 0);
        membrane.normal = gemmi::Vec3(1, 0, 0);
        membrane.h = 5;
        membrane.type = Tmdet::Types::MembraneType::PLAIN;
        secStructVectors.numCross(membrane, numBoth, numUp, numDown);
        assertTrue("Case3: numBoth == 0", numBoth == 0 && numUp == 0 && numDown == 0);
    }

    return 0;
}

void assertTrue(std::string caseName, bool condition) {
    std::cout << caseName << ": "
        << (condition ? "SUCCESS" : "FAIL") << std::endl;
}
