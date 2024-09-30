#include <iostream>
#include <sstream>
#include <string>
#include <filesystem>
#include <gemmi/cif.hpp>
#include <gemmi/gz.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/model.hpp>
#include <System/Environment.hpp>
#include <ValueObjects/TmdetStruct.hpp>
#include <DTOs/TmdetStruct.hpp>
#include <Utils/Dssp.hpp>
#include <Utils/SecStrVec.hpp>

void assertTrue(std::string testDescription, bool condition, int lineNumber);
void setup(Tmdet::Utils::SecStrVec &secStructVectors);

std::string fileName;
Tmdet::System::Environment environment;

int main(int argc, char *argv[], char **envp) {

    environment.init(envp, ".env");
    fileName = std::filesystem::path(__FILE__).filename();

    Tmdet::Utils::SecStrVec secStructVectors;
    setup(secStructVectors);

    //
    // Test intersections of plane and vectors
    int numBoth = 0;
    int numDown = 0;
    int numUp = 0;
    std::string testDescription;

    Tmdet::ValueObjects::Membrane membrane;

    // Test case 1
    {
        testDescription = "two vectors intersects with two planes (numBoth: 2)";
        membrane.origo = gemmi::Vec3(0, 0, 0);
        membrane.normal = gemmi::Vec3(1, 0, 0);
        membrane.h = 7;
        membrane.type = Tmdet::Types::MembraneType::PLAIN;
        secStructVectors.numCross(membrane, numBoth, numUp, numDown);
        assertTrue(testDescription, numBoth == 2 && numUp == 0 && numDown == 0, __LINE__);
    }

    // Test case 2
    {
        testDescription = "no intersections (two vectors between two planes)";
        membrane.origo = gemmi::Vec3(0, 0, 0);
        membrane.normal = gemmi::Vec3(1, 0, 0);
        membrane.h = 20; // between two planes
        membrane.type = Tmdet::Types::MembraneType::PLAIN;
        secStructVectors.numCross(membrane, numBoth, numUp, numDown);
        assertTrue(testDescription, numBoth == 0 && numUp == 0 && numDown == 0, __LINE__);
    }

    // Test case 3
    {
        testDescription = "no intersections (planes are behind vector starting points)";
        // membrane planes are "under" helicies
        membrane.origo = gemmi::Vec3(-21, 0, 0);
        membrane.normal = gemmi::Vec3(1, 0, 0);
        membrane.h = 5;
        membrane.type = Tmdet::Types::MembraneType::PLAIN;
        secStructVectors.numCross(membrane, numBoth, numUp, numDown);
        assertTrue(testDescription, numBoth == 0 && numUp == 0 && numDown == 0, __LINE__);
    }

    // Test case 4
    {
        testDescription = "no intersections with inner sphere cap (numUp: 2 and numDown: 0)";
        // curved membrane
        membrane.origo = gemmi::Vec3(0, 0, 0);
        membrane.normal = gemmi::Vec3(1, 0, 0);
        membrane.h = 2;
        membrane.curver = 5;
        membrane.type = Tmdet::Types::MembraneType::CURVED;
        secStructVectors.numCross(membrane, numBoth, numUp, numDown);
        // no intersection with inner sphere
        assertTrue(testDescription, numBoth == 0 && numUp == 2 && numDown == 0, __LINE__);
    }

    // Test case 5
    {
        testDescription = "no intersections with inner and outer caps (curves 'behind' vectors' start)";
        // curved membrane with no intersections
        membrane.origo = gemmi::Vec3(-20.196, -3.826, -1.197);
        membrane.normal = gemmi::Vec3(1, 0, 0);
        membrane.h = 2;
        membrane.curver = 3;
        membrane.type = Tmdet::Types::MembraneType::CURVED;
        secStructVectors.numCross(membrane, numBoth, numUp, numDown);
        // no intersection with inner sphere
        assertTrue(testDescription, numBoth == 0 && numUp == 0 && numDown == 0, __LINE__);
    }

    // Test case 6
    {
        testDescription = "only one vector has an intersection with outer sphere cap (numUp: 1)";
        // curved membrane with UP intersection
        membrane.origo = gemmi::Vec3(-20.196, -3.826, -1.197);
        membrane.normal = gemmi::Vec3(1, 0, 0);
        membrane.h = 2;
        membrane.curver = 3.4;
        membrane.type = Tmdet::Types::MembraneType::CURVED;
        secStructVectors.numCross(membrane, numBoth, numUp, numDown);
        // no intersection with inner sphere
        assertTrue(testDescription, numBoth == 0 && numUp == 1 && numDown == 0, __LINE__);
    }

    return 0;
}

void assertTrue(std::string testDescription, bool condition, int lineNumber) {
    std::cout << (condition ? "Passed: " : "Failed: ") << testDescription;
    if (!condition) {
        std::cout << " (at line " << fileName << ":" << lineNumber << ")";
    }
    std::cout << std::endl;
}

void setup(Tmdet::Utils::SecStrVec &secStructVectors) {
    auto inputPath = environment.get("PDB_CIF_DIR");

    inputPath += "/af/1afo_updated.cif.gz";

    gemmi::Structure pdb; 
    gemmi::cif::Document document;
    auto tmdetVO = Tmdet::ValueObjects::get(inputPath, pdb, document);

    Tmdet::Utils::Dssp dssp = Tmdet::Utils::Dssp(tmdetVO);
    dssp.calcDsspOnStructure();
    dssp.writeDsspOnStructure();

    secStructVectors.define(tmdetVO);
}