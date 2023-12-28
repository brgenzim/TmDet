#include <fstream>
#include <iostream>
#include <gemmi/cif.hpp>
#include <gemmi/gz.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/model.hpp>
#include <pdbDssp.hpp>
#include <pdbDist.hpp>

namespace cif = gemmi::cif;

__attribute__((noreturn)) void usage(const char *fileName);

using namespace UniTmp::PdbLib::Utils;

int main(int argc, char *argv[]) {

    auto input = (char *)nullptr;
    for(int i=1; i<argc-1; i++) {
        if (!strcmp(argv[i],"-i")) {
            input = argv[i+1];
        }
    }
    if (input == (char *)nullptr) {
        usage(argv[0]);
    }
    //cif::Document doc = cif::read(gemmi::MaybeGzipped(input));
    gemmi::Structure pdb = gemmi::make_structure(cif::read(gemmi::MaybeGzipped(input)));
    auto distHelper = pdbDistance(pdb, 9.0);
    distHelper.listBoxElements();
    distHelper.setNeighboursForResiduesByCa();
    distHelper.listCaNeighbours();
    UniTmp::PdbLib::Utils::pdbCalcDsspOnStructure(pdb);
    UniTmp::PdbLib::Utils::pdbWriteDsspOnStructure(pdb);
}

void usage(const char *fileName) {
    std::cerr << "Usage: " << fileName << R"(
    -i : input cif file
)";
    exit(EXIT_FAILURE);
}
