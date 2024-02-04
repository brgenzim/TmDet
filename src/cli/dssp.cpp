#include <fstream>
#include <iostream>
#include <gemmi/cif.hpp>
#include <gemmi/gz.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/model.hpp>
#include <pdbDssp.hpp>
#include <pdbDist.hpp>
#include <pdbArgs.hpp>

namespace cif = gemmi::cif;

__attribute__((noreturn)) void usage(const char *fileName);

using namespace UniTmp::PdbLib::Utils;

pdbArgs setArguments(int argc, char *argv[]);

int main(int argc, char *argv[]) {

    pdbArgs args = setArguments(argc,argv);
    auto input = args.getValueAsString("i");
    //cif::Document doc = cif::read(gemmi::MaybeGzipped(input));
    gemmi::Structure pdb = gemmi::make_structure(cif::read(gemmi::MaybeGzipped(input)));
    auto distHelper = pdbDistance(pdb, 9.0);
    distHelper.listBoxElements();
    distHelper.setNeighboursForResiduesByCa();
    distHelper.listCaNeighbours();
    UniTmp::PdbLib::Utils::pdbCalcDsspOnStructure(pdb);
    UniTmp::PdbLib::Utils::pdbWriteDsspOnStructure(pdb);
}

pdbArgs setArguments(int argc, char *argv[]) {
    pdbArgs args;
    args.define(true,"i","input","Input cif file path","string","");
    args.set(argc,argv);
    args.check();
    return args;
}