#include <iostream>
#include "../Types/PdbAtomTypes.hpp"
#include "../Types/PdbChainTypes.hpp"
#include "../Types/PdbMethodTypes.hpp"
#include "../Types/PdbMoleculeTypes.hpp"
#include <gemmi/cif.hpp>       // file -> cif::Document
#include <gemmi/gz.hpp>        // uncompressing on the fly
#include <gemmi/mmcif.hpp>     // cif::Document -> Structure

namespace cif = gemmi::cif;

int main(int argc, char *argv[]) {
    //std::cout << UniTmp::PdbLib::Types::PdbAtomTypes::all.at("AL").name << std::endl;
    //std::cout << UniTmp::PdbLib::Types::PdbAtomTypes::FE.vdw << std::endl;
    UniTmp::PdbLib::Types::PdbChainType c = UniTmp::PdbLib::Types::PdbChainTypes::NUC;
    //std::cout << c.description << std::endl;
    //std::cout << UniTmp::PdbLib::Types::PdbMethodTypes::Xray.name << std::endl;
    //std::cout << UniTmp::PdbLib::Types::PdbMethodTypes::fromName.at("Solution scattering").code << std::endl;
    //std::cout << UniTmp::PdbLib::Types::PdbMoleculeTypes::ALA.name << std::endl;

    cif::Document doc = cif::read(gemmi::MaybeGzipped("/home/tusi/works/pdbtm_4.0/PdbLib/data/Components-pub.cif"));
    for (cif::Block& block : doc.blocks) {
        std::cout << block.name << std::endl;
    }
    //gemmi::Structure structure = gemmi::make_structure(doc);
}