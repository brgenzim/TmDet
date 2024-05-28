#include <iostream>
#include <string>
#include <Utils/Args.hpp>
#include <Utils/Dssp.hpp>
#include <Utils/Surface.hpp>
#include <ValueObjects/TmdetStruct.hpp>
#include <DTOs/TmdetStruct.hpp>
#include <gemmi/cif.hpp>
#include <gemmi/cifdoc.hpp>
#include <gemmi/gz.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/model.hpp>

using namespace std;

Tmdet::Utils::Args setArguments(int argc, char *argv[]);
void notTransmembrane(string x, Tmdet::ValueObjects::TmdetStruct& tmdetVO);

int main(int argc, char *argv[]) {
    gemmi::Structure pdb;
    Tmdet::Utils::Args args = setArguments(argc,argv);
    string inputPath = args.getValueAsString("i");
    string xmlPath = args.getValueAsString("x");
    string outputPdbPath = args.getValueAsString("p");
    bool n = args.getValueAsBool("n");
    bool nd = args.getValueAsBool("nd");
    bool tm = args.getValueAsBool("tm");
    
    //input is mandatory
    gemmi::cif::Document document = gemmi::cif::read(gemmi::MaybeGzipped(inputPath));
    pdb = gemmi::make_structure(std::move(document));
    Tmdet::ValueObjects::TmdetStruct tmdetVO = Tmdet::ValueObjects::TmdetStruct(pdb, document);
    tmdetVO.inputPath = inputPath;
    Tmdet::DTOS::TmdetStruct::parse(tmdetVO);
    Tmdet::Utils::Dssp dssp = Tmdet::Utils::Dssp(tmdetVO);
    dssp.writeDsspOnStructure();
    Tmdet::Utils::Surface surf = Tmdet::Utils::Surface(tmdetVO);
    surf.main();
    surf.setOutsideSurface();
    Tmdet::DTOS::TmdetStruct::out(tmdetVO);
    if (n) {
        notTransmembrane(xmlPath, tmdetVO);
    }

}

Tmdet::Utils::Args setArguments(int argc, char *argv[]) {
    Tmdet::Utils::Args args;
    args.define(true,"i","input","Input PDB file path (in cif format)","string","");
    args.define(true,"x","xml","Input/output xml file path","string","");
    args.define(false,"p","pdb_out","Output pdb file path","string","");
    args.define(false,"n","not","Set transmembrane='not' in the xml file","bool","false");
    args.define(false,"nd","force_nodel","Force not to delete not connected chains","bool","false");
    args.define(false,"tm","force_transmembrane","Set type to transmembrane without making decision using Q value","bool","false");
    args.set(argc,argv);
    args.check();
    return args;
}

void notTransmembrane(string xmlPath, Tmdet::ValueObjects::TmdetStruct& tmdetVO) {
    Tmdet::DTOS::TmdetStruct::readXml(tmdetVO, xmlPath);
    tmdetVO.tmp = false;
    Tmdet::DTOS::TmdetStruct::writeXml(tmdetVO, xmlPath);
    exit(EXIT_SUCCESS);
}
