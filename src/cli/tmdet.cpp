#include <iostream>
#include <string>
#include <Utils/Args.hpp>
#include <ValueObjects/Struct.hpp>
#include <DTOs/Struct.hpp>
#include <gemmi/cif.hpp>
#include <gemmi/gz.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/model.hpp>


using namespace std;

Tmdet::Utils::Args setArguments(int argc, char *argv[]);
void notTransmembrane(string x);

int main(int argc, char *argv[]) {
    gemmi::Structure pdb;
    Tmdet::Utils::Args args = setArguments(argc,argv);
    string inputPath = args.getValueAsString("i");
    string xmlPath = args.getValueAsString("x");
    string outputPdbPath = args.getValueAsString("p");
    bool n = args.getValueAsBool("n");
    bool nd = args.getValueAsBool("nd");
    bool tm = args.getValueAsBool("tm");
    if (n) {
        notTransmembrane(xmlPath);
    }
    if (inputPath != "") {
        pdb = gemmi::make_structure(cif::read(gemmi::MaybeGzipped(inputPath)));
    }

}

Tmdet::Utils::Args setArguments(int argc, char *argv[]) {
    Tmdet::Utils::Args args;
    args.define(false,"i","input","Input PDB file path (in cif format)","string","");
    args.define(true,"x","xml","Input/output xml file path","string","");
    args.define(false,"p","pdb_out","Output pdb file path","string","");
    args.define(false,"n","not","Set transmembrane='not' in the xml file","bool","false");
    args.define(false,"nd","force_nodel","Force not to delete not connected chains","bool","false");
    args.define(false,"tm","force_transmembrane","Set type to transmembrane without making decision using Q value","bool","false");
    args.set(argc,argv);
    args.check();
    return args;
}

void notTransmembrane(string xmlPath) {
    Tmdet::ValueObjects::Struct tmdetVO;
    Tmdet::DTOS::Struct::read(tmdetVO, xmlPath);
    tmdetVO.tmp = false;
    Tmdet::DTOS::Struct::write(tmdetVO, xmlPath);
    exit(EXIT_SUCCESS);
}