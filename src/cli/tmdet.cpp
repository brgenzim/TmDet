#include <iostream>
#include <string>
#include <args.hpp>
#include <TmdetStruct.hpp>
#include <TmdetXml.hpp>

using namespace std;

UniTmp::Utils::Args setArguments(int argc, char *argv[]);
void notTransmembrane(string x);

int main(int argc, char *argv[]) {
    UniTmp::Utils::Args args = setArguments(argc,argv);
    string inputPath = args.getValueAsString("i");
    string xmlPath = args.getValueAsString("x");
    string outputPdbPath = args.getValueAsString("p");
    bool n = args.getValueAsBool("n");
    bool nd = args.getValueAsBool("nd");
    bool tm = args.getValueAsBool("tm");
    if (n) {
        notTransmembrane(xmlPath);
    }

}

UniTmp::Utils::Args setArguments(int argc, char *argv[]) {
    UniTmp::Utils::Args args;
    args.define(false,"i","input","Input PDB file path (either cif or ent)","string","");
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
    UniTmp::TmdetLib::TmdetStruct tmdet;
    tmdet.read(xmlPath);
    tmdet.tmp(false);
    tmdet.write(xmlPath);
    exit(EXIT_SUCCESS);
}