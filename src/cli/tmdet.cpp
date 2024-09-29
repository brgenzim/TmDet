#include <iostream>
#include <fstream>
#include <string>
#include <Services/ChemicalComponentDirectoryService.hpp>
#include <System/Arguments.hpp>
#include <System/Environment.hpp>
#include <Utils/Dssp.hpp>
#include <Utils/Surface.hpp>
#include <Utils/Symmetry.hpp>
#include <Utils/Fragment.hpp>
#include <ValueObjects/TmdetStruct.hpp>
#include <DTOs/TmdetStruct.hpp>

using namespace std;

Tmdet::System::Arguments getArguments(int argc, char *argv[]);
void notTransmembrane(string x, Tmdet::ValueObjects::TmdetStruct& tmdetVO);
Tmdet::System::Environment environment;

int main(int argc, char *argv[], char **envp) {

    //get and check command line arguments
    Tmdet::System::Arguments args = getArguments(argc,argv);

    //get environment file content and shell environment variables
    environment.init(envp,args.getValueAsString("e"));

    //check ccd and fetch it if missing
    if (!Tmdet::Services::ChemicalComponentDirectoryService::isBuilt()) {
        std::cerr << "Chemical component directory is not set, please wait while installing it." << std::endl;
        Tmdet::Services::ChemicalComponentDirectoryService::fetch();
        Tmdet::Services::ChemicalComponentDirectoryService::build();
    }
    string inputPath = args.getValueAsString("i");
    string xmlPath = args.getValueAsString("x");
    string outputPdbPath = args.getValueAsString("p");
    bool n = args.getValueAsBool("n");
    bool nd = args.getValueAsBool("nd");
    bool tm = args.getValueAsBool("tm");

    //input is mandatory, TmdetVO can not be created without gemmi structure and document
    gemmi::Structure pdb; 
    gemmi::cif::Document document;
    auto tmdetVO = Tmdet::ValueObjects::get(inputPath, pdb, document);

    //change xml file to TMP="no" if the protein is not transmembrane and exit
    if (n) {
        notTransmembrane(xmlPath, tmdetVO);
    }

    //do clustering on the whole structure
    auto fragmentEngine = Tmdet::Utils::Fragment(tmdetVO);
    fragmentEngine.run();
    
    //Tmdet::Utils::Symmetry symmetry;
    //auto result = symmetry.CheckSymmetry(tmdetVO);

    //Tmdet::Utils::Dssp dssp = Tmdet::Utils::Dssp(tmdetVO);
    //dssp.calcDsspOnStructure();
    //dssp.writeDsspOnStructure();
    //Tmdet::Utils::Surface surf = Tmdet::Utils::Surface(tmdetVO);
    //surf.main();
    //surf.setOutsideSurface();
    //Tmdet::DTOS::TmdetStruct::out(tmdetVO);
    

}

Tmdet::System::Arguments getArguments(int argc, char *argv[]) {
    Tmdet::System::Arguments args;
    args.define(false,"e","env","Path for environment variable file","string",".env");
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

void notTransmembrane(string xmlPath, Tmdet::ValueObjects::TmdetStruct& tmdetVO) {
    Tmdet::DTOS::TmdetStruct::readXml(tmdetVO, xmlPath);
    tmdetVO.tmp = false;
    Tmdet::DTOS::TmdetStruct::writeXml(tmdetVO, xmlPath);
    exit(EXIT_SUCCESS);
}
