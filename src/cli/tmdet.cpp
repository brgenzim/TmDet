#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <Services/ChemicalComponentDirectoryService.hpp>
#include <System/Arguments.hpp>
#include <System/Environment.hpp>
#include <System/Config.hpp>
#include <System/FilePaths.hpp>
#include <ValueObjects/Protein.hpp>
#include <DTOs/Protein.hpp>
#include <DTOs/Xml.hpp>
#include <Engine/Organizer.hpp>

using namespace std;

Tmdet::System::Environment environment;

Tmdet::System::Arguments getArguments(int argc, char *argv[]) {
    Tmdet::System::Arguments args;
    args.define(false,"e","env","Path for environment variable file","string",".env");
    args.define(false,"i","input","Input PDB file full path (in cif format)","string","");
    args.define(false,"c","code","Input PDB code (c or i is mandatory)","string","");
    args.define(false,"x","xml","Input/output xml file path","string","");
    args.define(false,"p","pdb_out","Output pdb file path","string","");
    args.define(false,"o","over_write","Over write xml file instead of update","bool","false");
    args.define(false,"n","not","Set transmembrane='not' in the xml file","bool","false");
    args.define(false,"nc","no_cache","Do not use cached data","bool","false");
    args.define(false,"nd","force_nodel","Force not to delete not connected chains","bool","false");
    args.define(false,"tm","force_transmembrane","Set type to transmembrane without making decision using Q value","bool","false");
    args.set(argc,argv);
    args.check();
    return args;
}

void notTransmembrane(const std::string& xmlPath) {
    Tmdet::ValueObjects::Protein protein;
    Tmdet::DTOs::Xml xml;
    xml.readXml(protein, xmlPath);
    protein.notTransmembrane();
    xml.writeXml(protein, xmlPath);
}

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

    //setting input, output paths
    //if code is given then system directories are used
    //else user should provide the full path of xml and cif files
    string code = args.getValueAsString("c");
    string xmlPath = (code != ""?Tmdet::System::FilePaths::xml(code):args.getValueAsString("x"));
    string inputPath = (code != ""?Tmdet::System::FilePaths::cif(code):args.getValueAsString("i"));
    string outputPdbPath = (code != ""?Tmdet::System::FilePaths::pdbOut(code):args.getValueAsString("p"));

    //check xmlPath
    if (xmlPath == "") {
        std::cerr << "Error: -x or -c is mandatory" << std::endl;
        args.list();
        exit(EXIT_FAILURE);
    }
    
    //change xml file to TMP="no" if the protein is set to not transmembrane and exit
    if (bool n = args.getValueAsBool("n"); n) {
        if ( !Tmdet::System::FilePaths::fileExists(xmlPath) ) {
            std::cerr << "Error: file not found: " << xmlPath << std::endl;
            exit(EXIT_FAILURE);
        }
        notTransmembrane(xmlPath);
        exit(EXIT_SUCCESS);
    }
    
    //if -n or --not is not set then input is mandatory
    if (inputPath == "") {
        std::cerr << "Error: -i or -c is mandatory" << std::endl;
        args.list();
        exit(EXIT_FAILURE);
    }
    if ( !Tmdet::System::FilePaths::fileExists(inputPath) ) {
        std::cerr << "Error: file not found: " << inputPath << std::endl;
        exit(EXIT_FAILURE);
    }
    auto protein = Tmdet::DTOs::Protein::get(inputPath);
    Tmdet::DTOs::Xml xml;

    //if xml file exists and overwirte is not set then read the content of the
    //xml file and adjust proteinVO accordingly
    if (bool o = args.getValueAsBool("o"); !o && Tmdet::System::FilePaths::fileExists(xmlPath) ) {
        std::cerr << "xml file exists, updating its content" << std::endl;
        xml.readXml(protein, xmlPath);
        //TODO align structure to tmdet data
    }
    else {
        if (Tmdet::System::FilePaths::fileExists(xmlPath)) {
            protein.clear();
            std::cerr << "overwriting xml file" << std::endl;
        }
    }
    
    
    //do the membrane region determination and annotation
    auto organizer = Tmdet::Engine::Organizer(protein, args);
}
