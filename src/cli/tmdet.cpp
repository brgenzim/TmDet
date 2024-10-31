// Copyright(c) 2003-present, Gabor E. Tusnady & tmdet contributors.
// Distributed under the MIT License (http://opensource.org/licenses/MIT)

#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <Config.hpp>
#include <Version.hpp>
#include <DTOs/Protein.hpp>
#include <DTOs/Xml/Reader3.hpp>
#include <DTOs/Xml/Writer.hpp>
#include <Engine/Organizer.hpp>
#include <Helpers/Pymol.hpp>
#include <Services/ChemicalComponentDirectoryService.hpp>
#include <System/Arguments.hpp>
#include <System/Environment.hpp>
#include <System/FilePaths.hpp>
#include <System/Logger.hpp>
#include <ValueObjects/Protein.hpp>

using namespace std;

Tmdet::System::Environment environment;
Tmdet::System::Logger logger;

Tmdet::System::Arguments getArguments(int argc, char *argv[]) {
    Tmdet::System::Arguments args;
    args.define(false,"e","env","Path for environment variable file","string",".env");
    args.define(false,"i","input","Input PDB file full path (in cif format)","string","");
    args.define(false,"c","code","Input PDB code (c or i is mandatory)","string","");
    args.define(false,"x","xml","Input/Output xml file path","string","");
    args.define(false,"po","pdb_out","Output pdb file path","string","");
    args.define(false,"n","not","Set transmembrane='not' in the xml file","bool","false");
    args.define(false,"nc","no_cache","Do not use cached data","bool","false");
    //args.define(false,"nd","force_nodel","Force not to delete not connected chains","bool","false");
    //args.define(false,"tm","force_transmembrane","Set type to transmembrane without making decision using Q value","bool","false");
    args.set(argc,argv);
    args.check();
    return args;
}

void notTransmembrane(const std::string& xmlPath) {
    Tmdet::ValueObjects::Protein protein;
    Tmdet::DTOs::Xml::Reader3 reader;
    reader.readXml(protein, xmlPath);
    protein.notTransmembrane();
    //xml.writeXml(protein, xmlPath);
}

int main(int argc, char *argv[], char **envp) {
    
    //get and check command line arguments
    Tmdet::System::Arguments args = getArguments(argc,argv);

    //get environment file content and shell environment variables
    environment.init(envp,args.getValueAsString("e"));

    //setting up logger
    std::ostream& coutRef = std::cout;
    std::ofstream logFile(environment.get("TMDET_LOG_FILE",DEFAULT_TMDET_LOG_FILE), std::ios_base::app);
    logger.addStream(coutRef);
    logger.addStream(logFile);
    logger.setLevel(Tmdet::System::level::TMDET_LOG_LEVEL);
    logger.info("Starting Tmdet(version: {})",Tmdet::version());
    std::string a;
    for(int i=1; i<argc; i++) {
        a += argv[i]; a += " ";
    }
    logger.info("command line arguments: {}",a);

    //check ccd and fetch it if missing
    if (!Tmdet::Services::ChemicalComponentDirectoryService::isBuilt()) {
        logger.warn("Chemical component directory is not set, please wait while installing it.");
        Tmdet::Services::ChemicalComponentDirectoryService::fetch();
        Tmdet::Services::ChemicalComponentDirectoryService::build();
    }

    //setting input, output paths
    //if code is given then system directories are used
    //else user should provide the full path of xml and cif files
    string code = args.getValueAsString("c");
    string xmlPath = (code != ""?Tmdet::System::FilePaths::xml(code,true):args.getValueAsString("x"));
    string inputPath = (code != ""?Tmdet::System::FilePaths::cif(code):args.getValueAsString("i"));
    string outputPdbPath = (code != ""?Tmdet::System::FilePaths::pdbOut(code):args.getValueAsString("p"));

    //check xmlPath
    if (xmlPath == "") {
        logger.error("argument -x or -c is mandatory");
        args.list();
        exit(EXIT_FAILURE);
    }
    
    //change xml file to TMP="no" if the protein is set to not transmembrane and exit
    if (bool n = args.getValueAsBool("n"); n) {
        if ( !Tmdet::System::FilePaths::fileExists(xmlPath) ) {
            logger.error("file not found: {}",xmlPath);
            exit(EXIT_FAILURE);
        }
        notTransmembrane(xmlPath);
        exit(EXIT_SUCCESS);
    }
    
    //if -n or --not is not set then input is mandatory
    if (inputPath == "") {
        logger.error("argument -i or -c is mandatory");
        args.list();
        exit(EXIT_FAILURE);
    }
    if ( !Tmdet::System::FilePaths::fileExists(inputPath) ) {
        logger.error("file not found: {}",inputPath);
        exit(EXIT_FAILURE);
    }
    auto protein = Tmdet::DTOs::Protein::get(inputPath);
    
    //if xml file exists and overwirte is not set then read the content of the
    //xml file and adjust proteinVO accordingly
    if (Tmdet::System::FilePaths::fileExists(xmlPath) ) {
        //Tmdet::DTOs::Xml::Reader3 xmlInput;
        //logger.warn("overwriting xml file");
        //xmlInput.readXml(protein, xmlPath);
    }
    
    //do the membrane region determination and annotation
    auto organizer = Tmdet::Engine::Organizer(protein, args);
    Tmdet::DTOs::Xml::Writer xmlOutput;
    xmlOutput.writeXml(protein, xmlPath);
    Tmdet::DTOs::Protein::writeCif(protein,outputPdbPath);
    auto pymol = Tmdet::Helpers::Pymol(protein);
    pymol.show();
}
