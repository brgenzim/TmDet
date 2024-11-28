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
#include <DTOs/Xml.hpp>
#include <Engine/Organizer.hpp>
#include <Helpers/Pymol.hpp>
#include <Services/ChemicalComponentDirectoryService.hpp>
#include <System/Arguments.hpp>
#include <System/Date.hpp>
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
    args.define(false,"r","run","Run the tmdet algorithm on the protein structure","bool","false");
    args.define(false,"pi","pdb_input","Input PDB file full path (in cif format)","string","");
    args.define(false,"c","code","Input PDB code (c or pi is mandatory)","string","");
    args.define(false,"x","xml","Input/Output xml file path","string","");
    args.define(false,"xi","xml_input","Input xml file path","string","");
    args.define(false,"xo","xml_output","Output xml file path","string","");
    args.define(false,"po","pdb_output","Output pdb file path","string","");
    args.define(false,"n","not","Set transmembrane='not' in the xml file","bool","false");
    args.define(false,"nc","no_cache","Do not use cached data","bool","false");
    args.define(false,"s","show","Show annotated structure by pymol","bool","false");
    args.define(false,"sp","show","Show parsed structure in the console","bool","false");
    args.define(false,"uc","unselect_chains","Unselect proteins chains","string","");
    args.define(false,"na","no_annotation","Do not make annotation","bool","false");
    args.define(false,"fa","force_nodel_antibody","Do not unselect antibodies in the structure","bool","false");
    args.define(false,"tm","force_transmembrane","Set type to transmembrane without making decision using Q value","bool","false");
    args.set(argc,argv);
    args.check();
    return args;
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
    INFO_LOG("Starting Tmdet(version: {})",Tmdet::version());
    INFO_LOG("Command line arguments: {}",args.getCommandLine());

    //check ccd and fetch it if missing
    if (!Tmdet::Services::ChemicalComponentDirectoryService::isBuilt()) {
        WARN_LOG("Chemical component directory is not set, please wait while installing it.");
        Tmdet::Services::ChemicalComponentDirectoryService::fetch();
        Tmdet::Services::ChemicalComponentDirectoryService::build();
    }

    //setting input, output paths
    //if code is given then system directories are used
    //else user should provide the full path of xml and cif files
    Tmdet::DTOs::Xml xml;
    string code = args.getValueAsString("c");
    string xmlInputPath = xml.setPath(code,args.getValueAsString("x"),args.getValueAsString("xi"));
    string xmlOutputPath = xml.setPath(code,args.getValueAsString("x"),args.getValueAsString("xo"));
    string pdbInputPath = (code != ""?Tmdet::System::FilePaths::cif(code):args.getValueAsString("pi"));
    string pdbOutputPath = (code != ""?Tmdet::System::FilePaths::pdbOut(code):args.getValueAsString("po"));

    //change xml file to TMP="no" if the protein is set to not transmembrane and exit
    if (bool n = args.getValueAsBool("n"); n) {
        xml.notTransmembrane(xmlInputPath,xmlOutputPath);
        exit(EXIT_SUCCESS);
    }
    
    //if -n or --not is not set then input is mandatory
    if (pdbInputPath == "") {
        ERROR_LOG("argument -pi or -c is mandatory");
        args.list();
        exit(EXIT_FAILURE);
    }

    //check existence of pdb input cif file
    if ( !Tmdet::System::FilePaths::fileExists(pdbInputPath) ) {
        logger.error("file not found: {}",pdbInputPath);
        exit(EXIT_FAILURE);
    }
    auto protein = Tmdet::DTOs::Protein::get(pdbInputPath);
    if (code != "") {
	protein.code = code;
    }
    bool na = args.getValueAsBool("na");

    //unselect antibodies if not prevented
    if (bool fa = args.getValueAsBool("fa"); !fa) {
        Tmdet::DTOs::Protein::unselectAntiBodyChains(protein);
    }

    //unselect chains given in the command line argument
    if (std::string uc = args.getValueAsString("uc"); uc != "") {
        Tmdet::DTOs::Protein::unselectChains(uc, protein);
    }

    //do the membrane region determination and annotation
    if (bool r = args.getValueAsBool("r"); r) {
        auto organizer = Tmdet::Engine::Organizer(protein, args);
        protein.version = Tmdet::version();
        protein.date = Tmdet::System::Date::get();

        //write xml output if required
        if (xmlOutputPath != "" && !na) {
            xml.write(xmlOutputPath, protein);
        }

        //write transformed pdb file if required and protein is tmp
        if (pdbOutputPath != "" && protein.tmp && !na) {
            Tmdet::DTOs::Protein::writeCif(protein,pdbOutputPath);
        }
    }
    else {
        //if xml input is given then read it's content
        if (xmlInputPath != "" && !na) {
            xml.read(xmlInputPath, protein);
        }    
    }

    //show protein by pymol if required
    if (bool s = args.getValueAsBool("s"); s && protein.tmp) {
        auto pymol = Tmdet::Helpers::Pymol(protein);
        pymol.show();
    }

    //print out protein structrure and properties
    if (bool sp = args.getValueAsBool("sp"); sp) {
        std::cout << Tmdet::DTOs::Protein::toString(protein) << std::endl;
    }
    
    //if no annotation flag is given, just witeout tm / not tm info
    if (na) {
       std::cout << protein.code << " " << (protein.tmp?"yes":"no") << " " << protein.qValue  << std::endl;
    }
}
