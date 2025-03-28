// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

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
#include <Engine/Fragmenter.hpp>
#include <Engine/Organizer.hpp>
#include <Helpers/Pymol.hpp>
#include <Services/ChemicalComponentDirectoryService.hpp>
#include <System/Arguments.hpp>
#include <System/Date.hpp>
#include <System/Environment.hpp>
#include <System/FilePaths.hpp>
#include <System/Logger.hpp>
#include <Utils/Dssp.hpp>
#include <Utils/Filter.hpp>
#include <Utils/MyDssp.hpp>
#include <Utils/NeighBors.hpp>
#include <Utils/SecStrVec.hpp>
#include <VOs/Protein.hpp>

using namespace std;

Tmdet::System::Environment environment;
Tmdet::System::Logger logger;

Tmdet::System::Arguments getArguments(int argc, char *argv[]) {
    Tmdet::System::Arguments args;
    
    //system
    args.define(false,false,"e","env","Path for environment variable file","string",".env");

    //path related
    args.define(false,true,"c","code","Input PDB code (c or pi is mandatory)","string","");
    args.define(false,false,"x","xml","Input/Output xml file path","string","");
    args.define(false,false,"pi","pdb_input","Input PDB file full path (in ent or cif format)","string","");
    args.define(false,false,"po","pdb_output","Output pdb file path","string","");
    args.define(false,false,"xi","xml_input","Input xml file path","string","");
    args.define(false,false,"xo","xml_output","Output xml file path","string","");
    args.define(false,true,"a","assembly","Set assembly id","int","1");
    args.define(false,true,"m","model","Set model id","int","0");
    
    //work 
    args.define(false,false,"r","run","Run the tmdet algorithm on the protein structure","bool","false");
    args.define(false,true,"n","not","Set transmembrane='not' in the xml file","bool","false");
    args.define(false,true,"cm","curved_membrane","Search for curved membrane","bool","false");
    args.define(false,true,"dm","duble_membrane","Enable duble membrane mode","bool","false");
    args.define(false,true,"fr","fragment_analysis","Investigate protein domains/fragments separately","bool","false");
    args.define(false,true,"bi","barrel_inside","Indicate chains those are within a barrel (but not part of barrel, like 5iv8)","string","");
    args.define(false,true,"ns","no_symmetry","Do not use symmetry axes as membrane normal","bool","false");
    args.define(false,true,"s","show","Show annotated structure by pymol","bool","false");
    args.define(false,true,"sp","show","Show parsed structure in the console","bool","false");
    args.define(false,true,"uc","unselect_chains","Unselect proteins chains","string","");
    args.define(false,true,"na","no_annotation","Do not make annotation","bool","false");
    args.define(false,true,"fa","force_nodel_antibody","Do not unselect antibodies in the structure","bool","false");
    args.define(false,true,"nc","no_cache","Do not use cached data","bool","false");
    args.define(false,true,"pf","pre_filter","Apply TmFilter","bool","false");
    args.define(false,true,"xf3","xml_out_fmt3","Set xml output format to v3","bool","false");
    
    //parameters
    args.define(false,true,"lq","lower_qvalue","Lower qValue, above it is membrane","float","30");
    args.define(false,true,"hq","higher_qvalue","Higher qValue, limit for transmembrane type","float","36");
    args.define(false,true,"hq2","higher_qvalue2","Higher qValue2, limit for second membrane","float","55");
    args.define(false,true,"minht","minimum_of_half_thickness","Minimum value of half thickness","float","10.0");
    args.define(false,true,"maxht","maximum_of_half_thickness","Maximum value of half thickness","float","20.0");
    args.define(false,true,"maxcht","maximum_of_curved_half_thickness","Maximum value of half thickness for curved membrane detection","float","14");
    args.define(false,true,"ihml","ifh_hydrph_limit","Hydrophobicity momentum limit for ifh detection","float","1.6");
    args.define(false,true,"ias","ifh_avg_surface","Average free solvent accessible surface limit for ifh detection","float","40");
    args.define(false,true,"ian","ifh_angle","Maximum angle between membrane plane and ifh","float","15");
    args.define(false,true,"iml","ifh_min_length","Minimum length of ifhs","int","6");
    args.define(false,true,"ba","boost_angle","Boost secondary structure element angle in optimization","float","0.7");
    args.define(false,true,"bba","boost_beta_angle","Boost beta sheet angle in optimization","float","0.55");
    args.define(false,true,"bp","boost_polarity","Boost polarity calculation in optimization","float","0.55");
    args.define(false,true,"lmhp","loop_min_helix_part","Minimum of a helix be part as re-entrant loop","float","0.25");
    args.define(false,true,"lmd","loop_min_depth","Minimum depth of a re-entrant loop in angstrom","float","3.0");
    args.define(false,true,"lmnss","loop_min_no_sec_str","Minimum number of residues in a re-entrant loop that has no secondary structure","int","1");
    args.define(false,true,"mums","max_unannotated_memb_segm","Maximum number of residues that can not be annotated","int","10");
    args.define(false,true,"sm","shift_membrane","Shift membrane with the given distance","float","0");
    args.define(false,true,"mltmh","min_length_of_tmh","Minimum length of transmembrane helix","int","12");
    args.define(false,true,"spen","straigth_penalty","Additional value for normalizing straigth","float","0.0");
    args.define(false,true,"mcbs","min_contacts_between_sheets","Minimum of contacts between sheets for barrel detection","int","5");
    args.define(false,true,"bh","broken_helix","Type of broken helix (loop or transmembrane helix)","string","L");
    args.define(false,true,"minbs","min_number_of_beta_sheets","Minimum number of beta sheets in beta barrel","int","8");
    
    args.set(argc,argv);
    args.check();
    return args;
}

int main(int argc, char *argv[], char **envp) {

    //get and check command line arguments
    Tmdet::System::Arguments args = getArguments(argc,argv);

    //get environment file content and shell environment variables
    environment.init(envp,args.getValueAsString("e"));
    args.setCommandLine();

    //setting up logger
    std::ostream& coutRef = std::cout;
    std::ofstream logFile(environment.get("TMDET_LOG_FILE",DEFAULT_TMDET_LOG_FILE), std::ios_base::app);
    logger.addStream(coutRef);
    logger.addStream(logFile);
    logger.setLevel(Tmdet::System::level::TMDET_LOG_LEVEL);

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
    if (bool xf3 = args.getValueAsBool("xf3"); xf3) {
        xml.setV3Fmt();
    }
    string code = args.getValueAsString("c");
    string xmlInputPath = xml.setPath(code,args.getValueAsString("x"),args.getValueAsString("xi"));
    string xmlOutputPath = xml.setPath(code,args.getValueAsString("x"),args.getValueAsString("xo"));
    string pdbInputPath = (code != ""?Tmdet::System::FilePaths::cif(code,args.getValueAsInt("a")):args.getValueAsString("pi"));
    string pdbOutputPath = (code != ""?Tmdet::System::FilePaths::pdbOut(code):args.getValueAsString("po"));

    //change xml file to TMP="no" if the protein is set to not transmembrane and exit
    if (bool n = args.getValueAsBool("n"); n) {
        xml.notTransmembrane(xmlInputPath, xmlOutputPath, args);
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
    auto protein = Tmdet::DTOs::Protein::get(pdbInputPath, args.getValueAsInt("m"));
    protein.forceSingleMembrane = !args.getValueAsBool("dm");

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

    if (int nr = protein.numberOfSelectedResidues(); nr > 20000) {
        WARN_LOG("Protein is too large (number of residues:{} > 20.000). Exiting.",nr);
        exit(EXIT_FAILURE);
    }

    //do the membrane region determination and annotation
    if (bool r = args.getValueAsBool("r"); r) {
        bool pf = args.getValueAsBool("pf");
        auto filter = Tmdet::Utils::Filter(protein,pf);
        auto dssp = Tmdet::Utils::Dssp(protein);
        auto mydssp = Tmdet::Utils::MyDssp(protein);
        auto ssVec = Tmdet::Utils::SecStrVec(protein);
        Tmdet::Utils::NeighBors::store(protein);

        if (bool fr = args.getValueAsBool("fr"); fr) {
            protein.forceSingleMembrane = true;
            auto fragmenter = Tmdet::Engine::Fragmenter(protein,args);
        }
        else {
            auto organizer = Tmdet::Engine::Organizer(protein, args);
        }
        protein.version = Tmdet::version();
        protein.date = Tmdet::System::Date::get();

        //write xml output if required
        if (xmlOutputPath != "" && !na) {
            xml.write(xmlOutputPath, protein, args);
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
        pymol.show(pdbOutputPath);
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
