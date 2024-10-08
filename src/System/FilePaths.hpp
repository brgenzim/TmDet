#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <System/Config.hpp>


namespace Tmdet::System {

    struct FilePaths {

        static bool fileExists (const std::string& name) {
            struct stat buffer;   
            return (stat (name.c_str(), &buffer) == 0); 
        }

        static std::string xml(const std::string& code) {
            return environment.get("TMDET_DATA_ROOT",DEFAULT_TMDET_DATA_ROOT)
                    + "/" + code.substr(1,2)
                    + "/" + code 
                    + ".xml";
        }

        static std::string cif(const std::string& code) {
            return environment.get("PDB_CIF_DIR",DEFAULT_PDB_CIF_DIR)
                    + "/" + code.substr(1,2)
                    + "/" + code 
                    + environment.get("PDB_CIF_EXT",DEFAULT_PDB_CIF_EXT);
        }

        static std::string ent(const std::string& code, const std::string& ext) {
            return environment.get("PDB_ENT_DIR",DEFAULT_PDB_ENT_DIR)
                    + "/" + code.substr(1,2)
                    + "/pdb" + code 
                    + ".ent.gz";
        }

        static std::string pdbOut(const std::string& code) {
            return environment.get("TMDET_DATA_ROOT",DEFAULT_TMDET_DATA_ROOT)
                    + "/" + code.substr(1,2)
                    + "/" + code 
                    + "_updated_tr.cif.gz";
        }

        static std::string cache(const std::string& hash) {
            return environment.get("TMDET_CACHE_ROOT",DEFAULT_TMDET_CACHE_ROOT)
                    + "/" + hash.substr(0,2)
                    + "/" + hash.substr(2,2)
                    + "/" + hash.substr(4,2);
        }
    };
}
