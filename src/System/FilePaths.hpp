// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <format>
#include <filesystem>
#include <sys/stat.h>
#include <unistd.h>
#include <Config.hpp>

/**
 * @brief namespace for tmdet system
 *
 * @namespace Tmdet
 * @namespace System
 */
namespace Tmdet::System {

    struct FilePaths {

        /**
         * @brief check if a file exists
         * 
         * @param name 
         * @return true 
         * @return false 
         */
        static bool fileExists (const std::string& name) {
            struct stat buffer;   
            return (stat (name.c_str(), &buffer) == 0); 
        }

        /**
         * @brief generate the path for xml file given the pdb code
         * 
         * @param code 
         * @param createDir 
         * @return std::string 
         */
        static std::string xml(const std::string& code, const bool createDir = false) {
            std::string path = std::format("{}/{}",
                    environment.get("PDBTM_DATA_ROOT",DEFAULT_TMDET_DATA_ROOT),
                    code.substr(1,2));
            if (createDir) {
                std::filesystem::create_directories(path);
            }
            return std::format("{}/{}.xml",path,code);
        }

        /**
         * @brief generate the path for a cif file given the pdb code
         * 
         * @param code 
         * @return std::string 
         */
        static std::string cif(const std::string& code, const int assemblyId = 0) {
            return std::format("{}/{}/{}{}",
                    environment.get("PDB_CIF_DIR",DEFAULT_PDB_CIF_DIR),
                    code.substr(1,2),code,
                    (assemblyId>0?std::format("-assembly{}.cif.gz",assemblyId):".cif.gz")
            );
        }

        /**
         * @brief generate the path for an ent file given the pdb code
         * 
         * @param code 
         * @param ext 
         * @return std::string 
         */
        static std::string ent(const std::string& code, const std::string& ext) {
            return std::format("{}/{}/pdb{}.ent.gz",
                    environment.get("PDB_ENT_DIR",DEFAULT_PDB_ENT_DIR),
                    code.substr(1,2),code);
        }

        /**
         * @brief generate the path for the transformed pdb file given the pdb code
         * 
         * @param code 
         * @return std::string 
         */
        static std::string pdbOut(const std::string& code) {
            return std::format("{}/{}/{}_updated_tr.cif.gz",
                    environment.get("PDBTM_DATA_ROOT",DEFAULT_TMDET_DATA_ROOT),
                    code.substr(1,2),code);
        }

        /**
         * @brief generate the cache directory name given the hash of the protein
         * 
         * @param hash 
         * @return std::string 
         */
        static std::string cache(const std::string& hash) {
            return std::format("{}/{}/{}/{}",
                    environment.get("TMDET_CACHE_ROOT",DEFAULT_TMDET_CACHE_ROOT),
                    hash.substr(0,2), hash.substr(2,2), hash.substr(4,2));
        }

        /**
         * @brief generate the temporary directory name given the hash 
         * 
         * @param hash 
         * @return std::string 
         */
        static std::string temp(const std::string& hash) {
            return std::format("{}/{}/{}/{}",
                    environment.get("TMDET_TEMP_ROOT",DEFAULT_TMDET_TEMP_ROOT),
                    hash.substr(0,2), hash.substr(2,2), hash.substr(4,2));
        }

        /**
         * @brief check if a file contains cif document
         *
         * @param path
         * @return bool
         */
        static bool isCif(const std::string& path) {
            return (path.substr(path.size()-4) == ".cif" || path.substr(path.size()-7) == ".cif.gz");
        }
    };
}
