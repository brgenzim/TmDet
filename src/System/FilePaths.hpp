#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <format>
#include <filesystem>
#include <sys/stat.h>
#include <unistd.h>
#include <Config.hpp>

namespace Tmdet::System {

    struct FilePaths {

        static bool fileExists (const std::string& name) {
            struct stat buffer;   
            return (stat (name.c_str(), &buffer) == 0); 
        }

        static std::string xml(const std::string& code, const bool createDir = false) {
            std::string path = std::format("{}/{}",
                    environment.get("TMDET_DATA_ROOT",DEFAULT_TMDET_DATA_ROOT),
                    code.substr(1,2));
            if (createDir) {
                std::filesystem::create_directories(path);
            }
            return std::format("{}/{}.xml",path,code);
        }

        static std::string cif(const std::string& code) {
            return std::format("{}/{}/{}{}",
                    environment.get("PDB_CIF_DIR",DEFAULT_PDB_CIF_DIR),
                    code.substr(1,2),code,
                    environment.get("PDB_CIF_EXT",DEFAULT_PDB_CIF_EXT));
        }

        static std::string ent(const std::string& code, const std::string& ext) {
            return std::format("{}/{}/pdb{}.ent.gz",
                    environment.get("PDB_ENT_DIR",DEFAULT_PDB_ENT_DIR),
                    code.substr(1,2),code);
        }

        static std::string pdbOut(const std::string& code) {
            return std::format("{}/{}/{}_updated_tr.cif.gz",
                    environment.get("TMDET_DATA_ROOT",DEFAULT_TMDET_DATA_ROOT),
                    code.substr(1,2),code);
        }

        static std::string cache(const std::string& hash) {
            return std::format("{}/{}/{}/{}",
                    environment.get("TMDET_CACHE_ROOT",DEFAULT_TMDET_CACHE_ROOT),
                    hash.substr(0,2), hash.substr(2,2), hash.substr(4,2));
        }
    };
}
