#ifndef __TMDET_SERVICES_CONFIGURATION__
#define __TMDET_SERVICES_CONFIGURATION__

#include <map>
#include <string>

namespace Tmdet::Services::ConfigurationService {

    namespace Keys {
        const std::string APP_NAME = "APP_NAME";
        const std::string APP_EXEC = "APP_EXEC";
        const std::string FRAGMENT_CIF_EXEC = "FRAGMENT_CIF_EXEC";
        const std::string TMDET_DIRECTORY = "TMDET_DIRECTORY";
        const std::string CHEMICAL_COMPONENT_DIRECTORY = "CHEMICAL_COMPONENT_DIRECTORY";
        const std::string CHEMICAL_COMPONENT_FILE = "CHEMICAL_COMPONENT_FILE";
        const std::string CHEMICAL_COMPONENT_DOWNLOAD_SCRIPT = "CHEMICAL_COMPONENT_DOWNLOAD_SCRIPT";
        const std::string RCSBROOT_DIRECTORY = "RCSBROOT";
        const std::string PDB_DIRECTORY = "PDB_DIRECTORY";

    }

    extern std::string AppName;

    extern void init();
    extern std::string getValue(std::string key);
    extern std::string getEnv(std::string key);
    extern bool hasKey(std::string& key);
    extern void dump();
}

#endif
