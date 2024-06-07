#include <iostream>
#include <filesystem>
#include <map>
#include <string>
#include <stdexcept>
#include <cstdlib>
#include <unistd.h>
#include <limits.h>
#include <Services/ConfigurationService.hpp>

namespace Tmdet::Services::ConfigurationService {

    std::map<std::string, std::string> Config;
    std::string AppName = "";

    void init() {

        // get exec paths
        char path[PATH_MAX];
        auto size = readlink("/proc/self/exe", path, PATH_MAX);
        if (size == -1) {
            throw std::runtime_error("readlink() call failed");
        }
        std::string appExec = std::string(path, size);
        std::filesystem::path fsAppExecPath(appExec);
        AppName = Config[Keys::APP_NAME] = fsAppExecPath.filename();
        Config[Keys::APP_EXEC] = appExec;
        Config[Keys::FRAGMENT_CIF_EXEC] = fsAppExecPath.replace_filename("fragment_cif");

        // get base directory
        std::string baseDir = fsAppExecPath.parent_path().string();
        if (baseDir.ends_with("build/src")) {
            baseDir += "/../..";
        }
        // set CCD and TMDET directories
        Config[Keys::TMDET_DIRECTORY] = baseDir;
        Config[Keys::CHEMICAL_COMPONENT_DIRECTORY] = baseDir + "/data/ccd";
        Config[Keys::CHEMICAL_COMPONENT_FILE] = Config[Keys::CHEMICAL_COMPONENT_DIRECTORY] + "/components.cif.gz";
        Config[Keys::CHEMICAL_COMPONENT_DOWNLOAD_SCRIPT] = Config[Keys::TMDET_DIRECTORY]
            + "/scripts/get-chemical-component-directory.sh";

        // set RCSBROOT
        // std::string value = getEnv("RCSBROOT");
        // if (value == "") {
        //     std::cerr << AppName << ": WARNING: RCSBROOT environment variable not set."
        //         << " Maxit-related PromotifService will fail"
        //         << std::endl;
        // }
        // Config[Keys::RCSBROOT_DIRECTORY] = value;

    }

    bool hasKey(std::string& key) {
        return Config.count(key) > 0;
    }

    std::string getValue(std::string key) {
        if (Config.count(key) == 0) {
            std::string app = Config[Keys::APP_NAME];
            throw std::runtime_error(app + ": key '" + key + "' not found in configuration");
        }
        return Config[key];
    }

    std::string getEnv(std::string key) {
        char* value = getenv(key.c_str());
        std::string result;
        if (value == nullptr) {
            return result;
        }
        result = std::string(value);
        return result;
    }

    void dump() {
        for (auto& [key, value] : Config) {
            std::cout << key << ": " << value << std::endl;
        }
    }
}
