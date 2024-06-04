#include <iostream>
#include <string>
#include <filesystem>
#include <stdexcept>
#include <cstdlib>
#include <Services/ConfigurationService.hpp>
#include <Services/ChemicalComponentDirectoryService.hpp>

namespace fs = std::filesystem;

namespace Tmdet::Services::ChemicalComponentDirectoryService {

    static void download();
    static void split();

    void build() {
        download();
        split();
    }

    bool isBuilt() {
        std::string lastFile(ConfigurationService::getValue(ConfigurationService::Keys::CHEMICAL_COMPONENT_DIRECTORY)
            + "/Z/Z/ZZZ.cif");

        return fs::exists(fs::path(lastFile));
    }

    std::string getComponentPath(std::string componentCode) {

    }

    void download() {
        std::string cmd("bash ");
        cmd += ConfigurationService::getValue(ConfigurationService::Keys::CHEMICAL_COMPONENT_DOWNLOAD_SCRIPT)
            + " " + ConfigurationService::getValue(ConfigurationService::Keys::CHEMICAL_COMPONENT_DIRECTORY);
        int exitCode = std::system(cmd.c_str());
        if (exitCode != 0) {
            std::string message(ConfigurationService::AppName);
            message += ": command failed: '" + cmd + "'";
            throw std::runtime_error(message);
        }
    }

    void split() {
        std::string cmd(ConfigurationService::getValue(ConfigurationService::Keys::FRAGMENT_CIF_EXEC));
        cmd += " -i " + ConfigurationService::getValue(ConfigurationService::Keys::CHEMICAL_COMPONENT_FILE)
            + " -d " + ConfigurationService::getValue(ConfigurationService::Keys::CHEMICAL_COMPONENT_DIRECTORY)
            + " -s > /dev/null 2>&1";
        int exitCode = std::system(cmd.c_str());
        if (exitCode != 0) {
            std::string message(ConfigurationService::AppName);
            message += ": command failed: '" + cmd + "'";
            throw std::runtime_error(message);
        }
    }

}
