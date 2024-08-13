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

    // These are default values and can be overwrite shell environment variables
    std::string TmdetDirectory{"/usr/local/share/tmdet"};
    std::string ChemicalComponentDirectory{TmdetDirectory + "/data/ccd"};
    std::string ChemicalComponentFile{ChemicalComponentDirectory + "/components.cif.gz"};
    std::string ChemicalComponentDownloadScript{TmdetDirectory + "/scripts/get-chemical-component-directory.sh"};
    std::string ChemicalComponentDirectoryUrl{"https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz"};
    std::string FragmentCifExec{TmdetDirectory + "/fragment_cif"};
    std::string PdbDataDirectory{"/zfs/databases/UniTmp/PDB/data/structures/divided/updated_mmcif/"};

    // set by init() call
    std::string AppName("not-set");


    namespace impl::Keys {
        const std::string TMDET_DIRECTORY = "TMDET_DIRECTORY";
        const std::string FRAGMENT_CIF_EXEC = "TMDET_FRAGMENT_CIF_EXEC";
        const std::string CHEMICAL_COMPONENT_DIRECTORY = "TMDET_CHEMICAL_COMPONENT_DIRECTORY";
        const std::string CHEMICAL_COMPONENT_DOWNLOAD_SCRIPT = "TMDET_CHEMICAL_COMPONENT_DOWNLOAD_SCRIPT";
        // used in test integration tests
        const std::string PDB_DATA_DIRECTORY = "TMDET_PDB_DATA_DIRECTORY";
    }

    void init() {

        // get exec paths
        char path[PATH_MAX];
        auto size = readlink("/proc/self/exe", path, PATH_MAX);
        if (size == -1) {
            throw std::runtime_error("readlink() call failed");
        }
        std::string appExec = std::string(path, size);
        std::filesystem::path fsAppExecPath(appExec);
        AppName = fsAppExecPath.filename();

        // Update default values if there is explicit value in environment variable
        std::string value;
        if ((value = getEnv(impl::Keys::TMDET_DIRECTORY)) != "") {
            TmdetDirectory = value;
            ChemicalComponentDirectory = TmdetDirectory + "/data/ccd";
            ChemicalComponentFile = ChemicalComponentDirectory + "/components.cif.gz";
            ChemicalComponentDownloadScript = TmdetDirectory + "/scripts/get-chemical-component-directory.sh";
        }
        if ((value = getEnv(impl::Keys::CHEMICAL_COMPONENT_DIRECTORY)) != "") {
            ChemicalComponentDirectory = value;
            ChemicalComponentFile = ChemicalComponentDirectory + "/components.cif.gz";
        }
        if ((value = getEnv(impl::Keys::FRAGMENT_CIF_EXEC)) != "") {
            FragmentCifExec = value;
        }
        if ((value = getEnv(impl::Keys::PDB_DATA_DIRECTORY)) != "") {
            PdbDataDirectory = value;
        }

        // Check required directories/files
        {
            std::filesystem::path path(TmdetDirectory);
            if (!std::filesystem::exists(path)) {
                std::cerr << TmdetDirectory << " does not exist. Update your "
                    << impl::Keys::TMDET_DIRECTORY << " variable." << std::endl;
                std::exit(1);
            }

            path = FragmentCifExec;
            if (!std::filesystem::exists(path)) {
                std::cerr << FragmentCifExec << " does not exist. Update your "
                    << impl::Keys::TMDET_DIRECTORY << " variable." << std::endl;
                std::exit(1);
            }

            path = PdbDataDirectory;
            if (!std::filesystem::exists(path) && AppName != "tmdet") {
                std::cerr << PdbDataDirectory << " does not exist. Update your "
                    << impl::Keys::PDB_DATA_DIRECTORY << " variable." << std::endl;
                std::exit(1);
            }

        }
    }

    std::string getEnv(std::string key) {
        char* value = getenv(key.c_str());
        if (value == nullptr) {
            return "";
        }
        return std::string(value);
    }

}
