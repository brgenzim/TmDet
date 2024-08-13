#ifndef __TMDET_SERVICES_CONFIGURATION__
#define __TMDET_SERVICES_CONFIGURATION__

#include <string>

namespace Tmdet::Services::ConfigurationService {

    // These are default values and can be overwrite shell environment variables.
    // Default values are set in cpp file.
    extern std::string TmdetDirectory;
    extern std::string ChemicalComponentDirectory;
    extern std::string ChemicalComponentFile;
    // TODO: This should be replaced by a CurlWrapperService function
    extern std::string ChemicalComponentDownloadScript;
    extern std::string FragmentCifExec;
    extern std::string PdbDataDirectory;

    // set by init() call
    extern std::string AppName;

    extern void init();
    extern std::string getEnv(std::string key);
}

#endif
