#ifndef __TMDET_SERVICES_CONFIGURATION__
#define __TMDET_SERVICES_CONFIGURATION__

#include <string>

namespace Tmdet::Services::ConfigurationService {

    // These are default values and can be overwrite shell environment variables.
    // Default values are set in cpp file - see them in parantheses below.
    // These can be overwritten by corresponding environment variables.
    //  - TMDET_DIRECTORY               (/usr/local/share/tmdet)
    //  - CHEMICAL_COMPONENT_DIRECTORY  (/usr/local/share/tmdet/data/ccd)
    //  - FRAGMENT_CIF_EXEC             (/usr/local/share/tmdet/fragment_cif)
    //  - PDB_DATA_DIRECTORY            (/zfs/databases/UniTmp/PDB/data/structures/divided/updated_mmcif)
    //  - UNITMP_SCHEMA                 (https://)
    //  - UNITMP_DOMAIN                 (unitmp.org)

    extern std::string TmdetDirectory;
    extern std::string ChemicalComponentDirectory;
    extern std::string ChemicalComponentFile;
    extern std::string ChemicalComponentDirectoryUrl;
    extern std::string UniTmpSchema;
    extern std::string UniTmpDomain;

    extern std::string FragmentCifExec;
    extern std::string PdbDataDirectory;

    // set by init() call
    extern std::string AppName;

    extern void init();
    extern std::string getEnv(std::string key);
}

#endif
