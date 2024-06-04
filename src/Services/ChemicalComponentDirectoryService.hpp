#ifndef __TMDET_SERVICES_CHEMICALCOMPONENTDIRECTORY__
#define __TMDET_SERVICES_CHEMICALCOMPONENTDIRECTORY__

#include <map>
#include <string>
#include <Types/Residue.hpp>

namespace Tmdet::Services::ChemicalComponentDirectoryService {

    extern void build();
    extern bool isBuilt();
    extern std::string getComponentPath(std::string componentCode);
}

#endif
