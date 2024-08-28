#ifndef __TMDET_SERVICES_CHEMICALCOMPONENTDIRECTORY__
#define __TMDET_SERVICES_CHEMICALCOMPONENTDIRECTORY__

#include <string>
#include <Types/Residue.hpp>

namespace Tmdet::Services::ChemicalComponentDirectoryService {

    extern void build();
    extern bool isBuilt();
    Types::Residue getComponentAsResidue(const std::string& threeLetterCode);
}

#endif
