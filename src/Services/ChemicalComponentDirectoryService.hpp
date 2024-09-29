#ifndef __TMDET_SERVICES_CHEMICAL_COMPONENT_DIRECTORY__
#define __TMDET_SERVICES_CHEMICAL_COMPONENT_DIRECTORY__

#include <string>
#include <Types/Residue.hpp>
#include <gemmi/cif.hpp>

namespace Tmdet::Services {

    class ChemicalComponentDirectoryService {
        private:
            static std::string createDir(const std::string& destDir, const std::string& cifName);
            static void writeCif(const char *cifPath, const gemmi::cif::Block& block);

        public:
            static void build();
            static void fetch();
            static bool isBuilt();
            static Types::Residue getComponentAsResidue(const std::string& threeLetterCode);

    };
}

#endif
