// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <string>
#include <Types/Residue.hpp>
#include <gemmi/cif.hpp>

/**
 * @brief namespace for tmdet services
 *
 * @namespace Tmdet
 * @namespace Services
 */
namespace Tmdet::Services {

    /**
     * @brief Service for acquiring CCD
     */
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
