// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <string>
#include <VOs/Protein.hpp>

namespace Tmdet::Utils {

    class CifUtil {
        public:
            static constexpr std::string TMDET_MEMBRANE_ASYM_ID = "TM_";

            /**
             * @brief Util function to get last part of a tag name.
             */
            static std::string getSuffix(const std::string& tag);

            /**
             * @brief Update _atom_site loop and add membrane-representation atoms.
             */
            static void prepareDocumentBlock(Tmdet::VOs::Protein& protein);

            /**
             * @brief Update the first data block name, if _entry.id exists and differs from it.
             */
            static void updateDataBlockNameIfNeeded(gemmi::cif::Document& document);

    };
}
