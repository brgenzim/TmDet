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
            static constexpr std::string ENTRY_PREFIX = "JOB_";

            /**
             * @brief Util function to get last part of a tag name.
             */
            static std::string getSuffix(const std::string& tag);

            /**
             * @brief Util function to get first part of a tag name.
             */
            static std::string getPrefix(const std::string& tag);


            /**
             * @brief Update _atom_site loop and add membrane-representation atoms.
             */
            static void prepareDocumentBlock(Tmdet::VOs::Protein& protein);

            /**
             * @brief Update the first data block name, if _entry.id exists and differs from it.
             */
            static void updateDataBlockNameIfNeeded(gemmi::cif::Document& document);

            /**
             * @brief Set _entry.id, if it has a too long name. It is important for TmdetWeb's job management.
             */
            static void setEntryIdFromFilePath(gemmi::cif::Document& documemt, const std::string& filePath);

            static gemmi::Vec3 multiply(const gemmi::Mat33& mx, const gemmi::Vec3& vec);
            static gemmi::Vec3 multiply(const gemmi::Vec3& vec, const double scalar);
            static gemmi::Vec3 add(const gemmi::Vec3& v1, const gemmi::Vec3& v2);
            static gemmi::Mat33 transpose(const gemmi::Mat33& mx);
    };
}
