#pragma once

#include <string>
#include <VOs/Protein.hpp>

namespace Tmdet::Utils {

    class CifUtil {
        public:
            static constexpr std::string TMDET_MEMBRANE_ASYM_ID = "TM_";

            static std::string getSuffix(const std::string& tag);

            /**
             * Update _atom_site loop and add membrane-representation atoms
             */
            static void prepareDocumentBlock(Tmdet::VOs::Protein& protein);

    };
}
