#ifndef __UNITMP_TMDETLIB_TMDET_PROTEIN_TYPES__
#define __UNITMP_TMDETLIB_TMDET_PROTEIN_TYPES__

#include <unordered_map>

namespace UniTmp::TmdetLib {

    struct TmdetProteinType {
        char *name;
        char *description;
        bool globular;
        bool membrane;
        bool alpha;
    };

    namespace TmdetProteinTypes {
        const TmdetProteinType ALPHA {
            "Tm_Alpha", 
            "Alpha helical transmembrane protein",
            false, true, true
        };
        const TmdetProteinType BETA {
            "Tm_Beta", 
            "Beta barrel transmembrane protein",
            false, true, false
        };
        const TmdetProteinType GLOB {
            "Soluble",
            "Globular, water soluble, not transmembrane protein",
            true, false, false
        };
        const TmdetProteinType CA_TM {
            "Ca_Tm", 
            "Transmembrane protein in low resolution, containing mostly Calpha or backbone atoms",
            false, true, true 
        };
        const TmdetProteinType CA_GLOB {
            "Ca_Globular", 
            "Globular, not transmembrane protein in low resolution, containing mostly Calpha or backbone atoms",
            true, false, false 
        };
        const TmdetProteinType NUCL {
            "Nucleotide", 
            "Nucleotide polimer chain",
            false, false, false 
        };
        const TmdetProteinType MIXED {
            "Tm_Mixed", 
            "Transmembrane protein containing both alpha helical membrane regions and beta barrel",
            false, true, true 
        };

        const std::unordered_map<const char*, TmdetProteinType> all {
            {"Tm_Alpha", ALPHA},
            {"Tm_Beta", BETA},
            {"Soluble", GLOB},
            {"Ca_Tm", CA_TM},
            {"Ca_Globular", CA_GLOB},
            {"Nucleotide", NUCL},
            {"Tm_Mixed", MIXED}
        };
    }

}

#endif
