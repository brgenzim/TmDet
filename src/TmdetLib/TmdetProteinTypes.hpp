#ifndef __UNITMP_TMDETLIB_TMDET_PROTEIN_TYPES__
#define __UNITMP_TMDETLIB_TMDET_PROTEIN_TYPES__

#include <unordered_map>

namespace UniTmp::TmdetLib {

    struct TmdetProteinType {
        char *description;
        bool globular;
        bool membrane;
        bool alpha;
    };

    const unordered_map<const char*, TmdetProteinType> TmdetProteinTypes = {
        { "Tm_Alpha", { 
                "Alpha helical transmembrane protein",
                false, true, true
            }
        },
        { "Tm_Beta", { 
                "Beta barrel transmembrane protein",
                false, true, false
            }
        },
        { "Soluble", {
                "Globular, water soluble, not transmembrane protein",
                true, false, false
            }
        },
        { "Ca_Tm", { 
                "Transmembrane protein in low resolution, containing mostly Calpha or backbone atoms",
                false, true, true 
            }
        },
        { "Ca_Globular", { 
                "Globular, not transmembrane protein in low resolution, containing mostly Calpha or backbone atoms",
                true, false, false 
            }
        },
        { "Nucleotide", { 
                "Nucleotide polimer chain",
                false, false, false 
            }
        },
        { "Tm_Mixed", { 
                "Transmembrane protein containing both alpha helical membrane regions and beta barrel",
                false, true, true 
            }
        }
    }

}

#endif
