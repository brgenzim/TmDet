#ifndef __TMDET_TYPES_PROTEIN__
#define __TMDET_TYPES_PROTEIN__

#include <unordered_map>
#include <string>

using namespace std;

namespace Tmdet::Types {

    struct Protein {
        string name;
        string description;
        bool globular;
        bool membrane;
        bool alpha;
    };
    
    const unordered_map<string, Protein> Proteins = {
        { "Tm_Alpha", { 
                "Tm_Alpha",
                "Alpha helical transmembrane protein",
                false, true, true
            }
        },
        { "Tm_Beta", { 
                "Tm_Beta",
                "Beta barrel transmembrane protein",
                false, true, false
            }
        },
        { "Soluble", {
                "Soluble",
                "Globular, water soluble, not transmembrane protein",
                true, false, false
            }
        },
        { "Ca_Tm", { 
                "Ca_Tm",
                "Transmembrane protein in low resolution, containing mostly Calpha or backbone atoms",
                false, true, true 
            }
        },
        { "Ca_Globular", { 
                "Ca_Globular",
                "Globular, not transmembrane protein in low resolution, containing mostly Calpha or backbone atoms",
                true, false, false 
            }
        },
        { "Nucleotide", { 
                "Nucleotide",
                "Nucleotide polimer chain",
                false, false, false 
            }
        },
        { "Tm_Mixed", { 
                "Tm_Mixed",
                "Transmembrane protein containing both alpha helical membrane regions and beta barrel",
                false, true, true 
            }
        }
    };
}

#endif
