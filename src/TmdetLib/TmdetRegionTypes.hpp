#ifndef __UNITMP_TMDETLIB_TMDET_REGION_TYPES__
#define __UNITMP_TMDETLIB_TMDET_REGION_TYPES__

#include <unordered_map>

namespace UniTmp::TmdetLib {

    struct TmdetRegionType {
        char *name;
        char code;
        char *description;
    };

    namespace TmdetRegionTypes {
        const TmdetRegionType TM {
            "transmembrane", 
            'M',
            "Region crossing the membrane"
        };
        const TmdetRegionType SIDE1 {
            "side1",
            '1',
            "Side one of the membrane"
        };
        const TmdetRegionType SIDE2 {
            "side2",
            '2',
            "Other side of the membrane"
        };
        const TmdetRegionType LOOP {
            "loop",
            'L',
            "Re-entrant membrane loop"
        };
        const TmdetRegionType IFH {
            "interfacial helix",
            'F',
            "Interfacial helix on membrane surface"
        };
        const TmdetRegionType INSIDE {
            "inside",
            'I',
            "Beta barrel inside element"
        };
        const TmdetRegionType CBEG {
            "chain begin",
            'N',
            "N-terminal chain segment"
        };
        const TmdetRegionType CEND {
            "chain end",
            'C',
            "C-terminal chain segment"
        };
        const TmdetRegionType UNK {
            "unkown",
            'U',
            "Unknow localization"
        };

        const std::unordered_map<const char*, TmdetRegionType> all {
            {"transmembrane", TM},
            {"side1", SIDE1},
            {"side2", SIDE2},
            {"loop", LOOP},
            {"interfacial helix", IFH},
            {"inside", INSIDE},
            {"chain begin", CBEG},
            {"chain end", CEND},
            {"unknown", UNK},
        };

        /*const std::unordered_map<const char, TmdetRegionType> allByCode {
            {'M', TM},
            {'1', SIDE1},
            {'2', SIDE2},
            {'L', LOOP},
            {'F', IFH},
            {'I', INSIDE},
            {'N', CBEG},
            {'C', CEND},
            {'U', UNK},
        };*/
    }

}

#endif
