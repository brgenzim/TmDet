#ifndef __UNITMP_TMDETLIB_TMDET_MEMBRANE_TYPES__
#define __UNITMP_TMDETLIB_TMDET_MEMBRANE_TYPES__

#include <unordered_map>

namespace UniTmp::TmdetLib {

    struct TmdetMembraneType {
        char *name;
        char *description;
    };

    namespace TmdetMembraneTypes {
        const TmdetMembraneType PLAIN {"Plain", "Simple plain membrane used most of membrane proteins"};
        const TmdetMembraneType CURVED {"Curved", "Curved membrane represented by a radius"};
        const TmdetMembraneType DOUBLE {"Double", "Double membrane for bacterial proteins"};

        const std::unordered_map<const char*, TmdetMembraneType> all {
            {"Plain", PLAIN},
            {"Curved", CURVED},
            {"Double", DOUBLE}
        };
    };

}

#endif
