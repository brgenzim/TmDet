#ifndef __UNITMP_TMDETLIB_TMDET_MEMBRANE__
#define __UNITMP_TMDETLIB_TMDET_MEMBRANE__

#include <string>
#include <vector>
#include <gemmi/unitcell.hpp>
#include <gemmi/math.hpp>
#include <TmdetMembraneTypes.hpp>

using namespace std;
using namespace gemmi;

namespace UniTmp::TmdetLib {

    struct _tmdetMembrane {
        Mat33 rot;
        Vec3 trans;
        double h;
        double curver;
        double sizer;
        TmdetMembraneType type;
    };

    class TmdetMembrane {
        public:
            vector<_tmdetMembrane> membranes;

            TmdetMembrane();
            ~TmdetMembrane();
            
    };
}

#endif