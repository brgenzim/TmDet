#ifndef __TMDET_VALUE_OBJECTS_TMATRIX__
#define __TMDET_VALUE_OBJECTS_TMATRIX__

#include <string>
#include <vector>
#include <gemmi/unitcell.hpp>
#include <gemmi/math.hpp>

using namespace std;
using namespace gemmi;

namespace Tmdet::ValueObjects {

    struct TMatrix {
        gemmi::Mat33 rot;
        gemmi::Vec3 trans;
    };

}

#endif