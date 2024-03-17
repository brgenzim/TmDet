#ifndef __TMDET_VALUE_OBJECTS_MEMBRANE__
#define __TMDET_VALUE_OBJECTS_MEMBRANE__

#include <string>
#include <vector>
#include <gemmi/unitcell.hpp>
#include <gemmi/math.hpp>
#include <Types/Membrane.hpp>
#include <ValueObjects/TMatrix.hpp>

using namespace std;
using namespace gemmi;

namespace Tmdet::ValueObjects {

    struct Membrane {
        TMatrix tmatrix;
        double h;
        double curver;
        double sizer;
        Tmdet::Types::Membrane type;
    };

}

#endif