#ifndef __TMDET_VALUE_OBJECTS_MEMBRANE__
#define __TMDET_VALUE_OBJECTS_MEMBRANE__

#include <gemmi/unitcell.hpp>
#include <gemmi/math.hpp>
#include <Types/Membrane.hpp>
#include <ValueObjects/TMatrix.hpp>

namespace Tmdet::ValueObjects {

    struct Membrane {
        TMatrix tmatrix;
        gemmi::Vec3 origo;
        gemmi::Vec3 normal;
        double h;
        double curver;
        double sizer;
        Tmdet::Types::Membrane type;
    };

}

#endif