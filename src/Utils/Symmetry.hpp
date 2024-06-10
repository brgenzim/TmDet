#ifndef __TMDET_UTILS_SYMMETRY__
#define __TMDET_UTILS_SYMMETRY__

#include <array>
#include <string>
#include <any>
#include <gemmi/model.hpp>
#include <ValueObjects/TmdetStruct.hpp>

using namespace std;

namespace Tmdet::Utils {

    struct _symmetryData {
        int id;
        int cont;
        int cluster;
        gemmi::Vec3 origo;
        gemmi::Vec3 axis;
        double rotAngle;
    };

    class Symmetry {
        
    };
}

#endif