#ifndef __TMDET_UTILS_SYMMETRY__
#define __TMDET_UTILS_SYMMETRY__

#include <vector>
#include <string>
#include <any>
#include <gemmi/model.hpp>
#include <ValueObjects/TmdetStruct.hpp>

// #define __SYM_DBG 1 // to enable debug messages of this feature

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
        public:
        std::vector<std::vector<_symmetryData>> CheckSymmetry(Tmdet::ValueObjects::TmdetStruct &tmdetVO);
    };
}

#endif