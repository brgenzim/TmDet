#pragma once

#include <vector>
#include <string>
#include <any>
#include <gemmi/model.hpp>
#include <ValueObjects/Protein.hpp>

// #define __SYM_DBG 1 // to enable debug messages of this feature

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
            std::vector<std::vector<_symmetryData>> CheckSymmetry(Tmdet::ValueObjects::Protein &proteinVO);
    };
}
