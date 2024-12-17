// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#include <string>
#include <vector>
#include <any>
#include <Config.hpp>
#include <Engine/PlaneSideDetector.hpp>
#include <System/Logger.hpp>
#include <Types/Region.hpp>
#include <VOs/Membrane.hpp>

namespace Tmdet::Engine {

    double PlaneSideDetector::getDistance(const gemmi::Vec3& vec) {
        return vec.z;
    }

    void PlaneSideDetector::setZs(const std::vector<Tmdet::VOs::Membrane>& membranes) {
        if (membranes.size() == 2) {
            if (membranes[1].origo < 0) {
                z1 = membranes[0].halfThickness;
                z2 = -membranes[0].halfThickness;
                z3 = membranes[1].origo+membranes[1].halfThickness;
                z4 = membranes[1].origo-membranes[1].halfThickness;
                o1 = membranes[0].origo;
                o2 = membranes[1].origo;
            }
            else {
                z1 = membranes[1].origo+membranes[1].halfThickness;
                z2 = membranes[1].origo-membranes[1].halfThickness;
                z3 = membranes[0].halfThickness;
                z4 = -membranes[0].halfThickness;
                o1 = membranes[1].origo;
                o2 = membranes[0].origo;
            }
        }
        else {
            z1 = membranes[0].halfThickness;
            z4 = -membranes[0].halfThickness;
            o1 = o2 = membranes[0].origo;
        }
    }
}
