// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#include <string>
#include <vector>
#include <any>
#include <Config.hpp>
#include <Engine/CurvedSideDetector.hpp>
#include <System/Logger.hpp>
#include <Types/Region.hpp>
#include <VOs/Membrane.hpp>

namespace Tmdet::Engine {

    double CurvedSideDetector::getDistance(const gemmi::Vec3& vec) {
        return vec.dist(gemmi::Vec3(0,0,protein.membranes[0].origo));
    }

    void CurvedSideDetector::setZs(const std::vector<Tmdet::VOs::Membrane>& membranes) {
        if (membranes.size() == 2) {
            int i1;
            int i2;
            if (membranes[1].origo + membranes[1].sphereRadius < membranes[0].origo + membranes[0].sphereRadius) {
                i1 = 0; i2 = 1;
            }
            else {
                i1 = 1; i2 = 0;
            }
            z1 = membranes[i1].origo + membranes[i1].sphereRadius + membranes[i1].halfThickness;
            z2 = membranes[i1].origo + membranes[i1].sphereRadius - membranes[i1].halfThickness;
            z3 = membranes[i2].origo + membranes[i2].sphereRadius + membranes[i2].halfThickness;
            z4 = membranes[i2].origo + membranes[i2].sphereRadius - membranes[i2].halfThickness;
            o1 = (z1 + z2) / 2;
            o2 = (z3 + z4) / 2;
        }
        else {
            z1 = membranes[0].sphereRadius + membranes[0].halfThickness;
            z4 = membranes[0].sphereRadius - membranes[0].halfThickness;
            o1 = o2 = membranes[0].sphereRadius;
        }
    }

}
