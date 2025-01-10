// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#include <gemmi/model.hpp>
#include <Engine/PlaneOptimizer.hpp>
#include <Engine/Rotator.hpp>
#include <Helpers/Vector.hpp>

namespace Tmdet::Engine {

    double PlaneOptimizer::distance(gemmi::Vec3& vec) {
        return normal.x * (vec.x - massCentre.x)
                + normal.y * (vec.y - massCentre.y)
                + normal.z * (vec.z - massCentre.z);
    }

    double PlaneOptimizer::getAngle(Tmdet::VOs::SecStrVec& vector) {
        return std::abs(Tmdet::Helpers::Vector::cosAngle(normal,vector.end - vector.begin));
    }

    void PlaneOptimizer::testMembraneNormal() {
        testMembraneNormalOne();
    }

    void PlaneOptimizer::testMembraneNormalFinal() {
        testMembraneNormalOne();
    }

    void PlaneOptimizer::setBestOrigo(double minz, double maxz) {
        //
    }

    void PlaneOptimizer::setMembraneOrigo(Tmdet::VOs::Membrane& membrane, double minz, double maxz) {
        double o = (minz+maxz) / 2 + minZ;
        if (protein.membranes.empty()) {
            massCentre += o * bestNormal;
            membrane.origo = 0;
            lastO = o;
        }
        else {
            membrane.origo = o - lastO;
        }
    }
}
