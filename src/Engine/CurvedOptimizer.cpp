// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#include <gemmi/model.hpp>
#include <Engine/CurvedOptimizer.hpp>
#include <Engine/Rotator.hpp>
#include <Helpers/Vector.hpp>

namespace Tmdet::Engine {

    const std::vector<double> origos = {
        //5, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 22, 25, 30, 40, 50, 60, 70,
        50, 55, 60, 65, 70, 75, 80, 85, 90, 100, 120, 140, 150, 160, 180, 200, 
        225, 250, 275, 300, 325, 350, 375, 400, 450, 500, 750, 1000, 2000
    };

    double CurvedOptimizer::distance(gemmi::Vec3& vec) {
        return origoVec3.dist(vec);
    }

    double CurvedOptimizer::getAngle(Tmdet::VOs::SecStrVec& vector) {
        gemmi::Vec3 vec = (vector.begin + vector.end) / 2;
        vec -= origoVec3;
        return std::abs(Tmdet::Helpers::Vector::cosAngle(vec,vector.end - vector.begin));
    }

    void CurvedOptimizer::testMembraneNormal() {
        for (auto& o: origos) {
            origo = -1.0 * o;
            origoVec3 = massCentre + origo * normal;
            DEBUG_LOG("Test Radius: {}",o);
            testMembraneNormalOne(); 
        }
    }

    void CurvedOptimizer::testMembraneNormalFinal() {
        origo = bestOrigo;
        origoVec3 = massCentre + bestOrigo * normal;
        DEBUG_LOG("Test Final Radius: {}",bestOrigo);
        testMembraneNormalOne(); 
    }

    void CurvedOptimizer::setBestOrigo(double minz, double maxz) {
        bestOrigo = origo;
        double mcz = distance(massCentre) - minZ;
        bestSphereRadius = -origo - mcz + (maxz + minz) / 2;
        DEBUG_LOG("setBestOrigo: origo:{:5.2f} mcz:{:5.2f} minZ:{:5.2f} maxZ:{:5.2f} minz:{:5.2f} maxz:{:5.2f} R:{:5.2f}",
            origo,mcz,minZ,maxZ,minz,maxz,bestSphereRadius);
    }

    void CurvedOptimizer::setMembraneOrigo(Tmdet::VOs::Membrane& membrane, double minz, double maxz) {
        membrane.type = Tmdet::Types::MembraneType::CURVED;
        membrane.sphereRadius = bestSphereRadius;
        membrane.origo = bestOrigo;
        DEBUG_LOG("setMembraneOrigo: origo:{} r:{}",membrane.origo,membrane.sphereRadius);
    }
    
}
