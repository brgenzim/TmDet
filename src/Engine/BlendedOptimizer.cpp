#include <gemmi/model.hpp>
#include <Engine/BlendedOptimizer.hpp>
#include <Engine/Rotator.hpp>

namespace Tmdet::Engine {

    const std::vector<double> origos = {
        10, 11, 12, 13, 14, 15, 17, 20, 25, 30, 40, 50, 60, 80, 100, 120, 140, 150, 160, 200, 250, 300, 350, 400
    };

    double BlendedOptimizer::distance(gemmi::Vec3& vec) {
        return origoVec3.dist(vec);
    }

    void BlendedOptimizer::testMembraneNormal() {
        for (auto& o: origos) {
            origo = -1.0 * o;
            origoVec3 = massCentre + origo * normal;
            testMembraneNormalOne();  
        }
    }

    void BlendedOptimizer::setBestOrigo(double minz, double maxz) {
        bestOrigo = origo;
        double mcz = distance(massCentre) - minZ;
        bestSphereRadius = -origo - mcz + (maxz + minz) / 2;
        DEBUG_LOG("setBestOrigo: origo:{:5.2f} mcz:{:5.2f} minZ:{:5.2f} maxZ:{:5.2f} minz:{:5.2f} maxz:{:5.2f} R:{:5.2f}",
            origo,mcz,minZ,maxZ,minz,maxz,bestSphereRadius);
    }

    void BlendedOptimizer::setMembraneOrigo(Tmdet::VOs::Membrane& membrane, double minz, double maxz) {
        membrane.type = Tmdet::Types::MembraneType::BLENDED;
        membrane.sphereRadius = bestSphereRadius;
        membrane.origo = bestOrigo;
        DEBUG_LOG("setMembraneOrigo: origo:{} r:{}",membrane.origo,membrane.sphereRadius);
    }
    
}
