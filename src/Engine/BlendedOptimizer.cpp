#include <gemmi/model.hpp>
#include <Engine/BlendedOptimizer.hpp>
#include <Engine/Rotator.hpp>

namespace Tmdet::Engine {

    const std::vector<double> origos = {
        10, 11, 12, 13, 14, 15, 17, 20, 25, 30, 40, 50, 60, 80, 100, 120, 150, 200, 250, 300, 350, 400
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
        bestSphereRadius = (maxz -minz) / 2 + minZ -5;
    }

    void BlendedOptimizer::setMembraneOrigo(Tmdet::VOs::Membrane& membrane, double minz, double maxz) {
        membrane.type = Tmdet::Types::MembraneType::BLENDED;
        membrane.sphereRadius = bestSphereRadius;
        membrane.origo = bestOrigo;
        DEBUG_LOG("setMembraneOrigo: origo:{} r:{}",membrane.origo,membrane.sphereRadius);
    }
    
}
