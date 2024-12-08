#pragma once

#include <Engine/SideDetector.hpp>
#include <Types/Region.hpp>
#include <VOs/Membrane.hpp>
#include <VOs/Protein.hpp>
#include <VOs/Residue.hpp>

/**
 * @brief namespace for tmdet engine
 */
namespace Tmdet::Engine {

    class PlaneSideDetector : public SideDetector {
        protected:
            double getDistance(const gemmi::Vec3& vec);
            void setZs(const std::vector<Tmdet::VOs::Membrane>& membranes);
            
        public:
            explicit PlaneSideDetector(Tmdet::VOs::Protein& protein) :
                SideDetector(protein) {
                    type="Plane";
                    run();
            }
            
    };
}