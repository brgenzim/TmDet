#pragma once

#include <gemmi/model.hpp>
#include <Engine/Optimizer.hpp>
#include <VOs/Protein.hpp>

/**
 * @brief namespace for tmdet engine
 */
namespace Tmdet::Engine {

    /**
     * @brief class for searching for plane membrane
     */
    class BlendedOptimizer : public Optimizer {
        private: 
            gemmi::Vec3 origoVec3;
            double origo;
            double bestSphereRadius;

        protected:
            
            double distance(gemmi::Vec3& vec);

            void setBestOrigo(double minz, double maxz);

            /**
             * @brief set membrane width using the best membrane definition
             *        and the qValues of slices
             * @return bool
             */
            bool getMembrane(Tmdet::VOs::Membrane& membrane, int count);

            /**
             * @brief set membrane origo
             */
            void setMembraneOrigo(Tmdet::VOs::Membrane& membrane, double minz, double maxz);

        public:
            BlendedOptimizer(Tmdet::VOs::Protein& protein) :
                Optimizer(protein) {
                    type = "Blended";
                }
            
            void testMembraneNormal();
    };
}
