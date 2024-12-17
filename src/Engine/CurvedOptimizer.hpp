// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <gemmi/model.hpp>
#include <Engine/Optimizer.hpp>
#include <VOs/Protein.hpp>

/**
 * @brief namespace for tmdet engine
 *
 * @namespace Tmdet
 * @namespace Engine
 */
namespace Tmdet::Engine {

    /**
     * @brief class for searching for curved membrane
     */
    class CurvedOptimizer : public Optimizer {
        private: 
            /**
             * @brief the vector that contains the origo
             */
            gemmi::Vec3 origoVec3;

            /**
             * @brief centre of the sphere
             */
            double origo;

            /**
             * @brief best sphere radius
             */
            double bestSphereRadius;

        protected:
            
            /**
             * @brief distance between a given coordinate from the origo
             * 
             * @param vec 
             * @return double 
             */
            double distance(gemmi::Vec3& vec);

            /**
             * @brief searcvh for the best sphere centre
             * 
             * @param minz 
             * @param maxz 
             */
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
            /**
             * @brief Construct a new Curved Optimizer object
             * 
             * @param protein 
             */
            CurvedOptimizer(Tmdet::VOs::Protein& protein) :
                Optimizer(protein) {
                    type = "Curved";
                }
            
            /**
             * @brief calculate Q value for a given normal
             */
            void testMembraneNormal();
    };
}
