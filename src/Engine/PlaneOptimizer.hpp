// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <gemmi/model.hpp>
#include <Engine/Optimizer.hpp>

/**
 * @brief namespace for tmdet engine
 *
 * @namespace Tmdet
 * @namespace Engine
 */
namespace Tmdet::Engine {

    /**
     * @brief class for searching for plane membrane
     */
    class PlaneOptimizer : public Optimizer {

        protected:
            /**
             * @brief coordinate of the best origo 
             */
            gemmi::Vec3 bestOrigo;

            /**
             * @brief distance of the atom from the membrane plane
             * 
             * @param vec 
             * @return double 
             */
            double distance(gemmi::Vec3& vec);

            /**
             * @brief Set place of the best origo
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
             * @brief Construct a new Plane Optimizer object
             * 
             * @param protein 
             */
            PlaneOptimizer(Tmdet::VOs::Protein& protein) :
                Optimizer(protein) {
                    type = "Plain";
                }
            
            /**
             * @brief calculate Q value for a given normal
             */
            void testMembraneNormal();
    };
}
