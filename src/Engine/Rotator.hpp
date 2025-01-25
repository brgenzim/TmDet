// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <gemmi/math.hpp>

/**
 * @brief namespace for tmdet engine
 *
 * @namespace Tmdet
 * @namespace Engine
 */
namespace Tmdet::Engine {

    /**
     * @brief rotate a normal vector around the 4 PI
     */
    class Rotator {
        private:

            /**
             * @brief current alpha angle of the normal vector
             *        vertical directtion
             */
            double alpha=0;

            /**
             * @brief step of alpha angle
             */
            double alpha_step;

            /**
             * @brief end of alpha angle rotation
             */
            double alpha_end;

            /**
             * @brief current beta angle of the normal vector
             *        horizontal angle
             */
            double beta=0;

            /**
             * @brief step of beta angle (depends on alpha)
             */
            double beta_step=2*M_PI;

            /**
             * @brief helper internal variables
             * 
             */
            double q=0;
            double qq=1;

            /**
             * @brief calculate next alpha value and set beta_step
             */
            void nextAlpha();

        public:

            /**
             * @brief Construct a new Rotator object
             * 
             */
            Rotator();

            /**
             * @brief Destroy the Rotator object
             * 
             */
            ~Rotator()=default;

            /**
             * @brief get the next normal vector and return false when there is no more
             * 
             * @param normal 
             * @return true 
             * @return false 
             */
            bool next(gemmi::Vec3& normal);

            void end90() {
                alpha_end = M_PI / 2;
            }

            void end180() {
                alpha_end = M_PI;
            }
    };
}
