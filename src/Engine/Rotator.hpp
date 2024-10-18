#pragma once

#include <gemmi/math.hpp>

/**
 * @brief namespace for tmdet engine
 */
namespace Tmdet::Engine {

    /**
     * @brief rotator class is to rotate a normal vector around the 4 PI
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
    };
}
