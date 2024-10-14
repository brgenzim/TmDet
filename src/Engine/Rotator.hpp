#pragma once

#include <vector>
#include <map>
#include <set>
#include <Config.hpp>

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
             * @brief distance for next step in vertical direction
             */
            double alpha_dist;

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
            double beta_step;

            /**
             * @brief helper internal variables
             * 
             */
            double q;
            double qq;

            /**
             * @brief calculate next alpha value and set beta_step
             */
            void nextAlpha();

        public:

            /**
             * @brief Construct a new Rotator object
             * 
             */
            explicit Rotator() {
                alpha_dist = std::stof(environment.get("TMDET_BALL_DIST",DEFAULT_TMDET_BALL_DIST)) / 2;
            }

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
