#pragma once

#include <string>
#include <vector>
#include <gemmi/unitcell.hpp>
#include <gemmi/math.hpp>

/**
 * @brief namespace for value objects
 */
namespace Tmdet::ValueObjects {

    /**
     * @brief description of transformation matrix to
     *        transform the protein structure that membrane
     *        plane will be in the x-y plane and z axes is
     *        the membrane normal, the origo is in the middle
     *        of the membrane
     * 
     */
    struct TMatrix {
        /**
         * @brief rotation
         */
        gemmi::Mat33 rot = {1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0};

        /**
         * @brief translation
         */
        gemmi::Vec3 trans = {0.0,0.0,0.0};
    };
}
