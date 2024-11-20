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

        void transform(gemmi::Vec3& vec) {
            vec.x += trans.x;
            vec.y += trans.y;
            vec.z += trans.z;
            double x = vec.x * rot[0][0]
                        + vec.y * rot[0][1]
                        + vec.z * rot[0][2];
            double y = vec.x * rot[1][0]
                        + vec.y * rot[1][1]
                        + vec.z * rot[1][2];
            double z = vec.x * rot[2][0]
                        + vec.y * rot[2][1]
                        + vec.z * rot[2][2];
            vec.x = x;
            vec.y = y;
            vec.z = z;

        }
    };
}
