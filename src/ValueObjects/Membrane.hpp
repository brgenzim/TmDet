#pragma once

#include <gemmi/unitcell.hpp>
#include <gemmi/math.hpp>
#include <Types/Membrane.hpp>
#include <ValueObjects/TMatrix.hpp>

/**
 * @brief namespace for value objects
 */
namespace Tmdet::ValueObjects {

    /**
     * @brief representation of a membrane
     */
    struct Membrane {

        /**
         * @brief transformation matrix for the protein
         *        to be the membrane plane in the xy plane
         *        and the membrane normal parallel to the z axes
         */
        TMatrix tmatrix;

        /**
         * @brief the last column in the transformation matrix
         */
        gemmi::Vec3 origo;
        gemmi::Vec3 normal;
        double h;
        double curver;
        double sizer;
        Tmdet::Types::Membrane type;
    };

}
