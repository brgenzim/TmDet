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

        /**
         * @brief the membrane normal vector after transformation
         *        i.e. x and y should be zero, while z is the half
         *        width of the membrane
         */
        gemmi::Vec3 normal;

        /**
         * @brief membrane width in angstrom
         */
        double h;

        /**
         * @brief the radius of the membrane curvature 
         *        zero for plain membrane
         * 
         */
        double curver = 0.0;

        /**
         * @brief the radius of membrane in the x-y plane that
         *        contains the protein (for Mol* representation)
         */
        double sizer;

        /**
         * @brief type of the membrane (plain or curved)
         */
        Tmdet::Types::Membrane type;
    };
}
