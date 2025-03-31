// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <gemmi/unitcell.hpp>
#include <gemmi/math.hpp>
#include <Types/Membrane.hpp>
#include <VOs/TMatrix.hpp>

/**
 * @brief namespace for value objects
 */
namespace Tmdet::VOs {

    /**
     * @brief representation of a membrane
     */
    struct Membrane {

        /**
         * @brief centre of the membrane after transformation 
         *        for the firs membrane it is always zero
         *
         *        For curved membrane the origo is the centre of the
         *        sphere and the membrane is between sphereRadius -
         *        halfThickness and sphereRadius + halThickness
         */
        double origo = 0.0;

        /**
         * @brief membrane half thickness in angstrom
         */
        double halfThickness = 8.0;

        /**
         * @brief the radius of the membrane curvature 
         *        zero for plain membrane
         * 
         */
        double sphereRadius = 0.0;

        /**
         * @brief the radius of membrane in the x-y plane that
         *        contains the protein (for Mol* representation)
         */
        double membraneRadius = 10.0;

        /**
         * @brief type of the membrane (plain or curved)
         */
        Tmdet::Types::Membrane type = Tmdet::Types::MembraneType::PLANE;
    };
}
