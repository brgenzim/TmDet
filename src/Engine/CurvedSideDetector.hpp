// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <Engine/SideDetector.hpp>
#include <Types/Region.hpp>
#include <VOs/Membrane.hpp>
#include <VOs/Protein.hpp>
#include <VOs/Residue.hpp>

/**
 * @brief namespace for tmdet engine
 *
 * @namespace Tmdet
 * @namespace Engine
 */
namespace Tmdet::Engine {

    /**
     * @brief set region types for residues of protein in curved membrane
     */
    class CurvedSideDetector : public SideDetector {
        protected:
            /**
             * @brief calculate distance of an atom from membrane plane
             * 
             * @param vec 
             * @return double 
             */
            double getDistance(const gemmi::Vec3& vec);

            /**
             * @brief Set the z coordinates of membrane boundaries
             * 
             * @param membranes 
             */
            void setZs(const std::vector<Tmdet::VOs::Membrane>& membranes);
            
        public:
            /**
             * @brief Construct a new Curved Side Detector object
             * 
             * @param protein 
             */
            explicit CurvedSideDetector(Tmdet::VOs::Protein& protein) :
                SideDetector(protein) {
                    type="Curved";
                    run();
            }
            
    };
}