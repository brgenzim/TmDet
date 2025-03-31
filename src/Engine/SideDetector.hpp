// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

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
     * @brief set region types for residues of protein either in plane or curved membrane
     */
    class SideDetector {
        protected:
            /**
             * @brief the protein value object
             */
            Tmdet::VOs::Protein& protein;

            /**
             * @brief membrane planes
             *
             * for one membrane:
             *              side1
             *        z1 -------------   +
             *             membrane    -----
             *        z4 -------------   -
             *              side2
             *
             * for two membranes: z1 > z2 > z3 > z4
             *              side1
             *        z1 ------------- top                        +
             *             membrane                             -----
             *        z2 ------------- bottom of upper membrane   -
             *
             *           intermembrane
             *           (or periplasm)
             *
             *        z3 ------------- top                        +
             *             membrane                             -----
             *        z4 ------------- bottom of lower membrane   -
             *               side2
             */
            double z1;
            double z2;
            double z3;
            double z4;
            double o1;
            double o2;

            /**
             * @brief type of the side detector (plabe or curved)
             */
            std::string type="";
            
            /**
             * @brief main entry point of side detection
             */
            void run();

            /**
             * @brief clear temporary data 
             * 
             */
            void end();

            /**
             * @brief Set type of residues according to their z coordinate
             * 
             * @param typeName 
             * @param membranes 
             */
            void setType(std::string typeName, const std::vector<Tmdet::VOs::Membrane>& membranes);

            /**
             * @brief Get the distance of the given atom
             * 
             * @param vec 
             * @return double 
             */
            virtual double getDistance(const gemmi::Vec3& vec) = 0;

            /**
             * @brief Set z1, z2, z3, z4 values
             * 
             * @param membranes 
             */
            virtual void setZs(const std::vector<Tmdet::VOs::Membrane>& membranes) = 0;

            /**
             * @brief Get side according the z coordinate
             * 
             * @param residue 
             * @param z 
             * @return Tmdet::Types::Region 
             */
            Tmdet::Types::Region getSideByZ(Tmdet::VOs::Residue& residue, double z) const;

            /**
             * @brief Set direction value of residues
             */
            void setDirection();

            /**
             * @brief get z coordinate for direction calculation
             * 
             * @param chain 
             * @param pos 
             * @return double 
             */
            double getZForDirection(Tmdet::VOs::Chain& chain, int pos);

            /**
             * @brief Get direction of the given residue
             * 
             * @param chain 
             * @param pos 
             * @return double 
             */
            double getResidueDirection(Tmdet::VOs::Chain& chain, int pos);
            
        public:
            /**
            * @brief Construct a new Side Detector object
            * 
            * @param protein 
            */
            explicit SideDetector(Tmdet::VOs::Protein& protein) :
                protein(protein) {
            }

            /**
             * @brief Destroy the Side Detector object
             */
            ~SideDetector() {
                end();
            }
    };
}