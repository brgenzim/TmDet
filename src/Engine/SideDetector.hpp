#pragma once

#include <Types/Region.hpp>
#include <VOs/Membrane.hpp>
#include <VOs/Protein.hpp>
#include <VOs/Residue.hpp>

/**
 * @brief namespace for tmdet engine
 */
namespace Tmdet::Engine {

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
            std::string type="";
            
            void run();
            void end();
            void setType(std::string typeName, const std::vector<Tmdet::VOs::Membrane>& membranes);
            virtual double getDistance(const gemmi::Vec3& vec) = 0;
            virtual void setZs(const std::vector<Tmdet::VOs::Membrane>& membranes) = 0;
            Tmdet::Types::Region getSideByZ(Tmdet::VOs::Residue& residue, double z) const;
            void setDirection();
            double getZForDirection(Tmdet::VOs::Chain& chain, int pos);
            double getResidueDirection(Tmdet::VOs::Chain& chain, int pos);
            
        public:
            explicit SideDetector(Tmdet::VOs::Protein& protein) :
                protein(protein) {
            }
            ~SideDetector() {
                end();
            }
    };
}