#pragma once

#include <Types/Region.hpp>
#include <ValueObjects/Membrane.hpp>
#include <ValueObjects/Protein.hpp>
#include <ValueObjects/Residue.hpp>

/**
 * @brief namespace for tmdet engine
 */
namespace Tmdet::Engine {

    class SideDetector {
        private:
            /**
             * @brief the protein value object
             */
            Tmdet::ValueObjects::Protein& protein;

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
            
            void run();
            void end();
            void setType(std::string typeName, const std::vector<Tmdet::ValueObjects::Membrane>& membranes);
            void setZs(const std::vector<Tmdet::ValueObjects::Membrane>& membranes);
            Tmdet::Types::Region getSideByZ(Tmdet::ValueObjects::Residue& residue, double z) const;
            void setDirection();
            double getResidueDirection(Tmdet::ValueObjects::Chain& chain, int pos);
            
        public:
            explicit SideDetector(Tmdet::ValueObjects::Protein& protein) :
                protein(protein) {
                    run();
            }
            ~SideDetector() {
                end();
            }
    };
}