#pragma once

#include <array>
#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include <Types/Region.hpp>
#include <Utils/SecStrVec.hpp>
#include <ValueObjects/Protein.hpp>
#include <ValueObjects/Residue.hpp>
#include <ValueObjects/Membrane.hpp>

/**
 * @brief namespace for tmdet engine
 */
namespace Tmdet::Engine {

    class Annotator {
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
             *        z1 -------------
             *             membrane
             *        z4 ------------- 
             *              side2
             *
             * for two membranes:
             *              side1
             *        z1 ------------- top
             *             membrane
             *        z2 ------------- bottom of upper membrane
             *
             *           intermembrane
             *           (or periplasm)
             *
             *        z3 ------------- top
             *             membrane
             *        z4 ------------- bottom of lower membrane
             *               side2
             */
            double z1;
            double z2;
            double z3;
            double z4;
            bool doubleMembrane = false;
            Tmdet::Utils::SecStrVec ssVec = Tmdet::Utils::SecStrVec(protein);

            void setZs();
             Tmdet::Types::Region getSideByZ(double z) const;
            void storeRegion(Tmdet::ValueObjects::Chain& chain, unsigned int beg, unsigned int end) const;
            void replaceRegion(const Tmdet::ValueObjects::SecStrVec& vector, Tmdet::Types::Region regionType);
            bool getNextRegion(Tmdet::ValueObjects::Chain& chain, int& begin, int& end) const;
            bool getNextDefined(Tmdet::ValueObjects::Chain& chain, int& begin) const;
            bool getNextSame(Tmdet::ValueObjects::Chain& chain, const int& begin, int& end) const;


        public:
            void detectSides();
            void detectAlphaHelices();
            void detectBarrel();
            void detectLoops();
            void detectInterfacialHelices();
            void getRegions();

            explicit Annotator(Tmdet::ValueObjects::Protein& protein) :
                protein(protein) {}
            ~Annotator()=default;
    };
}