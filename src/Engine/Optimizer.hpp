#pragma once

#include <array>
#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include <ValueObjects/Protein.hpp>
#include <ValueObjects/Residue.hpp>
#include <ValueObjects/Membrane.hpp>

/**
 * @brief namespace for tmdet engine
 */
namespace Tmdet::Engine {

    /**
     * @brief description of a slice
     */
    struct _slice {
        
        /**
         * @brief number of crossing straight secondary structure elements
         */
        int numStraight = 0;

        /**
         * @brief number of turns in the slice
         * 
         */
        int numTurn = 0;

        /**
         * @brief number of all residues (C alpha atoms) in the slice
         */
        int numCA = 0;

        /**
         * @brief outside water accessible surface of the atoms in the slice
         */
        double surf = 0.0;

        /**
         * @brief sum up voromqa_v1_energy_means producted with atom surface
         */
        double voronota = 0.0;

        /**
         * @brief calculated qValue for the slice
         */
        double qValue = 0.0;
    };

    /**
     * @brief class for searching for membrane plane
     */
    class Optimizer {
        private:

            /**
             * @brief flag for running (i.e. temporary containers are set)
             */
            bool run;

            /**
             * @brief minimum on z axes for the given normal vector
             */
            double min;

            /**
             * @brief maximum on z axes for the given normal vector
             */
            double max;

            /**
             * @brief 1 Angstrom wide slices of the protein along the z axes
             */
            std::vector<_slice> slices;

            /**
             * @brief the actual membrane normal
             */
            gemmi::Vec3 normal;

            /**
             * @brief the mass centre of the protein
             */
            gemmi::Vec3 massCentre;

            /**
             * @brief the protein structure in Protein Value Object
             */
            Tmdet::ValueObjects::Protein& protein;

            /**
             * @brief best qValue
             */
            double bestQ = 0.0;

            /**
             * @brief slice index of highest qValue
             */
            unsigned long int bestSliceIndex;

            /**
             * @brief membrane normal belonging to the best qValue
             */
            gemmi::Vec3 bestNormal;

            double lastO;

            /**
             * @brief initialize the algorithm
             */
            void init();

            /**
             * @brief end and clean of the algorithm
             */
            void end();

            /**
             * @brief Set the distances from the centre of membrane plane
             */
            void setDistances();

            /**
             * @brief calculate the distance of the atom from the centre of membrane plane
             * 
             * @param residue 
             */
            void setAtomDistances(Tmdet::ValueObjects::Residue& residue) const;

            /**
             * @brief Set the box containing the protein
             */
            void setBoundaries();

            /**
             * @brief sumup slice properties
             */
            void sumupSlices();

            /**
             * @brief add contribution of the residue to the appropriate slice
             * 
             * @param residue 
             */
            void residueToSlice(Tmdet::ValueObjects::Residue& residue);

            /**
             * @brief calculate the value of the objective function for one slice
             * 
             * @param s 
             * @return double 
             */
            double getQValueForSlice(const _slice& s) const;

            /**
             * @brief calculate the value of the objective function for each slice
             * 
             * @return std::vector<double> 
             */
            std::vector<double> getQValueForSlices();

            /**
             * @brief smooth Q values and return the largest one
             * 
             * @param qs 
             * @return double 
             */
            double smoothQValues(std::vector<double> qs);

            /**
             * @brief check if tests resulted valid membrane definition
             *
             * @return bool
             */
            bool isTransmembrane() const;

            /**
             * @brief set membrane width using the best membrane definition
             *        and the qValues of slices
             * @return bool
             */
            bool getMembrane(Tmdet::ValueObjects::Membrane& membrane) ;

        public:

            /**
             * @brief Construct a new Optim object
             * 
             * @param protein
             */
            explicit Optimizer(Tmdet::ValueObjects::Protein& protein) : 
                protein(protein) {}

            /**
             * @brief Destroy the Optim object
             * 
             */
            ~Optimizer() {
                end();
            };

            /**
             * @brief set the membrane normal
             * 
             * @param _normal
             */
            void setNormal(const gemmi::Vec3 _normal);

            /**
             * @brief calculate Q value for a given normal
             */
            void testMembraneNormal();

            /**
             * @brief copy membrane definition to the protein value object
             */
            void setMembranesToProtein();

            /**
             * @brief search for best membrane normal by rotating
             *        membrane normal around 4pi
             */
            void searchForMembraneNormal();

            /**
             * @brief 
             */
            void setProteinTMatrix(gemmi::Vec3& origo, gemmi::Vec3& normal) const;
    };
}
