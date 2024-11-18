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

#define RES(res,a) (protein.chains[res.chainIdx].residues[res.idx + a])

    /**
     * @brief description of a slice
     */
    struct _slice {
        
        /**
         * @brief apolar surface
         */
        double apol = 0.0;

        /**
         * @brief outside water accessible surface of the atoms in the slice
         */
        double surf = 0.0;

        /**
         * @brief turn ratio in the slice
         * 
         */
        double turn = 0.0;

        /**
         * @brief surface of residues being in turn
         */
        double tSurf = 0.0;

        /**
         * @brief ratio of straight element in the slice
         */
        double straight = 0.0;

        /**
         * @brief hidrophob ratio in the slice
         */
        double hydrophobicity = 0.0;

        /**
         * @brief sum up voromqa_v1_energy_means producted with atom surface
         */
        double voronota = 0.0;

        /**
         * @brief ratio of chain ends in the slice
         */
        double chainEnd = 0.0;

        /**
         * @brief number of C alpha atoms in the slice
         */
        int numCa = 0;

        /**
         * @brief number of all residues in the slice
         */
        int numAtom = 0;
        
        /**
         * @brief number of hydrophob residues in the slice
         */
        int numHyd = 0;

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
            bool run = false;

            /**
             * @brief minimum on z axes for the given normal vector
             */
            double min = 1e30;

            /**
             * @brief maximum on z axes for the given normal vector
             */
            double max = -1e30;

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

            /**
             * @brief best slices
             */
            std::vector<_slice> bestSlices;


            double lastO = 0;

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
             * @brief Set the straightness of a residue
             */
            void setStraight();

            /**
             * @brief Set the box containing the protein
             */
            void setBoundaries();

            /**
             * @brief smoothin apolar surface values
             */
            void smoothSurf();

            /**
             * @brief sumup slice properties
             */
            void sumupSlices();

            /**
             * @brief helper function to divide two numbers
             * 
             * @param numerator 
             * @param denominator 
             * @return double 
             */
            double divide(double numerator, double denominator);

            /**
             * @brief calculate the value of the objective function for one slice
             * 
             * @param s 
             * @return double 
             */
            double getQValueForSlice(const _slice& s);

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

            void clear();
    };
}
