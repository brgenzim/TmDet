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
        int numStraight;

        /**
         * @brief number of turns in the slice
         * 
         */
        int numTurn;

        /**
         * @brief number of all residues (C alpha atoms) in the slice
         */
        int numCA;

        /**
         * @brief outside water accessible surface of the atoms in the slice
         */
        double surf;

        /**
         * @brief sum up voromqa_v1_energy_means producted with atom surface
         */
        double voronota;
    };

    /**
     * @brief class for searching for membrane plane
     */
    class Optim {
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
             * @brief properties of the actual membrane definition
             */
            Tmdet::ValueObjects::Membrane membraneVO;

            /**
             * @brief the protein structure in Protein Value Object
             */
            Tmdet::ValueObjects::Protein& proteinVO;

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
            void setAtomDistances(Tmdet::ValueObjects::Residue& residue);

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

        public:

            /**
             * @brief Construct a new Optim object
             * 
             * @param proteinVO 
             */
            explicit Optim(Tmdet::ValueObjects::Protein& proteinVO) : 
                proteinVO(proteinVO) {}

            /**
             * @brief Destroy the Optim object
             * 
             */
            ~Optim()=default;

            /**
             * @brief get Q value for a given membrane normal
             * 
             * @return double 
             */
            double getQValue();
    };
}
