#pragma once

#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include <VOs/Protein.hpp>
#include <VOs/Residue.hpp>
#include <VOs/Membrane.hpp>
#include <VOs/Slice.hpp>


/**
 * @brief namespace for tmdet engine
 */
namespace Tmdet::Engine {

#define RES(res,a) (protein.chains[res.chainIdx].residues[res.idx + a])

    

    /**
     * @brief class for searching for membrane plane
     */
    class Optimizer {
        protected:

            std::string type = "";

            /**
             * @brief 1 Angstrom wide slices of the protein along the z axes
             */
            std::vector<Tmdet::VOs::Slice> slices;

            /**
             * @brief the actual membrane normal
             */
            gemmi::Vec3 normal;

            /**
             * @brief the protein structure in Protein Value Object
             */
            Tmdet::VOs::Protein& protein;

            /**
             * @brief best qValue
             */
            double bestQ = 0.0;

            /**
             * @brief membrane normal belonging to the best qValue
             */
            gemmi::Vec3 bestNormal;

            /**
             * @brief best slices
             */
            std::vector<Tmdet::VOs::Slice> bestSlices;

            /**
             * @brief minimum on z axes for the given normal vector
             */
            double minZ = 1e30;

            /**
             * @brief maximum on z axes for the given normal vector
             */
            double maxZ = -1e30;

            /**
             * @brief the mass centre of the protein
             */
            gemmi::Vec3 massCentre;

            double lastO = 0;

            double bestOrigo = 0;


            /**
             * @brief initialize the algorithm
             */
            void init();

            /**
             * @brief end and clean of the algorithm
             */
            void end();

            virtual double distance(gemmi::Vec3& vec) = 0;

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
             * @brief smooth Q values 
             */
            void smoothQValues();

            /**
             * @brief Get the best qValue that has the apropriate membrane width
             */
            void checkBestSlice();

            double getWidth(const int z, int& minz, int& maxz);

            void testMembraneNormalOne();

            /**
             * @brief check if tests resulted valid membrane definition
             *
             * @return bool
             */
            bool isTransmembrane() const;

            virtual void setBestOrigo(double minz, double maxz) = 0;

            /**
             * @brief set membrane width using the best membrane definition
             *        and the qValues of slices
             * @return bool
             */
            bool getMembrane(Tmdet::VOs::Membrane& membrane, int count);

            /**
             * @brief set membrane origo
             */
            virtual void setMembraneOrigo(Tmdet::VOs::Membrane& membrane, double minz, double maxz) = 0;

        public:

            /**
             * @brief Construct a new Optim object
             * 
             * @param protein
             */
            explicit Optimizer(Tmdet::VOs::Protein& protein) : 
                protein(protein) {
                    init();
                }

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
             * @brief get the normal vector for the best qValue
             */
            gemmi::Vec3 getBestNormal() const {
                return bestNormal;
            }

            /**
             * @brief calculate Q value for a given normal
             */
            virtual void testMembraneNormal() = 0;

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
            void setProteinTMatrix(gemmi::Vec3& origo) const;

            /**
             * @brief clear former results
             */
            void clear();

    };
}
