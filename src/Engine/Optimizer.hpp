// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include <System/Arguments.hpp>
#include <VOs/Protein.hpp>
#include <VOs/Residue.hpp>
#include <VOs/Membrane.hpp>
#include <VOs/Slice.hpp>


/**
 * @brief namespace for tmdet engine
 *
 * @namespace Tmdet
 * @namespace Engine
 */
namespace Tmdet::Engine {

#define RES(res,a) (protein.chains[res.chainIdx].residues[res.idx + a])

    /**
     * @brief class for searching for membrane plane
     */
    class Optimizer {
        protected:
            
            /**
             * @brief the protein structure in Protein Value Object
             */
            Tmdet::VOs::Protein& protein;

            /**
             * @brief command line arguments
             */
            Tmdet::System::Arguments& args;
            
            /**
             * @brief type of the optimizer (plane or curved)
             */
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

            /**
             * @brief last origo
             */
            double lastO = 0;

            /**
             * @brief best origo
             */
            double bestOrigo = 0;

            double bestMinZ;
            double bestMaxZ;
            double minHalfThickness;
            double maxHalfThickness;
            double maxCurvedHalfThickness;
            double lowerQ;
            double higherQ;
            double higherQ2;
            double ifhAngleLimit;
            int ifhResLimit;
            int ifhMinLength;
            double boostAngle;
            double boostBetaAngle;
            double boostPolarity;
            
            /**
             * @brief initialize the algorithm
             */
            void init();

            /**
             * @brief clean temporary data
             */
            void end();

            /**
             * @brief distance of the atom from the membrane plane
             *        or the centre of the sphere
             * 
             * @param vec 
             * @return double 
             */
            virtual double distance(gemmi::Vec3& vec) = 0;


            virtual double getAngle(Tmdet::VOs::SecStrVec& vector) = 0;

            virtual void testMembraneNormalFinal() = 0;

            /**
             * @brief Set the distances from the centre of membrane plane
             */
            void setDistances();

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

            /**
             * @brief Get the width of the membrane (region of slices those qValue
             *        is above TMDET_MEMBRANE_VALUE)
             * 
             * @param z 
             * @param minz 
             * @param maxz 
             * @return double 
             */
            double getWidth(const int z, int& minz, int& maxz);

            /**
             * @brief calculate qValue for the given membrane normal
             */
            void testMembraneNormalOne();

            /**
             * @brief check if tests resulted valid membrane definition
             *
             * @return bool
             */
            bool isTransmembrane() const;

            /**
             * @brief Set place of the best origo
             * 
             * @param minz 
             * @param maxz 
             */
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
            explicit Optimizer(Tmdet::VOs::Protein& protein, Tmdet::System::Arguments& args) : 
                protein(protein),
                args(args) {
                    init();
                }

            /**
             * @brief Destroy the Optim object
             * 
             */
            ~Optimizer();

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

            /**
             * @brief Get the type of the optimizer
             * 
             * @return std::string 
             */
            std::string getType() const {
                return type;
            }

    };
}
