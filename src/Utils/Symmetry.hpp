#pragma once

#include <vector>
#include <string>
#include <any>
#include <eigen3/Eigen/Dense>
#include <gemmi/model.hpp>
#include <ValueObjects/Protein.hpp>

// #define __SYM_DBG 1 // to enable debug messages of this feature

/**
 * @brief namespace for tmdet objects
 * @namespace Tmdet
 */
namespace Tmdet {
    
    /**
     * @brief namespace for various utilities
     * @namespace Utils
     */
    namespace Utils {

        /**
         * @brief temporary container for the results of symmetry operation
         */
        struct _symmetryData {
            int id = 0;
            int cidx1;
            int cidx2;
            std::string entityId;
            gemmi::Vec3 origo;
            gemmi::Vec3 axis;
            double rotAngle;

            double distance(struct _symmetryData& other) {
                return axis.dist(other.axis)
                        + abs(rotAngle - other.rotAngle);
            }

            bool same(struct _symmetryData& other) {
                return distance(other) < 7;
            }
        };

        /**
         * @brief class for symmetry related operations,
         *        it gives back the rotational axeses of 
         *        homooligomer parts (if any) of the structure
         */
        class Symmetry {
            private:
                /**
                 * @brief the structure in protein value object
                 */
                Tmdet::ValueObjects::Protein& protein;

                /**
                 * @brief calculated symmetry operation between
                 *        chain pairs
                 */
                std::vector<std::vector<_symmetryData>> sim;

                void run();
                void initSymetryContainer();
                std::vector<_symmetryData> getRotationalAxes();
                void searchForRotatedChains(int cidx1);
                void calculateRotationalOperation(int cidx1, int cidx2);
                void getCoordinates(int cidx1, int cidx2, std::vector<Eigen::Vector3d>& coord1, std::vector<Eigen::Vector3d>& coord2, Eigen::Vector3d& t1, Eigen::Vector3d& t2);
                void getSymmetryOperand(Eigen::Matrix4d& R, Eigen::Vector3d& t1, Eigen::Vector3d& t2, _symmetryData& simij);
                bool lsqFit(std::span<Eigen::Vector3d>& r1, std::span<Eigen::Vector3d>& r2, double& rmsd, Eigen::Matrix4d& Rot);
                std::vector<_symmetryData> clusterAxes(std::vector<_symmetryData>& axes);

            public:
                /**
                 * @brief Construct a new Symmetry object
                 * 
                 * @param protein 
                 */
                explicit Symmetry(Tmdet::ValueObjects::Protein& protein) :
                    protein(protein) {
                        run();
                }
                
                /**
                 * @brief Destroy the Symmetry object
                 */
                ~Symmetry()=default;

                /**
                 * @brief get definition of possible membrane planes
                 */
                std::vector<Tmdet::ValueObjects::Membrane> getMembraneAxes();
        };
    }
}
