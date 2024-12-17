// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <vector>
#include <string>
#include <any>
#include <format>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <gemmi/model.hpp>
#include <VOs/Protein.hpp>

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
            bool good = false;
            gemmi::Vec3 origo = {0.0,0.0,0.0};
            gemmi::Vec3 axis = {0.0,0.0,0.0};

            double distance(const struct _symmetryData& other) const {
                return axis.dist(other.axis);
            }

            bool same(const struct _symmetryData& other) const {
                return distance(other) < 7;
            }

            friend std::ostream& operator<<(std::ostream& os, const _symmetryData& other) {
                os << std::format("Origo: {:.4f} {:.4f} {:.4f}",other.origo.x,other.origo.y,other.origo.z) << std::endl;
                os << std::format(" Axes: {:.4f} {:.4f} {:.4f}",other.axis.x,other.axis.y,other.axis.z)  << std::endl;
                return os;
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
                Tmdet::VOs::Protein& protein;

                /**
                 * @brief calculated symmetry operation between
                 *        chain pairs
                 */
                std::vector<_symmetryData> sim;

                std::vector<Tmdet::VOs::Membrane> run();
                std::vector<_symmetryData> getRotationalAxes();
                bool searchForRotatedChains(const std::vector<std::string>& chainIds);
                int calculateRotationalOperation(int cidx1, int cidx2);
                void getCoordinates(int cidx1, int cidx2, std::vector<Eigen::Vector3d>& coord1, std::vector<Eigen::Vector3d>& coord2, Eigen::Vector3d& t1, Eigen::Vector3d& t2);
                void getSymmetryOperand(Eigen::Matrix4d& R, const Eigen::Vector3d& t1, const Eigen::Vector3d& t2, _symmetryData& simij) const;
                bool lsqFit(const std::span<Eigen::Vector3d>& r1, const std::span<Eigen::Vector3d>& r2, double& rmsd, Eigen::Matrix4d& Rot) const;
                Eigen::Matrix4d rotateZ(Eigen::Vector3d T) const;
                bool haveSameAxes() const;
                _symmetryData getAverageAxes() const;
                std::vector<gemmi::Vec3> clusterAxes(std::vector<_symmetryData> axes);
                
            public:
                /**
                 * @brief Construct a new Symmetry object
                 * 
                 * @param protein 
                 */
                explicit Symmetry(Tmdet::VOs::Protein& protein) :
                    protein(protein) {}
                
                /**
                 * @brief Destroy the Symmetry object
                 */
                ~Symmetry()=default;

                /**
                 * @brief get definition of possible membrane planes
                 */
                std::vector<gemmi::Vec3> getMembraneAxes();
        };
    }
}
