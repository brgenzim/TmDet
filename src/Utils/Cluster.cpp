#include <vector>
#include <any>
#include <mlpack/methods/emst/dtb.hpp>
#include <mlpack/methods/emst/emst.hpp>
#include <mlpack/methods/dbscan.hpp>
#include <mlpack/core.hpp>
#include <mlpack.hpp>
#include <gemmi/model.hpp>
#include <Utils/Cluster.hpp>

namespace Tmdet::Utils {

    void Cluster::run() {
        
        std::vector<gemmi::Position> ca_positions;

        // Extract Cα positions
        for (auto& chain : tmdetVO.chains) {
            for (auto& residue : chain.residues) {
                const gemmi::Atom* ca_atom = residue.gemmi.get_ca();
                if (ca_atom) {
                    ca_positions.push_back(ca_atom->pos);
                }
            }
        }
    
        // Create the distance matrix
        size_t n = ca_positions.size();
        arma::mat distance_matrix(n, n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i + 1; j < n; ++j) {
                double distance = ca_positions[i].dist(ca_positions[j]);
                distance_matrix(i, j) = distance_matrix(j, i) = distance;
            }
        }

        // Run hierarchical agglomerative clustering
        //arma::Row<size_t> assignments;
        //mlpack::tree::DBSCAN<> hdbscan(distance_matrix, assignments, 2);

        // Adjust distance matrix for DBSCAN (DBSCAN requires a set of points, not a distance matrix)
        arma::mat ca_positions_mat(3, n);
        for (size_t i = 0; i < n; ++i) {
            ca_positions_mat(0, i) = ca_positions[i].x;
            ca_positions_mat(1, i) = ca_positions[i].y;
            ca_positions_mat(2, i) = ca_positions[i].z;
        }

        // Set DBSCAN parameters
        double epsilon = 5.0; // Adjust as needed
        size_t minPoints = 3; // Minimum number of points in a cluster

        // Run DBSCAN
        arma::Row<size_t> assignments;
        mlpack::dbscan::DBSCAN<> dbscan(epsilon, minPoints);
        dbscan.Cluster(ca_positions_mat, assignments);

        //Write back the clustering results to residues any 
        size_t i = 0;
        for (auto& chain : tmdetVO.chains) {
            for (auto& residue : chain.residues) {
                const gemmi::Atom* ca_atom = residue.gemmi.get_ca();
                if (ca_atom) {
                     residue.temp.insert({"cluster",any_cast<int>(assignments[i++])});
                }
            }
        }
    }
}
