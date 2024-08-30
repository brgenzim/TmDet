#include <vector>
#include <any>
#include <iostream>
#include <cmath>
#include <limits>
//#include <mlpack/methods/emst/dtb.hpp>
//#include <mlpack/methods/emst/emst.hpp>
//#include <mlpack/methods/dbscan.hpp>
//#include <mlpack/core.hpp>
//#include <mlpack.hpp>
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
        size_t minPoints = 10; // Minimum number of points in a cluster

        /*
        // Run DBSCAN
        arma::Row<size_t> assignments;
        assignments.clear();
        mlpack::dbscan::DBSCAN<> dbscan(epsilon, minPoints);
        dbscan.Cluster(ca_positions_mat, assignments);

        //Write back the clustering results to residues any 
        size_t i = 0;
        for (auto& chain : tmdetVO.chains) {
            for (auto& residue : chain.residues) {
                const gemmi::Atom* ca_atom = residue.gemmi.get_ca();
                if (ca_atom) {
                     residue.temp.insert({"cluster",std::any_cast<size_t>(assignments[i++])});
                }
            }
        }*/


        

using namespace std;

struct Cluster {
    vector<int> points;
    double centroid;
};

double calculateCentroid(const vector<double>& data, const Cluster& cluster) {
    double sum = 0.0;
    for (int point : cluster.points) {
        sum += data[point];
    }
    return sum / cluster.points.size();
}

double calculateWardDistance(const Cluster& c1, const Cluster& c2, const vector<double>& data) {
    double centroid1 = calculateCentroid(data, c1);
    double centroid2 = calculateCentroid(data, c2);

    double withinClusterVariance1 = 0.0;
    for (int point : c1.points) {
        withinClusterVariance1 += pow(data[point] - centroid1, 2);
    }

    double withinClusterVariance2 = 0.0;
    for (int point : c2.points) {
        withinClusterVariance2 += pow(data[point] - centroid2, 2);
    }

    Cluster mergedCluster = c1;
    mergedCluster.points.insert(mergedCluster.points.end(), c2.points.begin(), c2.points.end());
    double mergedCentroid = calculateCentroid(data, mergedCluster);

    double withinClusterVarianceMerged = 0.0;
    for (int point : mergedCluster.points) {
        withinClusterVarianceMerged += pow(data[point] - mergedCentroid, 2);
    }

    return withinClusterVarianceMerged - (withinClusterVariance1 + withinClusterVariance2);
}

void hierarchicalClustering(const vector<double>& data) {
    int n = data.size();
    vector<Cluster> clusters(n);

    for (int i = 0; i < n; ++i) {
        clusters[i].points.push_back(i);
        clusters[i].centroid = data[i];
    }

    while (clusters.size() > 1) {
        double minDist = numeric_limits<double>::max();
        int c1Idx = -1, c2Idx = -1;

        for (int i = 0; i < clusters.size(); ++i) {
            for (int j = i + 1; j < clusters.size(); ++j) {
                double dist = calculateWardDistance(clusters[i], clusters[j], data);
                if (dist < minDist) {
                    minDist = dist;
                    c1Idx = i;
                    c2Idx = j;
                }
            }
        }

        Cluster mergedCluster = clusters[c1Idx];
        mergedCluster.points.insert(mergedCluster.points.end(), clusters[c2Idx].points.begin(), clusters[c2Idx].points.end());
        mergedCluster.centroid = calculateCentroid(data, mergedCluster);

        clusters.erase(clusters.begin() + c2Idx);
        clusters.erase(clusters.begin() + c1Idx);
        clusters.push_back(mergedCluster);

        cout << "Merged clusters " << c1Idx << " and " << c2Idx << " with distance " << minDist << endl;
    }

    cout << "Final cluster: ";
    for (int point : clusters[0].points) {
        cout << point << " ";
    }
    cout << endl;
}

int main() {
    vector<double> data = {1.0, 2.0, 5.0, 8.0, 10.0};
    hierarchicalClustering(data);
    return 0;
}

    }
}
