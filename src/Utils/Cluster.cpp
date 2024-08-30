#include <vector>
#include <any>
#include <iostream>
#include <cmath>
#include <limits>
#include <gemmi/model.hpp>
#include <Utils/Cluster.hpp>

namespace Tmdet::Utils {

    void Cluster::extractCAlphas() {
        for (auto& chain : tmdetVO.chains) {
            for (auto& residue : chain.residues) {
                const gemmi::Atom* ca_atom = residue.gemmi.get_ca();
                if (ca_atom) {
                    cas.push_back(ca_atom->pos);
                }
            }
        }
    }

    void Cluster::calculateCentroid(_cluster& cluster) {
        cluster.centroid = gemmi::Position(0.0, 0.0, 0.0);
        for (int point : cluster.points) {
            cluster.centroid += cas[point]->pos;
        }
        cluster.centroid /= cluster.points.size();
    }

    void Cluster::calculateVariance(_cluster& cluster) {
        cluster.variance = 0.0;
        for (int point : cluster.points) {
            double d = cas[point]->pos.dist(cluster.centroid);
            cluster.variance += d * d;
        }
        cluster.variance = sqrt(cluster.variance);
    }

    double Cluster::calculateWardDistance(const _cluster& c1, const _cluster& c2) {
        _cluster mergedCluster = c1;
        mergedCluster.points.insert(mergedCluster.points.end(), c2.points.begin(), c2.points.end());
        calculateCentroid(mergedCluster);
        calculateVariance(mergedCluster);
        return mergedCluster.variance - (c1.variance + c2.variance);
    }

    int Cluster::calculateChainBrakes(const _cluster& c1, const _cluster& c2) {
        int chainBrakes = 0;
        for (int point1 : c1.points) {
            for (int point2: c2.points) {
                chainBrakes += (cas[point1]->pos.dist(cas[point2]->pos) < CA_DIST);
            }
        }
        return chainBrakes;
    }

    void Cluster::run() {
        extractCAlphas();
        int n = ca_positions.size();
        std::vector<_cluster> clusters(n);

        for (int i = 0; i < n; ++i) {
            clusters[i].points.push_back(i);
            clusters[i].centroid = cas[i]->pos;
            clusters[i].variance = 0.0;
        }

        while (clusters.size() > 1) {
            double minDist = 1e30;
            int c1Idx = -1, c2Idx = -1;

            for (int i = 0; i < clusters.size(); ++i) {
                for (int j = i + 1; j < clusters.size(); ++j) {
                    double dist = calculateWardDistance(clusters[i], clusters[j]) 
                        + 100 * calculateChainBrakes(clusters[i], clusters[j]);
                    if (dist < minDist) {
                        minDist = dist;
                        c1Idx = i;
                        c2Idx = j;
                    }
                }
            }

            _cluster mergedCluster = clusters[c1Idx];
            mergedCluster.points.insert(mergedCluster.points.end(), clusters[c2Idx].points.begin(), clusters[c2Idx].points.end());
            calculateCentroid(mergedCluster);
            calculateVariance(mergedCluster);

            clusters.erase(clusters.begin() + c2Idx);
            clusters.erase(clusters.begin() + c1Idx);
            clusters.push_back(mergedCluster);

            std::cout << "Merged clusters " << cas[c1Idx].serial << " and " << cas[c2Idx].serial << " with distance " << minDist << std::endl;
        }

        std::cout << "Final cluster: ";
        for (int point : clusters[0].points) {
            std::cout << point << " ";
        }
        std::cout << std::endl;
    }

    
}
