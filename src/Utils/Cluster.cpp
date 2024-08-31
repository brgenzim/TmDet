#include <vector>
#include <any>
#include <iostream>
#include <cmath>
#include <limits>
#include <gemmi/model.hpp>
#include <Utils/Cluster.hpp>

namespace Tmdet::Utils {

    void Cluster::showData(_cluster& cluster) {
        for (int point : cluster.points) {
            std::cout << cas[point]->serial << ": " << cas[point]->pos.x << std::endl;
        }
        std::cout << "\t ===> " << cluster.centroid.x << std::endl;
    }


    void Cluster::extractCAlphas() {
        for (auto& chain : tmdetVO.chains) {
            for (auto& residue : chain.residues) {
                const gemmi::Atom* ca_atom = residue.gemmi.get_ca();
                if (ca_atom) {
                    cas.push_back(ca_atom);
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
        //cluster.variance = sqrt(cluster.variance);
    }

    double Cluster::calculateWardDistance(const _cluster& c1, const _cluster& c2) {
        //_cluster mergedCluster = c1;
        //mergedCluster.points.insert(mergedCluster.points.end(), c2.points.begin(), c2.points.end());
        //calculateCentroid(mergedCluster);
        //calculateVariance(mergedCluster);
        //return mergedCluster.variance - (c1.variance + c2.variance);
        return c1.centroid.dist(c2.centroid);
    }

    int Cluster::calculateChainBrakes(const _cluster& c1, const _cluster& c2) {
        int chainBrakes = 0;
        for (int point1 : c1.points) {
            for (int point2: c2.points) {
                chainBrakes += (cas[point1]->pos.dist(cas[point2]->pos) > CA_DIST);
            }
        }
        return chainBrakes;
    }

    void Cluster::run() {
        extractCAlphas();
        int n = cas.size();
        std::vector<_cluster> clusters(n);

        for (auto i = 0; i < n; ++i) {
            clusters[i].points.push_back(i);
            clusters[i].centroid = cas[i]->pos;
            clusters[i].variance = 0.0;
        }

        while (clusters.size() > 1) {
            double minDist = 1e30;
            double minwdist;
            double mincdist;
            int c1Idx = -1, c2Idx = -1;

            for (auto i = 0; i < clusters.size(); ++i) {
                for (auto j = i + 1; j < clusters.size(); ++j) {
                    double wdist = calculateWardDistance(clusters[i], clusters[j]);
                    double cdist = calculateChainBrakes(clusters[i], clusters[j]);
                    if (wdist + cdist  < minDist) {
                        minDist = wdist + cdist;
                        minwdist = wdist;
                        mincdist = cdist;
                        c1Idx = i;
                        c2Idx = j;
                    }
                }
            }

            _cluster mergedCluster = clusters[c1Idx];
            mergedCluster.points.insert(mergedCluster.points.end(), clusters[c2Idx].points.begin(), clusters[c2Idx].points.end());
            calculateCentroid(mergedCluster);
            //calculateVariance(mergedCluster);

            clusters.erase(clusters.begin() + c2Idx);
            clusters.erase(clusters.begin() + c1Idx);
            clusters.push_back(mergedCluster);

            std::cout << "Merged clusters " << cas[c1Idx]->serial << " and " << cas[c2Idx]->serial << " with distance " << minDist << "(" << minwdist << "+" << mincdist << ")" << std::endl;
            showData(mergedCluster);
        }

        std::cout << "Final cluster: ";
        for (int point : clusters[0].points) {
            std::cout << point << " ";
        }
        std::cout << std::endl;
    }

    
}
