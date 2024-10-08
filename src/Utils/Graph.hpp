#pragma once

#include <vector>
#include <map>
#include <set>

namespace Tmdet::Utils {

    class Graph {
        private:
            unsigned long V;
            std::vector<std::vector<unsigned int>> adj;
            std::vector<std::vector<bool>> edges;
            std::vector<unsigned int> clusters;
            unsigned int numClusters = 1;
            double ClusterCutLimit = 1.6;

            double graphValue(unsigned int clIdx, unsigned int cutPos);
            double graphValue2(unsigned int clIdx, unsigned int beg, unsigned int end);
            std::vector<std::pair<double,unsigned int>> scan(unsigned int clIdx);
            std::vector<std::pair<double,unsigned int>> smooth(std::vector<std::pair<double,unsigned int>> in);
            std::vector<unsigned int> minPositions(std::vector<std::pair<double,unsigned int>> in);
            void cut(unsigned int clIdx, double& minValue, unsigned int& bestBeg, unsigned int& bestEnd, int& bestCluster);

        public:
            explicit Graph(unsigned long V) : V(V), adj(V), edges(V, std::vector<bool>(V, false)), clusters(V,0) {}

            void addEdge(unsigned int u, unsigned int v);
            std::vector<unsigned int> optim();
            int getNumClusters() const { return numClusters;}
    };
}
