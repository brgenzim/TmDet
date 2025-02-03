// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <vector>
#include <map>
#include <set>
#include <VOs/Chain.hpp>

namespace Tmdet::Utils {

    class Graph {
        private:
            Tmdet::VOs::Chain& chain;
            unsigned long V;
            std::vector<std::vector<unsigned int>> adj;
            std::vector<std::vector<bool>> edges;
            std::vector<unsigned int> clusters;
            std::vector<bool> isSS;
            unsigned int numClusters = 1;
            double ClusterCutLimit = 1.6;
            const int manyContacts = 100;

            double graphValue(unsigned int clIdx, unsigned int cutPos);
            double graphValue2(unsigned int clIdx, unsigned int beg, unsigned int end);
            std::vector<std::pair<double,unsigned int>> scan(unsigned int clIdx);
            std::vector<std::pair<double,unsigned int>> smooth(std::vector<std::pair<double,unsigned int>> in);
            std::vector<unsigned int> minPositions(std::vector<std::pair<double,unsigned int>> in);
            void cut(unsigned int clIdx, double& minValue, unsigned int& bestBeg, unsigned int& bestEnd, int& bestCluster);

        public:
            explicit Graph(Tmdet::VOs::Chain& chain, unsigned long V) : 
                chain(chain),
                V(V), 
                adj(V), 
                edges(V, std::vector<bool>(V, false)), 
                clusters(V,0),
                isSS(V,false) {}

            void addEdge(unsigned int u, unsigned int v, bool uss, bool vss);
            std::vector<unsigned int> optim();
            int getNumClusters() const { return numClusters;}
    };
}
