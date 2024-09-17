#ifndef __TMDET_UTILS_GRAPH__
#define __TMDET_UTILS_GRAPH__

#include <vector>
#include <map>
#include <set>

using namespace std;

namespace Tmdet::Utils {

#define MIN_DIFFERENT_LIMIT 1.6

    class Graph {
        private:
            int V;
            vector<vector<int>> adj;
            vector<vector<bool>> edges;
            vector<int> clusters;
            int numClusters = 1;

            double graphValue(int clIdx, int cutPos);
            double graphValue2(int clIdx, int beg, int end);
            vector<pair<double,int>> scan(int clIdx);
            vector<pair<double,int>> smooth(vector<pair<double,int>> in);
            vector<int> minPositions(vector<pair<double,int>> in);
            void cut(int clIdx, double& minValue, int& bestBeg, int& bestEnd, int& bestCluster);

        public:
            explicit Graph(int V) : V(V), adj(V), edges(V, vector<bool>(V, false)), clusters(V,0) {}

            void addEdge(int u, int v);
            vector<int> optim();
            int getNumClusters() const { return numClusters;}

        
    };
}

#endif