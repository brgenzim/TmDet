#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <utility>
#include <Utils/Graph.hpp>

using namespace std;

namespace Tmdet::Utils {

    /**
     * add an edge to the graph
     * @param int u
     * @param int v
     * @return void
     */
    void Graph::addEdge(int u, int v) {
        if (u<0 || u>=V || v<0 || v>=V ) {
            return;
            //TODO: throw error
        }
        if (!edges[u][v]) {
            adj[u].push_back(v);
            adj[v].push_back(u);
            edges[u][v] = edges[v][u] =true;
        }
    }

    /**
     * Calculate the value for a given cluster
     * cutting in a given position
     * @param int clIdx
     * @param int cutPos
     * @return double
     */
    double Graph::graphValue(int clIdx, int cutPos) {
        double ret=0;
        int numLeft=0;
        int numRight=0;
        for(int v=0; v<V; v++) {
            if (clusters[v] == clIdx) {
                numLeft += (v<cutPos?1:0);
                numRight += (v<cutPos?0:1);
                for (const auto& u: adj[v]) {
                    if (clusters[u] == clIdx && v<cutPos && u>= cutPos) {
                        ret++;
                    }
                }
            }
        }
        if (numLeft>0||numRight>0) {
            ret /= (numLeft<numRight?numLeft:numRight);
        }
        ret += (numLeft<30||numRight<30?10000:0);
        return ret;
    }

    /**
     * Calculate the value for a given cluster
     * using two cutting positions
     * @param int clIdx
     * @param int beg
     * @param int end
     * @return double
     */
    double Graph::graphValue2(int clIdx, int beg, int end) {
        double ret=0;
        int numLeft=0;
        int numRight=0;
        for(int v=0; v<V; v++) {
            if (clusters[v] == clIdx) {
                numLeft += (v>=beg&&v<end?1:0);
                numRight += (v>=beg&&v<end?0:1);
                for (const auto& u: adj[v]) {
                    if (clusters[u] == clIdx && 
                            (
                                ((v<beg || v>=end) && u>=beg && u<end) ||
                                ((u<beg || u>=end) && v>=beg && v<end)
                            )
                        )
                        ret++;
                }
            }
        }
        if (numLeft>0&&numRight>0) {
            ret /= (numLeft<numRight?numLeft:numRight);
        }
        ret += (numLeft<30||numRight<30?10000:0);
        return ret;
    }

    /**
     * Scanning through polypeptid backbone
     * for determining best cutting point
     * @param int clIdx
     * @return std::vector<std::pair<double,int>>
     */
    std::vector<std::pair<double,int>> Graph::scan(int clIdx) {
        std::vector<std::pair<double,int>> ret(V);
        int j=0;
        for(int i=0; i<V; i++) {
            if (clusters[i] == clIdx) {
                ret[j] = std::pair<double,int>({graphValue(clIdx,i),i});
                j++;
            }
        }
        return ret;
    }

    /**
     * Smoothing scanned profile
     * @param std::vector<std::pair<double,int>> in
     * @return std::vector<std::pair<double,int>>
     */
    std::vector<std::pair<double,int>> Graph::smooth(vector<pair<double,int>> in) {
        std::vector<std::pair<double,int>> ret(V);
        for(unsigned int i=0; i<in.size(); i++) {
            int k=0;
            double q=0;
            for(int j=-10; j<11; j++) {
                if (j+(int)i>=0 && i+j<in.size()) {
                    q += in[i+j].first;
                    k++;
                }
            }
            q /= k;
            ret[i] = pair<double,int>({q,in[i].second});
        }
        return ret;
    }

    /**
     * Determining local minimum points in the smoothed profile
     * @param std::vector<std::pair<double,int>> in
     * @return std::vector<int>
     */
    vector<int> Graph::minPositions(vector<pair<double,int>> in) {
        std::vector<int> ret;
        ret.push_back(in[0].second);
        for(unsigned int i=1; i<in.size()-1; i++) {
            if (in[i].first<in[i-1].first && in[i].first<in[i+1].first) {
                ret.push_back(in[i].second);
            }
        }
        ret.push_back(in[in.size()-1].second);
        return ret;
    }

    /**
     * Determining the best cuts for a given cluster
     * @param int clIdx
     * @param double& minValue
     * @param int& bestBeg
     * @param int& bestEnd
     * @param int& bestCluster
     * @return void
     */
    void Graph::cut(int clIdx, double& minValue, int& bestBeg, int& bestEnd, int& bestCluster) {
        vector<pair<double,int>> p = scan(clIdx);
        vector<pair<double,int>> sp = smooth(p);
        vector<int> mps = minPositions(sp);

        for(unsigned int b=0; b<mps.size()-1; b++) {
            for(unsigned int e=b+1; e<mps.size(); e++) {
                double value = graphValue2(clIdx,mps[b],mps[e]);
                if (value < minValue) {
                    minValue = value;
                    bestBeg = mps[b];
                    bestEnd = mps[e];
                    bestCluster = clIdx;
                }
                //std::cout << b << ":" << e << "=" << value << std::endl;
            }
        }
        //std::cout << "Cut:" << clIdx << ":" << bestBeg << ":" << bestEnd << " = " << minValue << std::endl;
    }

    vector<int> Graph::optim() {
        while(true) {
            double minValue = MIN_DIFFERENT_LIMIT;
            int bestBeg = -1;
            int bestEnd = -1;
            int cl = -1;
            for(int i=0; i<numClusters; i++) {
                cut(i,minValue,bestBeg,bestEnd,cl);
            }
            //std::cout << cl << ":" << bestBeg << "-" << bestEnd << " =  "<< minValue << std::endl;
            if (cl == -1) {
                break;
            }
            if (minValue<MIN_DIFFERENT_LIMIT) {
                for (int i=bestBeg; i<bestEnd; i++) {
                    if (clusters[i] == cl) {
                        clusters[i] = numClusters;
                    }
                }
                numClusters++;
            }
        }
        return clusters;
    }
}

