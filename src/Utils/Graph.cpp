// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

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

namespace Tmdet::Utils {

    /**
     * @brief add an edge to the graph
     * @param unsigned int u
     * @param unsigned int v
     * @return void
     */
    void Graph::addEdge(unsigned int u, unsigned int v) {
        if (u>=V || v>=V ) {
            return;
        }
        if (!edges[u][v]) {
            adj[u].push_back(v);
            adj[v].push_back(u);
            edges[u][v] = edges[v][u] =true;
        }
    }

    /**
     * @brief Calculate the value for a given cluster
     *        cutting in a given position
     * @param unsigned int clIdx
     * @param unsigned int cutPos
     * @return double
     */
    double Graph::graphValue(unsigned int clIdx, unsigned int cutPos) {
        double ret=0;
        int numLeft=0;
        int numRight=0;
        for(unsigned long v=0; v<V; v++) {
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
        ret += (numLeft<20||numRight<20?10000:0);
        return ret;
    }

    /**
     * @brief Calculate the value for a given cluster
     *        using two cutting positions
     * @param unsigned int clIdx
     * @param unsigned int beg
     * @param unsigned int end
     * @return double
     */
    double Graph::graphValue2(unsigned int clIdx, unsigned int beg, unsigned int end) {
        double ret=0;
        unsigned int numLeft=0;
        unsigned int numRight=0;
        for(unsigned int v=0; v<V; v++) {
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
        ret += (numLeft<20||numRight<20?10000:0);
        return ret;
    }

    /**
     * @brief Scanning through polypeptid backbone
     *        for determining best cutting point
     * @param unsigned int clIdx
     * @return std::vector<std::pair<double,unsigned int>>
     */
    std::vector<std::pair<double,unsigned int>> Graph::scan(unsigned int clIdx) {
        std::vector<std::pair<double,unsigned int>> ret(V);
        unsigned int j=0;
        for(unsigned int i=0; i<V; i++) {
            if (clusters[i] == clIdx) {
                ret[j] = std::pair<double,unsigned int>({graphValue(clIdx,i),i});
                j++;
            }
        }
        return ret;
    }

    /**
     * @brief Smoothing scanned profile
     * @param std::vector<std::pair<double,unsigned int>> in
     * @return std::vector<std::pair<double,unsigned int>>
     */
    std::vector<std::pair<double,unsigned int>> Graph::smooth(std::vector<std::pair<double,unsigned int>> in) {
        std::vector<std::pair<double,unsigned int>> ret(V);
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
            ret[i] = std::pair<double,unsigned int>({q,in[i].second});
        }
        return ret;
    }

    /**
     * @brief Determining local minimum points in the smoothed profile
     * @param std::vector<std::pair<double,unsigned int>> in
     * @return std::vector<unsigned int>
     */
    std::vector<unsigned int> Graph::minPositions(std::vector<std::pair<double,unsigned int>> in) {
        std::vector<unsigned int> ret;
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
     * @brief Determining the best cuts for a given cluster
     * @param unsigned int clIdx
     * @param double& minValue
     * @param unsigned int& bestBeg
     * @param unsigned int& bestEnd
     * @param unsigned int& bestCluster
     * @return void
     */
    void Graph::cut(unsigned int clIdx, double& minValue, unsigned int& bestBeg, unsigned int& bestEnd, int& bestCluster) {
        std::vector<std::pair<double,unsigned int>> p = scan(clIdx);
        std::vector<std::pair<double,unsigned int>> sp = smooth(p);
        std::vector<unsigned int> mps = minPositions(sp);

        for(unsigned int b=0; b<mps.size()-1; b++) {
            for(unsigned int e=b+1; e<mps.size(); e++) {
                double value = graphValue2(clIdx,mps[b],mps[e]);
                if (value < minValue) {
                    minValue = value;
                    bestBeg = mps[b];
                    bestEnd = mps[e];
                    bestCluster = clIdx;
                }
            }
        }
    }

    /**
     * @brief generates fragments
     * 
     * @return vector<unsigned int> 
     */
    std::vector<unsigned int> Graph::optim() {
        while(true) {
            double minValue = ClusterCutLimit;
            unsigned int bestBeg = 0;
            unsigned int bestEnd = 0;
            int cl = -1;
            for(unsigned int i=0; i<numClusters; i++) {
                cut(i,minValue,bestBeg,bestEnd,cl);
            }
            if (cl == -1) {
                break;
            }
            if (minValue < ClusterCutLimit) {
                for (unsigned int i=bestBeg; i<bestEnd; i++) {
                    if (clusters[i] == (unsigned int)cl) {
                        clusters[i] = numClusters;
                    }
                }
                numClusters++;
            }
        }
        return clusters;
    }
}

