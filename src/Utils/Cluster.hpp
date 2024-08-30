#ifndef __TMDET_UTILS_CLUSTER__
#define __TMDET_UTILS_CLUSTER__

#include <gemmi/model.hpp>
#include <ValueObjects/TmdetStruct.hpp>

#define CA_DIST 4.0

namespace Tmdet::Utils {

    struct _cluster {
        std::vector<int> points;
        gemmi::Position centroid;
        double variance;
    };

    class Cluster {
        private:
            Tmdet::ValueObjects::TmdetStruct& tmdetVO;
            std::vector<gemmi::Atom *> cas;
            
            void extractCAlphas();
            void calculateCentroid(_cluster& cluster);
            void calculateVariance(_cluster& cluster);
            double calculateWardDistance(const _cluster& c1, const _cluster& c2);
            int calculateChainBrakes(const _cluster& c1, const _cluster& c2);
            
        public:
            Cluster(Tmdet::ValueObjects::TmdetStruct& _tmdetVO) : tmdetVO(_tmdetVO) {} ;
            ~Cluster() {};

            void run();
    };
}
#endif
