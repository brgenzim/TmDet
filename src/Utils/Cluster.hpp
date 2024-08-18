#ifndef __TMDET_UTILS_CLUSTER__
#define __TMDET_UTILS_CLUSTER__

#include <ValueObjects/TmdetStruct.hpp>

namespace Tmdet::Utils {

    class Cluster {
        private:
            Tmdet::ValueObjects::TmdetStruct& tmdetVO;

            
        public:
            Cluster(Tmdet::ValueObjects::TmdetStruct& _tmdetVO) : tmdetVO(_tmdetVO) {} ;
            ~Cluster() {};

            void run();
    };
}
#endif
