#ifndef __TMDET_UTILS_SECSTRVEC__
#define __TMDET_UTILS_SECSTRVEC__

#include <string>
#include <vector>
#include <any>
#include <unordered_map>
#include <Types/SecStruct.hpp>
#include <gemmi/math.hpp>

using namespace std;

namespace Tmdet::Utils {

    struct _secStrVec {
        Tmdet::Types::SecStruct type;
        gemmi::Vec3 begin;
        gemmi::Vec3 end;
    };

    class SecStrVec {
        private:
            vector<_secStrVec> vectors;
            bool ifCross(_secStrVec& vec, Tmdet::ValueObjects::Membrane& membraneVO);

        public:
            SecStrVec() {}
            ~SecStrVec() {}

            void define(Tmdet::ValueObjects::TmdetStruct& tmdetVO);
            void numCross(Tmdet::ValueObjects::Membrane& membraneVO, int &numBoth, int &numUp, int &numDown);
    };
}

#endif