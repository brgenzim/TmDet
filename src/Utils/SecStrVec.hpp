#ifndef __TMDET_UTILS_SECSTRVEC__
#define __TMDET_UTILS_SECSTRVEC__

#include <string>
#include <vector>
#include <any>
#include <unordered_map>
#include <Types/SecStruct.hpp>
#include <ValueObjects/TmdetStruct.hpp>
#include <gemmi/math.hpp>

using namespace std;

#define __SECSTRVEC_DBG 1

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
            bool getNextRegion(Tmdet::ValueObjects::Chain& chain, int& begin, int& end);
            bool getNextNotUnkown(Tmdet::ValueObjects::Chain& chain, int& begin);
            bool getNextSame(Tmdet::ValueObjects::Chain& chain, int& begin, int& end);
            _secStrVec getVector(Tmdet::ValueObjects::Chain& chain, int begin, int end);
            _secStrVec getAlphaVector(Tmdet::ValueObjects::Chain& chain, int begin, int end);
            gemmi::Vec3 getMeanPosition(Tmdet::ValueObjects::Chain& chain, int pos);
            _secStrVec getBetaVector(Tmdet::ValueObjects::Chain& chain, int begin, int end);


        public:
            SecStrVec() {}
            ~SecStrVec() {}

            void define(Tmdet::ValueObjects::TmdetStruct& tmdetVO);
            void numCross(Tmdet::ValueObjects::Membrane& membraneVO, int &numBoth, int &numUp, int &numDown);
    };
}

#endif
