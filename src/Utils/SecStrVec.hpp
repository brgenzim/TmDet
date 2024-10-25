#pragma once

#include <string>
#include <vector>
#include <any>
#include <unordered_map>
#include <Types/SecStruct.hpp>
#include <ValueObjects/Protein.hpp>
#include <gemmi/math.hpp>

namespace StructVO = Tmdet::ValueObjects;

//#define __SECSTRVEC_DBG 1

namespace Tmdet::Utils {

    struct _secStrVec {
        Tmdet::Types::SecStruct type;
        gemmi::Vec3 begin;
        gemmi::Vec3 end;
        int chainIdx;
        int begResIdx;
        int endResIdx;
    };

    class SecStrVec {
        private:
            std::vector<_secStrVec> vectors;
            Tmdet::ValueObjects::Protein& protein;

            void define();
            bool checkCross(_secStrVec& vec, Tmdet::ValueObjects::Membrane& membrane, int& numBoth, int& numUp, int& numDown);
            bool getNextRegion(StructVO::Chain& chain, int& begin, int& end);
            bool getNextNotUnkown(StructVO::Chain& chain, int& begin);
            bool getNextSame(StructVO::Chain& chain, int& begin, int& end);
            _secStrVec getVector(StructVO::Chain& chain, int begin, int end);
            _secStrVec getAlphaVector(StructVO::Chain& chain, int begin, int end);
            gemmi::Vec3 getMeanPosition(StructVO::Chain& chain, int pos);
            _secStrVec getBetaVector(StructVO::Chain& chain, int begin, int end);


        public:
            explicit SecStrVec(Tmdet::ValueObjects::Protein& protein) :
                protein(protein) {
                    define();
            }
            
            void numCrossingAlpha(Tmdet::ValueObjects::Membrane& membrane, int &numBoth, int &numUp, int &numDown);
            void numCrossingBeta(Tmdet::ValueObjects::Membrane& membrane, int &numBoth, int &numUp, int &numDown);
            

    };
}
