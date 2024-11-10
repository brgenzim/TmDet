#pragma once

#include <string>
#include <vector>
#include <any>
#include <unordered_map>
#include <gemmi/math.hpp>
#include <Types/SecStruct.hpp>
#include <ValueObjects/SecStrVec.hpp>
#include <ValueObjects/Protein.hpp>

namespace Tmdet::Utils {

    class SecStrVec {
        private:
            Tmdet::ValueObjects::Protein& protein;

            void define();            
            bool getNextRegion(Tmdet::ValueObjects::Chain& chain, int& begin, int& end) const;
            bool getNextNotUnkown(Tmdet::ValueObjects::Chain& chain, int& begin) const;
            bool getNextSame(Tmdet::ValueObjects::Chain& chain, const int& begin, int& end) const;
            Tmdet::ValueObjects::SecStrVec getVector(Tmdet::ValueObjects::Chain& chain, int begin, int end) const;
            Tmdet::ValueObjects::SecStrVec getAlphaVector(Tmdet::ValueObjects::Chain& chain, int begin, int end) const;
            gemmi::Vec3 getMeanPosition(Tmdet::ValueObjects::Chain& chain, int pos) const;
            Tmdet::ValueObjects::SecStrVec getBetaVector(Tmdet::ValueObjects::Chain& chain, int begin, int end) const;
            const gemmi::Atom* getCa(Tmdet::ValueObjects::Chain& chain, int begin, int end) const;
            void checkAlphaVectorsForSplitting();
            bool checkAlphaVectorForSplitting(const Tmdet::ValueObjects::SecStrVec& vector);
            std::vector<Tmdet::ValueObjects::SecStrVec> splitAlphaVector(const Tmdet::ValueObjects::SecStrVec& vector);
            bool getStraightVector(int chainIdx, int begResIdx, int endResIdxAll, Tmdet::ValueObjects::SecStrVec& vec);
            double getCaDist(int chainIdx, int resIdx);
            void checkAlphaVectorsForMerging();
            bool checkAlphaVectorForMerging(const Tmdet::ValueObjects::SecStrVec& v1, const Tmdet::ValueObjects::SecStrVec& v2) const;
            Tmdet::ValueObjects::SecStrVec mergeVectors(const Tmdet::ValueObjects::SecStrVec& v1, const Tmdet::ValueObjects::SecStrVec& v2) const;

        public:
            explicit SecStrVec(Tmdet::ValueObjects::Protein& protein) :
                protein(protein) {
                    define();
            }
            ~SecStrVec()=default;

            
    };
}
