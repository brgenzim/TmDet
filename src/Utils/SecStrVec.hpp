// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <string>
#include <vector>
#include <any>
#include <unordered_map>
#include <gemmi/math.hpp>
#include <Types/SecStruct.hpp>
#include <VOs/SecStrVec.hpp>
#include <VOs/Protein.hpp>

namespace Tmdet::Utils {

    class SecStrVec {
        private:
            Tmdet::VOs::Protein& protein;

            void define();            
            bool getNextRegion(Tmdet::VOs::Chain& chain, int& begin, int& end) const;
            bool getNextNotUnkown(Tmdet::VOs::Chain& chain, int& begin) const;
            bool getNextSame(Tmdet::VOs::Chain& chain, const int& begin, int& end) const;
            Tmdet::VOs::SecStrVec getVector(const Tmdet::VOs::Chain& chain, int begin, int end) const;
            Tmdet::VOs::SecStrVec getAlphaVector(const Tmdet::VOs::Chain& chain, const int begin, const int end) const;
            gemmi::Vec3 getMeanPosition(const Tmdet::VOs::Chain& chain, const int pos) const;
            Tmdet::VOs::SecStrVec getBetaVector(const Tmdet::VOs::Chain& chain, int begin, int end) const;
            const gemmi::Atom* getCa(const Tmdet::VOs::Chain& chain, int begin, int end) const;
            void checkVectorsForSplitting();
            bool isStraight(const Tmdet::VOs::Chain& chain, const int beg, const int end, const Tmdet::Types::SecStruct& type);
            std::vector<Tmdet::VOs::SecStrVec> splitVector(const Tmdet::VOs::SecStrVec& vector);
            void getStraightVector(Tmdet::VOs::Chain& chain, int& beg, int& end,const Tmdet::Types::SecStruct& type);
            void checkVectorsForMerging();
            bool checkVectorForMerging(const Tmdet::VOs::SecStrVec& v1, const Tmdet::VOs::SecStrVec& v2);
            Tmdet::VOs::SecStrVec mergeVectors(const Tmdet::VOs::SecStrVec& v1, const Tmdet::VOs::SecStrVec& v2) const;

        public:
            explicit SecStrVec(Tmdet::VOs::Protein& protein) :
                protein(protein) {
                    define();
            }
            ~SecStrVec()=default;

            
    };
}
