// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#include <iostream>
#include <array>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <sstream>
#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <Config.hpp>
#include <DTOs/SecStrVec.hpp>
#include <Helpers/Vector.hpp>
#include <System/Logger.hpp>
#include <Types/Residue.hpp>
#include <Utils/Surface.hpp>
#include <Utils/SecStrVec.hpp>
#include <VOs/Protein.hpp>
#include <VOs/Chain.hpp>
#include <VOs/SecStrVec.hpp>

using namespace std;

namespace Tmdet::Utils {
    
    void SecStrVec::define() {
        protein.secStrVecs.clear();
        for(auto& chain: protein.chains) {
            int begin = 0;
            int end = 0;
            while(getNextRegion(chain, begin, end)) {
                if (end - begin > 3) {
                    protein.secStrVecs.push_back(getVector(chain, begin, end - 1));
                }
                begin = end;
            }
        }
        checkVectorsForSplitting();
        if (protein.secStrVecs.size()>1) {
            checkVectorsForMerging();
        }
        for(unsigned long int i=0; auto& vector: protein.secStrVecs) {
            for (int j=vector.begResIdx; j<=vector.endResIdx; j++) {
                protein.chains[vector.chainIdx].residues[j].secStrVecIdx = (int)i;
            }
            i++;
        }
    }

    
    bool SecStrVec::getNextRegion(Tmdet::VOs::Chain& chain, int& begin, int& end) const {
        return (getNextNotUnkown(chain, begin) && getNextSame(chain, begin, end));
    }

    bool SecStrVec::getNextNotUnkown(Tmdet::VOs::Chain& chain, int& begin) const {
        while(begin < (int)chain.residues.size() 
            && !chain.residues[begin].selected
            && chain.residues[begin].ss == Tmdet::Types::SecStructType::U) {
            begin++;
        }
        return (begin < (int)chain.residues.size());
    }

    bool SecStrVec::getNextSame(Tmdet::VOs::Chain& chain, const int& begin, int& end) const {
        end = begin+1;
        while(end < (int)chain.residues.size() 
            && chain.residues[begin].selected
            && chain.residues[begin].ss.same(chain.residues[end].ss)
            && chain.orderDistance(end-1,end) == 1) {
            end++;
        }
        return true;
    }

    Tmdet::VOs::SecStrVec SecStrVec::getVector(const Tmdet::VOs::Chain& chain, int begin, int end) const {
        return (chain.residues[begin].ss.isAlpha()?
                    getAlphaVector(chain, begin, end):
                    getBetaVector(chain,begin,end));
    }

    Tmdet::VOs::SecStrVec SecStrVec::getAlphaVector(const Tmdet::VOs::Chain& chain, const int begin, const int end) const {
        gemmi::Vec3 b;
        gemmi::Vec3 e;
        if (chain.type == Tmdet::Types::ChainType::LOW_RES) {
            b = getMeanPosition(chain,begin);
            e = getMeanPosition(chain,end-2);
            auto v = (e -b).normalized();
            b -= 1 * v;
            e += 2 * v;
        }
        else {
            auto co = gemmi::Vec3(0,0,0);
            auto centre = gemmi::Vec3(0,0,0);
            bool first = true;
            gemmi::Vec3 firstCa;
            gemmi::Vec3 lastCa;
            for (int i=begin; i<=end; i++) {
                if (chain.residues[i].temp.contains("co")) {
                    co += any_cast<gemmi::Vec3>(chain.residues[i].temp.at("co"));
                    centre += any_cast<gemmi::Vec3>(chain.residues[i].temp.at("ca"));
                    if (first) {
                        firstCa = centre;
                        first = false;
                    }
                    lastCa = any_cast<gemmi::Vec3>(chain.residues[i].temp.at("ca"));
                }
            }
            centre /= (end-begin+1);
            co /= co.length();
            b = centre + (firstCa.dot(co) - centre.dot(co)) * co;
            e = centre + (lastCa.dot(co) - centre.dot(co)) * co;
        }
        return Tmdet::VOs::SecStrVec({
            Tmdet::Types::SecStructType::H,
            b, e, chain.idx, begin, end
        });
    }

    gemmi::Vec3 SecStrVec::getMeanPosition(const Tmdet::VOs::Chain& chain, const int pos) const {
        gemmi::Vec3 vec(0,0,0);
        int i = 0;
        int j = 0;
        while (i<3) {
            if (auto ca = chain.residues[pos+i].gemmi.get_ca(); ca != nullptr) {
                vec += ca->pos;
                j++;
            }
            i++;
        }
        if (j) {
            vec /= j;
        }
        return vec;
    }

    Tmdet::VOs::SecStrVec SecStrVec::getBetaVector(const Tmdet::VOs::Chain& chain, int begin, int end) const {
        auto ca = getCa(chain,begin,end);
        auto b = gemmi::Vec3(ca->pos);
        ca = getCa(chain,end,begin);
        auto e = gemmi::Vec3(ca->pos);
        auto v = (e -b).normalized();
        b -= v;
        e += v;
        return Tmdet::VOs::SecStrVec({
            Tmdet::Types::SecStructType::E,
            b, e, chain.idx, begin, end
        });
    }

    const gemmi::Atom* SecStrVec::getCa(const Tmdet::VOs::Chain& chain, int begin, int end) const {
        int direction = (end-begin) / std::abs(end-begin);
        auto i = begin;
        while((direction == 1 && i<end) || (direction == -1 && i>end)) {
            auto a = chain.residues[i].gemmi.get_ca();
            if (a == nullptr) {
                a = chain.residues[i].gemmi.get_c();
            }
            if (a == nullptr) {
                a = chain.residues[i].gemmi.get_n();
            }
            if (a != nullptr) {
                return a;
            }
            i += direction;
        }
        return (gemmi::Atom*)nullptr;
    }

    void SecStrVec::checkVectorsForSplitting() {
        auto vectors = protein.secStrVecs;
        protein.secStrVecs.clear();
        for(auto& vector: vectors) {
            if (vector.type.isAlpha()) {
                if (!isStraight(protein.chains[vector.chainIdx],vector.begResIdx,vector.endResIdx,vector.type)) {
                    for (auto part: splitVector(vector)) {
                        protein.secStrVecs.emplace_back(part);
                    }
                }
                else {
                    protein.secStrVecs.emplace_back(vector);
                }
            }
            else {
                protein.secStrVecs.emplace_back(vector);
            }
        }
    }

    bool SecStrVec::isStraight(const Tmdet::VOs::Chain& chain, const int beg, const int end, const Tmdet::Types::SecStruct& type) {
        if (end-beg<3) {
            return true;
        }
        auto vec = (getVector(chain,beg,end));
        for (int i=beg; i<=end; i++) {
            auto ca = chain.residues[i].getCa();
            if (ca !=nullptr && Tmdet::Helpers::Vector::distanceFromLine(vec.begin, vec.end, ca->pos) > (type.isAlpha()?5.5:16.0)) {
                return false;
            }
        }
        return true;
    }

    std::vector<Tmdet::VOs::SecStrVec> SecStrVec::splitVector(const Tmdet::VOs::SecStrVec& vector) {
        std::vector<Tmdet::VOs::SecStrVec> ret;
        Tmdet::VOs::Chain& chain = protein.chains[vector.chainIdx];
        int beg = vector.begResIdx;
        int end = vector.endResIdx;
        getStraightVector(chain,beg,end,vector.type);
        ret.push_back(getVector(chain,beg,end));
        if(!isStraight(chain,vector.begResIdx,beg,vector.type)) {
            for(auto part: splitVector(getVector(chain,vector.begResIdx,beg))) {
                ret.push_back(part);
            }
        }
        else {
            if (beg-vector.begResIdx > 3) {
                ret.push_back(getVector(chain,vector.begResIdx,beg-1));
            }
        }
        if(!isStraight(chain,end+1,vector.endResIdx,vector.type)) {
            for(auto part: splitVector(getVector(chain,end+1,vector.endResIdx))) {
                ret.push_back(part);
            }
        }
        else {
            if (vector.endResIdx-end > 3)
            ret.push_back(getVector(chain,end+1,vector.endResIdx));
        }
        return ret;
    }

    void SecStrVec::getStraightVector(Tmdet::VOs::Chain& chain, int& beg, int& end,const Tmdet::Types::SecStruct& type) {
        int maxLength = 0;
        int newBegin = beg;
        int newEnd = end;
        for(int i=beg; i<=end-4; i++) {
            for (int j=i+4; j<end; j++) {
                int length = j-i+1;
                if (isStraight(chain,i,j,type) && length>maxLength) {
                    maxLength = length;
                    newBegin = i;
                    newEnd = j;
                }
            }
        }
        beg = newBegin;
        end = newEnd;
    }

    void SecStrVec::checkVectorsForMerging() {
        auto vectors = protein.secStrVecs;
        protein.secStrVecs.clear();
        unsigned long int i = 1;
        bool needStore = false;
        auto prevVector = vectors[0];
        while (i<vectors.size()) {
            if (prevVector.type == vectors[i].type 
                //&& vectors[i].type.isAlpha()
                && checkVectorForMerging(prevVector,vectors[i])) {
                    auto mergedVector = mergeVectors(prevVector,vectors[i]);
                    protein.secStrVecs.emplace_back(mergedVector);
                    prevVector = mergedVector;
                    needStore = false;
            }
            else {
                protein.secStrVecs.emplace_back(prevVector);
                prevVector = vectors[i];
                needStore = true;
            }
            i++;
        }
        if (needStore) {
            protein.secStrVecs.emplace_back(prevVector);
        }
    }

    bool SecStrVec::checkVectorForMerging(const Tmdet::VOs::SecStrVec& v1, const Tmdet::VOs::SecStrVec& v2) {
        double d = v1.end.dist(v2.begin);
        double angle = Tmdet::Helpers::Vector::angle((v1.end-v1.begin),(v2.end-v2.begin));

        return (d < TMDET_SECSTRVEC_MERGE_DIST
                && v1.chainIdx==v2.chainIdx
                && isStraight(protein.chains[v1.chainIdx],v1.begResIdx,v2.endResIdx,v1.type)
                && angle < (v1.type.isAlpha()?20:30)
        );
    }

    Tmdet::VOs::SecStrVec SecStrVec::mergeVectors(const Tmdet::VOs::SecStrVec& v1, const Tmdet::VOs::SecStrVec& v2) const {
        return Tmdet::VOs::SecStrVec({
            v1.type,v1.begin,v2.end,v1.chainIdx,v1.begResIdx,v2.endResIdx
        });
    }
}
