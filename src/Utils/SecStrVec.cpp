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
#include <Helpers/Pymol.hpp>
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
        DEBUG_LOG("Processing SecStrVec::define()");
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
        checkAlphaVectorsForSplitting();
        if (protein.secStrVecs.size()>1) {
            checkAlphaVectorsForMerging();
        }
        for(unsigned long int i=0; auto& vector: protein.secStrVecs) {
            DEBUG_LOG("{}",Tmdet::DTOs::SecStrVec::toString(protein,vector));
            for (int j=vector.begResIdx; j<=vector.endResIdx; j++) {
                protein.chains[vector.chainIdx].residues[j].secStrVecIdx = (int)i;
            }
            i++;
        }
        DEBUG_LOG(" Processed SecStrVec::define(#vectors: {})",protein.secStrVecs.size());
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

    Tmdet::VOs::SecStrVec SecStrVec::getVector(Tmdet::VOs::Chain& chain, int begin, int end) const {
        DEBUG_LOG("SecStrVec to define: {}-{}-{}::{}",
            chain.id,chain.residues[begin].authId,chain.residues[end].authId,
            chain.residues[begin].ss.code);
        return (chain.residues[begin].ss.isAlpha()?
                    getAlphaVector(chain, begin, end):
                    getBetaVector(chain,begin,end));
    }

    Tmdet::VOs::SecStrVec SecStrVec::getAlphaVector(const Tmdet::VOs::Chain& chain, const int begin, const int end) const {
        /*auto b = getMeanPosition(chain,begin);
        auto e = getMeanPosition(chain,end-2);
        auto v = (e -b).normalized();
        b -= 1 * v;
        e += 2 * v;*/
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
        auto b = centre + (firstCa.dot(co) - centre.dot(co)) * co;
        auto e = centre + (lastCa.dot(co) - centre.dot(co)) * co;
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

    Tmdet::VOs::SecStrVec SecStrVec::getBetaVector(Tmdet::VOs::Chain& chain, int begin, int end) const {
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

    const gemmi::Atom* SecStrVec::getCa(Tmdet::VOs::Chain& chain, int begin, int end) const {
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

    void SecStrVec::checkAlphaVectorsForSplitting() {
        auto vectors = protein.secStrVecs;
        protein.secStrVecs.clear();
        for(auto& vector: vectors) {
            if (vector.type.isAlpha() && !isStraight(protein.chains[vector.chainIdx],vector.begResIdx,vector.endResIdx)) {
                for (auto part: splitAlphaVector(vector)) {
                    protein.secStrVecs.emplace_back(part);
                }
            }
            else {
                protein.secStrVecs.emplace_back(vector);
            }
        }
    }

    bool SecStrVec::isStraight(const Tmdet::VOs::Chain& chain, const int beg, const int end) {
        if (end-beg<3) {
            return true;
        }
        auto vec = getAlphaVector(chain,beg,end);
        for (int i=beg; i<=end; i++) {
            auto ca = chain.residues[i].getCa();
            if (ca !=nullptr && Tmdet::Helpers::Vector::distanceFromLine(vec.begin, vec.end, ca->pos) > 2.9) {
                return false;
            }
        }
        return true;
    }

    std::vector<Tmdet::VOs::SecStrVec> SecStrVec::splitAlphaVector(const Tmdet::VOs::SecStrVec& vector) {
        std::vector<Tmdet::VOs::SecStrVec> ret;
        Tmdet::VOs::Chain& chain = protein.chains[vector.chainIdx];
        int beg = vector.begResIdx;
        int end = vector.endResIdx;
        getStraightVector(chain,beg,end);
        ret.push_back(getAlphaVector(chain,beg,end));
        if(!isStraight(chain,vector.begResIdx,beg)) {
            for(auto part: splitAlphaVector(getAlphaVector(chain,vector.begResIdx,beg))) {
                ret.push_back(part);
            }
        }
        if(!isStraight(chain,end+1,vector.endResIdx)) {
            for(auto part: splitAlphaVector(getAlphaVector(chain,end+1,vector.endResIdx))) {
                ret.push_back(part);
            }
        }
        return ret;
    }

    void SecStrVec::getStraightVector(Tmdet::VOs::Chain& chain, int& beg, int& end) {
        int maxLength = 0;
        int newBegin = beg;
        int newEnd = end;
        for(int i=beg; i<=end-4; i++) {
            for (int j=i+4; j<end; j++) {
                int length = j-i+1;
                if (isStraight(chain,i,j) && length>maxLength) {
                    maxLength = length;
                    newBegin = i;
                    newEnd = j;
                }
            }
        }
        DEBUG_LOG("getStraight: {}:{}-{} -> {}-{}",chain.id,
            chain.residues[beg].authId,chain.residues[end].authId,
            chain.residues[newBegin].authId,chain.residues[newEnd].authId);
        beg = newBegin;
        end = newEnd;
    }

    void SecStrVec::checkAlphaVectorsForMerging() {
        auto vectors = protein.secStrVecs;
        protein.secStrVecs.clear();
        unsigned long int i = 0;
        unsigned long step = 1;
        while (i<vectors.size()-1) {
            DEBUG_LOG("checkAlphaVectorsForMerging: {} and {}",
                Tmdet::DTOs::SecStrVec::toString(protein,vectors[i]),
                Tmdet::DTOs::SecStrVec::toString(protein,vectors[i+1]));
            if (vectors[i].type.isAlpha() && vectors[i+1].type.isAlpha() && checkAlphaVectorForMerging(vectors[i],vectors[i+1])) {
                protein.secStrVecs.emplace_back(mergeVectors(vectors[i],vectors[i+1]));
                step = 2;
            }
            else {
                protein.secStrVecs.emplace_back(vectors[i]);
                step = 1;
            }
            i+=step;
        }
        if (step == 1) {
            protein.secStrVecs.emplace_back(vectors[i]);
        }
    }

    bool SecStrVec::checkAlphaVectorForMerging(const Tmdet::VOs::SecStrVec& v1, const Tmdet::VOs::SecStrVec& v2) {
        double d = v1.end.dist(v2.begin);
        DEBUG_LOG("check merge {} {}-{}: {}",
            protein.chains[v1.chainIdx].id,v1.endResIdx,v2.begResIdx,d);
        return (d < TMDET_SECSTRVEC_MERGE_DIST
                && v1.chainIdx==v2.chainIdx
                && isStraight(protein.chains[v1.chainIdx],v1.begResIdx,v2.endResIdx)
        );
    }

    Tmdet::VOs::SecStrVec SecStrVec::mergeVectors(const Tmdet::VOs::SecStrVec& v1, const Tmdet::VOs::SecStrVec& v2) const {
        DEBUG_LOG("merge");
        return Tmdet::VOs::SecStrVec({
            v1.type,v1.begin,v2.end,v1.chainIdx,v1.begResIdx,v2.endResIdx
        });
    }
}
