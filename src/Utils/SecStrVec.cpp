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
                if (end - begin > 4) {
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
            DEBUG_LOG("{}",Tmdet::DTOs::SecStrVec::toString(vector));
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
            && chain.residues[begin].ss.same(chain.residues[end].ss)) {
            end++;
        }
        return true;
    }

    Tmdet::VOs::SecStrVec SecStrVec::getVector(Tmdet::VOs::Chain& chain, int begin, int end) const {
        DEBUG_LOG("SecStrVec defined: {}-{}-{}",chain.id,begin,end);
        return (chain.residues[begin].ss.isAlpha()?
                    getAlphaVector(chain, begin, end):
                    getBetaVector(chain,begin,end));
    }

    Tmdet::VOs::SecStrVec SecStrVec::getAlphaVector(Tmdet::VOs::Chain& chain, int begin, int end) const {
        auto b = getMeanPosition(chain,begin);
        auto e = getMeanPosition(chain,end-3);
        auto v = (e -b).normalized();
        b -= 2 * v;
        e += 4 * v;
        return Tmdet::VOs::SecStrVec({
            Tmdet::Types::SecStructType::H,
            b, e, chain.idx, begin, end
        });
    }

    gemmi::Vec3 SecStrVec::getMeanPosition(Tmdet::VOs::Chain& chain, int pos) const {
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
            if (vector.type.isAlpha() && !checkAlphaVectorForSplitting(vector)) {
                for (auto part: splitAlphaVector(vector)) {
                    protein.secStrVecs.emplace_back(part);
                }
            }
            else {
                protein.secStrVecs.emplace_back(vector);
            }
        }
    }

    bool SecStrVec::checkAlphaVectorForSplitting(const Tmdet::VOs::SecStrVec& vector) {
        if (vector.endResIdx - vector.begResIdx < 15) {
            return true;
        }
        for (int i = vector.begResIdx; i<= vector.endResIdx; i++) {
            if (auto ca = protein.chains[vector.chainIdx].residues[i].gemmi.get_ca(); ca != nullptr 
                && Tmdet::Helpers::Vector::distanceFromLine(vector.begin, vector.end, ca->pos) > 7.0) {
                    return false;
            }
        }
        return true;
    }

    std::vector<Tmdet::VOs::SecStrVec> SecStrVec::splitAlphaVector(const Tmdet::VOs::SecStrVec& vector) {
        std::vector<Tmdet::VOs::SecStrVec> ret;
        int begResIdx = vector.begResIdx;
        Tmdet::VOs::SecStrVec straightVec;
        while(begResIdx < vector.endResIdx) {
            if (getStraightVector(vector.chainIdx,begResIdx, vector.endResIdx, straightVec)) {
                ret.emplace_back(straightVec);
            }
            begResIdx = straightVec.endResIdx + 1;
        }
        return ret;
    }

    bool SecStrVec::getStraightVector(int chainIdx, int begResIdx, int endResIdxAll, Tmdet::VOs::SecStrVec& vec) {
        int p1=begResIdx;
        bool ok = true;
        while (ok) {
            if (abs(getCaDist(chainIdx,p1) - 6.2) < 1 && p1+5<endResIdxAll) {
                p1++; 
            }
            else {
                ok = false;
            }
        }
        DEBUG_LOG("getStraightVector: {} {} {}",begResIdx,p1,endResIdxAll);
        if (begResIdx+1<p1) {
            vec = getAlphaVector(protein.chains[chainIdx],begResIdx,p1+4);
        }
        else {
            vec.endResIdx = p1;
        }
        return (p1-begResIdx>2);
    }

    double SecStrVec::getCaDist(int chainIdx, int resIdx) {
        if (resIdx+4>=(int)protein.chains[chainIdx].residues.size()) {
            return 1000;
        }
        auto a1 = protein.chains[chainIdx].residues[resIdx].gemmi.get_ca();
        auto a2 = protein.chains[chainIdx].residues[resIdx+4].gemmi.get_ca();
        return (a1!=nullptr&&a2!=nullptr?a1->pos.dist(a2->pos):-1000);
    }

    void SecStrVec::checkAlphaVectorsForMerging() {
        auto vectors = protein.secStrVecs;
        protein.secStrVecs.clear();
        unsigned long int i = 0;
        unsigned long step = 1;
        while (i<vectors.size()-1) {
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

    bool SecStrVec::checkAlphaVectorForMerging(const Tmdet::VOs::SecStrVec& v1, const Tmdet::VOs::SecStrVec& v2) const {
        double d = v1.end.dist(v2.begin);
        auto vv1 = v1.end - v1.begin;
        auto vv2 = v2.end - v2.begin;
        double a = Tmdet::Helpers::Vector::angle(vv1,vv2);
        DEBUG_LOG("check merge {} {}-{}: {} {}",
            protein.chains[v1.chainIdx].id,v1.endResIdx,v2.begResIdx,d,a);
        return (d < TMDET_SECSTRVEC_MERGE_DIST
                && a < TMDET_SECSTRVEC_MERGE_ANGLE
                && v1.chainIdx==v2.chainIdx
        );
    }

    Tmdet::VOs::SecStrVec SecStrVec::mergeVectors(const Tmdet::VOs::SecStrVec& v1, const Tmdet::VOs::SecStrVec& v2) const {
        DEBUG_LOG("merge");
        return Tmdet::VOs::SecStrVec({
            v1.type,v1.begin,v2.end,v1.chainIdx,v1.begResIdx,v2.endResIdx
        });
    }
}
