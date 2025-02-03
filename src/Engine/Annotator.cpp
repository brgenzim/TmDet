// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#include <any>
#include <memory>
#include <gemmi/model.hpp>
#include <Config.hpp>
#include <Helpers/Vector.hpp>
#include <Engine/Annotator.hpp>
#include <Engine/RegionHandler.hpp>
#include <Engine/CurvedSideDetector.hpp>
#include <Engine/PlaneSideDetector.hpp>
#include <Engine/BetaAnnotator.hpp>
#include <Helpers/String.hpp>
#include <System/Logger.hpp>
#include <Types/Chain.hpp>
#include <Types/Region.hpp>
#include <Utils/SecStrVec.hpp>
#include <Types/Membrane.hpp>

#define REGTYPE(r) any_cast<Tmdet::Types::Region>(r.temp.at("type"))
//#define REGTTYPE(r) any_cast<Tmdet::Types::Region>(r.temp.at("ttype"))
#define REGZTYPE(r) any_cast<Tmdet::Types::Region>(r.temp.at("ztype"))
#define REGTYPEW(r,w) any_cast<Tmdet::Types::Region>(r.temp.at(w))
#define REGZ(r) any_cast<double>(r.temp.at("z"))
#define REGHZ(r) any_cast<double>(r.temp.at("hz"))
#define REGDIR(r) any_cast<double>(r.temp.at("direction"))

namespace Tmdet::Engine {

    void Annotator::run() {
        DEBUG_LOG("Processing Annotator::run()");
        std::unique_ptr<Tmdet::Engine::SideDetector> sideDetector;
        protein.eachChain(
            [&](Tmdet::VOs::Chain& chain) -> void {
                chain.regions.clear();
            }
        );
        if (protein.membranes[0].type.isPlane()) {
            DEBUG_LOG("Plane sidedetector");
            sideDetector = std::make_unique<Tmdet::Engine::PlaneSideDetector>(protein);
        }
        else {
            DEBUG_LOG("Curved sidedetector");
            sideDetector = std::make_unique<Tmdet::Engine::CurvedSideDetector>(protein);
        }
        smoothRegions("type");
        detectLoops();        
        auto betaAnnotator = Tmdet::Engine::BetaAnnotator(protein,regionHandler);
        
        setChainsType();
        annotateChains();
        protein.eachSelectedChain(
            [&](Tmdet::VOs::Chain& chain) -> void {
                if (chain.type.isBeta()) {
                    betaAnnotator.detectBarrelInside(chain);
                }
            }
        );
        detectInterfacialHelices();
        int limit = protein.type.isBeta()?50:10;
        if (args.getValueAsBool("fr")) {
            limit = 0;
        }
        if (int nm = regionHandler.finalize<Tmdet::Types::Region>(); nm>limit) {
            DEBUG_LOG("Too many unhandled membrane segments: {}",nm);
            protein.notTransmembrane();
        }
        else {
            regionHandler.store<Tmdet::Types::Region>();
            finalCheck();
        }
        if (protein.tmp) {
            setMembraneSize();
        }
        DEBUG_LOG(" Processed Annotator::run({})",(protein.tmp?regionHandler.toString("type"):"not tmp"));
    }

    void Annotator::setChainsType() {
        DEBUG_LOG("Processing Annotator::setChainsType()");
        protein.eachChain(
            [&](Tmdet::VOs::Chain& chain) -> void {
                chain.type = (chain.selected?
                    Tmdet::Types::ChainType::NON_TM:
                    Tmdet::Types::ChainType::NOT_SELECTED);
            }
        );
        for(const auto& ssVec: protein.secStrVecs) {
            if (ssVec.barrelIdx != -1) {
                protein.chains[ssVec.chainIdx].type = Tmdet::Types::ChainType::BETA;
                protein.type = Tmdet::Types::ProteinType::TM_BETA;
                DEBUG_LOG("ChainType Beta: {}",protein.chains[ssVec.chainIdx].id);
            }
        }
        if (std::string bi = args.getValueAsString("bi"); bi != "" && protein.type.isBeta()) {
            for(auto chainId: Tmdet::Helpers::String::explode(",",bi)) {
                if (int chainIdx = protein.searchChainById(chainId); chainIdx != -1) {
                    protein.chains[chainIdx].type = Tmdet::Types::ChainType::BETA;
                }
                else {
                    WARN_LOG("Could not find chain: {}",bi);
                }
            }
        }
        protein.eachChain(
            [&](Tmdet::VOs::Chain& chain) -> void {
                int alpha=0;
                if (chain.type == Tmdet::Types::ChainType::NON_TM && (chain.isTmp || chain.length < 50)) {
                    for(const auto& residue: chain.residues) {
                        if (residue.selected) {
                            if (REGTYPE(residue) == Tmdet::Types::RegionType::MEMB) {
                                if (residue.ss.isAlpha()) {
                                    alpha++;
                                }
                            }
                        }
                    }
                    if (alpha>6) {
                        chain.type = Tmdet::Types::ChainType::ALPHA;
                        protein.type = (protein.type.isBeta()?
                                Tmdet::Types::ProteinType::TM_MIXED:
                                Tmdet::Types::ProteinType::TM_ALPHA);
                    }
                }
                DEBUG_LOG("ChainTypeAlpha: {} ({})",chain.id,alpha);
            }
        );
        DEBUG_LOG("Processed Annotator::setChainsType()");
    }

    void Annotator::annotateChains() {
        protein.eachSelectedChain(
            [&](Tmdet::VOs::Chain& chain) -> void {
                if (chain.type.isAlpha()) {
                    detectReEntrantLoops(chain);
                    detectTransmembraneHelices(chain);
                }
            }
        );
        DEBUG_LOG(" Processed: Annotator::annotateChains:\n {}",regionHandler.toString("type"));
    }

    void Annotator::smoothRegions(std::string what) {
        DEBUG_LOG("Processing: Annotator::smoothRegions");
        protein.eachSelectedChain(
            [&](Tmdet::VOs::Chain& chain) -> void {
                int beg = 0;
                int end = 0;
                while(regionHandler.getNext<Tmdet::Types::Region>(chain,beg,end,what)) {
                    if (REGTYPEW(chain.residues[beg],what) == Tmdet::Types::RegionType::MEMB && end-beg < 3) {
                        regionHandler.replace(chain,beg,end-1,REGZTYPE(chain.residues[beg]),what);
                    }
                    beg=end;
                }
                auto regions = regionHandler.getAll<Tmdet::Types::Region>(chain,what);
                for (unsigned int i=0; i<regions.size(); i++) {
                    if ( regions[i].type.isNotMembrane()
                            && regions[i].end-regions[i].beg < 3
                            && ((i>0 && regions[i-1].type.isNotAnnotatedMembrane() && regions[i-1].end - regions[i-1].beg < 4)
                            || (i<regions.size()-1 && regions[i+1].type.isNotAnnotatedMembrane() && regions[i+1].end - regions[i+1].beg < 4))) {
                        regionHandler.replace(chain,regions[i].beg,regions[i].end,Tmdet::Types::RegionType::MEMB,what);
                    }
                }
            }
        );
        DEBUG_LOG("Processed: Annotator::smoothRegions ()");
    }

    void Annotator::detectLoops() {
        protein.eachSelectedChain(
            [&](Tmdet::VOs::Chain& chain) -> void {
                int beg=0;
                int end=0;
                while(regionHandler.getNext<Tmdet::Types::Region>(chain,beg,end,"type")) {
                    if (REGTYPE(chain.residues[beg]) == Tmdet::Types::RegionType::MEMB) {
                        detectLoop(chain,beg,end-1);
                    }
                    beg=end;
                }
            }
        );
    }

    void Annotator::detectLoop(Tmdet::VOs::Chain& chain, int beg, int end) {
        DEBUG_LOG("detectLoop: {}:{}-{}",chain.id, chain.residues[beg].authId,chain.residues[end].authId);
        for(int i=beg+5; i<=end-5; i++) {
            if ( (std::abs(REGZ(chain.residues[i])) > 5 || REGHZ(chain.residues[i]) < 9.0)
                    && ((REGZ(chain.residues[i]) < REGZ(chain.residues[i-1])
                        && REGZ(chain.residues[i]) < REGZ(chain.residues[i+1]))
                        ||
                        (REGZ(chain.residues[i]) > REGZ(chain.residues[i-1])
                        && REGZ(chain.residues[i]) > REGZ(chain.residues[i+1]))
                    )
                    && (maxDist(chain,i,i-2,beg,-1) > 4 && maxDist(chain,i,i+2,end,1) > 4)
             ) {
                chain.residues[i].temp.at("type") = chain.residues[i].temp.at("ztype");
                DEBUG_LOG("detectLoop: isLoop: {}:{}",chain.id, chain.residues[i].authId);
            }
        }
    }

    double Annotator::maxDist(Tmdet::VOs::Chain& chain, int pos, int beg, int end, int dir) {
        double maxDist = 0;
        int i;
        for (i=beg; (i!=(end+dir) && REGDIR(chain.residues[beg]) * REGDIR(chain.residues[i]) > -10); i+=dir) {
            double dist = (REGZ(chain.residues[pos]) > 0?
                            REGZ(chain.residues[pos]) - REGZ(chain.residues[i]) :
                            REGZ(chain.residues[i]) - REGZ(chain.residues[pos]));
            if (dist>maxDist) {
                maxDist = dist;
            }
        }
        DEBUG_LOG("maxDist: {}:{} {}-{} ==> {} {}",
            chain.id,chain.residues[pos].authId,chain.residues[beg].authId,chain.residues[end].authId,i,maxDist);
        return (std::abs(pos-i)>3?maxDist:0);
    }

    bool Annotator::hasOtherSide(Tmdet::VOs::Chain& chain, int pos, int beg, int end) {
        for (int i=beg; i<=end; i++) {
            if (REGZ(chain.residues[pos])*REGZ(chain.residues[i])<0) {
                return true;
            }
        }
        return false;
    }

    void Annotator::detectInterfacialHelices() {
        
        float hfLimit = args.getValueAsFloat("hml");
        float avgSurface = args.getValueAsFloat("ias");
        ifhAngleLimit = args.getValueAsFloat("ian");
        DEBUG_LOG("Processing Annotator::detectInterfacialHelices(hfm:{}, as:{})",hfLimit,avgSurface);
        for(auto& membrane: protein.membranes) {
            auto alphaVecs = getParallelAlphas(membrane);
            if (!alphaVecs.empty() ) {
                for(const auto& vector: alphaVecs) {
                    if (vector.endResIdx -vector.begResIdx>2
                        && protein.chains[vector.chainIdx].residues[vector.begResIdx].selected
                        && protein.chains[vector.chainIdx].residues[vector.endResIdx].selected
                        && averageSurface(protein.chains[vector.chainIdx],vector.begResIdx,vector.endResIdx) > avgSurface
                        && sameSide(protein.chains[vector.chainIdx],vector.begResIdx,vector.endResIdx)
                        && hydrophocityMomentum(protein.chains[vector.chainIdx],vector.begResIdx,vector.endResIdx) > hfLimit) {
                            for (int i=vector.begResIdx; i<=vector.endResIdx; i++) {
                                if (protein.chains[vector.chainIdx].residues[i].selected
                                    && !REGTYPE(protein.chains[vector.chainIdx].residues[i]).isAnnotatedMembraneType()) {
                                    protein.chains[vector.chainIdx].residues[i].temp["type"] = std::any(Tmdet::Types::RegionType::IFH);
                                }
                            }
                    }
                }
            }
            DEBUG_LOG(" #alphaVecs: {}",alphaVecs.size());
        }
        DEBUG_LOG(" Processed Annotator::detectInterfacialHelices()");
    }

    void Annotator::detectReEntrantLoops(Tmdet::VOs::Chain& chain) {
        DEBUG_LOG("Processing Annotator::detectReEntrantLoops()");
        int begin = 0;
        int end = 0;
        while(regionHandler.getNext<Tmdet::Types::Region>(chain,begin,end,"type")) {
            end--;
            int numHelix;
            if (REGTYPE(chain.residues[begin]) == Tmdet::Types::RegionType::MEMB) {
                DEBUG_LOG("detectReEntrantLoops: {} {} {} zBeg:{} zEnd:{} sameSide:{}",
                    chain.id,chain.residues[begin].authId,chain.residues[end].authId,
                    REGZ(chain.residues[begin]),REGZ(chain.residues[end]),
                    sameSide(chain,begin,end)
                );
                if (sameSide(chain,begin,end)
                    && REGHZ(chain.residues[begin]) < 10.0 
                    && REGHZ(chain.residues[end]) < 10.0
                    && hasHelixTurnLoop(chain,begin,end,numHelix)) {
                        DEBUG_LOG("Loop found: {} {} {} {}",chain.id,chain.residues[begin].authId,
                                chain.residues[end].authId,numHelix);
                        if (numHelix==1) {
                            regionHandler.replace(chain,begin,end,Tmdet::Types::RegionType::LOOP);
                        }
                        else {
                            regionHandler.replace(chain,begin,end,Tmdet::Types::RegionType::TWO_H_LOOP);
                        }
                }
                else if (begin>0
                        && end<chain.length-1
                        && (chain.residues[begin].labelId - chain.residues[begin-1].labelId > 1
                            || chain.residues[end+1].labelId - chain.residues[end].labelId > 1)
                        && (end-begin > 6)
                        && (end-begin < 20)
                        && sameSide(chain,begin,end)
                    ) {
                        DEBUG_LOG("Broken loop found: {} {} {}",chain.id,chain.residues[begin].authId,
                                chain.residues[end].authId);
                        regionHandler.replace(chain,begin,end,Tmdet::Types::RegionType::LOOP);
                }
                else {
                    DEBUG_LOG("Not loop: {} {} {}",chain.id,begin,end);
                }
            }
            begin = end+1;
        }
        DEBUG_LOG(" Processed Annotator::detectReEntrantLoops()");
    }
        
    bool Annotator::hasHelixTurnLoop(Tmdet::VOs::Chain& chain, int begin, int end, int& numHelix) {
        int numNoSS = 0;
        std::unordered_map<int,int> helixCounts;
        double maxHz = 0;
        double hzBeg = REGHZ(chain.residues[begin]);
        double hzEnd = REGHZ(chain.residues[end]);
        double hz = (hzBeg>hzEnd?hzBeg:hzEnd);
        for(int i=begin; i<=end; i++) {
            int vecIdx = chain.residues[i].secStrVecIdx;
            if (vecIdx ==-1) {
                numNoSS++;
            }
            else if (protein.secStrVecs[chain.residues[i].secStrVecIdx].type.isAlpha()) {
                if (helixCounts.contains(vecIdx)) {
                    helixCounts[vecIdx]++;
                }
                else {
                    helixCounts.try_emplace(vecIdx,0);
                }

            }
            if (REGHZ(chain.residues[i]) > maxHz) {
                maxHz = REGHZ(chain.residues[i]);
            }
        }
        int maxCount=0;
        double maxPercent=0;
        numHelix = 0;
        for(auto& [idx,count]: helixCounts) {
            double percent = 100.0 * count / (protein.secStrVecs[idx].endResIdx - protein.secStrVecs[idx].begResIdx + 1);
            if (count > maxCount) {
                maxCount = count;
                maxPercent = percent;
            }
            numHelix += (percent>25);
        }
        DEBUG_LOG("hasOneHelix: {} {} {}: {}::{}% {} hzDiff:{} numHelix:{} return:{}",chain.id,
            chain.residues[begin].authId,chain.residues[end].authId,
            maxCount,maxPercent,numNoSS,maxHz-hz,numHelix,
            (maxPercent>25 &&  maxHz-hz > 3));
        return (maxPercent>25 &&  maxHz-hz > 3);
    }

    void Annotator::detectTransmembraneHelices(Tmdet::VOs::Chain& chain) {
        DEBUG_LOG("Processing Annotator::detectTransmembraneHelices()");
        int begin = 0;
        int end = 0;
        int lengthLimit = (args.getValueAsBool("fr")?14:10);
        while(regionHandler.getNext<Tmdet::Types::Region>(chain,begin,end,"type")) {
            DEBUG_LOG("detectTMH: {}:{}-{}",chain.id,chain.residues[begin].authId,chain.residues[end-1].authId);
            if (REGTYPE(chain.residues[begin]).isNotAnnotatedMembrane()
                && helixContent(chain,begin,end) > 0.3
                && (end-begin > lengthLimit || ((begin < 2 || end > chain.length -2) && end-begin > 6))
                && REGZ(chain.residues[begin]) * REGZ(chain.residues[end-1]) < 0) {
                    regionHandler.replace(chain,begin,end-1,Tmdet::Types::RegionType::HELIX);
            }
            begin = end;
        }
        DEBUG_LOG(" Processed Annotator::detectTransmembraneHelices()");
    }

    double Annotator::helixContent(Tmdet::VOs::Chain& chain, int beg, int end) {
        double ret = 0.0;
        for (int i=beg; i<=end; i++) {
            if (chain.residues[i].ss.isAlpha()) {
                ret++;
            }
        }
        ret /= (end-beg+1);
        return ret;
    }
    bool Annotator::sameSide(Tmdet::VOs::Chain& chain, int beg, int end) {
        return ((std::abs(REGZ(chain.residues[beg]))<4.0?0:REGZ(chain.residues[beg])) * 
                (std::abs(REGZ(chain.residues[end]))<4.0?0:REGZ(chain.residues[end])) > -0.1);
    }
    std::vector<Tmdet::VOs::SecStrVec> Annotator::getParallelAlphas(Tmdet::VOs::Membrane& membrane) {
        std::vector<Tmdet::VOs::SecStrVec> ret;
        for (auto& vector : protein.secStrVecs) {
            if (vector.type.isStrictAlpha() 
                && checkParallel(vector, membrane)) {
                DEBUG_LOG("\t===>parallel: {}",vector.type.code);
                ret.emplace_back(vector);
            }
        }
        return ret;
    }

    bool Annotator::checkParallel(Tmdet::VOs::SecStrVec& vec, Tmdet::VOs::Membrane& membrane) const {    
        auto begRes = protein.chains[vec.chainIdx].residues[vec.begResIdx];
        auto endRes = protein.chains[vec.chainIdx].residues[vec.endResIdx];
        DEBUG_LOG("checkParallel: {}:{}-{} {}",protein.chains[vec.chainIdx].id, 
                begRes.authId,endRes.authId,vec.type.name);
        if (protein.chains[vec.chainIdx].type.isAlpha()
            && protein.chains[vec.chainIdx].selected 
            && begRes.selected 
            && endRes.selected
            && (REGHZ(begRes) < (REGTYPE(begRes).isNotMembrane()?10.0:5.0) 
                    || REGHZ(endRes) < (REGTYPE(endRes).isNotMembrane()?10.0:5.0))
            //&& (REGTYPE(begRes).isNotAnnotatedMembrane() || REGTYPE(begRes).isNotMembrane())
            //&& (REGTYPE(endRes).isNotAnnotatedMembrane() || REGTYPE(endRes).isNotMembrane())
         ) {
            double angle;
            if (membrane.type.isPlane()) {
                angle = std::abs(90.0 - Tmdet::Helpers::Vector::angle((vec.end-vec.begin),gemmi::Vec3(0,0,1)));
            }
            else {
                auto normal = (vec.begin+vec.end) / 2 - gemmi::Vec3(0,0,membrane.origo);
                angle = std::abs(90.0 - Tmdet::Helpers::Vector::angle((vec.end-vec.begin),normal));
            }
            DEBUG_LOG("is Parallel Alpha: {}:{}-{} {}",protein.chains[vec.chainIdx].id, 
                begRes.authId,endRes.authId,angle);
            return ( angle < ifhAngleLimit);
        }
        return false;
    }

    double Annotator::averageSurface(Tmdet::VOs::Chain& chain, int beg, int end) {
        double surf = 0.0;
        for(int i=beg; i<=end; i++) {
            surf += chain.residues[i].surface;
        }
        surf /= (end-beg+1);
        DEBUG_LOG("Annotator::averageSurface: {} {} {}: {}",chain.id,chain.residues[beg].authId,chain.residues[end].authId,surf);
        return surf;
    }

    double Annotator::hydrophocityMomentum(Tmdet::VOs::Chain& chain, int beg, int end) {
        auto vec = gemmi::Vec3(1.0,0.0,0.0);
        auto sum = gemmi::Vec3(0.0,0.0,0.0);
        double ca = -0.173648178;
        double sa = 0.984807753;
        for(int i=beg; i<=end; i++) {
            sum += vec * (chain.residues[i].type.hsc + 12.3 ) / 16.0;
            double x = ca * vec.x - sa * vec.y;
            double y = sa * vec.x + ca * vec.y;
            vec.x = x;
            vec.y = y;
        }
        double hm = sqrt(sum.x * sum.x + sum.y * sum.y);
        DEBUG_LOG("Annotator::hydophobicityMomentum: {} {} {}: {}",chain.id,chain.residues[beg].authId,chain.residues[end].authId,hm);
        return std::abs(hm);
    }

    void Annotator::finalCheck() {
        int nA=0;
        int nB=0;
        protein.eachChain(
            [&](Tmdet::VOs::Chain &chain) {
                for(const auto& region: chain.regions) {
                    nA += region.type.isAlpha();
                    nB += region.type.isBeta();
                }
            }
        );
        DEBUG_LOG("Annotator::finalCheck: type: {} na: {} nb: {} tmp: {}",
            protein.type.name,nA,nB,protein.tmp?"yes":"no");
        if ((!protein.type.isAlpha() && nB < 8)
            || (nA == 0 && nB == 0)
            || (!protein.type.isBeta() && nA == 0)) {
            protein.notTransmembrane();
        }
        DEBUG_LOG("Annotator::finalCheck: type: {} na: {} nb: {} tmp: {}",
            protein.type.name,nA,nB,protein.tmp?"yes":"no");
    }

    void Annotator::setMembraneSize() {
        double minX = 10000;
        double maxX = -10000;
        double minY = 10000;
        double maxY = -10000;
        protein.eachResidue(
            [&](Tmdet::VOs::Residue& residue) {
                for(const auto& a: residue.atoms) {
                    minX = (a.gemmi.pos.x<minX?a.gemmi.pos.x:minX);
                    maxX = (a.gemmi.pos.x>maxX?a.gemmi.pos.x:maxX);
                    minY = (a.gemmi.pos.y<minY?a.gemmi.pos.y:minY);
                    maxY = (a.gemmi.pos.y>maxY?a.gemmi.pos.y:maxY);
                }
            }
        );
        double r = (maxX-minX>maxY-minY?maxX-minX:maxY-minY);
        r/=2;
        for (auto& membrane: protein.membranes) {
            membrane.membraneRadius = r;
        }
    }
}