#include <any>
#include <memory>
#include <gemmi/model.hpp>
#include <Config.hpp>
#include <Helpers/Vector.hpp>
#include <Engine/Annotator.hpp>
#include <Engine/RegionHandler.hpp>
#include <Engine/BlendedSideDetector.hpp>
#include <Engine/PlaneSideDetector.hpp>
#include <Engine/BetaAnnotator.hpp>
#include <System/Logger.hpp>
#include <Types/Chain.hpp>
#include <Types/Region.hpp>
#include <Utils/SecStrVec.hpp>
#include <Types/Membrane.hpp>

#define REGTYPE(r) any_cast<Tmdet::Types::Region>(r.temp.at("type"))
#define REGTTYPE(r) any_cast<Tmdet::Types::Region>(r.temp.at("ttype"))
#define REGZTYPE(r) any_cast<Tmdet::Types::Region>(r.temp.at("ztype"))
#define REGTYPEW(r,w) any_cast<Tmdet::Types::Region>(r.temp.at(w))
#define REGZ(r) any_cast<double>(r.temp.at("z"))
#define REGHZ(r) any_cast<double>(r.temp.at("hz"))

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
            DEBUG_LOG("Blended sidedetector");
            sideDetector = std::make_unique<Tmdet::Engine::BlendedSideDetector>(protein);
        }
            
        setChainsType();
        if (protein.membranes[0].type.isPlane()) {
            detectInterfacialHelices();
        }
        annotateChains();
        if (int nm = regionHandler.finalize(); nm>10) {
            DEBUG_LOG("Too many unhandled membrane segments: {}",nm);
            protein.notTransmembrane();
        }
        else {
            regionHandler.store();
            finalCheck();
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
                if (chain.selected) {
                    int alpha=0;
                    int beta=0;
                    for(const auto& residue: chain.residues) {
                        if (residue.selected) {
                            if (REGTTYPE(residue) == Tmdet::Types::RegionType::MEMB) {
                                if (residue.ss.isAlpha()) {
                                    alpha++;
                                }
                                if (residue.ss.isBeta() && !residue.isInside()) {
                                    beta++;
                                }
                            }
                        }
                    }
                    if (alpha>3*beta && alpha>0) {
                        chain.type = Tmdet::Types::ChainType::ALPHA;
                        protein.type = (protein.type==Tmdet::Types::ProteinType::TM_BETA?
                                Tmdet::Types::ProteinType::TM_MIXED:
                                Tmdet::Types::ProteinType::TM_ALPHA);
                    }
                    else if (beta>0) {
                        chain.type = Tmdet::Types::ChainType::BETA;
                        protein.type = (protein.type==Tmdet::Types::ProteinType::TM_ALPHA?
                                Tmdet::Types::ProteinType::TM_MIXED:
                                Tmdet::Types::ProteinType::TM_BETA);
                    }
                }
                DEBUG_LOG("ChainType: {}-{}",chain.id,chain.type.name);
            }
        );
        DEBUG_LOG("Processed Annotator::setChainsType()");
    }

    void Annotator::annotateChains() {
        protein.eachSelectedChain(
            [&](Tmdet::VOs::Chain& chain) -> void {
                smoothRegions(chain,"type");
                smoothRegions(chain,"ttype");
                detectLoops(chain);
                if (chain.type.isAlpha()) {
                    detectReEntrantLoops(chain);
                    detectTransmembraneHelices(chain);
                }
                if (chain.type == Tmdet::Types::ChainType::BETA || protein.type.isBeta()) {
                    auto annotator = Tmdet::Engine::BetaAnnotator(chain,regionHandler);
                }
            }
        );
        DEBUG_LOG(" Processed: Annotator::annotateChains:\n {}",regionHandler.toString("type"));
    }

    void Annotator::smoothRegions(Tmdet::VOs::Chain& chain, std::string what) {
        int beg=0;
        int end=0;
        while(regionHandler.getNext(chain,beg,end,what)) {
            if (REGTYPEW(chain.residues[beg],what) == Tmdet::Types::RegionType::MEMB && end-beg < 3) {
                regionHandler.replace(chain,beg,end-1,REGZTYPE(chain.residues[beg]),what);
            }
            beg=end;
        }
        beg=0;
        end=0;
        while(regionHandler.getNext(chain,beg,end,what)) {
            if (REGTYPEW(chain.residues[beg],what) != Tmdet::Types::RegionType::MEMB && end-beg < 3) {
                regionHandler.replace(chain,beg,end-1,Tmdet::Types::RegionType::MEMB,what);
            }
            beg=end;
        }
        
    }

    void Annotator::detectLoops(Tmdet::VOs::Chain& chain) {
        int beg=0;
        int end=0;
        while(regionHandler.getNext(chain,beg,end,"type")) {
            if (REGTYPE(chain.residues[beg]) == Tmdet::Types::RegionType::MEMB) {
                detectLoop(chain,beg,end);
            }
            beg=end;
        }
    }

    void Annotator::detectLoop(Tmdet::VOs::Chain& chain, int beg, int end) {
        for(int i=beg+5; i<end-5; i++) {
            if ( (!sameSide(chain,beg,i) || !sameSide(chain,i,end-1))
                && std::abs(REGZ(chain.residues[i])) > 5.0
                && regionHandler.notSameDirection(chain,i-1,i+1)
                && regionHandler.notSameDirection(chain,i-2,i+2)
                && regionHandler.notSameDirection(chain,i-3,i+3)) {
                    chain.residues[i].temp.at("type") = chain.residues[i].temp.at("ztype");
                    DEBUG_LOG("detectLoop: {}:{}",chain.id, chain.residues[i].authId);
            }
        }
    }

    void Annotator::detectInterfacialHelices() {
        DEBUG_LOG("Processing Annotator::detectInterfacialHelices()");
        for(auto& membrane: protein.membranes) {
            auto alphaVecs = getParallelAlphas(membrane);
            if (!alphaVecs.empty() ) {
                for(const auto& vector: alphaVecs) {
                    if (vector.endResIdx -vector.begResIdx>7
                        && protein.chains[vector.chainIdx].residues[vector.begResIdx].selected
                        && protein.chains[vector.chainIdx].residues[vector.endResIdx].selected) {
                        regionHandler.replace(protein.chains[vector.chainIdx],vector.begResIdx,vector.endResIdx,Tmdet::Types::RegionType::IFH);
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
        while(regionHandler.getNext(chain,begin,end,"type")) {
            end--;
            if (REGTYPE(chain.residues[begin]) == Tmdet::Types::RegionType::MEMB) {
                DEBUG_LOG("detectReEntrantLoops: {} {} {}",chain.id,begin,end);
                DEBUG_LOG("\tz coords: {} {}",REGZ(chain.residues[begin]),REGZ(chain.residues[end]));
                DEBUG_LOG("\tsameside: {}",sameSide(chain,begin,end));
                if (std::abs(REGZ(chain.residues[begin]) - REGZ(chain.residues[end])) < 5.0
                    && sameSide(chain,begin,end)
                    && hasHelixLoop(chain,begin,end)) {
                        DEBUG_LOG("Loop found: {} {} {}",chain.id,begin,end);
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
        
    bool Annotator::hasHelixLoop(Tmdet::VOs::Chain& chain, int begin, int end) {
        int numHelix = 0;
        int numNoSS = 0;
        int vecIdx = -1;
        bool twoHelices = false;
        for(int i=begin; i<=end; i++) {
            if (chain.residues[i].secStrVecIdx ==-1) {
                numNoSS++;
            }
            else if (protein.secStrVecs[chain.residues[i].secStrVecIdx].type == Tmdet::Types::SecStructType::H) {
                numHelix++;
                if (vecIdx == -1) {
                    vecIdx = chain.residues[i].secStrVecIdx;
                }
                else if (vecIdx != chain.residues[i].secStrVecIdx) {
                    twoHelices=true;
                }
            }
        }
        DEBUG_LOG("hasOneHelix: {} {} {}: {} {}",chain.id,begin,end,numHelix,numNoSS);
        return (numHelix > 5 && (numNoSS > 2||twoHelices));
    }

    void Annotator::detectTransmembraneHelices(Tmdet::VOs::Chain& chain) {
        DEBUG_LOG("Processing Annotator::detectTransmembraneHelices()");
        int begin = 0;
        int end = 0;
        while(regionHandler.getNext(chain,begin,end,"type")) {
            if (REGTYPE(chain.residues[begin]) == Tmdet::Types::RegionType::MEMB 
                && end-begin > 10
                && REGZ(chain.residues[begin]) * REGZ(chain.residues[end-1]) < 0) {
                    regionHandler.replace(chain,begin,end-1,Tmdet::Types::RegionType::HELIX);
            }
            begin = end;
        }
        DEBUG_LOG(" Processed Annotator::detectTransmembraneHelices()");
    }

   bool Annotator::sameSide(Tmdet::VOs::Chain& chain, int beg, int end) {
        return REGZ(chain.residues[beg]) * REGZ(chain.residues[end]) > 0;
    }
    std::vector<Tmdet::VOs::SecStrVec> Annotator::getParallelAlphas(Tmdet::VOs::Membrane& membrane) {
        std::vector<Tmdet::VOs::SecStrVec> ret;
        for (auto& vector : protein.secStrVecs) {
            if (vector.type.isStrictAlpha() && checkParallel(vector, membrane)) {
                DEBUG_LOG("\t===>parallel");
                ret.emplace_back(vector);
            }
        }
        return ret;
    }

    bool Annotator::checkParallel(Tmdet::VOs::SecStrVec& vec, Tmdet::VOs::Membrane& membrane) const {    
        auto begRes = protein.chains[vec.chainIdx].residues[vec.begResIdx];
        auto endRes = protein.chains[vec.chainIdx].residues[vec.endResIdx];
        if (begRes.selected && endRes.selected) {
            DEBUG_LOG("checkParallel: {}:{} {}:{} {}",begRes.authId,REGHZ(begRes),endRes.authId,REGHZ(endRes),
                Tmdet::Helpers::Vector::cosAngle((vec.end-vec.begin),gemmi::Vec3(0,0,1)));
            return ((REGHZ(begRes) < 5.0 || REGHZ(endRes) < 5.0)
                    && std::abs(Tmdet::Helpers::Vector::cosAngle((vec.end-vec.begin),gemmi::Vec3(0,0,1))) < 0.3
            );
        }
        return false;
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
        if ((protein.type.isBeta() && nB < 8)
            || (nA == 0 && nB == 0)
            || (protein.type.isAlpha() && nA == 0)) {
            protein.notTransmembrane();
        }
        DEBUG_LOG("Annotator::finalCheck: type: {} na: {} nb: {} tmp: {}",
            protein.type.name,nA,nB,protein.tmp?"yes":"no");
    }

}