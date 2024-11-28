#include <any>
#include <gemmi/model.hpp>
#include <Config.hpp>
#include <Helpers/Vector.hpp>
#include <Engine/Annotator.hpp>
#include <Engine/RegionHandler.hpp>
#include <Engine/SideDetector.hpp>
#include <Engine/BetaAnnotator.hpp>
#include <System/Logger.hpp>
#include <Types/Chain.hpp>
#include <Types/Region.hpp>
#include <Utils/SecStrVec.hpp>

#define REGTYPE(r) any_cast<Tmdet::Types::Region>(r.temp.at("type"))
#define REGTTYPE(r) any_cast<Tmdet::Types::Region>(r.temp.at("ttype"))
#define REGZTYPE(r) any_cast<Tmdet::Types::Region>(r.temp.at("ztype"))
#define REGTYPEW(r,w) any_cast<Tmdet::Types::Region>(r.temp.at(w))
#define REGZ(r) any_cast<double>(r.temp.at("z"))
#define REGHZ(r) any_cast<double>(r.temp.at("hz"))

namespace Tmdet::Engine {

    void Annotator::run() {
        DEBUG_LOG("Processing Annotator::run()");
        setChainsType();
        detectInterfacialHelices();
        annotateChains();
        //finalize(); //if needed
        regionHandler.store();
        DEBUG_LOG(" Processed Annotator::run({})",regionHandler.toString("type"));
    }

    void Annotator::setChainsType() {
        DEBUG_LOG("Processing Annotator::setChainsType()");
        protein.eachChain(
            [&](Tmdet::ValueObjects::Chain& chain) -> void {
                chain.type = (chain.selected?
                    Tmdet::Types::ChainType::NON_TM:
                    Tmdet::Types::ChainType::NOT_SELECTED);
                if (chain.selected) {
                    int alpha=0;
                    int beta=0;
                    for(const auto& residue: chain.residues) {
                        if (REGTTYPE(residue) == Tmdet::Types::RegionType::MEMB) {
                            if (residue.ss.isAlpha()) {
                                alpha++;
                            }
                            if (residue.ss.isBeta()) {
                                beta++;
                            }
                        }
                    }
                    if (alpha>beta) {
                        chain.type = Tmdet::Types::ChainType::ALPHA;
                        protein.type = (protein.type==Tmdet::Types::ProteinType::TM_BETA?
                                Tmdet::Types::ProteinType::TM_MIXED:
                                Tmdet::Types::ProteinType::TM_ALPHA);
                    }
                    else {
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
            [&](Tmdet::ValueObjects::Chain& chain) -> void {
                smoothRegions(chain,"type");
                smoothRegions(chain,"ttype");
                detectLoops(chain);
                if (chain.type == Tmdet::Types::ChainType::ALPHA) {
                    detectReEntrantLoops(chain);
                    detectTransmembraneHelices(chain);
                }
                else if (chain.type == Tmdet::Types::ChainType::BETA) {
                    auto annotator = Tmdet::Engine::BetaAnnotator(chain,regionHandler);
                }
            }
        );
    }

    void Annotator::smoothRegions(Tmdet::ValueObjects::Chain& chain, std::string what) {
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

    void Annotator::detectLoops(Tmdet::ValueObjects::Chain& chain) {
        int beg=0;
        int end=0;
        while(regionHandler.getNext(chain,beg,end,"type")) {
            if (REGTYPE(chain.residues[beg]) == Tmdet::Types::RegionType::MEMB) {
                detectLoop(chain,beg,end);
            }
            beg=end;
        }
    }

    void Annotator::detectLoop(Tmdet::ValueObjects::Chain& chain, int beg, int end) {
        for(int i=beg+5; i<end-5; i++) {
            if (any_cast<double>(chain.residues[i].temp.at("hz")) < 4.0
                && any_cast<double>(chain.residues[i-1].temp.at("direction")) * any_cast<double>(chain.residues[i].temp.at("direction")) < 0) {
                    chain.residues[i].temp.at("type") = chain.residues[i].temp.at("ztype");
            }
        }
    }

    void Annotator::detectInterfacialHelices() {
        DEBUG_LOG("Processing Annotator::detectInterfacialHelices()");
        for(auto& membrane: protein.membranes) {
            auto alphaVecs = getParallelAlphas(membrane);
            if (!alphaVecs.empty() ) {
                for(const auto& vector: alphaVecs) {
                    if (protein.chains[vector.chainIdx].type == Tmdet::Types::ChainType::ALPHA) {
                        regionHandler.replace(protein.chains[vector.chainIdx],vector.begResIdx,vector.endResIdx,Tmdet::Types::RegionType::IFH);
                    }
                }
            }
            DEBUG_LOG(" #alphaVecs: {}",alphaVecs.size());
        }
        DEBUG_LOG(" Processed Annotator::detectInterfacialHelices()");
    }

    void Annotator::detectReEntrantLoops(Tmdet::ValueObjects::Chain& chain) {
        DEBUG_LOG("Processing Annotator::detectReEntrantLoops()");
        int begin = 0;
        int end = 0;
        while(regionHandler.getNext(chain,begin,end,"type")) {
            end--;
            if (REGTYPE(chain.residues[begin]) == Tmdet::Types::RegionType::MEMB) {
                DEBUG_LOG("detectReEntrantLoops: {} {} {}",chain.id,begin,end);
                if (end - begin >= TMDET_REENTRANT_LOOP_MIN_LENGTH 
                    && begin > 1
                    && end < chain.length -1
                    && sameSide(chain,begin,end)
                    && hasHelixLoop(chain,begin,end)
                    //&& isOnOneSide(chain,begin,end)
                    ) {
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
        
    bool Annotator::hasHelixLoop(Tmdet::ValueObjects::Chain& chain, int begin, int end) {
        bool hasHelix = false;
        bool hasNoSS = false;
        for(int i=begin; i<=end; i++) {
            if (chain.residues[i].secStrVecIdx ==-1) {
                hasNoSS = true;
            }
            else if (protein.secStrVecs[chain.residues[i].secStrVecIdx].type == Tmdet::Types::SecStructType::H) {
                hasHelix = true;
            }
        }
        DEBUG_LOG("hasOneHelix: {} {} {}: {} {}",chain.id,begin,end,(hasHelix?"Yes":"No"),(hasNoSS?"Yes":"No"));
        return hasHelix && hasNoSS;
    }

    void Annotator::detectTransmembraneHelices(Tmdet::ValueObjects::Chain& chain) {
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

   /* bool Annotator::isOnOneSide(Tmdet::ValueObjects::Chain& chain, int beg, int end) {
        for (int i=beg+1; i<=end; i++) {
            if (!sameSide(chain,beg,i)) {
                return false;
            }
        }
        return true;
    }

    void Annotator::extendRegions(Tmdet::Engine::RegionHandler& regionHandler) {
        protein.eachSelectedChain(
            [&](Tmdet::ValueObjects::Chain& chain) -> void {
                int beg = 0;
                int end = 0;
                while(regionHandler.getNext(chain,beg,end,"type")) {
                    if (REGTYPE(chain.residues[beg]) == Tmdet::Types::RegionType::MEMB) {
                        auto numCross = getNumCross(chain,beg,end-1);
                        setTurnResidues(chain,numCross);
                        for( auto cr: numCross) {
                            extendRegion(chain,cr[0],cr[1],regionHandler);
                        }
                    }
                    beg = end;
                }
            }
        );
    }
*/
    bool Annotator::sameSide(Tmdet::ValueObjects::Chain& chain, int beg, int end) {
        DEBUG_LOG("sameSide: {} {} {} {}",chain.id,beg,end,(REGZ(chain.residues[beg]) * REGZ(chain.residues[end])>0));
        return REGZ(chain.residues[beg]) * REGZ(chain.residues[end]) > 0;
    }

 /*   
    std::vector<std::array<int,2>> Annotator::getNumCross(Tmdet::ValueObjects::Chain& chain, int beg, int end) {
        std::vector<std::array<int,2>> ret;
        for (int i=beg+1; i<=end; i++) {
            int b;
            int e;
            if (REGTTYPE(chain.residues[i-1]) != Tmdet::Types::RegionType::MEMB
                && REGTTYPE(chain.residues[i]) == Tmdet::Types::RegionType::MEMB) {
                    b=i;
            }
            if (REGTTYPE(chain.residues[i-1]) == Tmdet::Types::RegionType::MEMB
                && REGTTYPE(chain.residues[i]) != Tmdet::Types::RegionType::MEMB) {
                    e=i-1;
                    if (!sameSide(chain,b,e)) {
                        ret.push_back(std::array<int,2>{b,e});
                    }
                    else {
                        if (int t = getTurnResidue(chain,b,e); t!=-1) {
                            ret.push_back(std::array<int,2>{b,t-1});
                            ret.push_back(std::array<int,2>{t+1,e});
                        }
                    }
            }
        }
        return ret;
    }

    int Annotator::getTurnResidue(Tmdet::ValueObjects::Chain& chain, int beg, int end) {
        double max = -1e30;
        int ret = -1;
        double zbeg = REGZ(chain.residues[beg]);
        double zend = REGZ(chain.residues[end]);
        for (int i=beg; i<=end; i++) {
            double zi = REGZ(chain.residues[i]);
            double q = std::abs(zbeg-zi) + std::abs(zend-zi);
            if (max<q) {
                max = q;
                ret = i;
            }
        }
        return ret;
    }

    void Annotator::setTurnResidues(Tmdet::ValueObjects::Chain& chain, std::vector<std::array<int,2>>& numCross) {
        for (unsigned int i=1; i<numCross.size(); i++) {
            if (int t = getTurnResidue(chain,numCross[i-1][1],numCross[i][0]); t!=-1) {
                chain.residues[t].temp.at("type") = chain.residues[t].temp.at("ztype");
            }
        }
    }

    void Annotator::extendRegion(Tmdet::ValueObjects::Chain& chain, int beg, int end, Tmdet::Engine::RegionHandler& regionHandler) {
        int rbeg = beg;
        while(rbeg>=0 && REGTYPE(chain.residues[rbeg]) == Tmdet::Types::RegionType::MEMB) {
            rbeg--;
        }
        int rend = end;
        while(rend < chain.length && REGTYPE(chain.residues[rend]) == Tmdet::Types::RegionType::MEMB) {
            rend++;
        }
        regionHandler.replace(chain, rbeg+1, rend-1, 
                (chain.type==Tmdet::Types::ChainType::ALPHA?
                    Tmdet::Types::RegionType::HELIX:
                    Tmdet::Types::RegionType::BETA)
        );
        DEBUG_LOG("extendRegion: {} {} {} {}",chain.id,rbeg+1,rend-1,REGTYPE(chain.residues[beg]).name);
    }
*/
    std::vector<Tmdet::ValueObjects::SecStrVec> Annotator::getParallelAlphas(Tmdet::ValueObjects::Membrane& membrane) {
        std::vector<Tmdet::ValueObjects::SecStrVec> ret;
        for (auto& vector : protein.secStrVecs) {
            if (vector.type.isAlpha() && checkParallel(vector, membrane)) {
                ret.emplace_back(vector);
            }
        }
        return ret;
    }

    bool Annotator::checkParallel(Tmdet::ValueObjects::SecStrVec& vec, Tmdet::ValueObjects::Membrane& membrane) const {    
        double sign = vec.begin.z / std::abs(vec.begin.z);
        DEBUG_LOG("checkParallel: {} {} {}",sign,std::abs(sign * vec.begin.z - membrane.halfThickness),std::abs(sign * vec.end.z - membrane.halfThickness));
        return (std::abs(sign * vec.begin.z - membrane.halfThickness) < 10
                && std::abs(sign * vec.end.z - membrane.halfThickness) < 10
                && std::abs(Tmdet::Helpers::Vector::cosAngle((vec.end-vec.begin),gemmi::Vec3(0,0,1))) < 0.2
        );
    }


}