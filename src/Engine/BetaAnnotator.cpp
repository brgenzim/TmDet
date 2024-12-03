#include <any>
#include <gemmi/model.hpp>
#include <Config.hpp>
#include <Engine/BetaAnnotator.hpp>
#include <Engine/RegionHandler.hpp>
#include <System/Logger.hpp>
#include <Types/Region.hpp>
#include <ValueObjects/Chain.hpp>
#include <ValueObjects/HBond.hpp>
#include <ValueObjects/Residue.hpp>

namespace Tmdet::Engine {

    void BetaAnnotator::run() {
        DEBUG_LOG("Processing BetaAnnotator::run()");
        init();
        detectBarrel();
        DEBUG_LOG("Barrel end: {}",regionHandler.toString("type"));
        detectLoops();
        DEBUG_LOG("Loops end: {}",regionHandler.toString("type"));
        detectBarrelInside();
        end();
        DEBUG_LOG(" Processed BetaAnnotator::run()");
    }

    void BetaAnnotator::init() {
        DEBUG_LOG("Processing BetaAnnotator::init()");
        int max=0;
        chain.eachResidue(
            [&](Tmdet::ValueObjects::Residue& residue) -> void {
                if (max<residue.idx) {
                    max = residue.idx;
                }
            }
        );
        max++;
        reIndex = std::vector<int>(max,-1);
        int i=0;
        chain.eachResidue(
            [&](Tmdet::ValueObjects::Residue& residue) -> void {
                residue.temp.try_emplace("cluster",std::any(-1));
                residue.temp.try_emplace("from",std::any(-1));
                residue.temp.try_emplace("to",std::any(-1));
                reIndex[residue.idx] = i++;
            }
        );
        i=0;
        chain.eachResidue(
            [&](Tmdet::ValueObjects::Residue& residue) -> void {
                auto hbond = any_cast<Tmdet::ValueObjects::HBond>(residue.temp.at("hbond1"));
                
                if (hbond.toChainIdx == chain.idx && reIndex[hbond.toResIdx] != -1) {
                    chain.residues[reIndex[hbond.toResIdx]].temp["from"] = std::any(i);
                    residue.temp["to"] = std::any(reIndex[hbond.toResIdx]);
                }
                i++;
            }
        );
        DEBUG_LOG(" Processed BetaAnnotator::init()");
    }

    void BetaAnnotator::end() {
        DEBUG_LOG("Processing BetaAnnotator::end()");
        chain.eachResidue(
            [&](Tmdet::ValueObjects::Residue& residue) -> void {
                residue.temp.erase("cluster");
                residue.temp.erase("from");
                residue.temp.erase("to");
            }
        );
        DEBUG_LOG(" Processed BetaAnnotator::end()");
    }

    void BetaAnnotator::detectBarrel() {
        int cluster = 0;
        int maxCount = -1;
        int maxCluster = -2;
        for(int i=0; i<chain.length; i++) {
            if (any_cast<Tmdet::Types::Region>(chain.residues[i].temp.at("type")) == Tmdet::Types::RegionType::MEMB 
                && chain.residues[i].ss.isBeta()
                && any_cast<int>(chain.residues[i].temp.at("cluster")) == -1) {
                    auto count = setCluster(i,cluster,0);
                    if (count > maxCount) {
                        maxCount = count;
                        maxCluster = cluster;
                    }
                    cluster++;
                }
        }
        for(int i=0; i<chain.length; i++) {
            if (any_cast<int>(chain.residues[i].temp.at("cluster")) == maxCluster) {
                    chain.residues[i].temp.at("type") = std::any(Tmdet::Types::RegionType::BETA);
                }
        }
        //secondary structure was not set because of inappropriate data
        DEBUG_LOG("Max count: {}",maxCount);
        if (maxCount<40) {
            int beg=0;
            int end=0;
            while(regionHandler.getNext(chain,beg,end,"type")) {
                if (any_cast<Tmdet::Types::Region>(chain.residues[beg].temp.at("type")).isNotAnnotatedMembrane() 
                    && end-beg < 13
                    && end-beg > 5
                    && std::abs(averageDirection(beg,end-1)) > 10
                    && (averageOutSurface(beg,end-1) > 40 || averageBeta(beg,end-1) > 0.8)) {
                    regionHandler.replace(chain,beg,end-1,Tmdet::Types::RegionType::BETA,"type");
                }
                beg=end;
            }
        }
    }

    double BetaAnnotator::averageOutSurface(int beg, int end) {
        double surf = 0.0;
        for(int i=beg; i<=end; i++) {
            surf += chain.residues[i].outSurface;
        }
        surf /= (end-beg+1);
        DEBUG_LOG("BetaAnnotator::averageOutSurface: {} {} {}: {}",chain.id,chain.residues[beg].authId,chain.residues[end].authId,surf);
        return surf;
    }

    double BetaAnnotator::averageBeta(int beg, int end) {
        double beta = 0.0;
        for(int i=beg; i<=end; i++) {
            beta += (chain.residues[i].ss.isBeta()?1.0:0.0);
        }
        beta /= (end-beg+1);
        DEBUG_LOG("BetaAnnotator::averageBeta: {} {} {}: {}",chain.id,chain.residues[beg].authId,chain.residues[end].authId,beta);
        return beta;
    }

    double BetaAnnotator::averageDirection(int beg, int end) {
        double d = 0.0;
        for(int i=beg; i<=end; i++) {
            d += any_cast<double>(chain.residues[i].temp.at("direction"));
        }
        d /= (end-beg+1);
        DEBUG_LOG("BetaAnnotator::averageDirection: {} {} {}: {}",chain.id,chain.residues[beg].authId,chain.residues[end].authId,d);
        return d;
    }

    int BetaAnnotator::setCluster(int pos, int cluster, int count) {
        if (!chain.residues[pos].temp.contains("cluster") || any_cast<int>(chain.residues[pos].temp.at("cluster")) != -1
            || any_cast<Tmdet::Types::Region>(chain.residues[pos].temp.at("type")) != Tmdet::Types::RegionType::MEMB
            || chain.residues[pos].ss.isAlpha()) {
            return --count;
        }
        chain.residues[pos].temp.at("cluster") = std::any(cluster);
        if (pos>0 && (chain.residues[pos-1].ss.isBeta() || regionHandler.sameDirection(chain,pos-1,pos))) {
            count = setCluster(pos-1,cluster,++count);
        }
        if (pos<chain.length-1 && (chain.residues[pos+1].ss.isBeta() || regionHandler.sameDirection(chain,pos,pos+1))) {
            count = setCluster(pos+1,cluster,++count);
        }
        if (chain.residues[pos].ss.isBeta()) {
            int other = any_cast<int>(chain.residues[pos].temp["to"]);
            if (other != -1 ) {
                count = setCluster(other,cluster,++count);
            }
            other = any_cast<int>(chain.residues[pos].temp["from"]);
            if ( other != -1) {
                count = setCluster(other,cluster,++count);
            }
        }
        return count;
    }

    void BetaAnnotator::detectLoops() {
        int beg=0;
        int end=0;
        while(regionHandler.getNext(chain,beg,end,"type")) {
            if (any_cast<Tmdet::Types::Region>(chain.residues[beg].temp.at("type")).isBeta() 
                && end-beg < 3) {
                regionHandler.replace(chain,beg,end-1,Tmdet::Types::RegionType::MEMB,"type");
            }
            beg=end;
        }
        beg=0;
        end=0;
        while(regionHandler.getNext(chain,beg,end,"type")) {
            if (end <= chain.length-1
                && beg > 0
                && any_cast<Tmdet::Types::Region>(chain.residues[beg].temp.at("type")).isNotAnnotatedMembrane() 
                && (any_cast<Tmdet::Types::Region>(chain.residues[beg-1].temp.at("type")).isAnnotatedTransMembraneType()
                    || any_cast<Tmdet::Types::Region>(chain.residues[end].temp.at("type")).isAnnotatedTransMembraneType())
                && end-beg < 8) {
                regionHandler.replace(chain,beg,end-1,any_cast<Tmdet::Types::Region>(chain.residues[beg].temp.at("ztype")),"type");
            }
            beg=end;
        }
    }

    void BetaAnnotator::detectBarrelInside() {
        chain.eachResidue(
            [&](Tmdet::ValueObjects::Residue& residue) -> void {
                if (any_cast<Tmdet::Types::Region>(residue.temp.at("type")).isNotAnnotatedMembrane() 
                    && residue.isInside()) {
                    residue.temp.at("type") = (any_cast<double>(residue.temp["hz"]) > 3 ? 
                        std::any(Tmdet::Types::RegionType::MEMBINS) :
                        std::any(residue.temp.at("ztype")));
                }
            }
        );
        DEBUG_LOG("Barrel inside1: {}",regionHandler.toString("type"));
        int beg=0;
        int end=0;
        while(regionHandler.getNext(chain,beg,end,"type")) {
            if (any_cast<Tmdet::Types::Region>(chain.residues[beg].temp.at("type")).isNotAnnotatedMembrane()
                && end <= chain.length-1
                && beg > 0
                && (any_cast<Tmdet::Types::Region>(chain.residues[beg-1].temp.at("type")).isMembraneInside()
                    || any_cast<Tmdet::Types::Region>(chain.residues[end].temp.at("type")).isMembraneInside())) {
                regionHandler.replace(chain,beg,end-1,Tmdet::Types::RegionType::MEMBINS,"type");
            }
            beg=end;
        }
        DEBUG_LOG("Barrel inside2: {}",regionHandler.toString("type"));
        beg=0;
        end=0;
        while(regionHandler.getNext(chain,beg,end,"type")) {
            if (any_cast<Tmdet::Types::Region>(chain.residues[beg].temp.at("type")).isNotMembrane()
                && end <= chain.length-1
                && beg > 0
                && end-beg < 10
                && any_cast<Tmdet::Types::Region>(chain.residues[beg-1].temp.at("type")).isMembraneInside()
                && any_cast<Tmdet::Types::Region>(chain.residues[end].temp.at("type")).isMembraneInside()) {
                regionHandler.replace(chain,beg,end-1,Tmdet::Types::RegionType::MEMBINS,"type");
            }
            beg=end;
        }
        DEBUG_LOG("Barrel inside3: {}",regionHandler.toString("type"));
    }
}