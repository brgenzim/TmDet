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
        init();
        detectBarrel();
        DEBUG_LOG("Barrel end: {}",regionHandler.toString("type"));
        end();
        detectBarrelInside();
    }

    void BetaAnnotator::init() {
        chain.eachResidue(
            [&](Tmdet::ValueObjects::Residue& residue) -> void {
                residue.temp.try_emplace("cluster",std::any(-1));
            }
        );
    }

    void BetaAnnotator::end() {
        chain.eachResidue(
            [&](Tmdet::ValueObjects::Residue& residue) -> void {
                residue.temp.erase("cluster");
            }
        );
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
                    DEBUG_LOG("detectBarrel: {} {}",cluster,count);
                    cluster++;
                }
        }
        for(int i=0; i<chain.length; i++) {
            if (any_cast<int>(chain.residues[i].temp.at("cluster")) == maxCluster) {
                    chain.residues[i].temp.at("type") = std::any(Tmdet::Types::RegionType::BETA);
                    DEBUG_LOG("Beta:: {} {}",chain.id,chain.residues[i].authId);
                }
        }
    }

    int BetaAnnotator::setCluster(int pos, int cluster, int count) {
        if (any_cast<int>(chain.residues[pos].temp.at("cluster")) != -1 
            || any_cast<Tmdet::Types::Region>(chain.residues[pos].temp.at("type")) != Tmdet::Types::RegionType::MEMB) {
            return --count;
        }
        DEBUG_LOG("setCluster: checking chain:{} pos:{} cluster:{}",chain.id,chain.residues[pos].authId,cluster);
        chain.residues[pos].temp.at("cluster") = std::any(cluster);
        if (pos>0 && (chain.residues[pos-1].ss.isBeta() || sameDirection(pos-1,pos))) {
            DEBUG_LOG("setCluster-: from {} to {}",chain.residues[pos].authId,chain.residues[pos-1].authId);
            count = setCluster(pos-1,cluster,++count);
        }
        if (pos<chain.length-1 && (chain.residues[pos+1].ss.isBeta() || sameDirection(pos,pos+1))) {
            DEBUG_LOG("setCluster+: from {} to {}",chain.residues[pos].authId,chain.residues[pos+1].authId);
            count = setCluster(pos+1,cluster,++count);
        }
        if (chain.residues[pos].temp.contains("hbond1")) {
            auto to = any_cast<Tmdet::ValueObjects::HBond>(chain.residues[pos].temp.at("hbond1"));
            if (to.toChainIdx == chain.idx && to.toResIdx != -1 && chain.residues[pos].ss.isBeta()) {
                DEBUG_LOG("setCluster*: from {} to {}",chain.residues[pos].authId,chain.residues[to.toResIdx].authId);
                count = setCluster(to.toResIdx,cluster,++count);
            }
        }
        return count;
    }

    bool BetaAnnotator::sameDirection(int p1, int p2) {
        double d1 = any_cast<double>(chain.residues[p1].temp.at("direction"));
        double d2 = any_cast<double>(chain.residues[p2].temp.at("direction"));
        DEBUG_LOG("sameDirection: {} {} {} {}",chain.residues[p1].authId,chain.residues[p2].authId,d1,d2);
        return (d1*d2>0 && std::abs(d1)>8 && std::abs(d2)>8);
    }

    void BetaAnnotator::detectBarrelInside() {
        chain.eachResidue(
            [&](Tmdet::ValueObjects::Residue& residue) -> void {
                if (any_cast<Tmdet::Types::Region>(residue.temp.at("type")) == Tmdet::Types::RegionType::MEMB) {
                    residue.temp.at("type") = std::any(Tmdet::Types::RegionType::MEMBINS);
                }
            }
        );
        int beg=0;
        int end=0;
        while(regionHandler.getNext(chain,beg,end,"type")) {
            if (any_cast<Tmdet::Types::Region>(chain.residues[beg].temp.at("type")).isNotMembrane() 
                && end <= chain.length-1
                && beg > 0
                && end-beg < 4
                && (any_cast<Tmdet::Types::Region>(chain.residues[beg-1].temp.at("type")) == Tmdet::Types::RegionType::MEMBINS
                || any_cast<Tmdet::Types::Region>(chain.residues[end].temp.at("type")) == Tmdet::Types::RegionType::MEMBINS)) {
                regionHandler.replace(chain,beg,end-1,Tmdet::Types::RegionType::MEMBINS,"type");
            }
            beg=end;
        }
        beg=0;
        end=0;
        while(regionHandler.getNext(chain,beg,end,"type")) {
            if (any_cast<Tmdet::Types::Region>(chain.residues[beg].temp.at("type")) == Tmdet::Types::RegionType::MEMBINS
                && end - beg < 10) {
                regionHandler.replace(chain,beg,end-1,any_cast<Tmdet::Types::Region>(chain.residues[beg].temp.at("ztype")),"type");
            }
            beg=end;
        }
    }
}