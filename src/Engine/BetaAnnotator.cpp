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
        if (chain.type.isBeta()) {
            //detectOutside();
            detectBarrel();
        }
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
                residue.temp.try_emplace("out",std::any(false));
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
                residue.temp.erase("out");
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
                && /*isOut(chain.residues[i])*/ chain.residues[i].ss.isBeta()
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
        int beg=0;
        int end=0;
        while(regionHandler.getNext(chain,beg,end,"type")) {
            if (any_cast<Tmdet::Types::Region>(chain.residues[beg].temp.at("type")).isBeta() 
                && end-beg < 5) {
                    regionHandler.replace(chain,beg,end-1,Tmdet::Types::RegionType::MEMB,"type");
            }
            beg=end;
        }
        beg=0;
        end=0;
        while(regionHandler.getNext(chain,beg,end,"type")) {
            if (any_cast<Tmdet::Types::Region>(chain.residues[beg].temp.at("type")).isNotAnnotatedMembrane() 
                && beg > 0
                && end < chain.length-1
                && end-beg < 3
                && any_cast<Tmdet::Types::Region>(chain.residues[beg-1].temp.at("type")).isBeta()
                && any_cast<Tmdet::Types::Region>(chain.residues[end].temp.at("type")).isBeta()
                && regionHandler.sameDirection(chain,beg-1,end)) {
                    regionHandler.replace(chain,beg,end-1,Tmdet::Types::RegionType::BETA,"type");
            }
            beg=end;
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
            || chain.residues[pos].ss.isStrictAlpha()) {
            return --count;
        }
        chain.residues[pos].temp.at("cluster") = std::any(cluster);
        if (pos>0 && (/*isOut(chain.residues[pos-1])*/ chain.residues[pos-1].ss.isBeta() || regionHandler.sameDirection(chain,pos-1,pos))) {
            count = setCluster(pos-1,cluster,++count);
        }
        if (pos<chain.length-1 && (/*isOut(chain.residues[pos+1])*/ chain.residues[pos+1].ss.isBeta() || regionHandler.sameDirection(chain,pos,pos+1))) {
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
                && any_cast<double>(chain.residues[beg].temp.at("hz")) < 4.0
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
                    && !isOut(residue) /*!residue.isInside()*/) {
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

    void BetaAnnotator::detectOutside() {
        double r=0;
        int zmin=10000;
        int zmax=-10000;
        chain.eachResidue(
            [&](Tmdet::ValueObjects::Residue& residue) -> void {
                if (any_cast<Tmdet::Types::Region>(residue.temp.at("type")).isNotAnnotatedMembrane()) {
                    const auto ca = residue.getCa();
                    double dist = sqrt(ca->pos.x * ca->pos.x + ca->pos.y * ca->pos.y);
                    if (r<dist) {
                        r=dist;
                    }
                    double z = any_cast<double>(residue.temp.at("z"));
                    zmin = (z<zmin?(int)z:zmin);
                    zmax = (z>zmax?(int)z:zmax);
                }
            }
        );
        r+=5;
        std::vector<std::vector<double>> closestd;
        std::vector<std::vector<int>> closesti;
        for (int z=zmin; z<=zmax; z++) {
            std::vector<double> d(30,10000);
            closestd.emplace_back(d);
            std::vector<int> i(30,-1);
            closesti.emplace_back(i);
        }
        chain.eachResidue(
            [&](Tmdet::ValueObjects::Residue& residue) -> void {
                if (any_cast<Tmdet::Types::Region>(residue.temp.at("type")).isNotAnnotatedMembrane()) {
                    const auto ca = residue.getCa();
                    double dist = r - sqrt(ca->pos.x * ca->pos.x + ca->pos.y * ca->pos.y);
                    double deg = 180 * atan(std::abs(ca->pos.x)/std::abs(ca->pos.y)) / M_PI;
                    if (ca->pos.x<0) {
                        deg = (ca->pos.y>0?360-deg:180+deg);
                    }
                    else if (ca->pos.y < 0) {
                        deg = 180 - deg;
                    }
                    int ideg = (int)(deg / 12);
                    int z = (int)(any_cast<double>(residue.temp.at("z")) - zmin);
                    if (dist < closestd[z][ideg]) {
                        closestd[z][ideg] = dist;
                        closesti[z][ideg] = residue.idx;
                    }
                   /* if (z>0 && dist < closestd[z-1][ideg]) {
                        closestd[z-1][ideg] = dist;
                        closesti[z-1][ideg] = residue.idx;
                    }
                    if (z<(zmax-zmin) && dist < closestd[z+1][ideg]) {
                        closestd[z+1][ideg] = dist;
                        closesti[z+1][ideg] = residue.idx;
                    }*/
                    //DEBUG_LOG("detectOutside: chain:{} res:{} x:{:5.3f} y:{:5.3f} z:{} deg:{:5.1f} ideg:{} dist:{:5.3f} closestd:{:5.3f} closesti:{}",
                    //    chain.id,residue.authId,ca->pos.x, ca->pos.y, z, deg, ideg, dist,
                    //    closestd[z][ideg], chain.residues[closesti[z][ideg]].authId);
                }
            }
        );
        for (int z=zmin; z<=zmax; z++) {
            for (int a=0; a<30; a++) {
                if (closesti[z-zmin][a] != -1) {
                    chain.residues[closesti[z-zmin][a]].temp.at("out") = std::any(true);
                    DEBUG_LOG("detectOutside: chain:{} res:{}",chain.id,closesti[z-zmin][a]);
                }
            }
        }
    }

    bool BetaAnnotator::isOut(Tmdet::ValueObjects::Residue& residue) {
        return any_cast<bool>(residue.temp.at("out"));
    }
}