#include <string>

#include <Config.hpp>
#include <DTOs/Protein.hpp>
#include <DTOs/Region.hpp>
#include <Engine/Fragmenter.hpp>
#include <Engine/Organizer.hpp>
#include <Engine/PlaneSideDetector.hpp>
#include <Engine/RegionHandler.hpp>
#include <System/Logger.hpp>
#include <Utils/Fragment.hpp>
#include <VOs/Protein.hpp>
#include <VOs/Region.hpp>

namespace Tmdet::Engine {

    void Fragmenter::run() {
        DEBUG_LOG("Processing Fragmenter:run()");
        if (protein.chains.size()>1) {
            WARN_LOG("Fragment analysis can be run on protein containing exactly one chain");
            return;
        }
        auto fragmentUtil = Tmdet::Utils::Fragment(protein);
        auto numFrags = fragmentUtil.run();
        //DEBUG_LOG("Fragment results: {}",Tmdet::DTOs::Protein::toString(protein));
        DEBUG_LOG("Number of fragments: {}",numFrags);
        runOnFragments(numFrags);
        findClusters();
        runOnBestCluster(findBestCluster());
        finalize();
        DEBUG_LOG("Processed Fragmenter:run()");
    }

    void Fragmenter::runOnFragments(int numFragments) {
        DEBUG_LOG("Processing Fragmenter:runOnFragments()");
        for(int i=0; i<numFragments; i++) {
            DEBUG_LOG("FragmentId: {}",i);
            toString();
            std::string members="";
            for(auto& residue: protein.chains[0].residues) {
                if (residue.temp.contains("fragment")) {
                    residue.selected = (any_cast<int>(residue.temp.at("fragment")) == i);
                    members += (residue.selected?"1":"0");
                }
                else {
                    residue.selected = false;
                    members += "0";
                }
            }
            DEBUG_LOG("Fragment members: {}",members);
            protein.clear();
            auto organizer = Tmdet::Engine::Organizer(protein, args);
            DEBUG_LOG("After searching: fragment {} tmp {}",i,(protein.tmp?"yes":"no"));
            auto d = _fragmentData();
            d.id = i;
            d.clusterId = i;
            d.tmp = protein.tmp;
            if (protein.tmp) {
                d.membrane = protein.membranes[0];
                d.normal = organizer.getBestNormal();
                d.origo = protein.tmatrix.trans;
                d.regions = protein.chains[0].regions;
                for(auto& r: d.regions) {
                    DEBUG_LOG("Fragment result: fr:{} {}",i,Tmdet::DTOs::Region::toString(r));
                }
            }
            data.push_back(d);
        }
        DEBUG_LOG("Processed Fragmenter:runOnFragments()");
    }

    void Fragmenter::findClusters() {
        DEBUG_LOG("Processing Fragmenter:findClusters()");
        for(unsigned int i=0; i<data.size(); i++) {
            DEBUG_LOG("Fragment {}: tmp: {} clusterId: {}",i,(data[i].tmp?"yes":"no"),data[i].clusterId);
        }
        for(unsigned int i=0; i<data.size(); i++) {
            if (data[i].tmp) {
                for(unsigned int j=0; j<i; j++) {
                    if (data[j].tmp && checkAngle(i,j) /*&& checkOrigo(i,j)*/) {
                        for (unsigned int k=0; k<=i; k++) {
                            if (data[k].clusterId == i) {
                                data[k].clusterId = j;
                                DEBUG_LOG("Merging fragments: {} {}",i,j);
                            }
                        }
                    }
                }
            }
        }
        DEBUG_LOG("Processed Fragmenter:findClusters()");
    }

    bool Fragmenter::checkAngle(int i, int j) {
        return data[i].normal.angle(data[j].normal) < 15;
    }

   /* bool Fragmenter::checkOrigo(int i, int j) {
        return data[i].origo.dist(data[j].origo) < 5;
    }*/

    int Fragmenter::findBestCluster() {
        DEBUG_LOG("Processing Fragmenter:findBestCluster()");
        int largestClusterId=0;
        for (auto& d: data) {
            if ((int)d.clusterId > largestClusterId) {
                largestClusterId = d.clusterId;
            }
        }
        std::vector<int> counts(largestClusterId+1,0);
        for (auto& d: data) {
            for (auto& r: d.regions) {
                if (r.type.isAnnotatedTransMembraneType()) {
                    counts[d.clusterId]++;
                }
            }
        }
        int max=0;
        int bestClusterId=-1;
        for (int i=0; i<=largestClusterId; i++) {
            if (max<counts[i]) {
                max = counts[i];
                bestClusterId = i;
            }
        }
        DEBUG_LOG("Processed Fragmenter:findBestCluster({})", bestClusterId);
        return bestClusterId;
    }

    void Fragmenter::runOnBestCluster(int bestClusterId) {
        DEBUG_LOG("Processing Fragmenter:runOnBestCluster()");
        protein.clear();
        toString();
        std::string members = "";
        for(auto& residue: protein.chains[0].residues) {
            if (residue.temp.contains("fragment")) {
                residue.selected = false;
                for (auto& d: data){
                    DEBUG_LOG("Final selection: res:{} resFr:{} d.id{} d.clusterId{}",
                        residue.authId,any_cast<int>(residue.temp.at("fragment")),d.id,d.clusterId);
                    if ((int)d.clusterId == bestClusterId 
                        && (int)d.id == any_cast<int>(residue.temp.at("fragment"))) {
                            residue.selected = true;
                    }
                }
                members += (residue.selected?"1":"0");
            }
            else {
                residue.selected = false;
            }
        }
        DEBUG_LOG("Final members: {}",members);
        auto organizer = Tmdet::Engine::Organizer(protein, args);
        DEBUG_LOG("Processed Fragmenter:runOnBestCluster()");
    }

    void Fragmenter::finalize() {
        protein.chains[0].eachResidue(
            [&](Tmdet::VOs::Residue& residue) -> void {
                residue.selected = !residue.selected;
            }
        );
        auto sideDetector = Tmdet::Engine::PlaneSideDetector(protein);
        for(auto& region: protein.chains[0].regions) {
            for (int i=region.beg.labelId-1; i<region.end.labelId; i++) {
                protein.chains[0].residues[i].temp.insert({"type",region.type});
            }
        }
        protein.chains[0].eachResidue(
            [&](Tmdet::VOs::Residue& residue) -> void {
                residue.selected = true;
                if (!residue.temp.contains("type")) {
                    residue.temp.insert({"type",Tmdet::Types::RegionType::UNK});
                }
                if (any_cast<Tmdet::Types::Region>(residue.temp.at("type")).isNotAnnotatedMembrane()) {
                    residue.temp.at("type") = Tmdet::Types::RegionType::ERROR;
                }
            }
        );
        auto regionHandler = Tmdet::Engine::RegionHandler(protein);
        DEBUG_LOG("Final regions: {}",regionHandler.toString("type"));
        protein.chains[0].regions.clear();
        regionHandler.store();
    }

    void Fragmenter::toString() {
        std::string cl="0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
        
        protein.eachSelectedChain(
            [&](Tmdet::VOs::Chain& chain) -> void {
                std::string cls="";
                chain.eachResidue(
                    [&](Tmdet::VOs::Residue& residue) -> void {
                        if (residue.temp.contains("fragment")) {
                            cls += cl[any_cast<int>(residue.temp.at("fragment"))];
                        }
                        else {
                            cls += "-";
                        }
                    }
                );
                DEBUG_LOG("Fragment results chain:{} {}",chain.id,cls);
            }
        );
    }
}
