// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#include <string>

#include <Config.hpp>
#include <DTOs/Protein.hpp>
#include <DTOs/Region.hpp>
#include <Engine/Fragmenter.hpp>
#include <Engine/Organizer.hpp>
#include <Engine/PlaneSideDetector.hpp>
#include <Engine/RegionHandler.hpp>
#include <Helpers/Pymol.hpp>
#include <System/Logger.hpp>
#include <Utils/Fragment.hpp>
#include <VOs/Protein.hpp>
#include <VOs/Region.hpp>

namespace Tmdet::Engine {

    void Fragmenter::run() {
        DEBUG_LOG("Processing Fragmenter:run()");
        if (getNumberOfSelectedChains()>1) {
            WARN_LOG("Fragment analysis can be run on protein containing exactly one chain");
            return;
        }
        saveState();
        auto fragmentUtil = Tmdet::Utils::Fragment(protein);
        auto numFrags = fragmentUtil.run();
        //toPymol(); return;
        DEBUG_LOG("Number of fragments: {}",numFrags);
        runOnFragments(numFrags);
        findClusters();
        runOnBestCluster(findBestCluster());
        finalize();
        DEBUG_LOG("Processed Fragmenter:run()");
    }

    int Fragmenter::getNumberOfSelectedChains() {
        int ret = 0;
        protein.eachSelectedChain(
            [&](Tmdet::VOs::Chain& chain) {
                ret++;
                chIdx = chain.idx;
            }
        );
        return ret;
    }

    void Fragmenter::runOnFragments(int numFragments) {
        DEBUG_LOG("Processing Fragmenter:runOnFragments()");
        for(int i=0; i<numFragments; i++) {
            DEBUG_LOG("FragmentId: {}",i);
            std::string members="";
            for(auto& residue: protein.chains[chIdx].residues) {
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
                d.regions = protein.chains[chIdx].regions;
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
            DEBUG_LOG("Fragment {}: tmp: {} clusterId: {} {:.2f} {:.2f} {:.2f}",
                i,(data[i].tmp?"yes":"no"),data[i].clusterId,
                data[i].normal.x,data[i].normal.y,data[i].normal.z);
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
        DEBUG_LOG("Processed Fragmenter:findBestCluster({},{})", bestClusterId,max);
        return bestClusterId;
    }

    void Fragmenter::runOnBestCluster(int bestClusterId) {
        DEBUG_LOG("Processing Fragmenter:runOnBestCluster()");
        protein.clear();
        restoreState();
        for (auto& d: data){
            d.final = ((int)d.clusterId == bestClusterId);
        }
        std::string res="";
        for(auto& residue: protein.chains[chIdx].residues) {
            if (residue.temp.contains("fragment")) {
                residue.selected = false;
                for (auto& d: data){
                    if (d.final  && (int)d.id == any_cast<int>(residue.temp.at("fragment"))) {
                        residue.selected = true;
                    }
                }
            }
            else {
                residue.selected = false;
            }
            res += (residue.selected?"1":"0");
        }
        DEBUG_LOG("Final selected residues: {}",res);
        auto organizer = Tmdet::Engine::Organizer(protein, args);
        DEBUG_LOG("Processed Fragmenter:runOnBestCluster()");
    }

    void Fragmenter::finalize() {
        protein.chains[chIdx].eachResidue(
            [&](Tmdet::VOs::Residue& residue) -> void {
                residue.selected = !residue.selected;
            }
        );
        auto sideDetector = Tmdet::Engine::PlaneSideDetector(protein);
        for(auto& region: protein.chains[chIdx].regions) {
            for (int i=region.beg.idx; i<=region.end.idx; i++) {
                protein.chains[chIdx].residues[i].temp.insert({"type",region.type});
            }
        }
        protein.chains[chIdx].eachResidue(
            [&](Tmdet::VOs::Residue& residue) -> void {
                residue.selected = true;
                if (!residue.temp.contains("type")) {
                    residue.temp.insert({"type",Tmdet::Types::RegionType::UNK});
                }
                if (any_cast<Tmdet::Types::Region>(residue.temp.at("type")).isNotAnnotatedMembrane()) {
                    residue.temp.at("type") = Tmdet::Types::RegionType::ERROR_FP;
                }
            }
        );
        for (auto& d: data) {
            if (d.tmp && !d.final) {
                for (auto& r: d.regions) {
                    if (r.type.isAnnotatedTransMembraneType()) {
                        for (int i=r.beg.idx; i<=r.end.idx; i++) {
                            protein.chains[chIdx].residues[i].temp.at("type") = Tmdet::Types::RegionType::ERROR_FN;
                        }
                    }
                }
            }
        }
        auto regionHandler = Tmdet::Engine::RegionHandler(protein);
        DEBUG_LOG("Final regions: {}",regionHandler.toString("type"));
        protein.chains[chIdx].regions.clear();
        regionHandler.store<Tmdet::Types::Region>();
    }

    void Fragmenter::toPymol() {
        auto& chain = protein.chains[chIdx];
        chain.regions.clear();
        chain.regions.reserve(100);
        auto l = chain.length -1;
        Tmdet::VOs::Region region = {
            {chain.residues[0].authId, chain.residues[0].authIcode,chain.residues[0].labelId, 0},
            {chain.residues[l].authId, chain.residues[l].authIcode,chain.residues[l].labelId, l},
            Tmdet::Types::RegionType::UNK
        };
        chain.regions.push_back(region);
        int beg=0;
        int end=0;
        RegionHandler regionHandler(protein);
        while(regionHandler.getNext<int>(chain, beg, end, "fragment")) {
            int fr = any_cast<int>(chain.residues[beg].temp.at("fragment"))%14;
            Tmdet::VOs::Region region = {
                {chain.residues[beg].authId, chain.residues[beg].authIcode,chain.residues[beg].labelId, beg},
                {chain.residues[end-1].authId, chain.residues[end-1].authIcode,chain.residues[end-1].labelId, end-1},
                getRegionType(fr)
            };
            DEBUG_LOG("Region stored: {}:{}-{}-{}",chain.id,region.beg.authId,region.end.authId,region.type.code);
            chain.regions.push_back(region);
            beg=end;
        }
        auto pymol = Tmdet::Helpers::Pymol(protein);
        pymol.show(protein.inputFile);
    }

    Tmdet::Types::Region Fragmenter::getRegionType(int fr) {
        for(auto& [key,reg]: Tmdet::Types::Regions) {
            if (reg.id == fr) {
                return reg;
            }
        }
        return Tmdet::Types::RegionType::UNK;
    }

    void Fragmenter::saveState() {
        protein.eachResidue(
            [&](Tmdet::VOs::Residue& residue) {
                for(const auto& atom: residue.atoms) {
                    depo.push_back(atom.gemmi.pos);
                }
            }
        );
        for (const auto& ssVec: protein.secStrVecs) {
            depo.push_back(ssVec.begin);
            depo.push_back(ssVec.end);
        }
    }

    void Fragmenter::restoreState() {
        int i = 0;
        protein.eachResidue(
            [&](Tmdet::VOs::Residue& residue) {
                for(const auto& atom: residue.atoms) {
                    atom.gemmi.pos.x = depo[i].x;
                    atom.gemmi.pos.y = depo[i].y;
                    atom.gemmi.pos.z = depo[i].z;
                    i++;
                }
            }
        );
        for (auto& ssVec: protein.secStrVecs) {
            ssVec.begin.x = depo[i].x;
            ssVec.begin.y = depo[i].y;
            ssVec.begin.z = depo[i].z;
            i++;
            ssVec.end.x = depo[i].x;
            ssVec.end.y = depo[i].y;
            ssVec.end.z = depo[i].z;
            i++;
        }
    }

}
