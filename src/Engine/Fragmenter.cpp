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
#include <System/Logger.hpp>
#include <Utils/Fragment.hpp>
#include <VOs/Protein.hpp>
#include <VOs/Region.hpp>

namespace Tmdet::Engine {

    void Fragmenter::run() {
        saveState();
        auto fragmentUtil = Tmdet::Utils::Fragment(protein);
        auto numFrags = fragmentUtil.run();
        runOnFragments(numFrags);
        findClusters();
        runOnBestCluster(findBestCluster());
        if (protein.tmp) {
            finalize();
        }
    }

    void Fragmenter::runOnFragments(int numFragments) {
        for(int i=0; i<numFragments; i++) {
            protein.eachResidue(
                [&](Tmdet::VOs::Residue& residue) -> void {
                    if (residue.temp.contains("fragment")) {
                        residue.selected = (any_cast<int>(residue.temp.at("fragment")) == i);
                    }
                    else {
                        residue.selected = false;
                    }
                }
            );
            protein.clear();
            auto organizer = Tmdet::Engine::Organizer(protein, args);
            auto d = _fragmentData();
            d.id = i;
            d.clusterId = i;
            d.tmp = protein.tmp;
            if (protein.tmp) {
                d.membrane = protein.membranes[0];
                d.normal = organizer.getBestNormal();
                d.origo = protein.tmatrix.trans;
                protein.eachSelectedChain(
                    [&](Tmdet::VOs::Chain& chain) -> void {
                        for(auto& r: chain.regions) {
                            d.regions.push_back(r);
                            d.regionChainIndexes.push_back(chain.idx);
                        }
                    }
                );
            }
            data.push_back(d);
        }
    }

    void Fragmenter::findClusters() {
        for(unsigned int i=0; i<data.size(); i++) {
            if (data[i].tmp) {
                for(unsigned int j=0; j<i; j++) {
                    if (data[j].tmp && checkAngle(i,j)) {
                        for (unsigned int k=0; k<=i; k++) {
                            if (data[k].clusterId == i) {
                                data[k].clusterId = j;
                            }
                        }
                    }
                }
            }
        }
    }

    bool Fragmenter::checkAngle(int i, int j) {
        return data[i].normal.angle(data[j].normal) < 15;
    }

    int Fragmenter::findBestCluster() {
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
        return bestClusterId;
    }

    void Fragmenter::runOnBestCluster(int bestClusterId) {
        protein.clear();
        restoreState();
        for (auto& d: data){
            d.final = ((int)d.clusterId == bestClusterId);
        }
        protein.eachResidue(
            [&](Tmdet::VOs::Residue& residue) -> void {
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
            }
        );
        auto organizer = Tmdet::Engine::Organizer(protein, args);
    }

    void Fragmenter::finalize() {
        //invert selection
        protein.eachResidue(
            [&](Tmdet::VOs::Residue& residue) -> void {
                residue.selected = !residue.selected;
            }
        );
        //except for signal peptides
        protein.eachSelectedChain(
            [&](Tmdet::VOs::Chain& chain) -> void {
                if (chain.signalP[1] > 0) {
                    for (int i=0; i<chain.signalP[1]; i++) {
                        chain.residues[i].selected = false;
                    }
                }
            }
        );
        auto sideDetector = Tmdet::Engine::PlaneSideDetector(protein);
        for(auto& chain: protein.chains) {
            for(auto& region: chain.regions) {
                for (int i=region.beg.idx; i<=region.end.idx; i++) {
                    chain.residues[i].temp.insert({"type",region.type});
                }
            }
            chain.eachResidue(
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
                    for(unsigned int r=0; r<d.regions.size(); r++) {
                        if (d.regions[r].type.isAnnotatedTransMembraneType()) {
                            for (int i=d.regions[r].beg.idx; i<=d.regions[r].end.idx; i++) {
                                protein.chains[d.regionChainIndexes[r]].residues[i].temp.at("type") = Tmdet::Types::RegionType::ERROR_FN;
                            }
                        }
                    }
                }
            }
        }
        auto regionHandler = Tmdet::Engine::RegionHandler(protein);
        protein.eachSelectedChain(
            [&](Tmdet::VOs::Chain& chain) -> void {
                chain.regions.clear();
            }
        );
        regionHandler.store<Tmdet::Types::Region>();
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
