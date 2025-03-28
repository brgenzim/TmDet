// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#include <any>
#include <gemmi/model.hpp>
#include <Config.hpp>
#include <Engine/BetaAnnotator.hpp>
#include <Engine/RegionHandler.hpp>
#include <System/Logger.hpp>
#include <Types/Region.hpp>
#include <Utils/NeighBors.hpp>
#include <VOs/Chain.hpp>
#include <VOs/CR.hpp>
#include <VOs/HBond.hpp>
#include <VOs/Residue.hpp>

namespace Tmdet::Engine {

    void BetaAnnotator::run() {
        protein.numBarrels = 0;
        getSheets();
        if (numSheets>7) {
            setConnections();
            detectBarrels();
            if (protein.numBarrels > 0) {
                setBarrel();
            }
        }
    }

    void BetaAnnotator::getSheets() {
        int vectorIndex = 0;
        for(auto& ssVec: protein.secStrVecs) {
            if (ssVec.type.isBeta()) {
                if (protein.chains[ssVec.chainIdx].residues[ssVec.begResIdx].selected 
                    && protein.chains[ssVec.chainIdx].residues[ssVec.endResIdx].selected ) {
                    bool inMembrane = false;
                    for (int i=ssVec.begResIdx; i<=ssVec.endResIdx; i++) {
                        if (protein.chains[ssVec.chainIdx].residues[i].selected 
                            && any_cast<Tmdet::Types::Region>(protein.chains[ssVec.chainIdx].residues[i].temp["type"]).isNotAnnotatedMembrane()) {
                            inMembrane = true;
                        }
                    }
                    double angle = std::abs(90 - Tmdet::Helpers::Vector::angle(gemmi::Vec3(0,0,1),ssVec.end - ssVec.begin));
                    if (inMembrane &&  angle > 10) {
                        sheetIndex.push_back(vectorIndex);
                        ssVec.sheetIdx = numSheets++;
                    }
                }
            }
            vectorIndex++;
        }
    }

    void BetaAnnotator::setConnections() {
        for (int i=0; i<numSheets; i++) {
            connectome.push_back(std::vector<int>(numSheets,0));
        }
        for(auto& ssVec: protein.secStrVecs) {
            if (ssVec.sheetIdx != -1) {
                for(int i=ssVec.begResIdx; i<=ssVec.endResIdx; i++) {
                    auto residue1 =  protein.chains[ssVec.chainIdx].residues[i];
                    if (residue1.selected 
                        && any_cast<Tmdet::Types::Region>(residue1.temp["type"]).isNotAnnotatedMembrane()) {
                            for(auto cr : Tmdet::Utils::NeighBors::get(residue1)) {
                                if (auto residue2 = protein.chains[cr.chainIdx].residues[cr.residueIdx]; residue2.selected) {
                                    int ssVecIdx = residue2.secStrVecIdx;
                                    if ( ssVecIdx != -1 && any_cast<Tmdet::Types::Region>(residue2.temp["type"]).isNotAnnotatedMembrane()
                                        && protein.secStrVecs[ssVecIdx].sheetIdx != -1
                                        && protein.secStrVecs[ssVecIdx].sheetIdx != ssVec.sheetIdx) {
                                            double co_angle = (residue1.temp.contains("co")?Tmdet::Helpers::Vector::angle(
                                                    any_cast<gemmi::Vec3>(residue1.temp["co"]),
                                                    any_cast<gemmi::Vec3>(residue1.temp["ca"]) - any_cast<gemmi::Vec3>(residue2.temp["ca"])
                                                ):0);
                                            if (co_angle < 50 || co_angle > 130) {
                                                connectome[ssVec.sheetIdx][protein.secStrVecs[ssVecIdx].sheetIdx]++;
                                                connectome[protein.secStrVecs[ssVecIdx].sheetIdx][ssVec.sheetIdx]++;
                                            }
                                    }
                                }
                            }
                    }
                }
            }
        }
        for (int i=0; i<numSheets; i++) {
            std::string res = "";
            for (int j=0; j<numSheets; j++) {
                res += std::format("{:2d} ",connectome[i][j]);
            }
        }
    }

    void BetaAnnotator::detectBarrels() {
        for (int i=0; i<numSheets; i++) {
            if (protein.secStrVecs[sheetIndex[i]].barrelIdx == -1) {
                std::vector<bool> elements(numSheets,false);
                elements[i] = true;
                int nb = detectBarrelSheets(i,-1,elements);
                nb = detectBarrelSheets(i,-1,elements); 
                if (nb >= args.getValueAsInt("minbs")) {
                    setIndex(elements);
                }
            }
        }
    }

    int BetaAnnotator::detectBarrelSheets(int sheetNum, int prevSheet, std::vector<bool>& elements) {
        int max = 0;
        double percent = 0.0;
        int maxSheet = -1;
        int minContactBetweenSheets = args.getValueAsInt("mcbs");
        for (int i=0; i<numSheets; i++) {
            if (max<connectome[sheetNum][i]
                    && protein.secStrVecs[sheetIndex[i]].barrelIdx == -1
                    && !elements[i] ) {
                max = connectome[sheetNum][i];
                percent = 50 * max / (protein.secStrVecs[sheetIndex[i]].endResIdx-protein.secStrVecs[sheetIndex[i]].begResIdx+1);
                maxSheet = i;
            }
        }
        if ((minContactBetweenSheets<=5 && (max>=minContactBetweenSheets || percent > 20))
            || max>=minContactBetweenSheets) {
            elements[maxSheet] = true;
            detectBarrelSheets(maxSheet, sheetNum, elements);
        }
        std::string s="";
        int ret = 0;
        for(int i=0; i<numSheets; i++) {
            if (elements[i]) {
                s+=std::format("-{}",i);
                ret++;
            }
        }
        return ret;
    }

    void BetaAnnotator::setIndex(std::vector<bool>& elements) {
        int nb = 0;
        for (int i=0; i<numSheets; i++) {
            if (elements[i]) {
                nb++;
                protein.secStrVecs[sheetIndex[i]].barrelIdx = protein.numBarrels;
            }
        }
        numSheetsInBarrels.push_back(nb);
        protein.numBarrels++;
    }

    void BetaAnnotator::setBarrel() {
        for(const auto& ssVec: protein.secStrVecs) {
            if (protein.chains[ssVec.chainIdx].residues[ssVec.begResIdx].selected
                && protein.chains[ssVec.chainIdx].residues[ssVec.endResIdx].selected
                && ssVec.barrelIdx != -1) {
                for(int i=ssVec.begResIdx; i<=ssVec.endResIdx; i++) {
                    if (protein.chains[ssVec.chainIdx].residues[i].selected
                        && any_cast<Tmdet::Types::Region>(protein.chains[ssVec.chainIdx].residues[i].temp.at("type")).isNotAnnotatedMembrane()
                        && numConnects(protein.chains[ssVec.chainIdx],i) > 0) {
                        protein.chains[ssVec.chainIdx].residues[i].temp.at("type") = std::any(Tmdet::Types::RegionType::BETA);
                    }
                }
                if (any_cast<Tmdet::Types::Region>(protein.chains[ssVec.chainIdx].residues[ssVec.begResIdx].temp.at("ztype")) !=
                        any_cast<Tmdet::Types::Region>(protein.chains[ssVec.chainIdx].residues[ssVec.endResIdx].temp.at("ztype"))) {
                            if (ssVec.begResIdx>0 
                                && protein.chains[ssVec.chainIdx].residues[ssVec.begResIdx-1].temp.contains("ztype")) {
                                protein.chains[ssVec.chainIdx].residues[ssVec.begResIdx-1].temp.at("type") = 
                                    protein.chains[ssVec.chainIdx].residues[ssVec.begResIdx-1].temp.at("ztype");
                            }
                            if (ssVec.endResIdx<protein.chains[ssVec.chainIdx].length-1
                                && protein.chains[ssVec.chainIdx].residues[ssVec.endResIdx+1].temp.contains("ztype")) {
                                protein.chains[ssVec.chainIdx].residues[ssVec.endResIdx+1].temp.at("type") = 
                                    protein.chains[ssVec.chainIdx].residues[ssVec.endResIdx+1].temp.at("ztype");
                            }
                }
            }
        }
        if (protein.numBarrels == 1 && numSheetsInBarrels[0] == 8) {
            protein.eachSelectedResidue(
                [&](Tmdet::VOs::Residue& residue) -> void {
                    if (any_cast<Tmdet::Types::Region>(residue.temp.at("type")).isNotAnnotatedMembrane()) {
                        residue.temp.at("type") = std::any(Tmdet::Types::RegionType::BETA);
                    }
                }
            );
        }
        protein.eachSelectedChain(
            [&](Tmdet::VOs::Chain& chain) -> void {
                int beg=0;
                int end=0;
                while(regionHandler.getNext<Tmdet::Types::Region>(chain,beg,end,"type")) {
                    if (chain.residues[beg].selected
                        && any_cast<Tmdet::Types::Region>(chain.residues[beg].temp.at("type")).isBeta() 
                        && beg > 0
                        && chain.residues[beg-1].selected
                        && end < chain.length-1
                        && end-beg < 3) {
                            regionHandler.replace(chain,beg,end-1,any_cast<Tmdet::Types::Region>(chain.residues[beg-1].temp.at("ztype")),"type");
                    }
                    beg=end;
                }
            }
        );
    }

    int BetaAnnotator::numConnects(Tmdet::VOs::Chain& chain, int pos) {
        int ret = 0;
        for(const auto& cr: Tmdet::Utils::NeighBors::get(chain.residues[pos])) {
            if (int v = protein.chains[cr.chainIdx].residues[cr.residueIdx].secStrVecIdx; v != -1) {
                if (protein.secStrVecs[v].barrelIdx != -1) {
                    ret++;
                }
            }
        }
        return ret;
    }

    void BetaAnnotator::detectBarrelInside(Tmdet::VOs::Chain& chain) {
        chain.eachSelectedResidue(
            [&](Tmdet::VOs::Residue& residue) -> void {
                if (any_cast<Tmdet::Types::Region>(residue.temp.at("type")).isNotAnnotatedMembrane()) {

                    residue.temp.at("type") = (any_cast<double>(residue.temp["hz"]) > 0 ? 
                        std::any(Tmdet::Types::RegionType::MEMBINS) :
                        std::any(residue.temp.at("ztype")));
                }
            }
        );
        int beg=0;
        int end=0;
        while(regionHandler.getNext<Tmdet::Types::Region>(chain,beg,end,"type")) {
            if (chain.residues[beg].selected
                && any_cast<Tmdet::Types::Region>(chain.residues[beg].temp.at("type")).isNotAnnotatedMembrane()
                && (beg == 0 || (beg >0 && chain.residues[beg-1].selected))
                && (end == chain.length -1 || (end < chain.length-1 && chain.residues[end].selected))
                && ((beg == 0 || (beg >0 && any_cast<Tmdet::Types::Region>(chain.residues[beg-1].temp.at("type")).isMembraneInside()))
                    || (end == chain.length-1 || (end < chain.length-1 && any_cast<Tmdet::Types::Region>(chain.residues[end].temp.at("type")).isMembraneInside())))
                ) {
                regionHandler.replace(chain,beg,end-1,Tmdet::Types::RegionType::MEMBINS,"type");
            }
            beg=end;
        }
        beg=0;
        end=0;
        while(regionHandler.getNext<Tmdet::Types::Region>(chain,beg,end,"type")) {
            if (chain.residues[beg].selected
                && any_cast<Tmdet::Types::Region>(chain.residues[beg].temp.at("type")).isNotMembrane()
                && end <= chain.length-1
                && beg > 0
                && chain.residues[beg-1].selected
                && chain.residues[end].selected
                && end-beg < 3
                && (any_cast<Tmdet::Types::Region>(chain.residues[beg-1].temp.at("type")).isMembraneInside()
                || any_cast<Tmdet::Types::Region>(chain.residues[end].temp.at("type")).isMembraneInside())) {
                regionHandler.replace(chain,beg,end-1,Tmdet::Types::RegionType::MEMBINS,"type");
            }
            beg=end;
        }
    }

    
}