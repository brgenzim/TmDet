#include <any>
#include <gemmi/model.hpp>
#include <Config.hpp>
#include <Engine/Annotator.hpp>
#include <System/Logger.hpp>
#include <Types/Region.hpp>
#include <Utils/SecStrVec.hpp>

namespace Tmdet::Engine {

    void Annotator::detectSides() {
        DEBUG_LOG("Processing: Annotator::detectSides()");
        doubleMembrane = (protein.membranes.size() == 2);
        setZs();
        for(auto& chain: protein.chains) {
            for(auto& residue: chain.residues) {
                if (auto atom = residue.gemmi.get_ca(); atom != nullptr) {
                    residue.temp.try_emplace("type",std::any(getSideByZ(atom->pos.z)));
                }
            }
        }
        DEBUG_LOG(" Processed: Annotator::detectSides()");
    }

    void Annotator::setZs() {
        if (doubleMembrane) {
            if (protein.membranes[1].origo < 0) {
                z1 = protein.membranes[0].halfThickness;
                z2 = -protein.membranes[0].halfThickness;
                z3 = protein.membranes[1].origo+protein.membranes[1].halfThickness;
                z4 = protein.membranes[1].origo-protein.membranes[1].halfThickness;
            }
            else {
                z1 = protein.membranes[1].origo+protein.membranes[1].halfThickness;
                z2 = protein.membranes[1].origo-protein.membranes[1].halfThickness;
                z3 = protein.membranes[0].halfThickness;
                z4 = -protein.membranes[0].halfThickness;
            }
        }
        else {
            z1 = protein.membranes[0].halfThickness;
            z4 = -protein.membranes[0].halfThickness;
        }
    }

    Tmdet::Types::Region Annotator::getSideByZ(double z) const {
        Tmdet::Types::Region r;
        if (z > z1) {
            r = Tmdet::Types::RegionType::SIDE1;
        }
        else if (z < z4 ) {
            r = Tmdet::Types::RegionType::SIDE2;
        }
        else {
            if (doubleMembrane) {
                if (z > z2 || z < z3) {
                    r = Tmdet::Types::RegionType::MEMB;
                }
                else {
                    r = Tmdet::Types::RegionType::INTERMEMB;
                }
            }
            else {
                r = Tmdet::Types::RegionType::MEMB;
            }
        }
        return r;
    }

    void Annotator::storeRegion(Tmdet::ValueObjects::Chain& chain,unsigned int beg, unsigned int end) const {
        
        chain.regions.emplace_back(
            (int)(beg+1),
            (int)chain.residues[beg].gemmi.seqid.num,
            chain.residues[beg].gemmi.seqid.icode,
            (int)chain.residues[beg].gemmi.label_seq,
            (int)(end+1),
            (int)chain.residues[end].gemmi.seqid.num,
            chain.residues[end].gemmi.seqid.icode,
            (int)chain.residues[end].gemmi.label_seq,
            std::any_cast<Tmdet::Types::Region>(chain.residues[beg].temp.at("type"))
        );
    }

    void Annotator::getRegions() {
        DEBUG_LOG("Processing: Annotator::getRegions()");
        for(auto& chain: protein.chains) {
            if (chain.selected) {
                int begin = 0;
                int end = 0;
                while(getNextRegion(chain,begin,end)) {
                    if (end - begin > 1) {
                        storeRegion(chain,begin,end-1);
                    }
                    begin = end;
                }
            }
        }
        DEBUG_LOG(" Processed: Annotator::getRegions()");
    }

    bool Annotator::getNextRegion(Tmdet::ValueObjects::Chain& chain, int& begin, int& end) const {
        DEBUG_LOG("getNextRegion: {} {} {}",chain.id,begin,end);
        return (getNextDefined(chain, begin) && getNextSame(chain, begin, end));
    }

    bool Annotator::getNextDefined(Tmdet::ValueObjects::Chain& chain, int& begin) const {
        while(begin < (int)chain.residues.size() && !chain.residues[begin].temp.contains("type")) {
            begin++;
        }
        DEBUG_LOG("getNextDefined: {} {}",chain.id,begin);
        return (begin < (int)chain.residues.size());
    }

    bool Annotator::getNextSame(Tmdet::ValueObjects::Chain& chain, const int& begin, int& end) const {
        end = begin;
        while(end < (int)chain.residues.size() && chain.residues[end].temp.contains("type")
            && std::any_cast<Tmdet::Types::Region>(chain.residues[begin].temp.at("type")) 
                        == std::any_cast<Tmdet::Types::Region>(chain.residues[end].temp.at("type"))) {
            end++;
        }
        DEBUG_LOG("getNextSame: {} {} {}",chain.id,begin,end);
        return true;
    }

    void Annotator::detectAlphaHelices() {
        DEBUG_LOG("Processing Annotator::detectAlphaHelices()");
        for(auto& membrane: protein.membranes) {
            auto alphaVecs = ssVec.getCrossingAlphas(membrane);
            if (!alphaVecs.empty() ) {
                protein.type = Tmdet::Types::ProteinType::TM_ALPHA;
                for(const auto& vector: alphaVecs) {
                    replaceRegion(vector,Tmdet::Types::RegionType::HELIX);
                }
            }
            DEBUG_LOG(" #alphaVecs: {}",alphaVecs.size());
        }
        DEBUG_LOG(" Processed Annotator::detectAlphaHelices()");
    }

    void Annotator::detectBarrel() {
        DEBUG_LOG("Processing Annotator::detectBarrel()");
        for(auto& membrane: protein.membranes) {
            auto betaVecs = ssVec.getCrossingBetas(membrane);
            if (betaVecs.size() > 7) {
                protein.type = (protein.type == Tmdet::Types::ProteinType::TM_ALPHA?Tmdet::Types::ProteinType::TM_MIXED:Tmdet::Types::ProteinType::TM_BETA);
                for(const auto& vector: betaVecs) {
                    replaceRegion(vector,Tmdet::Types::RegionType::BETA);
                }
            }
        }
        DEBUG_LOG(" Processed Annotator::detectBarrel()");
    }

    void Annotator::detectInterfacialHelices() {
        DEBUG_LOG("Processing Annotator::detectInterfacialHelices()");
        for(auto& membrane: protein.membranes) {
            auto alphaVecs = ssVec.getParallelAlphas(membrane);
            if (!alphaVecs.empty() ) {
                for(const auto& vector: alphaVecs) {
                    for (int i = vector.begResIdx; i<= vector.endResIdx; i++) {
                        protein.chains[vector.chainIdx].residues[i].temp.at("type") = std::any(Tmdet::Types::RegionType::IFH);
                    }
                }
            }
            DEBUG_LOG(" #alphaVecs: {}",alphaVecs.size());
        }
        DEBUG_LOG(" Processed Annotator::detectInterfacialHelices()");
    }

    void Annotator::replaceRegion(const Tmdet::ValueObjects::SecStrVec& vector, Tmdet::Types::Region regionType) {
        DEBUG_LOG("Processing Annotator::replaceRegion({} {} {} --> {})",
            protein.chains[vector.chainIdx].id,vector.begResIdx,vector.endResIdx,regionType.name);
        for (int i = vector.begResIdx; i<= vector.endResIdx; i++) {
            if (std::any_cast<Tmdet::Types::Region>(protein.chains[vector.chainIdx].residues[i].temp.at("type")) == Tmdet::Types::RegionType::MEMB) {
                protein.chains[vector.chainIdx].residues[i].temp.at("type") = std::any(regionType);
            }
        }
        DEBUG_LOG(" Processed Annotator::replaceRegion()");
    }

}