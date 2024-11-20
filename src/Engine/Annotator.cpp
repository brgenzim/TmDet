#include <any>
#include <gemmi/model.hpp>
#include <Config.hpp>
#include <Helpers/Vector.hpp>
#include <Engine/Annotator.hpp>
#include <System/Logger.hpp>
#include <Types/Region.hpp>
#include <Utils/SecStrVec.hpp>

namespace Tmdet::Engine {

    void Annotator::run() {
        DEBUG_LOG("Processing Annotator::run()");
        detectSides();
        detectAlphaHelices();
        detectBarrel();
        detectInterfacialHelices();
        detectReEntrantLoops();
        finalize();
        getRegions();
        DEBUG_LOG(" Processed Annotator::run()");
    }

    void Annotator::detectSides() {
        DEBUG_LOG("Processing: Annotator::detectSides()");
        doubleMembrane = (protein.membranes.size() == 2);
        setZs();
        for(auto& chain: protein.chains) {
            for(auto& residue: chain.residues) {
                if (auto atom = residue.getCa(); atom != nullptr) {
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

    void Annotator::getRegions() {
        DEBUG_LOG("Processing: Annotator::getRegions()");
        for(auto& chain: protein.chains) {
            if (chain.selected) {
                int begin = 0;
                int end = 0;
                while(getNextRegion(chain,begin,end)) {
                    if (end - begin > 1) {
                        chain.regions.emplace_back(begin,end-1,std::any_cast<Tmdet::Types::Region>(chain.residues[begin].temp.at("type")));
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
            auto alphaVecs = getCrossingAlphas(membrane);
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
            auto betaVecs = getCrossingBetas(membrane);
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
            auto alphaVecs = getParallelAlphas(membrane);
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

    void Annotator::detectReEntrantLoops() {
        protein.eachSelectedChain(
            [&](Tmdet::ValueObjects::Chain& chain) -> void {
                int begin = 0;
                int end = 0;
                while(getNextRegion(chain,begin,end)) {
                    DEBUG_LOG("detectReEntrantLoops: {} {} {}",chain.id,begin,end);
                    if (end - begin > TMDET_REENTRANT_LOOP_MIN_LENGTH 
                        && begin > 1
                        && end < chain.length -1
                        && any_cast<Tmdet::Types::Region>(chain.residues[begin].temp.at("type")) == Tmdet::Types::RegionType::MEMB
                        && checkLoopEnds(chain,begin,end)
                        && hasOneHelix(chain,begin,end)
                        ) {
                            DEBUG_LOG("Loop found: {} {} {}",chain.id,begin,end);
                            for (int i = begin; i<= end; i++) {
                                chain.residues[i].temp.at("type") = std::any(Tmdet::Types::RegionType::LOOP);
                            }
                    }
                    else {
                        DEBUG_LOG("Not loop: {} {} {}",chain.id,begin,end);
                    }
                    begin = end;
                }
            }
        );
    }

    bool Annotator::checkLoopEnds(Tmdet::ValueObjects::Chain& chain, int begin, int end) {
        auto ca1 = chain.residues[begin].getCa();
        auto ca2 = chain.residues[end].getCa();
        DEBUG_LOG("checkLoopEnds: {} {}",ca1->pos.z,ca2->pos.z);
        return (std::abs(ca1->pos.z-ca2->pos.z) < 6 
                    && ca1->pos.z * ca2->pos.z > 0
                    && std::abs(ca1->pos.z) > TMDET_MEMBRANE_MIN_HALFTHICKNESS);
    }
        
    bool Annotator::hasOneHelix(Tmdet::ValueObjects::Chain& chain, int begin, int end) {
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

    void Annotator::finalize() {

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

    std::vector<Tmdet::ValueObjects::SecStrVec> Annotator::getCrossingAlphas(Tmdet::ValueObjects::Membrane& membrane) {
        std::vector<Tmdet::ValueObjects::SecStrVec> ret;
        for (auto& vector : protein.secStrVecs) {
            if (vector.type.isAlpha() && checkCross(vector, membrane)) {
                ret.emplace_back(vector);
            }
        }
        return ret;
    }

    std::vector<Tmdet::ValueObjects::SecStrVec> Annotator::getParallelAlphas(Tmdet::ValueObjects::Membrane& membrane) {
        std::vector<Tmdet::ValueObjects::SecStrVec> ret;
        for (auto& vector : protein.secStrVecs) {
            if (vector.type.isAlpha() && checkParallel(vector, membrane)) {
                ret.emplace_back(vector);
            }
        }
        return ret;
    }

    std::vector<Tmdet::ValueObjects::SecStrVec> Annotator::getCrossingBetas(Tmdet::ValueObjects::Membrane& membrane) {
        std::vector<Tmdet::ValueObjects::SecStrVec> ret;
        for (auto& vector : protein.secStrVecs) {
            if (vector.type.isBeta() && checkCross(vector, membrane)) {
                ret.emplace_back(vector);
            }
        }
        return ret;
    }

    bool Annotator::checkCross(Tmdet::ValueObjects::SecStrVec& vec, Tmdet::ValueObjects::Membrane& membrane) const {
        bool resultUp = false;
        bool resultDown = false;

        if (membrane.type == Tmdet::Types::MembraneType::PLAIN) {
            resultUp = Tmdet::Helpers::Vector::simplifiedDoesVectorCrossPlane(
                    vec.begin.z,vec.end.z,membrane.origo + 3);
            resultDown = Tmdet::Helpers::Vector::simplifiedDoesVectorCrossPlane(
                    vec.begin.z,vec.end.z,membrane.origo - 3);
        } else if (membrane.type == Tmdet::Types::MembraneType::CURVED) {
            resultUp = Tmdet::Helpers::Vector::doesVectorCrossSphere(
                    vec.begin, vec.end, gemmi::Vec3(0,0,membrane.origo), membrane.sphereRadius + 5.0);
            resultDown = Tmdet::Helpers::Vector::doesVectorCrossSphere(
                    vec.begin, vec.end, gemmi::Vec3(0,0,membrane.origo), membrane.sphereRadius - 5.0);
        } 
        DEBUG_LOG("==>{}:{}",resultUp,resultDown);
        return resultUp && resultDown;
    }

    bool Annotator::checkParallel(Tmdet::ValueObjects::SecStrVec& vec, Tmdet::ValueObjects::Membrane& membrane) const {    
        double sign = vec.begin.z / std::abs(vec.begin.z);
        return (std::abs(sign * vec.begin.z - membrane.halfThickness) < 8
                && std::abs(sign * vec.end.z - membrane.halfThickness) < 8
        );
    }


}