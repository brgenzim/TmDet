#include <string>
#include <vector>
#include <any>
#include <Config.hpp>
#include <Engine/SideDetector.hpp>
#include <System/Logger.hpp>
#include <Types/Region.hpp>
#include <ValueObjects/Membrane.hpp>

namespace Tmdet::Engine {

    void SideDetector::run() {
        DEBUG_LOG("Processing: SideDetector::run()");
        auto membranes = protein.membranes;
        for(auto& membrane: membranes){
            membrane.halfThickness = 4.0;
        }
        setType("ttype",membranes);
        for(auto& membrane: membranes){
            membrane.halfThickness = 0.0;
        }
        setType("ztype",membranes);
        setType("type",protein.membranes);
        setDirection();
        DEBUG_LOG(" Processed: SideDetector::run()");
    }

    void SideDetector::end() {
        protein.eachSelectedResidue(
            [&](Tmdet::ValueObjects::Residue& residue) -> void {
                residue.temp.erase("type");
                residue.temp.erase("ttype");
                residue.temp.erase("ztype");
                residue.temp.erase("z");
                residue.temp.erase("hz");
            }
        );
    }

    void SideDetector::setType(std::string typeName, const std::vector<Tmdet::ValueObjects::Membrane>& membranes) {
        setZs(membranes);
        protein.eachSelectedResidue(
            [&](Tmdet::ValueObjects::Residue& residue) -> void {
                if (auto atom = residue.getCa(); atom != nullptr) {
                    residue.temp.try_emplace(typeName,std::any(getSideByZ(residue, atom->pos.z)));
                }
                else {
                    residue.temp.try_emplace(typeName,std::any(Tmdet::Types::RegionType::UNK));
                    if (!residue.temp.contains("z")) {
                        residue.temp.try_emplace("z",std::any(0.0));
                        residue.temp.try_emplace("hz",std::any(0.0));
                    }
                }
            }
        );
    }

    void SideDetector::setZs(const std::vector<Tmdet::ValueObjects::Membrane>& membranes) {
        if (membranes.size() == 2) {
            if (membranes[1].origo < 0) {
                z1 = membranes[0].halfThickness;
                z2 = -membranes[0].halfThickness;
                z3 = membranes[1].origo+membranes[1].halfThickness;
                z4 = membranes[1].origo-membranes[1].halfThickness;
                o1 = membranes[0].origo;
                o2 = membranes[1].origo;
            }
            else {
                z1 = membranes[1].origo+membranes[1].halfThickness;
                z2 = membranes[1].origo-membranes[1].halfThickness;
                z3 = membranes[0].halfThickness;
                z4 = -membranes[0].halfThickness;
                o1 = membranes[1].origo;
                o2 = membranes[0].origo;
            }
        }
        else {
            z1 = membranes[0].halfThickness;
            z4 = -membranes[0].halfThickness;
            o1 = o2 = membranes[0].origo;
        }
    }

    Tmdet::Types::Region SideDetector::getSideByZ(Tmdet::ValueObjects::Residue& residue, double z) const {
        Tmdet::Types::Region r;
        double rz; //relativ z coordinate from membrane central plane
        double hz=0; //distance from closest membrane surface
        if (z > z1) {
            r = Tmdet::Types::RegionType::SIDE1;
            rz = z - o1;
            hz = z - z1;
        }
        else if (z < z4 ) {
            r = Tmdet::Types::RegionType::SIDE2;
            rz = z - o2;
            hz = z4 - z;
        }
        else {
            if (protein.membranes.size() == 2) {
                if (z > z2 || z < z3) {
                    r = Tmdet::Types::RegionType::MEMB;
                    if (z > (z2+z3) / 2) {
                        hz = (z1-z<z-z2?z1-z:z-z2);
                    }
                    else {
                        hz = (z3-z<z-z4?z3-z:z-z4);
                    }
                }
                else {
                    r = Tmdet::Types::RegionType::INTERMEMB;
                }
                if (z < (z2+z3) / 2) {
                    rz = z - o2;
                    hz = z - z3;
                }
                else {
                    rz = z - o1;
                    hz = z2 - z;
                }
            }
            else {
                r = Tmdet::Types::RegionType::MEMB;
                rz = z - o1;
                hz = (z>0?z1-z:z-z4);
            }
        }
        if (!residue.temp.contains("z")) {
            residue.temp.try_emplace("z",std::any(rz));
        }
        if (!residue.temp.contains("hz")) {
            residue.temp.try_emplace("hz",std::any(hz));
        }
        else {
            residue.temp.at("hz") = std::any(hz);
        }
        return r;
    }

    void SideDetector::setDirection() {
        protein.eachSelectedChain(
            [&](Tmdet::ValueObjects::Chain& chain) -> void {
                chain.residues[0].temp.try_emplace("direction",std::any(0.0));
                chain.residues[1].temp.try_emplace("direction",std::any(0.0));
                chain.residues[2].temp.try_emplace("direction",std::any(0.0));
                chain.residues[chain.length-3].temp.try_emplace("direction",std::any(0.0));
                chain.residues[chain.length-2].temp.try_emplace("direction",std::any(0.0));
                chain.residues[chain.length-1].temp.try_emplace("direction",std::any(0.0));
                for(int i=3; i<chain.length-3; i++) {
                    chain.residues[i].temp.try_emplace("direction",std::any(getResidueDirection(chain,i)));
                    DEBUG_LOG("RES: chain:{} res:{} type:{} ss:{} z:{} z1:{} hz:{} z4:{} direction:{}",chain.id,chain.residues[i].authId,
                        any_cast<Tmdet::Types::Region>(chain.residues[i].temp.at("type")).code,
                        chain.residues[i].ss.code,
                        any_cast<double>(chain.residues[i].temp.at("z")),
                        z1,
                        any_cast<double>(chain.residues[i].temp.at("hz")),
                        z4,
                        any_cast<double>(chain.residues[i].temp.at("direction"))
                    );
                }
                
            }
        );
    }

    double SideDetector::getResidueDirection(Tmdet::ValueObjects::Chain& chain, int pos) {
        return any_cast<double>(chain.residues[pos+3].temp.at("z"))
                + any_cast<double>(chain.residues[pos+2].temp.at("z"))
                + any_cast<double>(chain.residues[pos+1].temp.at("z"))
                - any_cast<double>(chain.residues[pos-1].temp.at("z"))
                - any_cast<double>(chain.residues[pos-2].temp.at("z"))
                - any_cast<double>(chain.residues[pos-3].temp.at("z"));
    }
}
