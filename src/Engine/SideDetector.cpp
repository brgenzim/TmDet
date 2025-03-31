// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#include <string>
#include <vector>
#include <any>
#include <Config.hpp>
#include <Engine/SideDetector.hpp>
#include <System/Logger.hpp>
#include <Types/Region.hpp>
#include <VOs/Membrane.hpp>

namespace Tmdet::Engine {

    void SideDetector::run() {
        auto membranes = protein.membranes;
        for(auto& membrane: membranes){
            membrane.halfThickness = 0.0;
        }
        setType("ztype",membranes);
        setType("type",protein.membranes);
        setDirection();
    }

    void SideDetector::end() {
        protein.eachSelectedResidue(
            [&](Tmdet::VOs::Residue& residue) -> void {
                residue.temp.erase("type");
                residue.temp.erase("ztype");
                residue.temp.erase("z");
                residue.temp.erase("hz");
            }
        );
    }

    void SideDetector::setType(std::string typeName, const std::vector<Tmdet::VOs::Membrane>& membranes) {
        setZs(membranes);
        protein.eachSelectedResidue(
            [&](Tmdet::VOs::Residue& residue) -> void {
                if (auto atom = residue.getCa(); atom != nullptr) {
                    residue.temp.try_emplace(typeName,std::any(getSideByZ(residue, getDistance(atom->pos))));
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

    Tmdet::Types::Region SideDetector::getSideByZ(Tmdet::VOs::Residue& residue, double z) const {
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
                hz = (z>o1?z1-z:z-z4);
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
            [&](Tmdet::VOs::Chain& chain) -> void {
                for(int i=0; i<chain.length; i++) {
                    if (chain.residues[i].selected) {
                        chain.residues[i].temp.try_emplace("direction",std::any(getResidueDirection(chain,i)));
                    }
                }
            }
        );
    }

    double SideDetector::getZForDirection(Tmdet::VOs::Chain& chain, int pos) {
        double ret = 0.0;
        if (pos>=0 && pos<chain.length && chain.residues[pos].temp.contains("z")) {
            ret = any_cast<double>(chain.residues[pos].temp.at("z"));
        }
        return ret;
    }

    double SideDetector::getResidueDirection(Tmdet::VOs::Chain& chain, int pos) {
        return getZForDirection(chain,pos+3)
                + getZForDirection(chain,pos+2)
                + getZForDirection(chain,pos+1)
                - getZForDirection(chain,pos-1)
                - getZForDirection(chain,pos-2)
                - getZForDirection(chain,pos-3);
    }
}
