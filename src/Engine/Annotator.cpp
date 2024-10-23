#include <any>
#include <gemmi/model.hpp>
#include <Config.hpp>
#include <Engine/Annotator.hpp>
#include <System/Logger.hpp>
#include <Types/Region.hpp>

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

    Tmdet::Types::Region Annotator::getSideByZ(double z) {
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
        Tmdet::ValueObjects::Region region = {
            (int)(beg+1),
            (int)chain.residues[beg].gemmi.seqid.num,
            chain.residues[beg].gemmi.seqid.icode,
            (int)chain.residues[beg].gemmi.label_seq,
            (int)(end+1),
            (int)chain.residues[end].gemmi.seqid.num,
            chain.residues[end].gemmi.seqid.icode,
            (int)chain.residues[end].gemmi.label_seq,
            std::any_cast<Tmdet::Types::Region>(chain.residues[beg].temp.at("type"))
        };
        chain.regions.emplace_back(region);
    }

    void Annotator::getRegions() {
        DEBUG_LOG("Processing: Annotator::getRegions()");
        for(auto& chain: protein.chains) {
            if (chain.selected) {
                unsigned int start = 0;
                for(unsigned int i = 0; i< chain.residues.size(); i++) {
                    if (std::any_cast<Tmdet::Types::Region>(chain.residues[start].temp.at("type")) 
                        != std::any_cast<Tmdet::Types::Region>(chain.residues[i].temp.at("type"))) {
                        storeRegion(chain,start,i-1);
                        start = i;
                    }
                }
                storeRegion(chain,start,chain.residues.size()-1);
            }
        }
        DEBUG_LOG(" Processed: Annotator::detectgetRegions()");
    }

    void Annotator::detectBarrel() {
        
    }

}