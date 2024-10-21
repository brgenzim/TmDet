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
        DEBUG_LOG("Processing: Annotator::detectSides()");
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


}