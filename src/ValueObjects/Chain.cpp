#include <string>
#include <vector>
#include <Types/Chain.hpp>
#include <ValueObjects/Chain.hpp>
#include <ValueObjects/Residue.hpp>
#include <ValueObjects/Region.hpp>
#include <gemmi/model.hpp>

namespace Tmdet::ValueObjects {

        void Chain::addStructure(const gemmi::Chain& _gemmi) {
                gemmi = _gemmi;
                id = _gemmi.name;
                for (const auto &residue : gemmi.residues) {
                    if (!residue.atoms.empty()) {
                        this->entityId = residue.entity_id;
                        break;
                    }
                }
        }

        gemmi::Vec3 Chain::centre() {
            gemmi::Vec3 ret(0,0,0);
            int nr = 0;
            for(auto& residue: residues) {
                if (!residue.atoms.empty()) {
                    ret += residue.centre();
                    nr++;
                }
            }
            if (nr>0) {
                ret /= nr;
            }
            return ret;
        }
}
