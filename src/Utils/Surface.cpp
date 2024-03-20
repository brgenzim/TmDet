#include <iostream>
#include <array>
#include <math.h>
#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <ValueObjects/TmdetStruct.hpp>
#include <Utils/Surface.hpp>

using namespace std;

namespace Tmdet::Utils {


    void Surface::main() {
        for(auto& chain: tmdetVO.chains) {
            for(auto& residue: chain.residues) {
                residue.surface = 0.0;
                for(auto& atom: residue.atoms) {
                    calcSurfaceOfAtom(atom);
                    residue.surface += atom.surface;
                }
            }
        }
    }

    void Surface::createTempData() {
        for(auto& chain: tmdetVO.chains) {
            for(auto& residue: chain.residues) {
                for(auto& atom: residue.atoms) {
                    calcTempDataOfAtom(atom);
                }
            }
        }
    }

    void Surface::calcTempDataOfAtom(Tmdet::ValueObjects::Atom& atom) {
        atom.temp.insert({{"surfTemp",any_cast<surfTemp>(surfTemp())}});
        auto& gatom = atom.gemmi;
        for(auto m : tmdetVO.neighbors.find_neighbors(gatom, 0.1, 7.0)) {
            auto cra = m->to_cra(tmdetVO.gemmi.models[0]);
            auto gbtom = cra.atom;
            double dx=gatom.pos.x-gbtom->pos.x;
            double dy=gatom.pos.y-gbtom->pos.y;
            double d2 = (dx*dx)+(dy*dy);
			STEMP(atom,d2).emplace_back(d2);
			STEMP(atom,d).emplace_back(sqrt(d2));
            STEMP(atom,beta).emplace_back(atan2(dx,dy)+M_PI);
        }
    }

    void Surface::calcSurfaceOfAtom(Tmdet::ValueObjects::Atom& atom) {
        atom.surface = 0.0;
        auto& gatom = atom.gemmi;
        double vdwa = gatom.element.vdw_r() + SURF_PROBSIZE;
        for(double z=gatom.pos.z-vdwa+SURF_ZSLICE/2; z<gatom.pos.z+vdwa+SURF_ZSLICE/2; z+=SURF_ZSLICE) {
            double ra2 = vdwa * vdwa - (gatom.pos.z - z) * (gatom.pos.z - z);
            double ra = sqrt(ra2);
            bool ss = true;
            for(auto m : tmdetVO.neighbors.find_neighbors(gatom, 0.1, 7.0)) {
                if (ss) {
                    auto cra = m->to_cra(tmdetVO.gemmi.models[0]);
                    auto gbtom = cra.atom;
                    double vdwb = gbtom->element.vdw_r() + SURF_PROBSIZE;
                    double rb2 = vdwb * vdwb - (gbtom->pos.z - z) * (gbtom->pos.z - z);
                    double rb = sqrt(rb2);
                    double di = gatom.pos.dist(gbtom->pos);
                    double d2i = di * di;
                    if ( di < fabs(ra-rb)) {
                        if (ra<rb) {
                            ss = false;
                        }
                    }
                    else if (di<ra+rb && di>fabs(ra-rb)) {

                    }
                }
            }
        }
    }
}
