#include <iostream>
#include <array>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <Types/Residue.hpp>
#include <ValueObjects/TmdetStruct.hpp>
#include <Utils/Surface.hpp>

using namespace std;

namespace Tmdet::Utils {


    void Surface::main() {
        initTempData();
        setContacts();
    }

    void Surface::initTempData() {
        for(auto& chain: tmdetVO.chains) {
            for(auto& residue: chain.residues) {
                for(auto& atom: residue.atoms) {
                    double vdw = Tmdet::Types::Residues.at(residue.gemmi.name).atoms.at(atom.gemmi.name).atom.vdw + SURF_PROBSIZE;
                    atom.temp.insert({{"vdw",any_cast<double>(vdw)}});
                }
            }
        }
    }

    void Surface::setContacts() {
        for(auto& chain: tmdetVO.chains) {
            for(auto& residue: chain.residues) {
                residue.surface = 0.0;
                for(auto& atom: residue.atoms) {
                    setContactsOfAtom(atom);
                    residue.surface += atom.surface;
                }
            }
        }
    }

    void Surface::setContactsOfAtom(Tmdet::ValueObjects::Atom& a_atom) {
        surfTemp st;
        for(auto m : tmdetVO.neighbors.find_neighbors(a_atom.gemmi, 0.1, 7.0)) {
            auto& b_atom = tmdetVO.chains.at(m->chain_idx).residues.at(m->residue_idx).atoms.at(m->atom_idx);
            double dist = a_atom.gemmi.pos.dist(b_atom.gemmi.pos);
            if ( dist < VDW(a_atom) + VDW(b_atom)) {
                setNeighbor(a_atom,b_atom,st);
            }
        }
        calcSurfaceOfAtom(a_atom, st);
    }

    void Surface::setNeighbor(Tmdet::ValueObjects::Atom& a_atom, Tmdet::ValueObjects::Atom& b_atom, surfTemp& st) {
        double dx=a_atom.gemmi.pos.x-b_atom.gemmi.pos.x;
        double dy=a_atom.gemmi.pos.y-b_atom.gemmi.pos.y;
        double d2 = (dx*dx)+(dy*dy);
        double d = sqrt(d2);
        double beta = atan2(dx,dy)+M_PI;
        surfNeighbor sn = {b_atom, d, d2, beta};
        st.neighbors.emplace_back(sn);
    }

    void Surface::calcSurfaceOfAtom(Tmdet::ValueObjects::Atom& a_atom, surfTemp& st) {
        a_atom.surface = 0.0;
        auto& a_gatom = a_atom.gemmi;
        double& vdwa = VDW(a_atom);
        for(double z=a_gatom.pos.z-vdwa+SURF_ZSLICE/2; z<a_gatom.pos.z+vdwa/*+SURF_ZSLICE/2*/; z+=SURF_ZSLICE) {
            a_atom.surface += calcSumArcsOfAtom(a_atom,st,calcArcsOfAtom(a_atom,st,z)) * SURF_ZSLICE;
        }
        a_atom.surface *= vdwa;
    }

    bool Surface::calcArcsOfAtom(Tmdet::ValueObjects::Atom& a_atom, surfTemp& st, double z) {
        double& vdwa = VDW(a_atom);
        double ra2 = vdwa * vdwa - (a_atom.gemmi.pos.z - z) * (a_atom.gemmi.pos.z - z);
        double ra = sqrt(ra2);
        bool ss = true;
        st.arc1.clear();
        st.arc2.clear();
        st.sorted.clear();
        for(auto& neighbor: st.neighbors) {
            if (ss) {
                auto& b_atom = neighbor.atom;
                auto& b_gatom = b_atom.gemmi;
                double vdwb = VDW(b_atom);
                if (fabs(b_gatom.pos.z -z) < vdwb) {    
                    double rb2 = vdwb * vdwb - (b_gatom.pos.z - z) * (b_gatom.pos.z - z);
                    double rb = sqrt(rb2);
                    if ( neighbor.d < fabs(ra-rb)) {
                        if (ra<rb) {
                            ss = false;
                        }
                    }
                    else if ( neighbor.d < ra+rb && neighbor.d > fabs(ra-rb)) {
                        double q = (neighbor.d2+ra2-rb2) / (2*neighbor.d*ra);
                        q = (q>1?1.0:q);
                        q = (q<-1?-1.0:q);
                        double alpha = acos(q);
                        double arc1 = neighbor.beta - alpha;
                        double arc2 = neighbor.beta + alpha;
                        arc1 = (arc1<0?arc1+2*M_PI:arc1);
                        arc2 = (arc2>2*M_PI?arc2-2*M_PI:arc2);
                        if (arc1 < arc2) {
                            st.arc1.emplace_back(arc1);
                            st.arc2.emplace_back(arc2);
                        }
                        else {
                            st.arc1.emplace_back(arc1);
                            st.arc2.emplace_back(2*M_PI);
                            st.arc1.emplace_back(0);
                            st.arc2.emplace_back(arc2);
                        }
                    }
                }
            }
        }
        return ss;
    }

    double Surface::calcSumArcsOfAtom(Tmdet::ValueObjects::Atom& atom, surfTemp& st, bool ss) {
        double arcsum = (ss?2.0*M_PI:0.0);
        int n = (int)(st.arc1.size());
        if (ss && n>0) {
            st.sorted = vector<int>(n);
            iota(st.sorted.begin(),st.sorted.end(),0);
            sort( st.sorted.begin(),st.sorted.end(), [&](int i,int j) {
                return st.arc1[i]<st.arc1[j];
            });
            for(int i=0; i<n;) {
                double beg = st.arc1.at(st.sorted.at(i));
                double end = st.arc2.at(st.sorted.at(i));
                int j;
                for(j=i+1; (j<n && st.arc1.at(st.sorted.at(j))<end); j++) {
                    if (st.arc2.at(st.sorted.at(j)) > end) {
                        end = st.arc2.at(st.sorted.at(j));
                    }
                }
                arcsum-=(end-beg);
                i=j;
            }
        }
        return arcsum;
    }
}
