#include <iostream>
#include <array>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <filesystem>
#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <Config.hpp>
#include <System/FilePaths.hpp>
#include <Types/Residue.hpp>
#include <ValueObjects/Protein.hpp>
#include <Utils/Surface.hpp>

using namespace std;
namespace StructVO = Tmdet::ValueObjects;

namespace Tmdet::Utils {

    void surfaceCache::proteinFromCache(Tmdet::ValueObjects::Protein& protein) {
        unsigned int i =0;
        for(auto& c : protein.chains) {
            if (c.selected) {
                chainFromCache(c,i);
            }
        }
    }

    void surfaceCache::chainFromCache(Tmdet::ValueObjects::Chain& chain, unsigned int& index) {
        for (auto& r : chain.residues) {
            r.surface = 0.0;
            for(auto& a : r.atoms) {
                a.surface = cache[index];
                r.surface += a.surface;
                index++;
            }
        }
    }

    void surfaceCache::proteinToCache(const Tmdet::ValueObjects::Protein& protein) {
        cache.clear();
        for( const auto& c : protein.chains) {
            if (c.selected) {
                chainToCache(c);
            }
        }
    }

    void surfaceCache::chainToCache(const Tmdet::ValueObjects::Chain& chain) {
        for (const auto& r : chain.residues) {
            for(const auto& a : r.atoms) {
                cache.emplace_back(a.surface);
            }
        }
    }

    bool surfaceCache::read(Tmdet::ValueObjects::Protein& protein) {
        logger.debug("Processing surfaceCache::read()");
        std::string hash = protein.hash();
        std::string dir = Tmdet::System::FilePaths::cache(hash);
        std::string path = dir + "/" + hash + "_" + protein.code + ".bin";
        std::ifstream file(path, ios::binary);
        if (!file.is_open()) {
            logger.warn("No data in surface cache for {}",protein.code);
            return false;
        }
        auto size = cache.size();
        file.read(reinterpret_cast<char*>(&size), sizeof(size));
        cache.resize(size);
        file.read(reinterpret_cast<char*>(&cache[0]), sizeof(cache));
        file.close();
        proteinFromCache(protein);
        logger.debug(" Processed surfaceCache::read()");
        return true;
    }

    void surfaceCache::write(const Tmdet::ValueObjects::Protein& protein) {
        logger.debug("Processing surfaceCache::write()");
        proteinToCache(protein);
        std::string hash = protein.hash();
        std::string dir = Tmdet::System::FilePaths::cache(hash);
        std::filesystem::create_directories(dir);
        std::string path = dir + "/" + hash + "_" + protein.code + ".bin";
        std::ofstream file(path, ios::binary);
        if (!file.is_open()) {
            logger.warn("Could not write surface cache. Path: {}",path);
            return;
        }
        auto size = cache.size();
        file.write(reinterpret_cast<const char*>(&size), sizeof(size));
        file.write(reinterpret_cast<const char*>(&cache), sizeof(cache));
        std::ostream_iterator<double> out_itr(file);
        std::ranges::copy(cache.begin(), cache.end(), out_itr);
        file.close();
        logger.debug(" Processed surfaceCache::write()");
    }
    
    void Surface::run() {
        surfaceCache cache;
        if (noCache || !cache.read(protein) ) {
            initTempData();
            setContacts();
            cache.write(protein);
        }
    }

    void Surface::initTempData() {
        logger.debug("Processing Surface::initTempData()");
        const double probSize = std::stof(environment.get("TMDET_SURF_PROBSIZE",DEFAULT_TMDET_SURF_PROBSIZE));
        for(auto& chain: protein.chains) {
            for(auto& residue: chain.residues) {
                for(auto& atom: residue.atoms) {
                    Tmdet::Types::Residue residueType = Tmdet::Types::ResidueType::getResidue(residue.gemmi.name);
                    double vdw = probSize;
                    if (residueType.atoms.contains(atom.gemmi.name)) {
                        vdw += residueType.atoms.at(atom.gemmi.name).atom.vdw;
                    } else {
                        vdw += Types::AtomType::DEFAULT_VDW;
                    }
                    atom.temp.insert({"vdw",any_cast<double>(vdw)});
                }
            }
        }
        logger.debug(" Processed Surface::initTempData()");
    }

    void Surface::setContacts() {
        logger.debug("Processing Surface::setContacts()");
        for(auto& chain: protein.chains) {
            for(auto& residue: chain.residues) {
                residue.surface = 0.0;
                for(auto& atom: residue.atoms) {
                    setContactsOfAtom(atom);
                    residue.surface += atom.surface;
                }
            }
        }
        logger.debug(" Processed Surface::setContacts()");
    }

    void Surface::setContactsOfAtom(Tmdet::ValueObjects::Atom& a_atom) {
        surfTemp st;
        for(auto m : protein.neighbors.find_neighbors(a_atom.gemmi, 0.1, 7.0)) {
            auto& b_atom = protein.chains.at(m->chain_idx).residues.at(m->residue_idx).atoms.at(m->atom_idx);
            double dist = a_atom.gemmi.pos.dist(b_atom.gemmi.pos);
            if ( dist < VDW(a_atom) + VDW(b_atom)) {
                setNeighbor(a_atom,b_atom,st);
            }
        }
        calcSurfaceOfAtom(a_atom, st);
    }

    void Surface::setNeighbor(const Tmdet::ValueObjects::Atom& a_atom, const Tmdet::ValueObjects::Atom& b_atom, surfTemp& st) {
        double dx=a_atom.gemmi.pos.x-b_atom.gemmi.pos.x;
        double dy=a_atom.gemmi.pos.y-b_atom.gemmi.pos.y;
        double d2 = (dx*dx)+(dy*dy);
        double d = sqrt(d2);
        double beta = atan2(dx,dy)+M_PI;
        surfNeighbor sn = {b_atom, d, d2, beta};
        st.neighbors.emplace_back(sn);
    }

    void Surface::calcSurfaceOfAtom(Tmdet::ValueObjects::Atom& a_atom, surfTemp& st) {
        const double zSlice = std::stof(environment.get("TMDET_SURF_ZSLICE",DEFAULT_TMDET_SURF_ZSLICE));
        a_atom.surface = 0.0;
        const auto& a_gatom = a_atom.gemmi;
        double vdwa = VDW(a_atom);
        for(double z=a_gatom.pos.z-vdwa+zSlice/2; z<a_gatom.pos.z+vdwa/*+SURF_ZSLICE/2*/; z+=zSlice) {
            a_atom.surface += calcSumArcsOfAtom(a_atom,st,calcArcsOfAtom(a_atom,st,z)) * zSlice;
        }
        a_atom.surface *= vdwa;
    }

    bool Surface::calcArcsOfAtom(Tmdet::ValueObjects::Atom& a_atom, surfTemp& st, double z) {
        double vdwa = VDW(a_atom);
        double ra2 = vdwa * vdwa - (a_atom.gemmi.pos.z - z) * (a_atom.gemmi.pos.z - z);
        double ra = sqrt(ra2);
        bool ss = true;
        st.arc1.clear();
        st.arc2.clear();
        st.sorted.clear();
        for(auto& neighbor: st.neighbors) {
            if (ss) {
                auto& b_atom = neighbor.atom;
                const auto& b_gatom = b_atom.gemmi;
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

    double Surface::calcSumArcsOfAtom(const Tmdet::ValueObjects::Atom& atom, surfTemp& st, bool ss) const {
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

    void Surface::setOutsideSurface() {
        boundingBox box;
        setBoundingBox(box);
        initFrame(box);
        setFrame(box);
        smoothFrame(box);
        findClosestAtoms(box);
        for (int z=0; z<(box.zmax-box.zmin); z++) {
	        for (int i=0; i<box.op; i++) {
	            if (box.closestAtoms[z][i]!=NULL) {
                    box.closestAtoms[z][i]->temp.at("outside") = any_cast<double>(box.closestAtoms[z][i]->surface);
                }
            }
        }
	}

    void Surface::setBoundingBox(boundingBox& box) const {
        box.xmin=10000; box.xmax=-10000;        
	    box.ymin=10000; box.ymax=-10000;        
        box.zmin=10000; box.zmax=-10000;        
	    
        for(const auto& chain: protein.chains) {
            for(const auto& residue: chain.residues) {
                for(const auto& atom: residue.atoms) {
                    box.xmin = MIN(box.xmin,atom.gemmi.pos.x);
                    box.xmax = MAX(box.xmax,atom.gemmi.pos.x);
                    box.ymin = MIN(box.ymin,atom.gemmi.pos.y);
                    box.ymax = MAX(box.ymax,atom.gemmi.pos.y);
                    box.zmin = MIN(box.zmin,atom.gemmi.pos.z);
                    box.zmax = MAX(box.zmax,atom.gemmi.pos.z);
		        }
 	        }
        }
	    box.xmin=ceil(box.xmin-1)-2; box.xmax=ceil(box.xmax)+2;
	    box.ymin=ceil(box.ymin-1)-2; box.ymax=ceil(box.ymax)+2;
        box.zmin=ceil(box.zmin-1)-2; box.zmax=ceil(box.zmax)+2;

        box.l1 = (box.xmax - box.xmin) + 1; 
	    box.l2 = 2 * (box.xmax - box.xmin) + 2;
	    box.l3 = 2 * (box.xmax - box.xmin) + 2 + (box.ymax - box.ymin) + 1;	
	    box.op = 2 * (box.xmax - box.xmin) + 2 * (box.ymax - box.ymin) + 2;
    }

    void Surface::initFrame(boundingBox& box) {
        for (int i=0; i<(box.zmax-box.zmin)+1; i++) {
            box.closestAtoms.emplace_back(vector<Tmdet::ValueObjects::Atom *>(box.op));
            box.frame.emplace_back(vector<gemmi::Position>(box.op));
        }
        box.closestDist = vector<double>(box.zmax-box.zmin+1);
    
        for (int z=0; z < (box.zmax-box.zmin)+1; z++) {
            int i,k;
        	for (i=0, k=0; i<box.l1; i++,k++)  {
            	box.frame[z][i].x = k+box.ymin;
                box.frame[z][i].y = box.ymax;
                box.frame[z][i].z = 0;
            }
            for (i=box.l1, k=0; i<box.l2; i++,k++) {
            	box.frame[z][i].x = k+box.ymin;
                box.frame[z][i].y = box.ymin;
                box.frame[z][i].z = 0;
            }
            for (i=box.l2, k=0; i<box.l3; i++,k++) {
            	box.frame[z][i].x = k+box.xmin;
                box.frame[z][i].y = box.ymax;
                box.frame[z][i].z = 0;
            }
            for (i=box.l3, k=0; i<box.op; i++,k++) {
            	box.frame[z][i].x = k+box.xmin;
                box.frame[z][i].y = box.ymin;
                box.frame[z][i].z = 0;
            }
        }
    }

    void Surface::setFrame(boundingBox& box) {
        for(auto& chain: protein.chains) {
            for(auto& residue: chain.residues) {
                for(auto& atom: residue.atoms) {
                    if (atom.temp.find("outside") == atom.temp.end()) {
		                atom.temp.insert({"outside",any_cast<double>(0.0)});
                    }
                    else {
                        atom.temp.at("outside") = any_cast<double>(0.0);
                    }
                    double x = floor(atom.gemmi.pos.x);
                    double y = floor(atom.gemmi.pos.y);  	    
                    int z = floor(atom.gemmi.pos.z) - box.zmin;
                
                    int k= (int)(x-box.xmin);
                    box.frame[z][k].x = x;
                    box.frame[z][k].y = MIN(box.frame[z][k].y,y-2);
                    box.frame[z][k].z = (int)(atom.gemmi.pos.z); 
                    
                    box.frame[z][k+box.l1].x = x;
                    box.frame[z][k+box.l1].y = MAX(box.frame[z][k+box.l1].y,y+2);
                    box.frame[z][k+box.l1].z = (int)(atom.gemmi.pos.z);
                    
                    k=(int)(y-box.ymin);
                    box.frame[z][k+box.l2].x = MIN(box.frame[z][k+box.l2].x, x-2);
                    box.frame[z][k+box.l2].y = y;
                    box.frame[z][k+box.l2].z = (int)(atom.gemmi.pos.z);
                    
                    box.frame[z][k+box.l3].x = MAX(box.frame[z][k+box.l3].x, x+2);
                    box.frame[z][k+box.l3].y = y;
                    box.frame[z][k+box.l3].z = (int)(atom.gemmi.pos.z);
                }
            }
        }
    }

    void Surface::smoothFrame(boundingBox& box) {
        for (int z=0; z<(box.zmax-box.zmin)+1; z++) {
            int i,j;
        	for (i=0; i<box.l1; i++) {
                for (j=-2; j<2; j++) {
                    if ((i+j>=0) && (i+j<box.l1) && box.frame[z][i+j].z!=0) {
                        box.frame[z][i].y = MIN(box.frame[z][i+j].y, box.frame[z][i].y);
                    }
                }
            }
            for (i=box.l1; i<box.l2; i++) {
                for (j=-2; j<2; j++) {
                    if ((i+j>=box.l1) && (i+j<box.l2) && box.frame[z][i+j].z!=0) {
                        box.frame[z][i].y = MAX(box.frame[z][i+j].y, box.frame[z][i].y);
                    }
                }
            }
            for (i=box.l2; i<box.l3; i++) {
                for (j=-2; j<2; j++) {
                    if ((i+j>=box.l2) && (i+j<box.l3) && box.frame[z][i+j].z!=0) {
                        box.frame[z][i].x = MIN(box.frame[z][i+j].x, box.frame[z][i].x);
                    }
                }
            }
            for (i=box.l3; i<box.op; i++) {
                for (j=-2; j<2; j++) {
                    if ((i+j>=box.l3) && (i+j<box.op) && box.frame[z][i+j].z!=0) {
                        box.frame[z][i].x = MAX(box.frame[z][i+j].x, box.frame[z][i].x);
                    }
                }
            }
        }
    }

    void Surface::findClosestAtoms(boundingBox& box) {
        for (int i=0; i<box.op; i++) {
		    for (int j=0; j<box.zmax-box.zmin; j++) {
			    box.closestDist[j] = 10000000;
            }
            for(auto& chain: protein.chains) {
                for(auto& residue: chain.residues) {
                    for(auto& atom: residue.atoms) {
				        int z = (int)(atom.gemmi.pos.z-box.zmin);
				        double dist = 
                            (atom.gemmi.pos.x - box.frame[z][i].x) * 
                            (atom.gemmi.pos.x - box.frame[z][i].x) +
                            (atom.gemmi.pos.y - box.frame[z][i].y) *
                            (atom.gemmi.pos.y - box.frame[z][i].y);
				        if (dist < box.closestDist[z]) {
				        	box.closestAtoms[z][i] = &atom;
					        box.closestDist[z]=dist;
				        }
                    }
                }
			}	
		}
 	}
}
