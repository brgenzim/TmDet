// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#include <iostream>
#include <array>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <any>
#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <Config.hpp>
#include <Engine/Optimizer.hpp>
#include <Engine/Rotator.hpp>
#include <Helpers/Vector.hpp>
#include <Types/Residue.hpp>
#include <Types/Chain.hpp>
#include <VOs/Protein.hpp>
#include <VOs/Slice.hpp>
#include <Utils/Surface.hpp>


namespace Tmdet::Engine {

    Optimizer::~Optimizer() {
        DEBUG_LOG("Destroying Optimizer");
        end();
    }

    void Optimizer::init() {
        DEBUG_LOG("Processing: Optimizer::init()");
        massCentre = protein.centre();
        protein.eachSelectedResidue(
            [](Tmdet::VOs::Residue& residue) -> void {
                residue.temp.try_emplace("dist",std::any(0.0));
                residue.temp.try_emplace("straight",std::any(0.0));
                for(auto& atom: residue.atoms) {
                    atom.temp.try_emplace("dist",std::any(0.0));
                }
            }
        );
        DEBUG_LOG(" Processed: Optimizer::init()");
    }

    void Optimizer::end() {
        DEBUG_LOG("Processing: Optimizer::end()");
        protein.eachSelectedResidue(
            [](Tmdet::VOs::Residue& residue) -> void {
                residue.temp.erase("dist");
                residue.temp.erase("straight");
                for(auto& atom: residue.atoms) {
                    atom.temp.erase("dist");
                }
            }
        );
        DEBUG_LOG(" Processed: Optimizer::end()");
    }

    void Optimizer::clear() {
        bestQ = 0;
    }

    void Optimizer::setDistances() {
        DEBUG_LOG("Processing: Optimizer::setDistances()");
        protein.eachSelectedResidue(
            [&](Tmdet::VOs::Residue& residue) -> void {
                bool hasCA = false;
                double d = 0.0;
                gemmi::Vec3 ca;
                for(auto& atom: residue.atoms) {
                    d = distance(atom.gemmi.pos);
                    atom.temp.at("dist") = std::any(d);
                    ca = atom.gemmi.pos;
                    if (atom.gemmi.name == "CA") {
                        residue.temp.at("dist") = std::any(d);
                        hasCA = true;
                    }
                }
                if (!hasCA) {
                    residue.temp.at("dist") = std::any(d);
                }
            }
        );
        DEBUG_LOG(" Processed: Optimizer::setDistances()");
    }

    void Optimizer::setStraight() {
        DEBUG_LOG("Processing: Optimizer::setStraight()");
        protein.eachSelectedResidue(
            [&](Tmdet::VOs::Residue& residue) -> void {
                if (residue.idx>2 
                    && residue.idx<protein.chains[residue.chainIdx].length-3
                    && protein.chains[residue.chainIdx].residues[residue.idx-3].selected
                    && protein.chains[residue.chainIdx].residues[residue.idx+3].selected) {
                    double d1 = any_cast<double>(RES(residue,-3).temp.at("dist"));
                    double d2 = any_cast<double>(residue.temp.at("dist"));
                    double d3 = any_cast<double>(RES(residue,+3).temp.at("dist"));
                    auto ca1 = RES(residue,-3).getCa();
                    auto ca3 = RES(residue,+3).getCa();
                    double q1 = d2-d1;
                    double q2 = d3-d2;
                    double s = 1.0;
                    if (ca1 != nullptr && ca3 != nullptr and q1*q2>0) {
                        if (double d = ca1->pos.dist(ca3->pos); d >1.0) {
                            s = std::abs(d3-d1)/d;
                        }
                    }
                    residue.temp.at("straight") = std::any(s);
                }
            }
        );
        DEBUG_LOG(" Processed: Optimizer::setStraight()");
    }

    void Optimizer::setBoundaries() {
        DEBUG_LOG("Processing: Optimizer::setBoundaries()");
        minZ = 1e30;
        maxZ = -1e30;
        protein.eachSelectedResidue(
            [&](Tmdet::VOs::Residue& residue) -> void {
                auto dist = any_cast<double>(residue.temp.at("dist"));
                minZ = (dist < minZ ? dist : minZ);
                maxZ = (dist > maxZ ? dist : maxZ);
            }
        );
        minZ-=1; maxZ+=1;
        DEBUG_LOG("Box size: {} {}: {}",minZ,maxZ,(unsigned int)(maxZ-minZ));
        slices.clear();
        slices.resize((unsigned int)(maxZ-minZ));
        DEBUG_LOG(" Processed: Optimizer::setBoundaries()");
    }

    void Optimizer::sumupSlices() {
        DEBUG_LOG("Processing: Optimizer::sumupSlices({})",slices.size());
        protein.eachSelectedResidue(
            [&](Tmdet::VOs::Residue& residue) -> void {
                auto sliceIndex = (unsigned int)(any_cast<double>(residue.temp.at("dist")) - minZ);
                if (protein.secStrVecs.size()<10 || protein.chains[residue.chainIdx].type == Tmdet::Types::ChainType::LOW_RES) {
                    slices[sliceIndex].straight += any_cast<double>(residue.temp.at("straight"));
                    slices[sliceIndex].numCa++;
                }
                if (protein.chains[residue.chainIdx].type == Tmdet::Types::ChainType::LOW_RES) {
                    slices[sliceIndex].surf++;
                    slices[sliceIndex].apol += (residue.type.hsc + 12.3 ) / 16.0;
                }
                else {
                    for(const auto& atom: residue.atoms) {
                        slices[sliceIndex].surf += atom.outSurface * (residue.ss.isBeta()?1.1:0.9);
                        if (residue.type.atoms.contains(atom.gemmi.name)) {
                            slices[sliceIndex].apol += residue.type.apol * atom.outSurface 
                                * (residue.ss.isBeta()?1.1:0.9);
                        }
                    }
                }
            }
        );
        for (auto& vector: protein.secStrVecs) {
            if (protein.chains[vector.chainIdx].selected
                && protein.chains[vector.chainIdx].residues[vector.begResIdx].selected
                && protein.chains[vector.chainIdx].residues[vector.endResIdx].selected) {
                double cosAngle = (type=="Plane"?std::abs(Tmdet::Helpers::Vector::cosAngle(
                                    normal,vector.end - vector.begin)):1);
                cosAngle = (cosAngle<0.3?-3:cosAngle);
                int d1 = distance(vector.begin);
                int d2 = distance(vector.end);
                int dbeg = (d1<d2?d1:d2) - minZ; 
                int dend = (d1>d2?d1:d2) - minZ; 
                dbeg=(dbeg<0?0:dbeg);
                dend=(dend<0?0:dend);
                if (dbeg>=(int)slices.size()) {
                    dbeg=(int)slices.size()-1;
                }
                if (dend>=(int)slices.size()) {
                    dend=(int)slices.size()-1;
                }
                DEBUG_LOG("sumupSliceVectors: beg:{} end:{} size:{}",dbeg,dend,slices.size());
                for(int i=dbeg; i<= dend; i++) {
                    slices[i].straight+= (vector.type.isBeta()?1:cosAngle);
                    slices[i].ifh += ((vector.type.isAlpha() && cosAngle < 0.2)?1:0);
                    slices[i].numCa++;
                }
                slices[dbeg].ssEnd++;
                slices[dend].ssEnd++;
            }
        }
        DEBUG_LOG(" Processed: Optimizer::sumupSlices()");
    }

    double Optimizer::divide(double numerator, double denominator) {
        double q= ( denominator<1e-5?0.0:numerator/denominator);
        return (q>1?1.0:q);
    }

    void Optimizer::smoothQValues() {
        std::string n="***************************************************************************";
        std::string m="mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm";
        std::string b="bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb";
        auto s = slices.size();
        for(unsigned long int i = 0; i<s; i++) {
            int k=0;
            double smoothedStraight = 0.0;
            double smoothedApol = 0.0;
            double smoothedSsEnd = 0.0;
            double smoothedSurf = 0.0;
            double smoothedCa = 0.0;
            double smoothedIfh = 0.0;
            for (int j=-2; j<=2; j++) {
                if (j+(int)i>=0 && j+i<s) {
                    smoothedStraight += slices[j+i].straight;
                    smoothedApol += slices[j+i].apol;
                    smoothedSsEnd += slices[j+i].ssEnd;
                    smoothedSurf += slices[i+j].surf;
                    smoothedIfh += slices[i+j].ifh;
                    smoothedCa += slices[i+j].numCa;
                    k++;
                }
            }
            smoothedStraight /= k;
            smoothedApol /= k;
            smoothedSsEnd /= k;
            smoothedSurf /= k;
            smoothedIfh /= k;
            smoothedCa /= k;
            smoothedStraight = divide(smoothedStraight,smoothedCa);
            smoothedApol = divide(smoothedApol, smoothedSurf);
            smoothedSsEnd = divide(smoothedSsEnd, smoothedCa);
            slices[i].smoothedIfh = divide(smoothedIfh, smoothedCa);

            
            slices[i].rawQ = 130.0 * smoothedStraight * (1.0 - slices[i].smoothedIfh) * (1.0 - smoothedSsEnd) * smoothedApol;

             DEBUG_LOG("{} straight:{:5.2f} apol:{:5.2f} surf:{:5.2f} end:{:5.2f} ifh:{:5.2f} q:{:5.2f}",
                i,smoothedStraight,smoothedApol,slices[i].surf,smoothedSsEnd,slices[i].smoothedIfh,slices[i].rawQ);
        }
        for(unsigned long int i = 0; i<s; i++) {
            int k=0;
            double q = 0.0;
            for (int j=-8; j<=8; j++) {
                if (j+(int)i>=0 && j+i<s) {
                    q += slices[j+i].rawQ;
                    k++;
                }
            }
            q/=k;
            slices[i].qValue = q;
            
            if (q>TMDET_MINIMUM_QVALUE) {
                DEBUG_LOG("{}:{:.2f} {}{}{}",
                    i, q,
                    n.substr(0,TMDET_MEMBRANE_QVALUE),
                    m.substr(0,(int)(TMDET_MINIMUM_QVALUE-TMDET_MEMBRANE_QVALUE)),
                    b.substr(0,(int)(q-TMDET_MINIMUM_QVALUE+1)));
            }
            else if (q>TMDET_MEMBRANE_QVALUE) {
                DEBUG_LOG("{}:{:.2f} {}{}",
                    i, q,
                    n.substr(0,TMDET_MEMBRANE_QVALUE),
                    m.substr(0,(int)(q-TMDET_MEMBRANE_QVALUE+1)));
            }
            else if (q>0) {
                DEBUG_LOG("{}:{:.2f} {}",i,q,n.substr(0,(int)(q)));
            }
            else {
                DEBUG_LOG("{}:{:.2f}",i,q);
            }
        }
    }

    void Optimizer::checkBestSlice() {
        auto s = slices.size();
        int minHW = TMDET_MEMBRANE_MIN_HALFTHICKNESS;
        if ((int)s > 2 * minHW) {
            for(int i = minHW; i < (int)s - minHW; i++) {
                int minz;
                int maxz;
                double q = getWidth(i,minz,maxz);
                if (i-minz > minHW && maxz -i > minHW && q>TMDET_MINIMUM_QVALUE) {
                    if (q>bestQ) {
                        bestQ = q;
                        bestMinZ = minz;
                        bestMaxZ = maxz;
                        bestNormal = normal;
                        bestSlices = slices;
                        setBestOrigo(minz,maxz);
                        DEBUG_LOG("BestWidth: {} {} {} {} {}",minz,i,maxz,q,bestQ);
                    }
                }
            }
        }
    }

    double Optimizer::getWidth(const int z, int& minz, int& maxz) {
        minz = z;
        while(minz>0 
            && slices[minz].qValue > TMDET_MINIMUM_QVALUE 
            && slices[minz].smoothedIfh < 0.07
            && z-minz < TMDET_MEMBRANE_MAX_HALFTHICKNESS) {
            minz--;
        }
        maxz = z+1;
        while(maxz<(int)slices.size()
            && slices[maxz].qValue > TMDET_MINIMUM_QVALUE
            && slices[maxz].smoothedIfh < 0.07
            && maxz-z < TMDET_MEMBRANE_MAX_HALFTHICKNESS) {
            maxz++;
        }
        if (maxz-minz>5) {
                while(minz>0 
                    && slices[minz].qValue > TMDET_MEMBRANE_QVALUE  
                    && slices[minz].smoothedIfh < 0.07
                    && z-minz < TMDET_MEMBRANE_MAX_HALFTHICKNESS) {
                minz--;
            }
            while(maxz<(int)slices.size() 
                    && slices[maxz].qValue > TMDET_MEMBRANE_QVALUE
                    && slices[maxz].smoothedIfh < 0.07
                    && maxz-z < TMDET_MEMBRANE_MAX_HALFTHICKNESS) {
                maxz++;
            }
        }
        return slices[z].qValue;
    }

    void Optimizer::testMembraneNormalOne() {
        setDistances();
        setStraight();
        setBoundaries();
        sumupSlices();
        smoothQValues();
        checkBestSlice();
    }

    void Optimizer::setNormal(gemmi::Vec3 _normal) {
        normal = _normal;
    }

    void Optimizer::searchForMembraneNormal() {
        Tmdet::Engine::Rotator rotator;
        while(rotator.next(normal)) {
            testMembraneNormal();
        }
    }

    bool Optimizer::isTransmembrane() const {
        DEBUG_LOG("Processing Optimizer::isTransmembrane(): {}",bestQ);
        return bestQ > TMDET_MINIMUM_QVALUE;
    }

    void Optimizer::setMembranesToProtein() {
        DEBUG_LOG("Processing Optimizer::setMembranesToProtein()");
        protein.qValue = bestQ;
        if (!isTransmembrane()) {
            DEBUG_LOG("Processing Optimizer::setMembranesToProtein(): not transmembrane protein");
            protein.tmp = false;
            return;
        }
        normal = bestNormal;
        bestQ = 0;
        testMembraneNormal();
        Tmdet::VOs::Membrane membrane;
        protein.membranes.clear();
        int i=0;
        int maxM = (protein.forceSingleMembrane?1:2);
        while(getMembrane(membrane,i) && i < maxM) {
            protein.membranes.push_back(membrane);
            i++;
        }
        if (!protein.membranes.empty()) {
            protein.tmp = true;
            setProteinTMatrix(massCentre);
        }
        DEBUG_LOG(" Processed Optimizer::setMembranesToProtein()");
    }

    bool Optimizer::getMembrane(Tmdet::VOs::Membrane& membrane, int count) {
        DEBUG_LOG("Processing Optimizer::getMembrane()");
        DEBUG_LOG("bestQ: {}",bestQ);
        DEBUG_LOG("bestSlices: {}",bestSlices[10].qValue);
        DEBUG_LOG("bestnormal: [{}, {}, {}]",bestNormal.x, bestNormal.y, bestNormal.z);
        
        /*double q = -1e30;
        for (unsigned long int i = 0; i <bestSlices.size(); i++) {
            if (bestSlices[i].qValue > q) {
                q = bestSlices[i].qValue;
                bestZ = i;
            }
        }
        */
        bestQ = 0;
        slices = bestSlices;
        checkBestSlice();
        unsigned long int bestZ = (bestMinZ+bestMaxZ) / 2;

        DEBUG_LOG("\tLargest qValue: {}",bestQ);
        if (bestQ < TMDET_MINIMUM_QVALUE) {
            DEBUG_LOG(" Processed Optimizer::getMembrane(q:{}): no more membrane",bestQ);
            return false;
        }
        if (count && (bestZ<12 || bestZ>bestSlices.size()-12)) {
            DEBUG_LOG(" Processed Optimizer::getMembrane(z:{}): no more membrane",bestZ);
            return false;
        }
/*
        unsigned long int minz = bestZ; 
        while (minz>0 && bestSlices[minz].qValue > TMDET_MINIMUM_QVALUE && bestZ - minz < TMDET_MEMBRANE_MAX_HALFTHICKNESS) {
            minz--;
        }
        unsigned long int maxz = bestZ; 
        while (maxz<bestSlices.size() && bestSlices[maxz].qValue > TMDET_MINIMUM_QVALUE && maxz - bestZ < TMDET_MEMBRANE_MAX_HALFTHICKNESS) {
            maxz++;
        }
        
        if (maxz - minz < 5) {
            DEBUG_LOG(" Processed Optimizer::getMembrane(): not transmembrane, membrane minQValue is small: {} {} {}",minz,bestZ,maxz);
            return false;
        }
        
        while (minz>2 && bestSlices[minz].qValue > TMDET_MEMBRANE_QVALUE && bestZ - minz < TMDET_MEMBRANE_MAX_HALFTHICKNESS) {
            minz--;
        }
        while (maxz<bestSlices.size() && bestSlices[maxz].qValue > TMDET_MEMBRANE_QVALUE && maxz - bestZ < TMDET_MEMBRANE_MAX_HALFTHICKNESS) {
            maxz++;
        }*/
        membrane.halfThickness = (bestMaxZ-bestMinZ) / 2;

        if (membrane.halfThickness < TMDET_MEMBRANE_MIN_HALFTHICKNESS) {
            DEBUG_LOG(" Processed Optimizer::getMembrane(): not transmembrane, membrane thickness is small: {} {} {}",bestMinZ,bestZ,bestMaxZ);
            return false;
        }

        unsigned long int i = bestZ; 
        while (i>0 && bestSlices[i].qValue > TMDET_MEMBRANE_QVALUE) {
            bestSlices[i].qValue=0;
            i--;
        }
        /*int j=0;
        while (i>0 && j<60) {
            bestSlices[i].qValue=0;
            i--; j++;
        }*/
        i = bestZ+1; 
        while (i<bestSlices.size() && bestSlices[i].qValue > TMDET_MEMBRANE_QVALUE) {
            bestSlices[i].qValue=0;
            i++;
        }
        /*j=0;
        while (i<bestSlices.size() && j<60) {
            bestSlices[i].qValue=0;
            i++; j++;
        }*/
        setMembraneOrigo(membrane,bestMinZ,bestMaxZ);

        DEBUG_LOG(" Processed Optimizer::getMembrane({}-{}-{}) thickness: {} ",
            bestMinZ,bestZ,bestMaxZ,membrane.halfThickness);
        return true;
    }

    void Optimizer::setProteinTMatrix(gemmi::Vec3& origo) const {
        DEBUG_LOG("TMatrix - bestnormal: [{}, {}, {}]",bestNormal.x, bestNormal.y, bestNormal.z);
        DEBUG_LOG("TMatrix - origo: [{}, {}, {}]",origo.x, origo.y, origo.z);
        double x = bestNormal.x;
	    double y = bestNormal.y;
	    double z = bestNormal.z;
	    if (double d = sqrt(y*y+z*z); d>1e-5) {
		    double sa=z/d;
		    double ca=y/d;
		    protein.tmatrix.rot[0][0] = d;
            protein.tmatrix.rot[0][1] = -x*ca;
            protein.tmatrix.rot[0][2] = -x*sa;
            
		    protein.tmatrix.rot[1][0] = 0;
            protein.tmatrix.rot[1][1] = sa;
            protein.tmatrix.rot[1][2] = -ca;
            
		    protein.tmatrix.rot[2][0] = x;
            protein.tmatrix.rot[2][1] = d*ca;
            protein.tmatrix.rot[2][2] = d*sa;
        }
        else {
        	protein.tmatrix.rot[0][0]=0;
            protein.tmatrix.rot[0][1]=0;
            protein.tmatrix.rot[0][2]=1;

            protein.tmatrix.rot[1][0]=0;
            protein.tmatrix.rot[1][1]=1;
            protein.tmatrix.rot[1][2]=0;

            protein.tmatrix.rot[2][0]=1;
            protein.tmatrix.rot[2][1]=0;
            protein.tmatrix.rot[2][2]=0;
        }
        protein.tmatrix.trans = -origo;
    }
}
