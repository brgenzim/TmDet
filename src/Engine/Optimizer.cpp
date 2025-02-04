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
#include <System/Arguments.hpp>
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
        minHalfThickness = args.getValueAsFloat("minht");
        maxHalfThickness = args.getValueAsFloat("maxht");
        maxCurvedHalfThickness = args.getValueAsFloat("maxcht");
        lowerQ = args.getValueAsFloat("lq");
        higherQ = args.getValueAsFloat("hq");
        higherQ2 = args.getValueAsFloat("hq2");
        massCentre = protein.centre();
        ifhAngleLimit = args.getValueAsFloat("ian");
        ifhResLimit = args.getValueAsInt("irl");
        boostAngle = args.getValueAsFloat("ba");
        boostPolarity = args.getValueAsFloat("bp");

        protein.eachSelectedResidue(
            [](Tmdet::VOs::Residue& residue) -> void {
                residue.temp.try_emplace("dist",std::any(0.0));
                for(auto& atom: residue.atoms) {
                    atom.temp.try_emplace("dist",std::any(0.0));
                }
            }
        );
        protein.eachSelectedChain(
            [&](Tmdet::VOs::Chain& chain) -> void {
                for (int i=0; i<chain.length; i++) {
                    if (chain.residues[i].selected) {
                        int k=0;
                        int w = (chain.residues[i].ss.isBeta()?0:8);
                        double apol = 0;
                        for (int j=-w; j<=w; j++) {
                            if (i+j>=0&&i+j<chain.length&&chain.residues[i+j].selected) {
                                //apol += (chain.residues[i+j].type.hsc + 12.3 ) / 16.0;
                                apol += chain.residues[i+j].type.apol;
                                k++;
                            }
                        }
                        chain.residues[i].apol = apol/k;
                        double q = (chain.residues[i].type.hsc + 12.3 ) / 16.0;
                        DEBUG_LOG("APOL: {}:{}({}) w:{} {},{} = {}",
                            chain.id,chain.residues[i].authId,
                            chain.residues[i].type.name,
                            w,chain.residues[i].type.hsc,q,
                            chain.residues[i].apol);
                    }
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
                if (protein.chains[residue.chainIdx].type == Tmdet::Types::ChainType::LOW_RES) {
                    slices[sliceIndex].surf+=5;
                    slices[sliceIndex].apol += 5*(residue.type.hsc + 12.3 ) / 16.0;
                }
                else {
                    for(const auto& atom: residue.atoms) {
                        slices[sliceIndex].surf += atom.outSurface /* (residue.ss.isBeta()?1.1:0.9)*/;
                        if (residue.type.atoms.contains(atom.gemmi.name)) {
                            slices[sliceIndex].apol += residue.apol * atom.outSurface 
                                /* (residue.ss.isBeta()?1.1:0.9)*/;
                        }
                    }
                }
            }
        );
        for (auto& vector: protein.secStrVecs) {
            if (protein.chains[vector.chainIdx].selected
                && protein.chains[vector.chainIdx].residues[vector.begResIdx].selected
                && protein.chains[vector.chainIdx].residues[vector.endResIdx].selected) {
                double cosAngle = getAngle(vector);
                double angle = std::abs(90.0 - acos(cosAngle) * 180.0 / M_PI);
                cosAngle = (cosAngle<0.2?-3:cosAngle);
                if (protein.hasIdenticalChains) {
                    cosAngle = (cosAngle>boostAngle?1:cosAngle);
                }
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
                DEBUG_LOG("sumupSliceVectors: begz:{} endz:{} beg:{} end:{} size:{} cosAngle:{} angle:{}",
                    dbeg,dend,
                    protein.chains[vector.chainIdx].residues[vector.begResIdx].authId,
                    protein.chains[vector.chainIdx].residues[vector.endResIdx].authId,
                    slices.size(),cosAngle,angle);
                //for(int i=dbeg; i<= dend; i++) {
                if (vector.type.isAlpha()
                    && protein.chains[vector.chainIdx].residues[vector.begResIdx].selected
                    && protein.chains[vector.chainIdx].residues[vector.endResIdx].selected
                    && dend-dbeg<5 
                    && angle < ifhAngleLimit)  {
                    for (int j=vector.begResIdx; j<=vector.endResIdx; j++) {
                        if (protein.chains[vector.chainIdx].residues[j].selected) {
                            int i = distance(any_cast<gemmi::Vec3&>(protein.chains[vector.chainIdx].residues[j].temp["ca"])) - minZ;
                            if (i>=0 && i<(int)slices.size()) {
                                slices[i].ifh++;
                            }
                        }
                    }
                }
                for(int i=dbeg; i<= dend; i++) {
                    if (vector.type.isBeta()) {
                        slices[i].beta += (cosAngle>0.55?1:cosAngle);
                    }
                    else {
                        slices[i].alpha += cosAngle;
                    }
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
                    smoothedStraight += (slices[j+i].alpha>slices[j+i].beta?slices[i+j].alpha:(slices[j+i].beta>4?slices[j+i].beta:0));
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
            smoothedApol = (smoothedSurf>10?divide(smoothedApol, smoothedSurf):0); smoothedApol = (smoothedApol>boostPolarity?1:smoothedApol);
            smoothedSsEnd = divide(smoothedSsEnd, smoothedCa);
            //slices[i].smoothedIfh = divide(smoothedIfh, smoothedCa);
            smoothedIfh /= ifhResLimit; smoothedIfh = (smoothedIfh>1?1:smoothedIfh);

            
            slices[i].rawQ = 100.0 * smoothedStraight * (1.0 - smoothedSsEnd) * (1.0 - smoothedIfh) * smoothedApol;

             DEBUG_LOG("{} straight:{:5.2f} apol:{:5.2f} surf:{:5.2f} end:{:5.2f} ifh:{} q:{:5.2f}",
                i,smoothedStraight,smoothedApol,slices[i].surf,smoothedSsEnd,slices[i].ifh,slices[i].rawQ);
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
            
            if (q>higherQ) {
                DEBUG_LOG("{}:{:.2f} {}{}{}",
                    i, q,
                    n.substr(0,lowerQ),
                    m.substr(0,(int)(higherQ-lowerQ)),
                    b.substr(0,(int)(q-higherQ+1)));
            }
            else if (q>lowerQ) {
                DEBUG_LOG("{}:{:.2f} {}{}",
                    i, q,
                    n.substr(0,lowerQ),
                    m.substr(0,(int)(q-lowerQ+1)));
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
        
        if ((int)s > 2 * minHalfThickness) {
            for(int i = minHalfThickness/2; i < (int)s - minHalfThickness/2; i++) {
                int minz;
                int maxz;
                if (slices[i].qValue>higherQ) {
                    double q = getWidth(i,minz,maxz);
                    if (maxz-minz > 2 * minHalfThickness && q>higherQ) {
                        if (q>bestQ) {
                            bestQ = q;
                            bestMinZ = minz; bestMinZ -= 0.5;
                            bestMaxZ = maxz; bestMaxZ -= 0.5;
                            bestNormal = normal;
                            bestSlices = slices;
                            setBestOrigo(minz,maxz);
                            DEBUG_LOG("BestWidth: {} {} {} {} {}",minz,i,maxz,q,bestQ);
                        }
                    }
                }
            }
        }
    }

    double Optimizer::getWidth(const int z, int& minz, int& maxz) {
        double maxHT = (type=="Plane"?maxHalfThickness:maxCurvedHalfThickness);
        double ifhLimit = (type=="Plane"?ifhResLimit:10000);
        minz = z;
        maxz = z+1;
        bool ok = true;
        while(ok) {
            ok = false;
            if (minz > 0
                && slices[minz].qValue > higherQ 
                && slices[minz].smoothedIfh < ifhLimit
                && maxz-minz < 2*maxHT) {
                    minz--;
                    ok = true;
            }
            if (maxz < (int)slices.size()-1
                && slices[maxz].qValue > higherQ
                && slices[maxz].smoothedIfh < ifhLimit
                && maxz-minz < 2*maxHT) {
                    maxz++;
                    ok = true;
            }
        }
        
        if (maxz-minz>7) {
            ok=true;
            while(ok) {
                ok = false;
                if (minz > 0
                    && slices[minz].qValue > lowerQ 
                    && slices[minz].smoothedIfh < ifhResLimit
                    && maxz-minz < 2*maxHT) {
                        minz--;
                        ok = true;
                }
                if (maxz < (int)slices.size()-1
                    && slices[maxz].qValue > lowerQ
                    && slices[maxz].smoothedIfh < ifhResLimit
                    && maxz-minz < 2*maxHT) {
                        maxz++;
                        ok = true;
                }
            }
            minz++;
            maxz--;
        }
        DEBUG_LOG("getWidth(): {} {} {}",minz,maxz,slices[z].qValue);
        return slices[z].qValue;
    }

    void Optimizer::testMembraneNormalOne() {
        setDistances();
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
        if (type == "Curved") {
            rotator.end180();
        }
        while(rotator.next(normal)) {
            testMembraneNormal();
        }
    }

    bool Optimizer::isTransmembrane() const {
        DEBUG_LOG("Processing Optimizer::isTransmembrane(): {}",bestQ);
        return bestQ > higherQ;
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
        testMembraneNormalFinal();
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
        DEBUG_LOG(" Processed Optimizer::setMembranesToProtein({}:{})",protein.tmp,protein.membranes.size());
    }

    bool Optimizer::getMembrane(Tmdet::VOs::Membrane& membrane, int count) {
        DEBUG_LOG("Processing Optimizer::getMembrane({})",count+1);
        DEBUG_LOG("bestQ: {}",bestQ);
        DEBUG_LOG("bestSlices: {}",bestSlices[10].qValue);
        DEBUG_LOG("bestnormal: [{}, {}, {}]",bestNormal.x, bestNormal.y, bestNormal.z);
        
        bestQ = 0;
        slices = bestSlices;
        checkBestSlice();
        unsigned long int bestZ = (bestMinZ+bestMaxZ) / 2;

        DEBUG_LOG("\tLargest qValue: {}",bestQ);
        if (bestQ < (count?higherQ2:higherQ)) {
            DEBUG_LOG(" Processed Optimizer::getMembrane(q:{}): no more membrane",bestQ);
            return false;
        }
        if (count && (bestZ<12 || bestZ>bestSlices.size()-12)) {
            DEBUG_LOG(" Processed Optimizer::getMembrane(z:{}): no more membrane",bestZ);
            return false;
        }
        membrane.halfThickness = (bestMaxZ-bestMinZ) / 2;

        if (membrane.halfThickness < minHalfThickness) {
            DEBUG_LOG(" Processed Optimizer::getMembrane(): not transmembrane, membrane thickness is small: {} {} {}",bestMinZ,bestZ,bestMaxZ);
            return false;
        }

        unsigned long int i = bestZ; 
        while (i>0 && bestSlices[i].qValue > lowerQ) {
            bestSlices[i].qValue=0;
            i--;
        }
        i = bestZ+1; 
        while (i<bestSlices.size() && bestSlices[i].qValue > lowerQ) {
            bestSlices[i].qValue=0;
            i++;
        }
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
