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
        end();
    }

    void Optimizer::init() {
        minHalfThickness = args.getValueAsFloat("minht");
        maxHalfThickness = args.getValueAsFloat("maxht");
        maxCurvedHalfThickness = args.getValueAsFloat("maxcht");
        lowerQ = args.getValueAsFloat("lq");
        higherQ = args.getValueAsFloat("hq");
        higherQ2 = args.getValueAsFloat("hq2");
        massCentre = protein.centre();
        ifhAngleLimit = args.getValueAsFloat("ian");
        boostAngle = args.getValueAsFloat("ba");
        boostBetaAngle = args.getValueAsFloat("bba");
        boostPolarity = args.getValueAsFloat("bp");
        int numRes = 0;
        protein.eachSelectedResidue(
            [&](Tmdet::VOs::Residue& residue) -> void {
                residue.temp.try_emplace("dist",std::any(0.0));
                for(auto& atom: residue.atoms) {
                    atom.temp.try_emplace("dist",std::any(0.0));
                }
                numRes++;
            }
        );
        addStraigth = args.getValueAsFloat("spen");
        addStraigth +=  (numRes>500?0.5:0) ;
        addStraigth +=  (numRes>1000?0.5:0) ;
        protein.eachSelectedChain(
            [&](Tmdet::VOs::Chain& chain) -> void {
                for (int i=0; i<chain.length; i++) {
                    if (chain.residues[i].selected) {
                        int k=0;
                        int w = (chain.residues[i].ss.isBeta()?0:8);
                        double apol = 0;
                        for (int j=-w; j<=w; j++) {
                            if (i+j>=0&&i+j<chain.length&&chain.residues[i+j].selected) {
                                apol += chain.residues[i+j].type.apol;
                                k++;
                            }
                        }
                        chain.residues[i].apol = apol/k;
                    }
                }
            }
        );
    }

    void Optimizer::end() {
        protein.eachSelectedResidue(
            [](Tmdet::VOs::Residue& residue) -> void {
                residue.temp.erase("dist");
                for(auto& atom: residue.atoms) {
                    atom.temp.erase("dist");
                }
            }
        );
    }

    void Optimizer::clear() {
        bestQ = 0;
    }

    void Optimizer::setDistances() {
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
    }

    void Optimizer::setBoundaries() {
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
        slices.clear();
        slices.resize((unsigned int)(maxZ-minZ));
    }



    void Optimizer::sumupSlices() {
        protein.eachSelectedResidue(
            [&](Tmdet::VOs::Residue& residue) -> void {
                auto sliceIndex = (unsigned int)(any_cast<double>(residue.temp.at("dist")) - minZ);
                if (protein.chains[residue.chainIdx].type == Tmdet::Types::ChainType::LOW_RES) {
                    slices[sliceIndex].surf+=10;
                    slices[sliceIndex].apol += 10 * residue.apol;
                }
                else {
                    for(const auto& atom: residue.atoms) {
                        slices[sliceIndex].surf += atom.outSurface;
                        if (residue.type.atoms.contains(atom.gemmi.name)) {
                            slices[sliceIndex].apol += atom.outSurface * (residue.ss.isBeta()?1.1:1.0) * ( 1-
                                    (residue.type.atoms.at(atom.gemmi.name).mean - Tmdet::Types::voronotaMeanMin) / 
                                        (Tmdet::Types::voronotaMeanMax - Tmdet::Types::voronotaMeanMin));
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
                if (vector.type.isAlpha()) {
                    cosAngle = (cosAngle<0.25?-3:cosAngle);
                    if (protein.hasIdenticalChains) {
                        cosAngle = (cosAngle>boostAngle?1:cosAngle);
                    }
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
                if (protein.chains[vector.chainIdx].residues[vector.begResIdx].selected
                    && protein.chains[vector.chainIdx].residues[vector.endResIdx].selected)  {
                    for (int j=vector.begResIdx; j<=vector.endResIdx; j++) {
                        if (protein.chains[vector.chainIdx].residues[j].selected) {
                            int i = distance(any_cast<gemmi::Vec3&>(protein.chains[vector.chainIdx].residues[j].temp["ca"])) - minZ;
                            if (i>=0 && i<(int)slices.size()) {
                                
                                if (vector.type.isBeta()) {
                                    slices[i].beta+= (cosAngle>boostBetaAngle?1:cosAngle);
                                }
                                else {
                                    slices[i].alpha+= cosAngle;
                                    if (dend-dbeg<5 
                                        && angle < ifhAngleLimit
                                        && vector.endResIdx - vector.begResIdx + 1 >= ifhMinLength) {
                                        slices[i].ifh++;
                                    }
                                }
                                slices[i].numCa++;
                            }
                        }
                    }
                }
                for(int i=dbeg; i<= dend; i++) {
                    if (vector.type.isBeta()) {
                        slices[i].beta += (cosAngle>boostBetaAngle?1:cosAngle);
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
                    smoothedStraight += (slices[j+i].beta>7?slices[j+i].beta:slices[j+i].alpha);
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
            smoothedStraight = 1.2*divide(smoothedStraight,(smoothedCa + addStraigth));
            smoothedStraight = (smoothedStraight<0?0:smoothedStraight);
            smoothedApol = (smoothedSurf>5?divide(smoothedApol, smoothedSurf):0); 
            smoothedApol = (smoothedApol>boostPolarity?1:smoothedApol);
            smoothedSsEnd = divide(smoothedSsEnd, smoothedCa);
            smoothedIfh = 0;
            slices[i].smoothedIfh = smoothedIfh;

            slices[i].rawQ = 100.0 * smoothedStraight * (1.0 - smoothedSsEnd) * (1.0 - smoothedIfh) * smoothedApol;
        }
        for(unsigned long int i = 0; i<s; i++) {
            int k=0;
            double q = 0.0;
            for (int j=-8; j<=8; j++) {
                if (j+(int)i>=0 && j+i<s) {
                    int w = (9 - std::abs(j));
                    q += w * slices[j+i].rawQ;
                    k += w;
                }
            }
            q/=k;
            slices[i].qValue = q;
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
                    if (maxz-minz > 2 * minHalfThickness 
                        && ((q>bestQ && minz>5 && maxz<(int)s-5)
                            || q>bestQ + 6)) {
                        bestQ = q;
                        bestMinZ = minz; bestMinZ -= 0.5;
                        bestMaxZ = maxz; bestMaxZ -= 0.5;
                        bestNormal = normal;
                        bestSlices = slices;
                        setBestOrigo(minz,maxz);
                    }
                }
            }
        }
    }

    double Optimizer::getWidth(const int z, int& minz, int& maxz) {
        double maxHT = (type=="Plane"?maxHalfThickness:maxCurvedHalfThickness);
        double ifhLimit = (type=="Plane"?0.1:1.1);
        minz = z;
        maxz = z+1;
        double q = 0.0;
        bool ok = true;
        while(ok) {
            ok = false;
            if (minz > 0
                && slices[minz].qValue > higherQ 
                && slices[minz].smoothedIfh < ifhLimit
                && maxz-minz < 2*maxHT) {
                    q += slices[minz].qValue;
                    minz--;
                    ok = true;
            }
            if (maxz < (int)slices.size()-1
                && slices[maxz].qValue > higherQ
                && slices[maxz].smoothedIfh < ifhLimit
                && maxz-minz < 2*maxHT) {
                    q += slices[maxz].qValue;
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
                    && slices[minz].smoothedIfh < ifhLimit
                    && maxz-minz < 2*maxHT) {
                        q += slices[minz].qValue;
                        minz--;
                        ok = true;
                }
                if (maxz < (int)slices.size()-1
                    && slices[maxz].qValue > lowerQ
                    && slices[maxz].smoothedIfh < ifhLimit
                    && maxz-minz < 2*maxHT) {
                        q += slices[maxz].qValue;
                        maxz++;
                        ok = true;
                }
            }
            minz++;
            maxz--;
        }
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
        return bestQ > higherQ;
    }

    void Optimizer::setMembranesToProtein() {
        protein.qValue = bestQ;
        if (!isTransmembrane()) {
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
    }

    bool Optimizer::getMembrane(Tmdet::VOs::Membrane& membrane, int count) {
        
        bestQ = 0;
        slices = bestSlices;
        checkBestSlice();
        unsigned long int bestZ = (bestMinZ+bestMaxZ) / 2;

        if (bestQ < (count?higherQ2:higherQ)) {
            return false;
        }
        if (count && (bestZ<12 || bestZ>bestSlices.size()-12)) {
            return false;
        }
        membrane.halfThickness = (bestMaxZ-bestMinZ) / 2;

        if (membrane.halfThickness < minHalfThickness) {
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

        return true;
    }

    void Optimizer::setProteinTMatrix(gemmi::Vec3& origo) const {
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
