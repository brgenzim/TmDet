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
#include <ValueObjects/Protein.hpp>
#include <Utils/Surface.hpp>


namespace Tmdet::Engine {

    void Optimizer::init() {
        DEBUG_LOG("Processing: Optimizer::init()");
        run = true;
        massCentre = protein.centre();
        protein.eachSelectedResidue(
            [](Tmdet::ValueObjects::Residue& residue) -> void {
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
            [](Tmdet::ValueObjects::Residue& residue) -> void {
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

    double Optimizer::distance(gemmi::Vec3& vec) {
        return normal.x * (vec.x - massCentre.x)
                + normal.y * (vec.y - massCentre.y)
                + normal.z * (vec.z - massCentre.z);
    }

    void Optimizer::setDistances() {
        DEBUG_LOG("Processing: Optimizer::setDistances()");
        protein.eachSelectedResidue(
            [&](Tmdet::ValueObjects::Residue& residue) -> void {
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
            [&](Tmdet::ValueObjects::Residue& residue) -> void {
                if (residue.idx>2 && residue.idx<protein.chains[residue.chainIdx].length-3) {
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
        min = 1e30;
        max = -1e30;
        protein.eachSelectedResidue(
            [&](Tmdet::ValueObjects::Residue& residue) -> void {
                auto dist = any_cast<double>(residue.temp.at("dist"));
                min = (dist < min ? dist : min);
                max = (dist > max ? dist : max);
            }
        );
        min-=3; max+=3;
        DEBUG_LOG("Box size: {} {}: {}",min,max,(unsigned int)(max-min));
        slices.clear();
        slices.resize((unsigned int)(max-min));
        DEBUG_LOG(" Processed: Optimizer::setBoundaries()");
    }

    void Optimizer::sumupSlices() {
        DEBUG_LOG("Processing: Optimizer::sumupSlices()");
        protein.eachSelectedResidue(
            [&](Tmdet::ValueObjects::Residue& residue) -> void {
                auto sliceIndex = (unsigned int)(any_cast<double>(residue.temp.at("dist")) - min);
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
                        slices[sliceIndex].surf += atom.outSurface * (residue.ss.isBeta()?1.4:0.9);
                        if (residue.type.atoms.contains(atom.gemmi.name)) {
                            slices[sliceIndex].apol += residue.type.apol * atom.outSurface 
                                * (residue.ss.isBeta()?1.2:0.9);
                            //slices[sliceIndex].apol += atom.outSurface * (residue.type.hsc + 12.3 ) / 16.0;
                            /*slices[sliceIndex].apol += 
                                atom.outSurface * 
                                ( 1.2-
                                //"voronota frustration: it is small if residue does not like to be on surface"
                                    (residue.type.atoms.at(atom.gemmi.name).mean - Tmdet::Types::voronotaMeanMin) / 
                                        (Tmdet::Types::voronotaMeanMax - Tmdet::Types::voronotaMeanMin));*/
                        }
                    }
                }
            }
        );
        for (auto& vector: protein.secStrVecs) {
            if (protein.chains[vector.chainIdx].selected) {
                double cosAngle = std::abs(Tmdet::Helpers::Vector::cosAngle(
                                    normal,vector.end - vector.begin));
                int d1 = distance(vector.begin);
                int d2 = distance(vector.end);
                int dbeg = (d1<d2?d1:d2) - min; 
                dbeg=(dbeg<0?0:dbeg);
                int dend = (d1>d2?d1:d2) - min; 
                if (dend>=(int)slices.size()) {
                    dend=(int)slices.size()-1;
                }
                for(int i=dbeg; i<= dend; i++) {
                    slices[i].straight+=cosAngle * (vector.type.isBeta()?1.5:1.0);
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
        std::string b="bbbbbbbbbbbbbbbbbbbbbbbbbbbb";
        auto s = slices.size();
        for(unsigned long int i = 0; i<s; i++) {
            int k=0;
            double smoothedStraight = 0.0;
            double smoothedApol = 0.0;
            double smoothedSsEnd = 0.0;
            double smoothedSurf = 0.0;
            double smoothedCa = 0.0;
            for (int j=-2; j<=2; j++) {
                if (j+(int)i>=0 && j+i<s) {
                    smoothedStraight += slices[j+i].straight;
                    smoothedApol += slices[j+i].apol;
                    smoothedSsEnd += slices[j+i].ssEnd;
                    smoothedSurf += slices[i+j].surf;
                    smoothedCa += slices[i+j].numCa;
                    k++;
                }
            }
            smoothedStraight /= k;
            smoothedApol /= k;
            smoothedSsEnd /= k;
            smoothedSurf /= k;
            smoothedCa /= k;
            smoothedStraight = divide(smoothedStraight,smoothedCa);
            smoothedApol = divide(smoothedApol, smoothedSurf);
            smoothedSsEnd = divide(smoothedSsEnd, smoothedCa);
            
            slices[i].rawQ = 130.0 * smoothedStraight *(1.0 - smoothedSsEnd) * smoothedApol;

             DEBUG_LOG("{} straight:{:5.2f} apol:{:5.2f} surf:{:5.2f} q:{:5.2f}",
                i,smoothedStraight,smoothedApol,slices[i].surf,slices[i].qValue);
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
                DEBUG_LOG("{}:{} {}{}{}",
                    i, q,
                    n.substr(0,TMDET_MEMBRANE_QVALUE),
                    m.substr(0,(int)(TMDET_MINIMUM_QVALUE-TMDET_MEMBRANE_QVALUE)),
                    b.substr(0,(int)(q-TMDET_MINIMUM_QVALUE+1)));
            }
            else if (q>TMDET_MEMBRANE_QVALUE) {
                DEBUG_LOG("{}:{} {}{}",
                    i, q,
                    n.substr(0,TMDET_MEMBRANE_QVALUE),
                    m.substr(0,(int)(q-TMDET_MEMBRANE_QVALUE+1)));
            }
            else {
                DEBUG_LOG("{}:{} {}",i,q,n.substr(0,(int)(q)));
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
                        bestNormal = normal;
                        bestSlices = slices;
                        DEBUG_LOG("BestWidth: {} {} {} {} {}",minz,i,maxz,q,bestQ);
                    }
                }
            }
        }
    }

    double Optimizer::getWidth(const int z, int& minz, int& maxz) {
        minz = z;
        while(minz>0 && slices[minz].qValue > TMDET_MEMBRANE_QVALUE) {
            minz--;
        }
        maxz = z+1;
        while(maxz<(int)slices.size() && slices[maxz].qValue > TMDET_MEMBRANE_QVALUE) {
            maxz++;
        }
        return slices[z].qValue;
    }

    void Optimizer::setNormal(gemmi::Vec3 _normal) {
        normal = _normal;
    }

    void Optimizer::testMembraneNormal() {
        if (!run) {
            init();
        }
        setDistances();
        setStraight();
        setBoundaries();
        sumupSlices();
        smoothQValues();
        checkBestSlice();
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
        testMembraneNormal();
        Tmdet::ValueObjects::Membrane membrane;
        protein.membranes.clear();
        int i=0;
        while(getMembrane(membrane,i) && i < 2) {
            protein.membranes.push_back(membrane);
            i++;
        }
        if (!protein.membranes.empty()) {
            protein.tmp = true;
            setProteinTMatrix(massCentre);
        }
        DEBUG_LOG(" Processed Optimizer::setMembranesToProtein()");
    }

    void Optimizer::searchForMembraneNormal() {
        Tmdet::Engine::Rotator rotator;
        while(rotator.next(normal)) {
            testMembraneNormal();
        }
    }

    bool Optimizer::getMembrane(Tmdet::ValueObjects::Membrane& membrane, int count) {
        DEBUG_LOG("Processing Optimizer::getMembrane()");
        DEBUG_LOG("bestQ: {}",bestQ);
        DEBUG_LOG("bestSlices: {}",bestSlices[10].qValue);
        DEBUG_LOG("bestnormal: [{}, {}, {}]",bestNormal.x, bestNormal.y, bestNormal.z);
        unsigned long int bestZ = -1;
        double q = -1e30;
        for (unsigned long int i = 0; i <bestSlices.size(); i++) {
            if (bestSlices[i].qValue > q) {
                q = bestSlices[i].qValue;
                bestZ = i;
            }
        }
        DEBUG_LOG("\tLargest qValue: {}",q);
        if (q < TMDET_MINIMUM_QVALUE) {
            DEBUG_LOG(" Processed Optimizer::getMembrane(q:{}): no more membrane",q);
            return false;
        }
        if (count && (bestZ<12 || bestZ>bestSlices.size()-12)) {
            DEBUG_LOG(" Processed Optimizer::getMembrane(z:{}): no more membrane",bestZ);
            return false;
        }

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
        }
        membrane.halfThickness = (maxz-minz) / 2;

        if (membrane.halfThickness < TMDET_MEMBRANE_MIN_HALFTHICKNESS) {
            DEBUG_LOG(" Processed Optimizer::getMembrane(): not transmembrane, membrane thickness is small: {} {} {}",minz,bestZ,maxz);
            return false;
        }

        unsigned long int i = bestZ; 
        while (i>0 && bestSlices[i].qValue > TMDET_MEMBRANE_QVALUE) {
            bestSlices[i].qValue=0;
            i--;
        }
        int j=0;
        while (i>0 && j<60) {
            bestSlices[i].qValue=0;
            i--; j++;
        }
        i = bestZ+1; 
        while (i<bestSlices.size() && bestSlices[i].qValue > TMDET_MEMBRANE_QVALUE) {
            bestSlices[i].qValue=0;
            i++;
        }
        j=0;
        while (i<bestSlices.size() && j<60) {
            bestSlices[i].qValue=0;
            i++; j++;
        }
        double o = (minz+maxz) / 2 + min;
        if (protein.membranes.empty()) {
            massCentre += o * bestNormal;
            membrane.origo = 0;
            lastO = o;
        }
        else {
            membrane.origo = o - lastO;
        }

        DEBUG_LOG(" Processed Optimizer::getMembrane({}-{}-{}) thickness: {} ",
            minz,bestZ,maxz,membrane.halfThickness);
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
