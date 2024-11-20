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
                residue.temp.try_emplace("turn",std::any(0.0));
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
                residue.temp.erase("turn");
                for(auto& atom: residue.atoms) {
                    atom.temp.erase("dist");
                }
            }
        );
        DEBUG_LOG(" Processed: Optimizer::end()");
    }

    void Optimizer::setDistances() {
        DEBUG_LOG("Processing: Optimizer::setDistances()");
        protein.eachSelectedResidue(
            [&](Tmdet::ValueObjects::Residue& residue) -> void {
                bool hasCA = false;
                double d = 0.0;
                gemmi::Vec3 ca;
                for(auto& atom: residue.atoms) {
                    d = normal.x * (atom.gemmi.pos.x - massCentre.x);
                    d += normal.y * (atom.gemmi.pos.y - massCentre.y);
                    d += normal.z * (atom.gemmi.pos.z - massCentre.z);
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
                DEBUG_LOG("Dist: {} {}: {} ({})",protein.chains[residue.chainIdx].id,residue.idx,any_cast<double>(residue.temp.at("dist")),residue.atoms.size());
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
                    
                    if (ca1 != nullptr && ca3 != nullptr and q1*q2>0) {
                        if (double d = ca1->pos.dist(ca3->pos); d >1.0) {
                            residue.temp.at("straight") = std::abs(d3-d1)/d;
                        }
                    }    
                    if (q1*q2<0) {
                        residue.temp.at("turn") = 1.0;
                    }
                }
                if (residue.secStrVecIdx != -1) {
                    residue.temp.at("straight") =
                            std::abs(Tmdet::Helpers::Vector::cosAngle(normal,protein.secStrVecs[residue.secStrVecIdx].end - protein.secStrVecs[residue.secStrVecIdx].begin));
                }
                
                DEBUG_LOG("Straight: {} {}: {}",protein.chains[residue.chainIdx].id,residue.idx,any_cast<double>(residue.temp.at("straight")));
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
        min--; max++;
        DEBUG_LOG("Box size: {} {}: {}",min,max,(unsigned int)(max-min));
        slices.clear();
        slices.resize((unsigned int)(max-min));
        DEBUG_LOG(" Processed: Optimizer::setBoundaries()");
    }

    void Optimizer::sumupSlices() {
        DEBUG_LOG("Processing: Optimizer::sumupSlices()");
        protein.eachSelectedResidue(
            [&](Tmdet::ValueObjects::Residue& residue) -> void {
                double resSurf = (any_cast<double>(residue.temp.at("outside")));
                auto sliceIndex = (unsigned int)(any_cast<double>(residue.temp.at("dist")) - min);
                DEBUG_LOG("in sumupSlices: {} {} {} {}",min,any_cast<double>(residue.temp.at("dist")),max,sliceIndex);
                slices[sliceIndex].numCa++;
                slices[sliceIndex].straight += any_cast<double>(residue.temp.at("straight")) ;//* (residue.ss.isBeta()?1.3:1.0);
                if (residue.idx == 0 || residue.idx == protein.chains[residue.chainIdx].length -1) {
                    slices[sliceIndex].chainEnd += 1.0;
                }
                if (protein.chains[residue.chainIdx].type == Tmdet::Types::ChainType::LOW_RES) {
                    slices[sliceIndex].numHyd++;
                    slices[sliceIndex].hydrophobicity += (residue.type.hsc + 12.3 ) / 16.0;
                    slices[sliceIndex].turn += any_cast<double>(residue.temp.at("turn"));
                    slices[sliceIndex].tSurf++;
                }
                else {
                    slices[sliceIndex].turn += resSurf * any_cast<double>(residue.temp.at("turn"));
                    slices[sliceIndex].tSurf += resSurf;
                    for(const auto& atom: residue.atoms) {
                        double atomSurf =  any_cast<double>(atom.temp.at("outside"));
                        slices[sliceIndex].numAtom++;
                        slices[sliceIndex].surf += atomSurf;
                        if (residue.type.atoms.contains(atom.gemmi.name)) {
                            slices[sliceIndex].apol += residue.type.apol * atomSurf;
                            slices[sliceIndex].voronota += 
                                atomSurf * 
                                ( 1-
                                //"voronota frustration: it is small if residue does not like to be on surface"
                                    (residue.type.atoms.at(atom.gemmi.name).mean - Tmdet::Types::voronotaMeanMin) / 
                                        (Tmdet::Types::voronotaMeanMax - Tmdet::Types::voronotaMeanMin));
                        }
                    }
                }
            }
        );
        DEBUG_LOG(" Processed: Optimizer::sumupSlices()");
    }

    void Optimizer::smoothSurf() {
        auto s = slices.size();
        std::vector<double> apol(s,0.0);
        std::vector<double> surf(s,0.0);
        for(unsigned long int i = 0; i<s; i++) {
            int k=0;
            for (int j=-10; j<=10; j++) {
                if (j+(int)i>=0 && j+i<s) {
                    apol[i] += slices[j+i].apol;
                    surf[i] += slices[j+i].surf;
                    k++;
                }
            }
            apol[i] /= k;
            surf[i] /= k;
        }
        for(unsigned long int i = 0; i<s; i++) {
            slices[i].apol = apol[i];
            slices[i].surf = surf[i];
        }
    }

    double Optimizer::divide(double numerator, double denominator) {
        return ( denominator<1e-5?0.0:numerator/denominator);
    }

    double Optimizer::getQValueForSlice(const _slice& s) {
        double straight = divide(s.straight,s.numCa);
        double chainEnd = 1.0 - divide(s.chainEnd, s.numCa);
        double turn = 1.0 - divide(s.turn, s.tSurf);
        double hydrophobicity = divide(s.hydrophobicity, s.numHyd);
        //double apol = divide(s.apol, s.surf);
        double voronota = divide(s.voronota, s.surf);
        return 120.0
            * straight
            * (chainEnd==0.0?1.0:chainEnd)
            * (turn==0.0?1.0:turn)
            * (hydrophobicity==0.0?1.0:hydrophobicity)
            * (s.surf>1?voronota:1.0);
    }

    std::vector<double> Optimizer::getQValueForSlices() {
        auto s = slices.size();
        std::vector<double> qs(s,0);
        for(unsigned long int i = 0; i<s; i++) {
            qs[i] = getQValueForSlice(slices[i]);
        }
        return qs;
    }


    double Optimizer::smoothQValues(std::vector<double> qs) {
        std::string n="***************************************************************************";
        std::string m="mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm";
        std::string b="bbbbbbbbbbbbbbbbbbbbbbbbbbbb";
        double maxQ = -1e30;
        auto s = slices.size();
        for(unsigned long int i = 0; i<s; i++) {
            int k=0;
            double q=0;
            for (int j=-5; j<=5; j++) {
                if (j+(int)i>=0 && j+i<s) {
                    q += qs[j+i];
                    k++;
                }
            }
            q /= k;
            q=(q<1?1:q);
            slices[i].qValue = q;
            if (i>TMDET_MEMBRANE_MIN_HALFTHICKNESS && i< slices.size() - TMDET_MEMBRANE_MIN_HALFTHICKNESS && q > maxQ) {
                maxQ = q;
                bestSliceIndex = i;
            }
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
        return maxQ;
    }

    void Optimizer::setNormal(gemmi::Vec3 _normal) {
        normal = _normal;
        DEBUG_LOG("Optimizer::setNormal: {}:{}:{}",normal.x,normal.y,normal.z);
    }

    void Optimizer::clear() {
        bestQ = 0;
    }

    void Optimizer::testMembraneNormal() {
        if (!run) {
            init();
        }
        setDistances();
        setStraight();
        setBoundaries();
        sumupSlices();
        smoothSurf();
        double q = smoothQValues(getQValueForSlices());
        if (q>bestQ) {
            bestQ = q;
            bestNormal = normal;
            bestSlices = slices;
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
        Tmdet::ValueObjects::Membrane membrane;
        protein.membranes.clear();
        while(getMembrane(membrane) && protein.membranes.size() < 2) {
            protein.membranes.push_back(membrane);
        }
        if (!protein.membranes.empty()) {
            protein.tmp = true;
            setProteinTMatrix(massCentre,bestNormal);
        }
        DEBUG_LOG(" Processed Optimizer::setMembranesToProtein()");
    }

    void Optimizer::searchForMembraneNormal() {
        Tmdet::Engine::Rotator rotator;
        while(rotator.next(normal)) {
            testMembraneNormal();
        }
    }

    bool Optimizer::getMembrane(Tmdet::ValueObjects::Membrane& membrane) {
        DEBUG_LOG("Processing Optimizer::getMembrane()");
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
            DEBUG_LOG(" Processed Optimizer::getMembrane({}): no more membrane",q);
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
        i = bestZ+1; 
        while (i<bestSlices.size() && bestSlices[i].qValue > TMDET_MEMBRANE_QVALUE) {
            bestSlices[i].qValue=0;
            i++;
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

    void Optimizer::setProteinTMatrix(gemmi::Vec3& origo, gemmi::Vec3& normal) const {
        double x = normal.x;
	    double y = normal.y;
	    double z = normal.z;
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
