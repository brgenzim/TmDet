#include <iostream>
#include <array>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <any>
#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <Config.hpp>
#include <Types/Residue.hpp>
#include <Types/Chain.hpp>
#include <Engine/Optimizer.hpp>
#include <Engine/Rotator.hpp>
#include <Utils/Surface.hpp>
#include <ValueObjects/Protein.hpp>

using namespace std;

namespace Tmdet::Engine {

    void Optimizer::initChain(Tmdet::ValueObjects::Chain& chain) {
        for(auto& residue: chain.residues) {
            residue.temp.try_emplace("turn",any(0));
            residue.temp.try_emplace("straight",any(0));
            residue.temp.try_emplace("win",any(0));
            residue.temp.try_emplace("dist",any(0.0));
            for(auto& atom: residue.atoms) {
                atom.temp.try_emplace("dist",any(0.0));
            }
        }
    }

    void Optimizer::init() {
        DEBUG_LOG("Processing: Optimizer::init()");
        run = true;
        massCentre = protein.centre();
        protein.eachSelectedChain(this,&Optimizer::initChain);
        /*EACH_SELECTED_CHAIN(protein) {
            initChain(chain);
        }*/
        DEBUG_LOG(" Processed: Optimizer::init()");
    }

    void Optimizer::end() {
        DEBUG_LOG("Processing: Optimizer::end()");
        EACH_SELECTED_CHAIN(protein) {
            for(auto& residue: chain.residues) {
                residue.temp.erase("turn");
                residue.temp.erase("straight");
                residue.temp.erase("win");
                residue.temp.erase("dist");
                for(auto& atom: residue.atoms) {
                    atom.temp.erase("dist");
                }
            }
        }
        DEBUG_LOG(" Processed: Optimizer::end()");
    }

    void Optimizer::setDistances() {
        DEBUG_LOG("Processing: Optimizer::setDistances()");
        EACH_SELECTED_CHAIN(protein) {
            for(auto& residue: chain.residues) {
                setAtomDistances(residue);
            }
        }
        DEBUG_LOG(" Processed: Optimizer::setDistances()");
    }

    void Optimizer::setAtomDistances(Tmdet::ValueObjects::Residue& residue) const {
        bool hasCA = false;
        double d;
        for(auto& atom: residue.atoms) {
            d = normal.x * (atom.gemmi.pos.x - massCentre.x);
            d += normal.y * (atom.gemmi.pos.y - massCentre.y);
            d += normal.z * (atom.gemmi.pos.z - massCentre.z);
            atom.temp.at("dist") = any(d);
            if (atom.gemmi.name == "CA") {
                residue.temp.at("dist") = any(d);
                hasCA = true;
            }
        }
        if (!hasCA) {
            residue.temp.at("dist") = any(d);
        }
    }

    void Optimizer::setBoundaries() {
        DEBUG_LOG("Processing: Optimizer::setBoundaries()");
        min = 1e30;
        max = -1e30;
        EACH_SELECTED_CHAIN(protein) {
            for(auto& residue: chain.residues) {
                auto dist = any_cast<double>(residue.temp.at("dist"));
                min = (dist < min ? dist : min);
                max = (dist > max ? dist : max);
            }
        }
        min--; max++;
        DEBUG_LOG("Box size: {} {}: {}",min,max,(unsigned int)(max-min));
        slices.clear();
        slices.resize((unsigned int)(max-min));
        DEBUG_LOG(" Processed: Optimizer::setBoundaries()");
    }

    void Optimizer::sumupSlices() {
        DEBUG_LOG("Processing: Optimizer::sumupSlices()");
        EACH_SELECTED_CHAIN(protein) {
            for(auto& residue: chain.residues) {
                residueToSlice(residue);
            }
        }
        DEBUG_LOG(" Processed: Optimizer::sumupSlices()");
    }

    void Optimizer::residueToSlice(Tmdet::ValueObjects::Residue& residue) {
        auto sliceIndex = (unsigned int)(any_cast<double>(residue.temp.at("dist")) - min);
        slices[sliceIndex].numCA++;
        if (residue.ss.isStrictAlpha()) {
            slices[sliceIndex].numAlpha++;
        }
        else if (residue.ss.isBeta()) {
            slices[sliceIndex].numBeta++;
        }
        else if (residue.ss.isStrictTurn()) {
            slices[sliceIndex].numTurn++;
        }
        if (protein.chains[residue.chainIdx].type == Tmdet::Types::ChainType::LOW_RES) {
            slices[sliceIndex].surf += (any_cast<double>(residue.temp.at("outside")));
            slices[sliceIndex].voronota += (any_cast<double>(residue.temp.at("outside"))) *
                                            (residue.type.hsc + 12.3 ) / 16;
        }
        else {
            for(const auto& atom: residue.atoms) {
                slices[sliceIndex].surf +=  any_cast<double>(atom.temp.at("outside"));
                if (residue.type.atoms.contains(atom.gemmi.name)) {
                    slices[sliceIndex].voronota += 
                        (any_cast<double>(atom.temp.at("outside")) * 
                        residue.type.apol *
                        ( 1-
                        //"voronota frustration: it is small if residue does not like to be on surface"
                            (residue.type.atoms.at(atom.gemmi.name).mean - Tmdet::Types::voronotaMeanMin) / 
                                (Tmdet::Types::voronotaMeanMax - Tmdet::Types::voronotaMeanMin)));
                }
            }
        }
    }

    double Optimizer::getQValueForSlice(const _slice& s) const {
        double q = 0.0;
        if (s.numCA > 0 && s.surf >0 ) {
            q += 35 * ((double)s.numAlpha + (double)s.numBeta * 1.3) / s.numCA;
            q -= 30 * (double)s.numTurn / s.numCA;
            q += 65 * s.voronota / s.surf;
        }
        return q;
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
            for (int j=-7; j<=7; j++) {
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
    }

    void Optimizer::clear() {
        bestQ = 0;
    }

    void Optimizer::testMembraneNormal() {
        if (!run) {
            init();
        }
        setDistances();
        setBoundaries();
        sumupSlices();
        double q = smoothQValues(getQValueForSlices());
        if (q>bestQ) {
            bestQ = q;
            bestNormal = normal;
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
        setNormal(bestNormal);
        testMembraneNormal();
        Tmdet::ValueObjects::Membrane membrane;
        protein.membranes.clear();
        while(getMembrane(membrane) && protein.membranes.size() < 2) {
            protein.membranes.push_back(membrane);
        }
        if (!protein.membranes.empty()) {
            protein.tmp = true;
            setProteinTMatrix(massCentre,normal);
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
        for (unsigned long int i = TMDET_MEMBRANE_MIN_HALFTHICKNESS; i <slices.size()-TMDET_MEMBRANE_MIN_HALFTHICKNESS; i++) {
            if (slices[i].qValue > q) {
                q = slices[i].qValue;
                bestZ = i;
            }
        }
        DEBUG_LOG("\tLargest qValue: {}",q);
        if (q < TMDET_MINIMUM_QVALUE) {
            DEBUG_LOG(" Processed Optimizer::getMembrane({}): no more membrane",q);
            return false;
        }

        unsigned long int minz = bestZ; 
        while (minz>2 && slices[minz].qValue > TMDET_MINIMUM_QVALUE && bestZ - minz < TMDET_MEMBRANE_MAX_HALFTHICKNESS) {
            minz--;
        }
        unsigned long int maxz = bestZ; 
        while (maxz<slices.size()-2 && slices[maxz].qValue > TMDET_MINIMUM_QVALUE && maxz - bestZ < TMDET_MEMBRANE_MAX_HALFTHICKNESS) {
            maxz++;
        }
        
        if (maxz - minz < 5) {
            DEBUG_LOG(" Processed Optimizer::getMembrane(): not transmembrane, membrane minQValue is small: {} {} {}",minz,bestZ,maxz);
            return false;
        }
        
        while (minz>2 && slices[minz].qValue > TMDET_MEMBRANE_QVALUE && bestZ - minz < TMDET_MEMBRANE_MAX_HALFTHICKNESS) {
            minz--;
        }
        while (maxz<slices.size()-2 && slices[maxz].qValue > TMDET_MEMBRANE_QVALUE && maxz - bestZ < TMDET_MEMBRANE_MAX_HALFTHICKNESS) {
            maxz++;
        }
        membrane.halfThickness = (maxz-minz) / 2;

        if (membrane.halfThickness < TMDET_MEMBRANE_MIN_HALFTHICKNESS) {
            DEBUG_LOG(" Processed Optimizer::getMembrane(): not transmembrane, membrane thickness is small: {} {} {}",minz,bestZ,maxz);
            return false;
        }

        unsigned long int i = bestZ; 
        while (i>2 && slices[i].qValue > TMDET_MEMBRANE_QVALUE) {
            slices[i].qValue=0;
            i--;
        }
        i = bestZ+1; 
        while (i<slices.size()-2 && slices[i].qValue > TMDET_MEMBRANE_QVALUE) {
            slices[i].qValue=0;
            i++;
        }

        double o = (minz+maxz) / 2 + min;
        if (protein.membranes.empty()) {
            massCentre += o * normal;
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
