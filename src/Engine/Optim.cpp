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
#include <Engine/Optim.hpp>
#include <Engine/Rotator.hpp>
#include <Utils/Surface.hpp>

using namespace std;

namespace Tmdet::Engine {

    void Optim::init() {
        logger.debug("Processing: Optim::init()");
        run = true;
        massCentre = protein.centre();
        for(auto& chain: protein.chains) {
            for(auto& residue: chain.residues) {
                residue.temp.try_emplace("turn",0);
                residue.temp.try_emplace("straight",0);
                residue.temp.try_emplace("win",any(0));
                residue.temp.try_emplace("dist",any(0.0));
                for(auto& atom: residue.atoms) {
                    atom.temp.try_emplace("dist",any(0.0));
                }
            }
        }
        logger.debug(" Processed: Optim::init()");
    }

    void Optim::end() {
        logger.debug("Processing: Optim::end()");
        for(auto& chain: protein.chains) {
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
        logger.debug(" Processed: Optim::end()");
    }

    void Optim::setDistances() {
        logger.debug("Processing: Optim::setDistances()");
        for(auto& chain: protein.chains) {
            for(auto& residue: chain.residues) {
                setAtomDistances(residue);
            }
        }
        logger.debug(" Processed: Optim::setDistances()");
    }

    void Optim::setAtomDistances(Tmdet::ValueObjects::Residue& residue) const {
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

    void Optim::setBoundaries() {
        logger.debug("Processing: Optim::setBoundaries()");
        min = 10000;
        max = -10000;
        for(auto& chain: protein.chains) {
            for(auto& residue: chain.residues) {
                for(auto& atom: residue.atoms) {
                    auto dist = any_cast<double>(atom.temp.at("dist"));
                    min = (dist < min ? dist : min);
                    max = (dist > max ? dist : max);
                }
            }
        }
        logger.debug("Box size: {} {}: {}",min,max,(unsigned int)(max-min+2));
        slices.clear();
        slices.resize((unsigned int)(max-min+2));
        logger.debug(" Processed: Optim::setBoundaries()");
    }

    void Optim::sumupSlices() {
        logger.debug("Processing: Optim::sumupSlices()");
        for(auto& chain: protein.chains) {
            for(auto& residue: chain.residues) {
                residueToSlice(residue);
            }
        }
        logger.debug(" Processed: Optim::sumupSlices()");
    }

    void Optim::residueToSlice(Tmdet::ValueObjects::Residue& residue) {
        auto sliceIndex = (unsigned int)(any_cast<double>(residue.temp.at("dist")) - min);
        slices[sliceIndex].numCA++;
        if (residue.ss.isAlpha() || residue.ss.isBeta()) {
            slices[sliceIndex].numStraight++;
        }
        else if (residue.ss.isTurn()) {
            slices[sliceIndex].numTurn++;
        }
        for(const auto& atom: residue.atoms) {
            slices[sliceIndex].surf +=  any_cast<double>(atom.temp.at("outside"));
            if (residue.type.atoms.contains(atom.gemmi.name)) {
                slices[sliceIndex].voronota += 
                    (any_cast<double>(atom.temp.at("outside")) * ( 1-
                    //"voronota frustration: it is small if residue does not like to be on surface"
                        (residue.type.atoms.at(atom.gemmi.name).mean - Tmdet::Types::voronotaMeanMin) / 
                            (Tmdet::Types::voronotaMeanMax - Tmdet::Types::voronotaMeanMin)));
            }
        }
    }

    double Optim::getQValueForSlice(const _slice& s) const {
        double q = 0.0;
        if (s.numCA > 0 && s.surf >0 ) {
            q += 40 * (double)s.numStraight / s.numCA;
            q += 10 * (1.0 - (double)s.numTurn / s.numCA);
            q += 60 * s.voronota / s.surf;
            logger.debug("Slice: ca: {} nStr: {} nT: {} voronota: {} surf: {} q: {}",s.numCA,s.numStraight,s.numTurn,s.voronota,s.surf,q);
        }
        return q;
    }

    std::vector<double> Optim::getQValueForSlices() {
        auto s = slices.size();
        std::vector<double> qs(s,0);
        for(unsigned long int i = 0; i<s; i++) {
            qs[i] = getQValueForSlice(slices[i]);
        }
        return qs;
    }


    double Optim::smoothQValues(std::vector<double> qs) {
        std::string m="**************************************************************************************************";
        std::string o="++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++";
        double maxQ = -1e30;
        auto s = slices.size();
        for(unsigned long int i = 0; i<s; i++) {
            int k=0;
            double q=0;
            for (int j=-3; j<=3; j++) {
                if (j+(int)i>=0 && j+i<s) {
                    q += qs[j+i];
                    k++;
                }
            }
            q /= k;
            q=(q<1?1:q);
            slices[i].qValue = q;
            if (q > maxQ) {
                maxQ = q;
                bestSliceIndex = i;
            }
            if (q>TMDET_MEMBRANE_QVALUE) {
                logger.debug("{}{}",
                    o.substr(0,TMDET_MEMBRANE_QVALUE),
                    m.substr(0,(int)(q-TMDET_MEMBRANE_QVALUE)));
            }
            else {
                logger.debug("{}",o.substr(0,(int)(q)));
            }
        }
        return maxQ;
    }

    void Optim::setNormal(gemmi::Vec3 _normal) {
        normal = _normal;
    }

    void Optim::testMembraneNormal() {
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

    bool Optim::isTransmembrane() const {
        return bestQ > TMDET_MINIMUM_QVALUE;
    }

    void Optim::setMembranesToProtein() {
        if (!isTransmembrane()) {
            return;
        }
        setNormal(bestNormal);
        testMembraneNormal();
        Tmdet::ValueObjects::Membrane membrane;
        protein.membranes.clear();
        while(getMembrane(membrane)) {
            protein.membranes.push_back(membrane);
        }
        if (!protein.membranes.empty()) {
            protein.tmp = true;
        }
    }

    void Optim::searchForMembraneNormal() {
        Tmdet::Engine::Rotator rotator;
        while(rotator.next(normal)) {
            testMembraneNormal();
        }
    }

    bool Optim::getMembrane(Tmdet::ValueObjects::Membrane& membrane) {
        logger.debug("Processing Optim::getMembrane()");
        unsigned long int bestZ = -1;
        double q = -1e30;
        for (unsigned long int i = 0; i <slices.size(); i++) {
            if (slices[i].qValue > q) {
                q = slices[i].qValue;
                bestZ = i;
            }
        }

        if (bestZ < TMDET_MEMBRANE_QVALUE) {
            logger.debug(" Processed Optim::getMembrane(): no more membrane");
            return false;
        }
        
        unsigned long int minz = bestZ; 
        while (minz>2 && slices[minz].qValue > TMDET_MEMBRANE_QVALUE) {
            minz--;
        }
        unsigned long int maxz = bestZ; 
        while (maxz<slices.size()-2 && slices[maxz].qValue > TMDET_MEMBRANE_QVALUE) {
            maxz++;
        }
        membrane.halfThickness = (maxz-minz) / 2;

        if (membrane.halfThickness < TMDET_MEMBRANE_MIN_HALFTHICKNESS) {
            logger.debug(" Processed Optim::getMembrane(): not transmembrane, membrane thickness is small");
            return false;
        }

        for (unsigned long int i=minz; i<=maxz; i++) {
            slices[i].qValue = 0;
        }

        double o = (minz+maxz) / 2;
        membrane.origo = massCentre + o * normal;
        membrane.normal = normal;
        setMembraneTMatrix(membrane);

        logger.debug(" Processed Optim::getMembrane() thickness: {} ",membrane.halfThickness);
        return true;
    }

    void Optim::setMembraneTMatrix(Tmdet::ValueObjects::Membrane& membrane) const {
        double x = membrane.normal.x;
	    double y = membrane.normal.y;
	    double z = membrane.normal.z;
	    if (double d = sqrt(y*y+z*z); d>1e-5) {
		    double sa=z/d;
		    double ca=y/d;
		    membrane.tmatrix.rot[0][0] = d;
            membrane.tmatrix.rot[0][1] = -x*ca;
            membrane.tmatrix.rot[0][2] = -x*sa;
            
		    membrane.tmatrix.rot[1][0] = 0;
            membrane.tmatrix.rot[1][1] = sa;
            membrane.tmatrix.rot[1][2] = -ca;
            
		    membrane.tmatrix.rot[2][0] = x;
            membrane.tmatrix.rot[2][1] = d*ca;
            membrane.tmatrix.rot[2][2] = d*sa;
        }
        else {
        	membrane.tmatrix.rot[0][0]=0;
            membrane.tmatrix.rot[0][1]=0;
            membrane.tmatrix.rot[0][2]=1;

            membrane.tmatrix.rot[1][0]=0;
            membrane.tmatrix.rot[1][1]=1;
            membrane.tmatrix.rot[1][2]=0;

            membrane.tmatrix.rot[2][0]=1;
            membrane.tmatrix.rot[2][1]=0;
            membrane.tmatrix.rot[2][2]=0;
        }
        membrane.tmatrix.trans = -membrane.origo;
    }
}
