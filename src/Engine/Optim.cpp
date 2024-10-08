#include <iostream>
#include <array>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <Types/Residue.hpp>
#include <Engine/Optim.hpp>
#include <Utils/Surface.hpp>

using namespace std;

namespace Tmdet::Engine {

    void Optim::init() {
        run = true;
        for(auto& chain: proteinVO.chains) {
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
    }

    void Optim::end() {
        for(auto& chain: proteinVO.chains) {
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
        run = false;
    }

    void Optim::setDistances() {
        for(auto& chain: proteinVO.chains) {
            for(auto& residue: chain.residues) {
                setAtomDistances(residue);
            }
        }
    }

    void Optim::setAtomDistances(Tmdet::ValueObjects::Residue& residue) {
        bool hasCA = false;
        double d;
        for(auto& atom: residue.atoms) {
            d = membraneVO.normal.x * (atom.gemmi.pos.x - membraneVO.origo.x);
            d += membraneVO.normal.y * (atom.gemmi.pos.y - membraneVO.origo.y);
            d += membraneVO.normal.z * (atom.gemmi.pos.z - membraneVO.origo.z);
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
        min = 10000;
        max = -10000;
        for(auto& chain: proteinVO.chains) {
            for(auto& residue: chain.residues) {
                for(auto& atom: residue.atoms) {
                    auto dist = any_cast<double>(atom.temp.at("dist"));
                    min = (dist < min ? dist : min);
                    max = (dist < max ? dist : max);
                }
            }
        }
    }

    void Optim::sumupSlices() {
        for(auto& chain: proteinVO.chains) {
            for(auto& residue: chain.residues) {
                residueToSlice(residue);
            }
        }
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
            //TODO check if it is better
            //sliceIndex = (unsigned int)(any_cast<double>(residue.temp.at("dist")) - min);
            slices[sliceIndex].surf += atom.surface;
            if (residue.type.atoms.contains(atom.gemmi.name)) {
                slices[sliceIndex].voronota += atom.surface * residue.type.atoms.at(atom.gemmi.name).mean;
            }
        }
    }

    double Optim::getQValueForSlice(const _slice& s) {
        double q;
        q = (double)s.numStraight / s.numCA;
        q -= (double)s.numTurn / s.numCA;
        q += s.voronota / s.surf;
        return q;
    }

    std::vector<double> Optim::getQValueForSlices() {
        unsigned int s = slices.size();
        std::vector<double> qs(s,0);
        for(unsigned int i = 0; i<s; i++) {
            qs[i] = getQValueForSlice(slices[i]);
        }
        return qs;
    }


    double Optim::smoothQValues(std::vector<double> qs) {
        double maxQ = -1e30;
        unsigned int s = slices.size();
        for(unsigned int i = 0; i<s; i++) {
            int k=0;
            double q=0;
            for (int j=-3; j<=3; j++) {
                if (j+(int)i>=0 && j+i<s) {
                    q += qs[j+i];
                    k++;
                }
            }
            q /= k;
            if (q > maxQ) {
                maxQ = q;
            }
        }
        return maxQ;
    }

    double Optim::getQValue() {
        setDistances();
        setBoundaries();
        sumupSlices();
        return smoothQValues(getQValueForSlices());
    }
}
