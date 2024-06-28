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

    void Optim::init() {
        run = true;
        for(auto& chain: tmdetVO.chains) {
            for(auto& residue: chain.residues) {
                residue.temp.insert({"turn",any_cast<int>(0)});
                residue.temp.insert({"straight",any_cast<int>(0)});
                residue.temp.insert({"win",any_cast<int>(0)});
                for(auto& atom: residue.atoms) {
                    atom.temp.insert({"dist",any_cast<float>(0.0)});
                }
            }
        }
    }

    void Optim::end() {
        for(auto& chain: tmdetVO.chains) {
            for(auto& residue: chain.residues) {
                residue.temp.erase("turn");
                residue.temp.erase("straight");
                residue.temp.erase("win");
                for(auto& atom: residue.atoms) {
                    atom.temp.erase("dist");
                }
            }
        }
        run = false;
    }

    void Optim::setDistances() {
        for(auto& chain: tmdetVO.chains) {
            for(auto& residue: chain.residues) {
                for(auto& atom: residue.atoms) {
                    float d = membraneVO.normal.x * (atom.pos.x - membraneVO.origo.x);
                    d += membraneVO.normal.y * (atom.pos.y - membraneVO.origo.y);
                    d += membraneVO.normal.z * (atom.pos.z - membraneVO.origo.z);
                    atom.temp.at("dist") = any_cast<float>(d);
                }
            }
        }
    }

    void Optim::setBoundaries() {
        min = 10000;
        max = -10000;
        for(auto& chain: tmdetVO.chains) {
            for(auto& residue: chain.residues) {
                for(auto& atom: residue.atoms) {
                    min = (atom.temp.at("dist")<min?atom.temp.at("dist"):min);
                    max = (atom.temp.at("dist")<max?atom.temp.at("dist"):max);
                }
            }
        }
    }
}