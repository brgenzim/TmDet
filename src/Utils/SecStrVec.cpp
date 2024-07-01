#include <iostream>
#include <array>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <Types/Residue.hpp>
#include <Utils/Surface.hpp>
#include <Utils/SecStrVec.hpp>

using namespace std;

namespace Tmdet::Utils {
    static void dumpVectorsForPyMOL(vector<_secStrVec> &vectors);

    void SecStrVec::define(Tmdet::ValueObjects::TmdetStruct& tmdetVO) {
        vectors.clear();
        for(auto& chain: tmdetVO.chains) {
            int begin, end;
            begin = end = 0;
            while(getNextRegion(chain, begin, end)) {
                vectors.push_back(getVector(chain, begin, end));
            }
        }
        dumpVectorsForPyMOL(vectors);
    }

    bool SecStrVec::ifCross(_secStrVec& vec, Tmdet::ValueObjects::Membrane& membraneVO) {
        // más az implementáció PLANE és CURVE esetén
    }

    bool SecStrVec::getNextRegion(Tmdet::ValueObjects::Chain& chain, int& begin, int& end) {
        begin = end + 1;
        return (getNextNotUnkown(chain, begin) && getNextSame(chain, begin, end));
    }

    bool SecStrVec::getNextNotUnkown(Tmdet::ValueObjects::Chain& chain, int& begin) {
        while(begin < (int)chain.residues.size() && chain.residues[begin].ss == Tmdet::Types::SecStructType::U) {
            begin++;
        }
        return (begin < (int)chain.residues.size());
    }

    bool SecStrVec::getNextSame(Tmdet::ValueObjects::Chain& chain, int& begin, int& end) {
        end = begin;
        while(end < (int)chain.residues.size() && chain.residues[begin].ss == chain.residues[end].ss) {
            end++;
        }
        return true;
    }

    _secStrVec SecStrVec::getVector(Tmdet::ValueObjects::Chain& chain, int& begin, int& end) {
        return (chain.residues[begin].ss.isAlpha()?
                    getAlphaVector(chain, begin, end):
                    getBetaVector(chain,begin,end));
    }

    _secStrVec SecStrVec::getAlphaVector(Tmdet::ValueObjects::Chain& chain, int& begin, int& end) {
        _secStrVec vec;
        vec.begin = getMeanPosition(chain,begin);
        vec.end = getMeanPosition(chain,end-3);
        vec.type = Tmdet::Types::SecStructType::H;

        return vec;
    }

    gemmi::Vec3 SecStrVec::getMeanPosition(Tmdet::ValueObjects::Chain& chain, int pos) {
        gemmi::Vec3 vec;
        int i = 0;
        for (; i<3; i++) {
            if (auto CA = chain.residues[pos+i].gemmi.get_ca()) {
                vec += CA->pos;
                i++;
            }
        }
        vec /= i;
        return vec;
    }

    _secStrVec SecStrVec::getBetaVector(Tmdet::ValueObjects::Chain& chain, int& begin, int& end) {
        _secStrVec vec;
        vec.begin = chain.residues[begin].gemmi.get_ca()->pos;
        vec.end = chain.residues[end].gemmi.get_ca()->pos;
        vec.type = Tmdet::Types::SecStructType::E;

        return vec;
    }

    void dumpVectorsForPyMOL(vector<_secStrVec> &vectors) {
        int counter = 1;
        for (auto& vector : vectors) {
            std::cout << "cgo_arrow [ " << vector.begin.x << ", " << vector.begin.y << ", " << vector.begin.z << " ], "
                << "[ " << vector.end.x << ", " << vector.end.y << ", " << vector.end.z << " ], "
                << "name=" << vector.type.name << counter++ << ", color=yellow" << std::endl;
        }
    }
}
