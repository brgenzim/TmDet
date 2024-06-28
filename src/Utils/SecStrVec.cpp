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

    void SecStrVec::define(Tmdet::ValueObjects::TmdetStruct& tmdetVO) {
        vectors.clear();
        for(auto& chain: tmdetVO.chains) {
            int begin, end;
            while(getNextRegion(chain, begin, end)) {
                vectors.push_back(getVector(chain, begin, end));
            }
        }
    }

    bool SecStrVec::ifCross(_secStrVec& vec, Tmdet::ValueObjects::Membrane& membraneVO) {

    }

    bool SecStrVev::getNextRegion(Tmdet::ValueObjects::Chain& chain, int& begin, int& end) {
        return (getNextNotUnkown(chain, begin) && getNextSame(chain, begin, end));
    }

    bool SecStrVec::getNextNotUnkown(Tmdet::ValueObjects::Chain& chain, int& begin) {
        while(begin<chain.residues.size() && chain.residues[begin].ss != Tmdet::Types::SecStructType::U) {
            begin++;
        }
        return (begin==chain.residues.size());
    }

    bool SecStrVec::getNextSame(Tmdet::ValueObjects::Chain& chain, int& begin, int& end) {
        end = begin;
        while(end<chain.residues.size() && chain.residues[begin].ss == chain.residues[end].ss) {
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
        vec.type = Tmdet::Types::SecStruct::H;
    }

    gemmi::Vec3 SecStrVec::getMeanPosition(Tmdet::ValueObjects::Chain& chain, int pos) {
        gemmi::Vec3 vec = 0;
        for (int i=0, int c=0; i<3; i++) {
            if (auto CA = chain.residues[pos+i].gemmi.get_ca()) {
                vec.begin += CA.pos;
                j++;
            }
        }
        vec /= j;
        return vec;
    }
    
    _secStrVec SecStrVec::getBetaVector(Tmdet::ValueObjects::Chain& chain, int& begin, int& end) {
        
    }
}
 