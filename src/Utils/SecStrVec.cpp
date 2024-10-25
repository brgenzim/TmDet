#include <iostream>
#include <array>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <sstream>
#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <Config.hpp>
#include <Helpers/Vector.hpp>
#include <Helpers/Pymol.hpp>
#include <System/Logger.hpp>
#include <Types/Residue.hpp>
#include <Utils/Surface.hpp>
#include <Utils/SecStrVec.hpp>
#include <ValueObjects/Protein.hpp>
#include <ValueObjects/Chain.hpp>

using namespace std;

namespace Tmdet::Utils {
    
    void SecStrVec::define() {
        DEBUG_LOG("Processing SecStrVec::define()");
        vectors.clear();
        for(auto& chain: protein.chains) {
            int begin, end;
            begin = end = 0;
            while(getNextRegion(chain, begin, end)) {
                if (end - begin > 1) {
                    vectors.push_back(getVector(chain, begin, end - 1));
                }
                begin = end;
            }
        }
        DEBUG_LOG("{}",Tmdet::Helpers::Pymol::dumpSecStrVec(vectors,"yellow"));
        DEBUG_LOG(" Processed SecStrVec::define(#vectors: {})",vectors.size());
    }

    void SecStrVec::numCrossingAlpha(Tmdet::ValueObjects::Membrane& membrane, int &numBoth, int &numUp, int &numDown) {
        numBoth = numUp = numDown = 0;
        for (auto& vector : vectors) {
            if (vector.type.isAlpha()) {
                checkCross(vector, membrane, numBoth, numUp, numDown);
            }
        }
    }

    void SecStrVec::numCrossingBeta(Tmdet::ValueObjects::Membrane& membrane, int &numBoth, int &numUp, int &numDown) {
        numBoth = numUp = numDown = 0;
        for (auto& vector : vectors) {
            if (vector.type.isBeta()) {
                checkCross(vector, membrane, numBoth, numUp, numDown);
            }
        }
    }

    bool SecStrVec::checkCross(_secStrVec& vec, Tmdet::ValueObjects::Membrane& membrane, int& numBoth, int& numUp, int& numDown) {
        bool resultUp = false;
        bool resultDown = false;

        if (membrane.type == Tmdet::Types::MembraneType::PLAIN) {
            resultUp = Tmdet::Helpers::Vector::doesVectorCrossPlane(
                    vec.begin, vec.end, gemmi::Vec3(0,0,1), gemmi::Vec3(0,0,membrane.origo + 5.0));
            resultDown = Tmdet::Helpers::Vector::doesVectorCrossPlane(
                    vec.begin, vec.end, gemmi::Vec3(0,0,1), gemmi::Vec3(0,0,membrane.origo - 5.0));
        } else if (membrane.type == Tmdet::Types::MembraneType::CURVED) {
            resultUp = Tmdet::Helpers::Vector::doesVectorCrossSphere(
                    vec.begin, vec.end, gemmi::Vec3(0,0,membrane.origo), sphereRadius + 5.0);
            resultDown = Tmdet::Helpers::Vector::doesVectorCrossSphere(
                    vec.begin, vec.end, gemmi::Vec3(0,0,membrane.), membrane.sphereRadius - 5.0);
        } else {
            throw runtime_error("Unexpected membrane type: " + membrane.type.name);
        }

        if (resultUp && resultDown) {
            numBoth++;
        } else if (resultUp) {
            numUp++;
        } else if (resultDown) {
            numDown++;
        }

        return resultUp || resultDown;
    }

    

    
    bool SecStrVec::getNextRegion(Tmdet::ValueObjects::Chain& chain, int& begin, int& end) {
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

    _secStrVec SecStrVec::getVector(Tmdet::ValueObjects::Chain& chain, int begin, int end) {
        return (chain.residues[begin].ss.isAlpha()?
                    getAlphaVector(chain, begin, end):
                    getBetaVector(chain,begin,end),
                    chain.idx, begin, end);
    }

    _secStrVec SecStrVec::getAlphaVector(Tmdet::ValueObjects::Chain& chain, int begin, int end) {
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
            auto CA = chain.residues[pos+i].gemmi.get_ca();
            vec += CA->pos;
        }
        vec /= i;
        return vec;
    }

    _secStrVec SecStrVec::getBetaVector(Tmdet::ValueObjects::Chain& chain, int begin, int end) {
        _secStrVec vec;
        vec.begin = chain.residues[begin].gemmi.get_ca()->pos;
        vec.end = chain.residues[end].gemmi.get_ca()->pos;
        vec.type = Tmdet::Types::SecStructType::E;

        return vec;
    }

}
