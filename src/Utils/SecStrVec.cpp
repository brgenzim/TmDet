#include <iostream>
#include <array>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <sstream>
#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <Types/Residue.hpp>
#include <Utils/Surface.hpp>
#include <Utils/SecStrVec.hpp>

using namespace std;

namespace Tmdet::Utils {
    static bool isVectorCrossingPlane(_secStrVec &vector, gemmi::Vec3 &planePoint, gemmi::Vec3 &planeNormal);

#ifdef __SECSTRVEC_DBG
    static string vec3ToString(gemmi::Vec3 &vector);
    static void dumpVectorsForPyMOL(vector<_secStrVec> &vectors);
#endif

    void SecStrVec::define(Tmdet::ValueObjects::TmdetStruct& tmdetVO) {
        vectors.clear();
        for(auto& chain: tmdetVO.chains) {
            int begin, end;
            begin = end = 0;
            while(getNextRegion(chain, begin, end)) {
                if (end - begin > 1) {
                    vectors.push_back(getVector(chain, begin, end - 1));
                }
                begin = end;
            }
        }
#ifdef __SECSTRVEC_DBG
        dumpVectorsForPyMOL(vectors);
#endif
    }

    void SecStrVec::numCross(Tmdet::ValueObjects::Membrane& membraneVO, int &numBoth, int &numUp, int &numDown) {
        // TODO: debug code
        auto result = ifCross(vectors[1], membraneVO);
        cout << "Crossing: " << result << endl;
    }


    bool SecStrVec::ifCross(_secStrVec& vec, Tmdet::ValueObjects::Membrane& membraneVO) {
        bool result = false;
        // más az implementáció PLANE és CURVE esetén
        if (membraneVO.type.name == Tmdet::Types::MembraneType::PLAIN.name) {
            // ATTENTION: membraneVO.normal must be normalized
            auto normal = membraneVO.normal;
            auto membranePoint = membraneVO.origo + normal * membraneVO.h;

            result = isVectorCrossingPlane(vec, membranePoint, normal);
            if (!result) {
                membranePoint = membraneVO.origo - normal * membraneVO.h;
                result = isVectorCrossingPlane(vec, membranePoint, normal);
            }
        } else {
            // TODO: CURVE
        }
        return result;
    }

    bool isVectorCrossingPlane(_secStrVec &vector, gemmi::Vec3 &planePoint, gemmi::Vec3 &planeNormal) {
        // https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection#Algebraic_form
        // plane points: (p - p0)*n = 0 (n is normal vector of plane)
        // line points:  p = l0 + l*d (l is unit vector)
        // ergo: d = ((p0 - l0)*n)/(l*n)

        auto l0 = vector.begin;
        auto l = (vector.end - vector.begin).normalized();
        auto n = planeNormal;
        auto p0 = planePoint;
        auto numerator = (p0 - l0).dot(n);
        auto denominator = l.dot(n);
        if (numerator != 0 && denominator == 0) {
            return false; // parallel
        } else if (numerator == 0 && denominator == 0) {
            return true; // vector lies in the plane
        } else { // denominator != 0
#ifdef __SECSTRVEC_DBG
            // intersection point - based on vector equation:
            // l0 + l*d
            auto intersectionPoint = l0 + l * (numerator / denominator);
            cout << endl << "Intersection Point: " << vec3ToString(intersectionPoint) << endl;
#endif
            return true; // intersects the plane
        }
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
                    getBetaVector(chain,begin,end));
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


    //
    // Util/Debug functions
    //

#ifdef __SECSTRVEC_DBG

    string vec3ToString(gemmi::Vec3 &vector) {
        std::stringstream stream;
        stream << "[ " << vector.x << ", " << vector.y << ", " << vector.z << " ]";
        return stream.str();
    }

    void dumpVectorsForPyMOL(vector<_secStrVec> &vectors) {
        int counter = 1;
        for (auto& vector : vectors) {
            std::cout << "cgo_arrow [ " << vector.begin.x << ", " << vector.begin.y << ", " << vector.begin.z << " ], "
                << "[ " << vector.end.x << ", " << vector.end.y << ", " << vector.end.z << " ], "
                << "name=" << vector.type.name << counter++ << ", color=yellow" << std::endl;
        }
    }

#endif

}
