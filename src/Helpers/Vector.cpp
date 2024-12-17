// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#include <string>
#include <format>
#include <Config.hpp>
#include <System/Logger.hpp>
#include <gemmi/model.hpp>
#include <Helpers/Vector.hpp>

namespace Tmdet::Helpers::Vector {

    bool isPointOnVector(gemmi::Vec3 &vector, gemmi::Vec3 &vectorBegin, gemmi::Vec3 &point) {
        // if same direction and between begin-end (length less than vector length)
        return ((point - vectorBegin).dot(vector) > 0
                && (point - vectorBegin).length_sq() < vector.length_sq());
    }

    std::string vec3ToString(const gemmi::Vec3& vec) {
        return std::format("[{}, {}, {}]",vec.x, vec.y, vec.z);
    }

    bool doesVectorCrossPlane(gemmi::Vec3& begin, gemmi::Vec3& end, gemmi::Vec3 planeNormal, gemmi::Vec3 planePoint) {
        bool result = false;

        auto vectorDiff = end - begin;
        auto l = vectorDiff.normalized();
        auto numerator = (planePoint - begin).dot(planeNormal);
        auto denominator = l.dot(planeNormal);
        if (numerator != 0 && denominator == 0) {
            result = false; // parallel
        } else if (numerator == 0 && denominator == 0) {
            result = true; // vector lies in the plane
        } else { // denominator != 0
            // intersection point - based on vector equation: begin + l*d
            auto intersectionPoint = begin + l * (numerator / denominator);
            if (Tmdet::Helpers::Vector::isPointOnVector(vectorDiff, begin, intersectionPoint)) {
                result = true; // intersects the plane and between BEGIN and END
                DEBUG_LOG("Intersection Point: {}",vec3ToString(intersectionPoint));
            }
        }

        return result;
    }

    bool simplifiedDoesVectorCrossPlane(double begin, double end, double plane) {
        DEBUG_LOG("simplifiedDoesVectorCrossPlane: {} {} {}",begin,end,plane);
        return ((begin<plane && plane<end) || (begin>plane && plane>end));
    }

    bool doesVectorCrossSphere(gemmi::Vec3& begin, gemmi::Vec3& end, gemmi::Vec3 origo, double radius) {
        bool result = false;

        auto vectorDiff = end - begin;
        auto l = vectorDiff.normalized();
        auto m = begin - origo;
        auto A = l.dot(l);
        auto B = 2 * m.dot(l);
        auto C = m.dot(m) - radius * radius;
        auto discriminant = B*B - 4*A*C;

        // one intersection point
        if (discriminant == 0) {
            auto d = -B / 2*A;
            // intersection point - based on vector equation: begin + l*d
            auto intersectionPoint = begin + l * d;
            // if same direction and between begin-end (length less than vector length)
            if (Tmdet::Helpers::Vector::isPointOnVector(vectorDiff, begin, intersectionPoint)) {
                result = true; // intersects the plane and between BEGIN and END
                DEBUG_LOG("Intersection Point: {}",Tmdet::Helpers::Vector::vec3ToString(intersectionPoint));
            }
        } else if (discriminant > 0) {
            // line and sphere have two intersection points
            auto sqrtDiscriminant = sqrt(discriminant);
            auto d1 = (-B + sqrtDiscriminant) / 2*A;
            auto d2 = (-B - sqrtDiscriminant) / 2*A;
            auto intersectionPoint1 = begin + l * d1;
            auto intersectionPoint2 = begin + l * d2;
            if (Tmdet::Helpers::Vector::isPointOnVector(vectorDiff, begin, intersectionPoint1)) {
                result = true;
                DEBUG_LOG("Intersection Point1: {}",Tmdet::Helpers::Vector::vec3ToString(intersectionPoint1));
            }
            if (Tmdet::Helpers::Vector::isPointOnVector(vectorDiff, begin, intersectionPoint2)) {
                result = true;
                DEBUG_LOG("Intersection Point2: {}",Tmdet::Helpers::Vector::vec3ToString(intersectionPoint2));
            }
        }

        return result;
    }

    double distanceFromLine(const gemmi::Vec3& a, const gemmi::Vec3& b, const gemmi::Vec3& c) {
        auto ab = b - a;
        auto ac = c - a;
        auto crossProduct = ab.cross(ac);
        return crossProduct.length() / ab.length();
    }

    double angle(const gemmi::Vec3& a, const gemmi::Vec3& b) {
        return acos(cosAngle(a,b)) * 180.0 / M_PI;
    }

    double cosAngle(const gemmi::Vec3& a, const gemmi::Vec3& b) {
        return (a.dot(b) / (a.length() * b.length()) );
    }
}