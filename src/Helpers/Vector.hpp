#pragma once

#include <string>
#include <format>
#include <System/Logger.hpp>
#include <gemmi/model.hpp>

namespace Tmdet::Helpers::Vector {

    /**
     * @brief check if a point is on a vector
     * 
     * @param vector 
     * @param vectorBegin 
     * @param point 
     * @return true 
     * @return false 
     */
    extern bool isPointOnVector(gemmi::Vec3 &vector, gemmi::Vec3 &vectorBegin, gemmi::Vec3 &point);

    /**
     * @brief convert gemmi::Vec3 into string
     * 
     * @param vec 
     * @return std::string 
     */
    extern std::string vec3ToString(const gemmi::Vec3& vec);

    /**
     * @brief check if a vector crosses a plane
     * 
     *        https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection#Algebraic_form
     *        plane points: (p - planePoint)*planeNormal = 0 (planeNormal is normal vector of plane)
     *        line points:  p = begin + l*d (l is unit vector)
     *        ergo: d = ((planePoint - begin)*planeNormal)/(l*planeNormal)
     *
     * @param begin 
     * @param end 
     * @param planeNormal 
     * @param planePoint 
     * @return true 
     * @return false 
     */
    extern bool doesVectorCrossPlane(gemmi::Vec3& begin, gemmi::Vec3& end, gemmi::Vec3 planeNormal, gemmi::Vec3 planePoint);

    extern bool simplifiedDoesVectorCrossPlane(double begin, double end, double plane);

    /**
     * @brief check if a vector crosses a sphere
     * 
     *        sphere points: |p-origo| = r (r is radii)
     *        line points:  p = begin + l*d (l is unit vector, begin is a fix point of the line)
     *        ergo: |begin + l*d - origo| = r
     *        ... after some arragements:
     *        d = (-B +/- SQRT(B^2 - 4AC)) / 2A,
     *        
     *        where:
     *        
     *          m = begin-origo
     *          A = l*l
     *          B = 2*m*l
     *          C = m*m - r^2
     * @param vector 
     * @param spherePoint 
     * @param origo 
     * @return true 
     * @return false 
     */
    extern bool doesVectorCrossSphere(gemmi::Vec3& begin, gemmi::Vec3& end, gemmi::Vec3 origo, double radius);

    extern double distanceFromLine(const gemmi::Vec3& a, const gemmi::Vec3& b, const gemmi::Vec3& c);

    extern double angle(const gemmi::Vec3& a, const gemmi::Vec3& b);

    extern double cosAngle(const gemmi::Vec3& a, const gemmi::Vec3& b);
}