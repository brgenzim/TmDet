#include <iostream>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <span>
#include <eigen3/Eigen/Dense>
#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <Config.hpp>
#include <Types/Residue.hpp>
#include <ValueObjects/Protein.hpp>
#include <Utils/Symmetry.hpp>
#include <Utils/Oligomer.hpp>
//#include <System/Logger.hpp>

using namespace std;

#define TINY   1.0e-10
#define RAD2DEG   57.29576

namespace Tmdet::Utils {

    static Eigen::Matrix4d Rotate_Z(Eigen::Vector3d T);
    static double cosineAngleOfVectors(Eigen::Vector3d u, gemmi::Vec3 v);

    void Symmetry::run() {
        logger.debug("Processing Symmetry::run()");
        initSymetryContainer();
        auto axes = getRotationalAxes();
        clusterAxes(axes);
        logger.debug(" Processed Symmetry::run()");
    }

    void Symmetry::initSymetryContainer() {
        logger.debug("Processing Symmetry::initSymmetryContainer()");
        const auto& nc = protein.chains.size();
        sim.resize(nc);
        for (auto& s: sim) {
            s.resize(nc);
        }
        logger.debug(" Processed Symmetry::initSymmetryContainer()");
    }

    std::vector<_symmetryData> Symmetry::getRotationalAxes() {
        logger.debug("Processing Symmetry::getRotationalAxes()");
        for (const auto& chain: protein.chains) {
            if (chain.selected && Oligomer::isEntityOligomerized(chain.entityId, protein.gemmi)) {
                searchForRotatedChains(chain.idx);
            }
        }
        std::vector<_symmetryData> axes;
        for (const auto& chain1: protein.chains ) {
            for (const auto& chain2: protein.chains) {
                auto& s = sim[chain1.idx][chain2.idx];
                if (s.id == 1) {
                    s.id = 0;
                    s.cidx1 = chain1.idx;
                    s.cidx2 = chain2.idx;
                    s.entityId = chain1.entityId;
                    axes.emplace_back(s);
                    logger.debug("Chains {} and {} are rotated", chain1.id, chain2.id);
                    logger.debug("\t Origo: {} {} {}", s.origo.x, s.origo.y, s.origo.z);
                    logger.debug("\tNormal: {} {} {}", s.axis.x, s.axis.y, s.axis.z);
                    logger.debug("\t Angle: {}", s.rotAngle);
                }
            }
        }
        logger.debug(" Processed Symmetry::getRotationalAxes()");
        return axes;
    }

    void Symmetry::searchForRotatedChains(int cidx1) {
        logger.debug("Processing Symmetry::searchForRotatedChains({})",protein.chains[cidx1].id);
        for(const auto& chain2: protein.chains) {
            if (chain2.idx > cidx1 && chain2.selected && protein.chains[cidx1].entityId == chain2.entityId) {
                calculateRotationalOperation(cidx1,chain2.idx);
            }
        }
        logger.debug(" Processed Symmetry::searchForRotatedChains({})",protein.chains[cidx1].id);
    }

    void Symmetry::calculateRotationalOperation(int cidx1, int cidx2) {
        logger.debug("Processing Symmetry::calculateRotationalOperation({}:{})",
            protein.chains[cidx1].id,protein.chains[cidx2].id);
        std::vector<Eigen::Vector3d> coord1;
        std::vector<Eigen::Vector3d> coord2;
        Eigen::Vector3d t1;
        Eigen::Vector3d t2;
        getCoordinates(cidx1,cidx2,coord1,coord2, t1, t2);
        double rmsd;
        Eigen::Matrix4d R;
        std::span<Eigen::Vector3d> coord1Slice(coord1.begin(), coord1.size());
        std::span<Eigen::Vector3d> coord2Slice(coord2.begin(), coord2.size());
        lsqFit(coord1Slice, coord2Slice, rmsd, R);
        getSymmetryOperand( R, t1, t2, sim[cidx1][cidx2]);
        double distance = (t2 - t1).squaredNorm();
        logger.debug("results {} {} distance: {} rmsd: {}",cidx1,cidx2,distance,rmsd);
        if (distance > 2.0 && rmsd < 5 /** orig: 3 */) {
            sim[cidx1][cidx2].id = sim[cidx1][cidx2].id = 1;
        }

        if (distance < 2.0 && R(0, 0) > 0.98 && R(1, 1) > 0.98 && R(2, 2) > 0.98) {
            logger.warn("Chain {} and {} are on top of each other", 
                protein.chains[cidx1].id, protein.chains[cidx2].id);
            sim[cidx1][cidx2].id = sim[cidx1][cidx2].id = -1;
        }
        else if (fabs(cosineAngleOfVectors((t1 - t2), sim[cidx1][cidx2].axis)) > 0.15) {
            logger.warn("Chain {} and {} are not simply rotated",
                protein.chains[cidx1].id, protein.chains[cidx2].id);
            sim[cidx1][cidx2].id = sim[cidx1][cidx2].id = -1;
        }
        logger.debug(" Processed Symmetry::calculateRotationalOperation({}:{})",
            protein.chains[cidx1].id,protein.chains[cidx2].id);
    }
    
    void Symmetry::getCoordinates(int cidx1, int cidx2, std::vector<Eigen::Vector3d>& coord1, std::vector<Eigen::Vector3d>& coord2, Eigen::Vector3d& t1, Eigen::Vector3d& t2) {
        logger.debug("Processing Symmetry::getCoordinates({}:{})",
        protein.chains[cidx1].id,protein.chains[cidx2].id);
        const auto& length = protein.chains[cidx1].length;
        gemmi::Vec3 centre1;
        gemmi::Vec3 centre2;
        unsigned int nca = 0;
        for (int idx = 0; idx<length; idx++) {
            auto ca1 = protein.chains[cidx1].residues[idx].gemmi.get_ca();
            auto ca2 = protein.chains[cidx2].residues[idx].gemmi.get_ca();
            if (ca1 && ca2) {
                centre1 += ca1->pos;
                centre2 += ca2->pos;
                nca++;
            }
        }
        centre1 /= nca;
        centre2 /= nca;
        t1 = Eigen::Vector3d(centre1.x, centre1.y, centre1.z);
        t2 = Eigen::Vector3d(centre2.x, centre2.y, centre2.z);
        coord1.resize(nca);
        coord2.resize(nca);
        nca = 0;
        for (int idx = 0; idx<length; idx++) {
            auto ca1 = protein.chains[cidx1].residues[idx].gemmi.get_ca();
            auto ca2 = protein.chains[cidx2].residues[idx].gemmi.get_ca();
            if (ca1 && ca2) {
                coord1[nca] = Eigen::Vector3d(ca1->pos.x - centre1.x, ca1->pos.y - centre1.y, ca1->pos.z - centre1.z);
                coord2[nca] = Eigen::Vector3d(ca2->pos.x - centre2.x, ca2->pos.y - centre2.y, ca2->pos.z - centre2.z);
                nca++;
            }
        }
        logger.debug(" Processed Symmetry::getCoordinates({}:{})",
            protein.chains[cidx1].id,protein.chains[cidx2].id);
    }
   
    void Symmetry::getSymmetryOperand(Eigen::Matrix4d& R, Eigen::Vector3d& t1, Eigen::Vector3d& t2, _symmetryData& simij) {
        logger.debug("Processing Symmetry::getSymmetryOperand()");
        Eigen::Matrix4d LR, Rot;
        Eigen::Vector3d av, ori, axis;
        double tg;

        if ((R(0, 0)<0.999)&&(R(1, 1)<0.999) /*&&(R(0, 0)>-0.999)*/)
        {
            /* Find vector which is invariant of the rotation
            lets fix its size by setting z to 1.0
            */
            axis.z() = 1.0;
            axis.y() = 
                ((1.0-R(0, 0))*R(1, 2)+R(1, 0)*R(0, 2))/
                ((1.0-R(1, 1))*(1.0-R(0, 0))-R(0, 1)*R(1, 0));
            if (R(0, 0)<1.0) {
                axis.x() = (R(0, 1)*axis.y() + R(0, 2))/
                    (1.0-R(0, 0));
            } else {
                axis.x() = 0;
            }
        } else {
            if (R(0, 0)>0.999) {
                axis.x() = 1;
            }
            axis.y() = 0;
            axis.z() = 0;
            if (R(1, 1)>0.999) {
                axis.x() = 0;
            }
            axis.y() = 1;
            axis.z() = 0;

        }
        axis = axis / axis.norm();
        LR=Rotate_Z(axis);

        Rot= LR.inverse() * (R * LR);

        if (Rot(0, 0)>1) Rot(0, 0)=1.0;
        if (Rot(0, 0)<-1) Rot(0, 0)=-1.0;

        if (Rot(0, 0)!=1.0)
            tg=Rot(0, 1)/(1-Rot(0, 0));
            else tg=1.0;

        Eigen::Vector3d midpoint = (t1 + t2) * 0.5;
        av = midpoint - t1;
        ori = av.cross(axis);
        if (tg!=0) ori = ori * tg;

        ori += midpoint;
        simij.origo = gemmi::Vec3(ori.x(), ori.y(), ori.z());

        simij.rotAngle = RAD2DEG*acos(Rot(0, 0));
        simij.axis = gemmi::Vec3(axis.x(), axis.y(), axis.z());
        logger.debug(" Processed Symmetry::getSymmetryOperand()");
    }

    bool Symmetry::lsqFit(std::span<Eigen::Vector3d>& r1, std::span<Eigen::Vector3d>& r2, double& rmsd, Eigen::Matrix4d& Rot) {
        logger.debug("Processing Symmetry::Lsq_fit");
        Eigen::Matrix3d U = Eigen::Matrix3d::Zero();
        auto nr_atoms = r1.size();

        for (int n = 0; n < nr_atoms; n++) {
            U += /*weight[n] * */ r1[n] * r2[n].transpose();
        }

        double det_U = U.determinant();
        if (std::abs(det_U) < TINY) {
            logger.error("{}.{}: determinant of U equals to zero", __FILE__, __LINE__);
            return false;
        }

        double sign_detU = (det_U > 0) ? 1.0 : -1.0;

        Eigen::MatrixXd omega(6, 6);
        omega.setZero();
        omega.block<3, 3>(0, 3) = U;
        omega.block<3, 3>(3, 0) = U.transpose();

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(omega);
        if (eigensolver.info() != Eigen::Success) {
            logger.error("{}.{}: Eigen decomposition failed", __FILE__, __LINE__);
            return false;
        }

        Eigen::VectorXd eva_omega = eigensolver.eigenvalues();
        Eigen::MatrixXd eve_omega = eigensolver.eigenvectors();

        if (det_U < 0.0) {
            if (std::abs(eva_omega(1) - eva_omega(2)) < TINY) {
                logger.error("{}.{}: determinant of U < 0 && degenerated eigenvalues", __FILE__, __LINE__);
                return false;
            }
        }
        
        Eigen::Matrix3d H, K;
        for (int i = 0; i < 3; i++) {
            H.col(i) = M_SQRT2 * eve_omega.col(i).head<3>();
            K.col(i) = M_SQRT2 * eve_omega.col(i).tail<3>();
        }

        double sigma = (H.col(1).cross(H.col(2))).dot(H.col(0));
        if (sigma <= 0.0) {
            H.col(2) = -H.col(2);
            K.col(2) = -K.col(2);
        }

        Eigen::Matrix3d R = K * H.transpose();
        R *= sign_detU;
        R *= -1;

        for (int n = 0; n < nr_atoms; n++) {
            r1[n] = R * r1[n];
        }

        rmsd = 0.0;
        for (int n = 0; n < nr_atoms; n++) {
            Eigen::Vector3d dr = r1[n] - r2[n];
            //rmsd += dr.squaredNorm();
            rmsd += dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        }
        rmsd /= nr_atoms;
        rmsd = std::sqrt(rmsd);

        Rot.block<3,3>(0,0) = R;
        Rot(3, 3) = 1.0;
        for (int i = 0; i < 3; i++) {
            Rot(i, 3) = 0.0;
            Rot(3, i) = 0.0;
        }
        logger.debug(" Processed Symmetry::Lsq_fit");
        return true;
    }

    Eigen::Matrix4d Rotate_Z(Eigen::Vector3d T) {
        double a, b, c, d, cosa, sina;
        Eigen::Matrix4d R1, R2, R;
        Eigen::Vector3d normalized;

        normalized = T / T.norm();

        a = normalized.x(); b = normalized.y(); c = normalized.z();

        d = sqrt(b*b + c*c);
        if (d != 0) {
            sina = c/d; cosa = b/d;

            R1 <<  1,      0,     0,  0,
                   0,   sina,  cosa,  0,
                   0,  -cosa,  sina,  0,
                   0,      0,     0,  1;

            R2 <<  d,      0,     a,  0,
                   0,      1,     0,  0,
                  -a,      0,     d,  0,
                   0,      0,     0,  1;

            R = R1 * R2;
        } else {
            R <<  0,  1,  0,  0,
                  0,  0,  1,  0,
                  1,  0,  0,  0,
                  0,  0,  0,  1;
        }

        return R;
    }

    double cosineAngleOfVectors(Eigen::Vector3d u, gemmi::Vec3 v) {
        Eigen::Vector3d v_(v.x, v.y, v.z);
        double dotProduct = u.dot(v_);
        return dotProduct / (u.norm() * v_.norm());
    }

    std::vector<_symmetryData> Symmetry::clusterAxes(std::vector<_symmetryData>& axes) {
        std::vector<_symmetryData> ret;
        int cl=0;
        for (unsigned long int i = 0; i<axes.size(); i++) {
            if (axes[i].id == 0) {
                for (unsigned long int j = i+1; j<axes.size(); j++) {
                    if (axes[i].entityId == axes[j].entityId) {
                        logger.debug("Sim dist: {}:{} - {}:{}= {}",
                            protein.chains[axes[i].cidx1].id,
                            protein.chains[axes[i].cidx2].id,
                            protein.chains[axes[j].cidx1].id,
                            protein.chains[axes[j].cidx2].id,
                            axes[i].distance(axes[j]));
                        if (axes[i].same(axes[j])) {
                            logger.debug("Axes {} and {} are the same",i,j);
                            if (axes[i].id==0) {
                                cl++;
                                axes[i].id=cl;
                            }
                            if (axes[j].id == 0) {
                                axes[j].id = axes[i].id;
                            }
                            else {
                                for (unsigned long int k = 0; k<axes.size(); k++) {
                                    if (axes[k].id == axes[j].id) {
                                        axes[k].id = axes[i].id;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        for (auto& a: axes) {
            logger.debug("Cluster: {} {} {}",protein.chains[a.cidx1].id, protein.chains[a.cidx2].id, a.id);
        }
        return ret;
    }

}
