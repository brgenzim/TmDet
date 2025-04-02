// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

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
#include <System/Logger.hpp>
#include <Types/Residue.hpp>
#include <VOs/Protein.hpp>
#include <Utils/Symmetry.hpp>
#include <Utils/Oligomer.hpp>

using namespace std;

namespace Tmdet::Utils {

    std::vector<gemmi::Vec3> Symmetry::getMembraneAxes() {
        return clusterAxes(getRotationalAxes());
    }

    std::vector<_symmetryData> Symmetry::getRotationalAxes() {
        std::vector<_symmetryData> axes;
        for (const auto& entity: Tmdet::Utils::Oligomer::getHomoOligomerEntities(protein.gemmi)) {
            if (searchForRotatedChains(entity.subchains) && haveSameAxes()) {
                    axes.emplace_back(getAverageAxes());
            }
        }
        return axes;
    }

    bool Symmetry::searchForRotatedChains(const std::vector<std::string>& chainIds) {
        std::string ids = "";
        int i=0;
        for(const auto& c: chainIds) {
            ids += c+":";
            i++;
        }
        if (i>0) {
            protein.hasIdenticalChains = true;
        }
        int numChains = 0 ;
        int numRotated = 0;
        sim.clear();
        for(const auto& chain1Id: chainIds) {
            if (auto cidx1 = protein.searchChainByLabId(chain1Id); cidx1 != -1) {
                for(const auto& chain2Id: chainIds) {
                    if (auto cidx2 = protein.searchChainByLabId(chain2Id); cidx2 != -1) {
                        if (cidx1 != cidx2) {
                            numRotated += calculateRotationalOperation(cidx1,cidx2);
                            numChains++;
                        }
                    }
                }
                break;
            }
        }
        return numRotated > 0;
    }

    int Symmetry::calculateRotationalOperation(int cidx1, int cidx2) {

        std::vector<Eigen::Vector3d> coord1;
        std::vector<Eigen::Vector3d> coord2;
        Eigen::Vector3d t1;
        Eigen::Vector3d t2;
        getCoordinates(cidx1,cidx2,coord1,coord2, t1, t2);
        double rmsd;
        Eigen::Matrix4d R;
        _symmetryData curSim;
        std::span<Eigen::Vector3d> coord1Slice(coord1.begin(), coord1.size());
        std::span<Eigen::Vector3d> coord2Slice(coord2.begin(), coord2.size());
        lsqFit(coord1Slice, coord2Slice, rmsd, R);
        getSymmetryOperand( R, t1, t2, curSim);
        double distance = (t2 - t1).squaredNorm();
        if (distance > 2.0 && rmsd < 12 ) {
            curSim.good = true;
            sim.emplace_back(curSim);
        }

        return curSim.good?1:0;
    }

    void Symmetry::getCoordinates(int cidx1, int cidx2, std::vector<Eigen::Vector3d>& coord1, std::vector<Eigen::Vector3d>& coord2, Eigen::Vector3d& t1, Eigen::Vector3d& t2) {

        gemmi::Vec3 centre1;
        gemmi::Vec3 centre2;
        unsigned int nca = 0;
        int idx1 = 0;
        int idx2 = 0;
        while (idx1 < protein.chains[cidx1].length && idx2 < protein.chains[cidx2].length) {
            while (idx1 < protein.chains[cidx1].length
                && protein.chains[cidx1].residues[idx1].authId < protein.chains[cidx2].residues[idx2].authId) {
                idx1++;
            }
            while (idx2 < protein.chains[cidx2].length
                && protein.chains[cidx1].residues[idx1].authId > protein.chains[cidx2].residues[idx2].authId) {
                idx2++;
            }
            if (idx1 < protein.chains[cidx1].length && idx2 < protein.chains[cidx2].length) {
                auto ca1 = protein.chains[cidx1].residues[idx1].gemmi.get_ca();
                auto ca2 = protein.chains[cidx2].residues[idx2].gemmi.get_ca();
                if (ca1 && ca2) {
                    centre1 += ca1->pos;
                    centre2 += ca2->pos;
                    nca++;
                }
                idx1++;
                idx2++;
            }
        }
        centre1 /= nca;
        centre2 /= nca;
        t1 = Eigen::Vector3d(centre1.x, centre1.y, centre1.z);
        t2 = Eigen::Vector3d(centre2.x, centre2.y, centre2.z);
        coord1.resize(nca);
        coord2.resize(nca);
        nca = 0;
        idx1=idx2=0;
        while (idx1 < protein.chains[cidx1].length && idx2 < protein.chains[cidx2].length) {
            while (idx1 < protein.chains[cidx1].length
                && protein.chains[cidx1].residues[idx1].authId < protein.chains[cidx2].residues[idx2].authId) {
                idx1++;
            }
            while (idx2 < protein.chains[cidx2].length
                && protein.chains[cidx1].residues[idx1].authId > protein.chains[cidx2].residues[idx2].authId) {
                idx2++;
            }
            if (idx1 < protein.chains[cidx1].length && idx2 < protein.chains[cidx2].length) {
                auto ca1 = protein.chains[cidx1].residues[idx1].gemmi.get_ca();
                auto ca2 = protein.chains[cidx2].residues[idx2].gemmi.get_ca();
                if (ca1 && ca2) {
                    coord1[nca] = Eigen::Vector3d(ca1->pos.x - centre1.x, ca1->pos.y - centre1.y, ca1->pos.z - centre1.z);
                    coord2[nca] = Eigen::Vector3d(ca2->pos.x - centre2.x, ca2->pos.y - centre2.y, ca2->pos.z - centre2.z);
                    nca++;
                }
                idx1++;
                idx2++;
            }
        }
    }

    void Symmetry::getSymmetryOperand(Eigen::Matrix4d& R, const Eigen::Vector3d& t1, const Eigen::Vector3d& t2, _symmetryData& simij) const {
        Eigen::Matrix4d LR;
        Eigen::Matrix4d Rot;
        Eigen::Vector3d av;
        Eigen::Vector3d ori;
        Eigen::Vector3d axis;
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
                axis.y() = 0;
                axis.z() = 0;
            }
            if (R(1, 1)>0.999) {
                axis.x() = 0;
                axis.y() = 1;
                axis.z() = 0;
            }
        }

        axis = axis / axis.norm();
        LR=rotateZ(axis);

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
        simij.axis = gemmi::Vec3(axis.x(), axis.y(), axis.z());
    }

    bool Symmetry::lsqFit(const std::span<Eigen::Vector3d>& r1, const std::span<Eigen::Vector3d>& r2, double& rmsd, Eigen::Matrix4d& Rot) const {
        Eigen::Matrix3d U = Eigen::Matrix3d::Zero();
        auto rsize = r1.size();

        for (unsigned int n = 0; n < rsize; n++) {
            U += r1[n] * r2[n].transpose();
        }

        double det_U = U.determinant();
        if (std::abs(det_U) < TMDET_TINY) {
            return false;
        }

        double sign_detU = (det_U > 0) ? 1.0 : -1.0;

        Eigen::MatrixXd omega(6, 6);
        omega.setZero();
        omega.block<3, 3>(0, 3) = U;
        omega.block<3, 3>(3, 0) = U.transpose();

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(omega);
        if (eigensolver.info() != Eigen::Success) {
            return false;
        }

        Eigen::VectorXd eva_omega = eigensolver.eigenvalues();
        Eigen::MatrixXd eve_omega = eigensolver.eigenvectors();

        if (det_U < 0.0 && std::abs(eva_omega(1) - eva_omega(2)) < TMDET_TINY) {
            return false;
        }
        
        Eigen::Matrix3d H;
        Eigen::Matrix3d K;
        for (int i = 0; i < 3; i++) {
            H.col(i) = M_SQRT2 * eve_omega.col(i).head<3>();
            K.col(i) = M_SQRT2 * eve_omega.col(i).tail<3>();
        }

        if ((H.col(1).cross(H.col(2))).dot(H.col(0)) <= 0.0) {
            H.col(2) = -H.col(2);
            K.col(2) = -K.col(2);
        }

        Eigen::Matrix3d R = K * H.transpose();
        R *= sign_detU;
        R *= -1;

        for (unsigned int n = 0; n < rsize; n++) {
            r1[n] = R * r1[n];
        }

        rmsd = 0.0;
        for (unsigned int n = 0; n < rsize; n++) {
            Eigen::Vector3d dr = r1[n] - r2[n];
            rmsd += dr.squaredNorm();
        }
        rmsd /= (double)rsize;
        rmsd = std::sqrt(rmsd);

        Rot.block<3,3>(0,0) = R;
        Rot(3, 3) = 1.0;
        for (int i = 0; i < 3; i++) {
            Rot(i, 3) = 0.0;
            Rot(3, i) = 0.0;
        }
        return true;
    }

    Eigen::Matrix4d Symmetry::rotateZ(Eigen::Vector3d T) const {
        double a;
        double b;
        double c;
        double d;
        double cosa;
        double sina;
        Eigen::Matrix4d R1;
        Eigen::Matrix4d R2;
        Eigen::Matrix4d R;
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

    bool Symmetry::haveSameAxes() const {
        bool first = true;
        _symmetryData sim1;
        for(const auto& s: sim) {
            if (s.good) {
                if (first) {
                    sim1 = s;
                }
                else {
                    if (!sim1.same(s)) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    _symmetryData Symmetry::getAverageAxes() const {
        _symmetryData ret;
        int n=0;
        gemmi::Vec3 f = sim[0].axis;
        for(const auto& s: sim) {
            if (s.good && f.dist(s.axis) < TMDET_TINY) {
                ret.origo += s.origo;
                ret.axis += s.axis;
                n++;
            }
        }
        if (n > 0) {
            ret.origo /= n;
            ret.axis /= n;
            ret.good = true;
        }
        return ret;
    }

    std::vector<gemmi::Vec3> Symmetry::clusterAxes(std::vector<_symmetryData> axes) {
        std::vector<gemmi::Vec3> ret;
        for(unsigned int i = 0; i<axes.size(); i++) {
            if (axes[i].good) {
                gemmi::Vec3 m;
                m = axes[i].axis;
                int k = 1;
                for(unsigned int j = i+1; j<axes.size(); j++) {
                    if (axes[j].good) {
                        axes[j].good = false;
                        m += axes[j].axis;
                        k++;
                    }
                }
                m /= k;
                ret.emplace_back(m);
            }
        }
        return ret;
    }
}
