#include <algorithm>
#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <Types/Residue.hpp>
#include <ValueObjects/TmdetStruct.hpp>
#include <Utils/Symmetry.hpp>

#include <iostream>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <span>

using namespace std;

#define TINY   1.0e-10
#define SMALL  1.0e-5
#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))
#define	ABS(x) ((x) < 0 ? -(x) : (x))
#define RAD2DEG   57.29576

namespace Tmdet::Utils {

    static int Match_Sequence(char *seq1,char *seq2,int *pos1,int *pos2);
    static void CM_Translate(int position, std::vector<Eigen::Vector3d> &koord, std::vector<double> &mass, Eigen::Vector3d &CM);
    static bool Lsq_fit(int nr_atoms, std::span<Eigen::Vector3d>& r1, std::span<Eigen::Vector3d>& r2, std::vector<double>& weight, double& rmsd, Eigen::Matrix4d& Rot);
    static void Get_Sim_Op(Eigen::Matrix4d& R, Eigen::Vector3d& t1, Eigen::Vector3d& t2, _symmetryData& simij);
    static Eigen::Matrix4d Rotate_Z(Eigen::Vector3d T);
    static double cosineAngleOfVectors(Eigen::Vector3d u, gemmi::Vec3 v);

    std::vector<std::vector<_symmetryData>> Symmetry::CheckSymmetry(Tmdet::ValueObjects::TmdetStruct &tmdetVO) {

        char *seq1, *seq2;
        std::vector<std::vector<_symmetryData>> sim;

        int nc = tmdetVO.chains.size();
        sim.resize(nc);
        for (int i=0;i<nc;i++) {
            sim[i].resize(nc);
            for (int j=0;j<nc;j++) {
                sim[i][j].id=0;
            }
        }

        auto ch1 = tmdetVO.chains.begin();
        for (int i = 0; ch1 != tmdetVO.chains.end(); ch1++, i++) {
            if (!ch1->selected) {
                continue;
            }
#ifdef __SYM_DBG
            std::cout << "Working on chain " << ch1->id << std::endl;
#endif
            int nall1 = ch1->residues.size();
            std::vector<double> w1;
            w1.resize(nall1);
            std::vector<Eigen::Vector3d> koord1;
            koord1.resize(nall1);

            seq1 = (char *)calloc((nall1 + 1), sizeof(char));

            int nca1 = 0;
            auto r1 = ch1->residues.begin();
            for (; r1 != ch1->residues.end(); r1++) {
                auto a1 = r1->gemmi.get_ca();
                if (a1 != NULL) {
                    koord1[nca1] = Eigen::Vector3d(a1->pos.x, a1->pos.y, a1->pos.z);
                    w1[nca1] = 1;
                }
                else
                    w1[nca1]=0;
                seq1[nca1] = Tmdet::Types::ResidueType::getResidue(r1->gemmi.name).a1;
                nca1++;
            }
            seq1[nca1] = '\0';
            if (nca1 < 4) {
                free(seq1);
                continue;
            }

            auto ch2 = tmdetVO.chains.begin();
            for (int j = 0; j < i; ch2++, j++) {
                if (!ch2->selected) {
                    continue;
                }
#ifdef __SYM_DBG
                std::cout << "Working on chains " << ch1->id << " " << ch2->id << std::endl;
#endif
                int nall2 = ch2->residues.size();
                std::vector<double> w2;
                w2.resize(nall2);
                std::vector<Eigen::Vector3d> koord2;
                koord2.resize(nall2);
                seq2 = (char *)calloc((nall2+1), sizeof(char));

                int nca2 = 0;
                auto r2 = ch2->residues.begin();
                for (; r2 != ch2->residues.end(); r2++) {
                    auto a2 = r2->gemmi.get_ca();
                    if (a2 != NULL) {
                        koord2[nca2] = Eigen::Vector3d(a2->pos.x, a2->pos.y, a2->pos.z);
                        w2[nca2]=1;
                    }
                    else
                        w2[nca2]=0;
                    seq2[nca2] = Tmdet::Types::ResidueType::getResidue(r2->gemmi.name).a1;
                    nca2++;
                }
                seq2[nca2] = '\0';
                if (nca2 < 4) {
                    free(seq2);
                    continue;
                }
                nca1 = 0;
                r1 = ch1->residues.begin();
                for (; r1 != ch1->residues.end(); r1++) {
                    auto a1 = r1->gemmi.get_ca();
                    if (a1 != NULL) {
                        koord1[nca1] = Eigen::Vector3d(a1->pos.x, a1->pos.y, a1->pos.z);
                        nca1++;
                    }
                }

                int pos1, pos2;
                if (Match_Sequence(seq1,seq2,&pos1,&pos2)==0) {
                    free(seq2);
                    continue;
                }
                int nall = std::min(nca1 - pos1, nca2 - pos2);
                std::vector<double> weight;
                weight.resize(nall);
                for (int k=0;k<nall;k++) {
                    if ((w1[k+pos1]!=0) && (w2[k+pos2]!=0))
                        weight[k]=1;
                    else
                        weight[k]=0;
                }

                Eigen::Vector3d t1, t2;
                CM_Translate(pos1, koord1, weight, t1);
                CM_Translate(pos2, koord2, weight, t2);

                double rmsd;
                Eigen::Matrix4d R; // rotation matrix
                std::span<Eigen::Vector3d> koord1Slice(koord1.begin()+pos1, nall);
                std::span<Eigen::Vector3d> koord2Slice(koord2.begin()+pos2, nall);
                Lsq_fit(nall, koord1Slice, koord2Slice, weight, rmsd, R);
                Get_Sim_Op( R, t1, t2, sim[i][j]);
                double distance = (t2 - t1).norm();

#ifdef __SYM_DBG
                std::cout << ch1->gemmi.name << "-" << ch2->gemmi.name << ":"
                    << " Distance: " << distance << ";    RMSD: " << rmsd << std::endl;
                std::cout << R << std::endl;
#endif

                if (distance > 2.0 && rmsd < 3) {
                    sim[i][j].id = sim[j][i].id = 1;
                }

                if (distance < 2.0 && R(0, 0) > 0.98 && R(1, 1) > 0.98 && R(2, 2) > 0.98) {
                    fprintf(stderr,"Chain %5d and %5d are on top of each other \n", i, j);
                    sim[i][j].id = sim[j][i].id = -1;
                }
                else if (fabs(cosineAngleOfVectors((t1 - t2), sim[i][j].axis)) > 0.15) {
                    fprintf(stdout,"Chain %5d and %5d are not simply rotated\n", i, j);
                    sim[i][j].id = sim[j][i].id = -1;
                }
            }
        }
#ifdef __SYM_DBG
        std::cout << "End of check symmetry" << std::endl;

        for (int i = 0; i < nc; i++) {
            for (int j = 0; j < nc; j++) {
                gemmi::Vec3 axis = sim[i][j].axis;
                gemmi::Vec3 origo = sim[i][j].origo;
                std::cout << "sim[" << i << "][" << j << "]:"
                    << " Origo: " << origo.x << " " << origo.y << " " << origo.z
                    << " Normal: " << axis.x << " " << axis.y << " " << axis.z
                    << " Angle: " << sim[i][j].rotAngle << std::endl;
            }
        }
#endif
        return sim;
    }

    int Match_Sequence(char *seq1,char *seq2,int *pos1,int *pos2) {

        char s1[6], s2[6];
        char *mp;

        sprintf(s1,"%.5s",seq1);s1[5]='\0';
        sprintf(s2,"%.5s",seq2);s2[5]='\0';

        mp=strstr(seq1,s2);
        if (mp!=NULL) *pos1=mp-seq1,*pos2=0;
        else
        {
            mp=strstr(seq2,s1);
            if (mp!=NULL) *pos1=0,*pos2=mp-seq2;
        }
        if (mp!=NULL)
        {
            /*printf(" %5d %5d\n",*pos1,*pos2);*/
        return 1;

        }
        return 0;
    }

    void CM_Translate(int position, std::vector<Eigen::Vector3d> &koord, std::vector<double> &mass, Eigen::Vector3d &CM) {
        /*
         * translates the coordinates of the molecules such
         * that the center of mass is in the origin.
         */
        CM = Eigen::Vector3d::Zero();
        double total_mass = 0.0;
        int length = mass.size();
        for (int i = 0; i < length; i++) if (mass[i]!=0) {
            total_mass += mass[i];
            CM += mass[i] * koord[position + i];
        }
        CM /= total_mass;

        /* ---- translating CM coordinates to origin  -------- */
        for (int i = 0; i < length; i++) if (mass[i]!=0) {
            koord[position + i] -= CM;
        }
    }

    void Get_Sim_Op(Eigen::Matrix4d& R, Eigen::Vector3d& t1, Eigen::Vector3d& t2, _symmetryData& simij)
    {
        Eigen::Matrix4d LR, Rot;
        Eigen::Vector3d av, ori, axis;
        double tg;

        /*Show_Matrix(R,"Rotori");*/

        //   Rinv= inverseMat(R);


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

        /*Show_Matrix(LR,"LR"); */

        Rot= LR.inverse() * (R * LR);

        /*Show_Matrix(Rot,"Rot"); */


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


        /*write_point(simij.ori,"ori");
        printf("%6.2f %6.2f %6.2f    ",
            simij.axis.x,simij.axis.y,simij.axis.z);
        printf("by %10.f %10.f%10.f%10.f    degree\n",
        RAD2DEG*acos(Rot(0, 0)),RAD2DEG*acos(Rot(1, 1)),
        RAD2DEG*asin(Rot(0, 1)),RAD2DEG*asin(-Rot(1, 0)));
        */

    }

    bool Lsq_fit(int nr_atoms, std::span<Eigen::Vector3d>& r1, std::span<Eigen::Vector3d>& r2, std::vector<double>& weight, double& rmsd, Eigen::Matrix4d& Rot) {
        Eigen::Matrix3d U = Eigen::Matrix3d::Zero();

        // Calculate matrix U
        for (int n = 0; n < nr_atoms; n++) {
            if (weight[n] > 0) {
                U += weight[n] * r1[n] * r2[n].transpose();
            }
        }

        double det_U = U.determinant();
        if (std::abs(det_U) < 1e-10) {
            std::cerr << "determinant of U equals to zero" << std::endl;
            return false;
        }

        double sign_detU = (det_U > 0) ? 1.0 : -1.0;

        // Construct omega
        Eigen::MatrixXd omega(6, 6);
        omega.setZero();
        omega.block<3, 3>(0, 3) = U;
        omega.block<3, 3>(3, 0) = U.transpose();

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(omega);
        if (eigensolver.info() != Eigen::Success) {
            std::cerr << "Eigen decomposition failed" << std::endl;
            return false;
        }

        Eigen::VectorXd eva_omega = eigensolver.eigenvalues();
        Eigen::MatrixXd eve_omega = eigensolver.eigenvectors();

        if (det_U < 0.0) {
            if (std::abs(eva_omega(1) - eva_omega(2)) < 1e-10) {
                std::cerr << "determinant of U < 0 && degenerated eigenvalues" << std::endl;
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
        // correction: something differs from the old function in the above code.
        R.transposeInPlace();
        R *= -1;

        // Rotate r1
        for (int n = 0; n < nr_atoms; n++) {
            if (weight[n] > 0) {
                r1[n] = R * r1[n];
            }
        }

        // Calculate RMSD
        rmsd = 0.0;
        for (int n = 0; n < nr_atoms; n++) {
            if (weight[n] > 0) {
                Eigen::Vector3d dr = r1[n] - r2[n];
                rmsd += dr.squaredNorm();
            }
        }
        rmsd /= nr_atoms;
        rmsd = std::sqrt(rmsd);

        // Set rotation matrix
        Rot.block<3,3>(0,0) = R;
        Rot(3, 3) = 1.0;
        for (int i = 0; i < 3; i++) {
            Rot(i, 3) = 0.0;
            Rot(3, i) = 0.0;
        }

        return true;
    }

    Eigen::Matrix4d Rotate_Z(Eigen::Vector3d T) {
        float  a, b, c, d, cosa, sina;
        Eigen::Matrix4d R1, R2, R;
        Eigen::Vector3d normalized;

        normalized = T / T.norm();

        a = normalized.x(); b = normalized.y(); c = normalized.z();

        d = sqrt(b*b + c*c);
        if (d != 0) {
            sina = c/d; cosa = b/d;

            R1(0, 0) = 1;   R1(0, 1) = 0;     R1(0, 2) = 0;    R1(0, 3) = 0;
            R1(1, 0) = 0;   R1(1, 1) = sina;  R1(1, 2) = cosa; R1(1, 3) = 0;
            R1(2, 0) = 0;   R1(2, 1) = -cosa; R1(2, 2) = sina; R1(2, 3) = 0;
            R1(3, 0) = 0;   R1(3, 1) = 0;     R1(3, 2) = 0;    R1(3, 3) = 1;


            R2(0, 0) = d;   R2(0, 1) = 0;     R2(0, 2) = a;    R2(0, 3) = 0;
            R2(1, 0) = 0;   R2(1, 1) = 1;     R2(1, 2) = 0;    R2(1, 3) = 0;
            R2(2, 0) = -a;  R2(2, 1) = 0;     R2(2, 2) = d;    R2(2, 3) = 0;
            R2(3, 0) = 0;   R2(3, 1) = 0;     R2(3, 2) = 0;    R2(3, 3) = 1;

            R = R1 * R2;
        } else {
            R(0, 0) = 0; R(0, 1) = 1; R(0, 2) = 0; R(0, 3) = 0;
            R(1, 0) = 0; R(1, 1) = 0; R(1, 2) = 1; R(1, 3) = 0;
            R(2, 0) = 1; R(2, 1) = 0; R(2, 2) = 0; R(2, 3) = 0;
            R(3, 0) = 0; R(3, 1) = 0; R(3, 2) = 0; R(3, 3) = 1;
        }

        return R;
    }

    double cosineAngleOfVectors(Eigen::Vector3d u, gemmi::Vec3 v) {
        Eigen::Vector3d v_(v.x, v.y, v.z);
        double dotProduct = u.dot(v_);
        return dotProduct / (u.norm() * v_.norm());
    }
}
