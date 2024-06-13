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

using namespace std;

#define TINY   1.0e-10
#define SMALL  1.0e-5
#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))
#define	ABS(x) ((x) < 0 ? -(x) : (x))
typedef struct {
    double el[4][4];
} matrix4;

namespace Tmdet::Utils {

    static int Match_Sequence(char *seq1,char *seq2,int *pos1,int *pos2);
    static void CM_Translate(int nr, std::vector<gemmi::Vec3> &koord, std::vector<double> &mass, gemmi::Vec3 &CM);
    static bool Lsq_fit(int nr_atoms, std::vector<gemmi::Vec3> &r1, std::vector<gemmi::Vec3> &r2, std::vector<double> &weigth, double &rmsd, matrix4 *Rot);
    static void eigen(double *a, double *r__, int *n, int *mv);

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
            w1.reserve(nall1);
            std::vector<gemmi::Vec3> koord1;
            koord1.reserve(nall1);
            seq1 = (char *)calloc((nall1 + 1), sizeof(char));

            int nca1 = 0;
            auto r1 = ch1->residues.begin();
            for (; r1 != ch1->residues.end(); r1++) {
                auto a1 = r1->gemmi.get_ca();
                if (a1 != NULL) {
                    koord1[nca1] = a1->pos;
                    w1[nca1] = 1;
                }
                else
                    w1[nca1]=0;
                seq1[nca1] = Tmdet::Types::ResidueType::getResidue(r1->gemmi.name).a1;
                nca1++;
            }
            seq1[nca1] = '\0';
            if (nca1 < 4) {
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
                w2.reserve(nall2);
                std::vector<gemmi::Vec3> koord2;
                koord2.reserve(nall2);
                seq2 = (char *)calloc((nall2+1), sizeof(char));

                int nca2 = 0;
                auto r2 = ch2->residues.begin();
                for (; r2 != ch2->residues.end(); r2++) {
                    auto a2 = r2->gemmi.get_ca();
                    if (a2 != NULL) {
                        koord2[nca2] = a2->pos;
                        w2[nca2]=1;
                    }
                    else
                        w2[nca2]=0;
                    seq2[nca2] = Tmdet::Types::ResidueType::getResidue(r2->gemmi.name).a1;
                    nca2++;
                }
                seq2[nca2] = '\0';
                if (nca2 < 4) {
                    continue;
                }
                {
                    // TODO: unneeded duplication? removable?
                    // nca1=0;
                    // for (r1=ch1->fres; r1!=NULL; r1=r1->next) {
                    //     a1=pdb_findAtomInRes(r1,"CA");
                    //     if (a1!=NULL) {
                    //         koord1[nca1].x=x(a1->coord);
                    //         koord1[nca1].y=y(a1->coord);
                    //         koord1[nca1].z=z(a1->coord);
                    //     }
                    //     nca1++;
                    // }
                }

                int pos1, pos2;
                if (Match_Sequence(seq1,seq2,&pos1,&pos2)==0) {
                    free(seq2);
                    continue;
                }
                int nall = std::min(nca1 - pos1, nca2 - pos2);
                std::vector<double> weight;
                weight.reserve(nall);
                for (int k=0;k<nall;k++) {
                    if ((w1[k+pos1]!=0) && (w2[k+pos2]!=0))
                        weight[k]=1;
                    else
                        weight[k]=0;
                }

                gemmi::Vec3 t1, t2;

                CM_Translate(pos1, koord1, weight, t1);
                CM_Translate(pos2, koord2, weight, t2);
                double rmsd;
                matrix4 R;
                //Lsq_fit(nall,&koord1[pos1],&koord2[pos2], weight, &rmsd, &R);
            //     Get_Sim_Op( R, t1, t2, &sim[i][j]);
            //     dist = vec_diff_2(t1,t2);

            //     if ((dist>2.0) && (rmsd <3))
            //         sim[i][j].id=sim[j][i].id=1;

            //     if ((dist<2.0)&&((R.el[0][0]>0.98)&&(R.el[1][1]>0.98)&&(R.el[2][2]>0.98))) {
            //         fprintf(stderr,"Chain %5d and %5d are on top of each other \n",i,j);
            //         sim[i][j].id=sim[j][i].id=-1;
            //     }
            //     else if (fabs(vec_cos(vec_sub(t1,t2),sim[i][j].norm))>0.15) {
            //         fprintf(stdout,"Chain %5d and %5d are not simply rotated\n",i,j);
            //         sim[i][j].id=sim[j][i].id=-1;
            //     }
            }
        }
#ifdef __SYM_DBG
        std::cout << "End of check symmetry" << std::endl;
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

    void CM_Translate(int position, std::vector<gemmi::Vec3> &koord, std::vector<double> &mass, gemmi::Vec3 &CM) {
        /*
         * translates the coordinates of the molecules such
         * that the center of mass is in the origin.
         */
        CM.x = 0.0;
        CM.y = 0.0;
        CM.z = 0.0;
        double total_mass = 0.0;
        int length = mass.size();
        for (int i = 0; i < length; i++) if (mass[i]!=0) {
            total_mass += mass[i];
            CM.x += mass[i] * koord[position + i].x;
            CM.y += mass[i] * koord[position + i].y;
            CM.z += mass[i] * koord[position + i].z;
        }
        CM.x /= total_mass;
        CM.y /= total_mass;
        CM.z /= total_mass;

        /* ---- translating CM coordinates to origin  -------- */
        for (int i = 0; i < length; i++) if (mass[i]!=0) {
            koord[position + i].x -= CM.x;
            koord[position + i].y -= CM.y;
            koord[position + i].z -= CM.z;
        }
    }

    bool Lsq_fit(int nr_atoms, std::vector<gemmi::Vec3> &r1, std::vector<gemmi::Vec3> &r2, std::vector<double> &weigth, double &rmsd, matrix4 *Rot) {
        /*
            Least square fit routine to superimpose coordinates r1 onto
            coordinates r2.

            omega is a symmetric matrix in symmetric storage mode:
            The element [i][j] is stored at position [i*(i+1)/2+j] with i>j

        */

        int    i, j, ii, jj, n;
        double  U[3][3], det_U, sign_detU, sigma, dr_sqrlength;
        double  H[3][3], K[3][3], R[3][3];
        double  omega[21], eve_omega[36], eva_omega[6];
        gemmi::Vec3  dr;


        /* ----- CALCULATE THE MATRIX U AND ITS DETERMINANT ----- */
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                U[i][j] = 0.0;

        for (n = 0; n < nr_atoms; n++)  if (weigth[n]>0)   
        {
            U[0][0] +=  weigth[n] * r1[n].x * r2[n].x;
            U[0][1] +=  weigth[n] * r1[n].x * r2[n].y;
            U[0][2] +=  weigth[n] * r1[n].x * r2[n].z;
            U[1][0] +=  weigth[n] * r1[n].y * r2[n].x;
            U[1][1] +=  weigth[n] * r1[n].y * r2[n].y;
            U[1][2] +=  weigth[n] * r1[n].y * r2[n].z;
            U[2][0] +=  weigth[n] * r1[n].z * r2[n].x;
            U[2][1] +=  weigth[n] * r1[n].z * r2[n].y;
            U[2][2] +=  weigth[n] * r1[n].z * r2[n].z;
        }

        det_U = U[0][0]*U[1][1]*U[2][2] + U[0][2]*U[1][0]*U[2][1] +
                U[0][1]*U[1][2]*U[2][0] - U[2][0]*U[1][1]*U[0][2] -
                U[2][2]*U[1][0]*U[0][1] - U[2][1]*U[1][2]*U[0][0];

        if (ABS(det_U) < TINY) {
            fprintf(stderr, "determinant of U equals to zero\n");
            return false;
        }

        sign_detU = det_U / ABS(det_U);  /* sign !!! */


    /* ----- CONSTRUCT OMEGA, DIAGONALIZE IT AND DETERMINE H AND K --- */

        for (i = 0; i < 6; i++)
            for (j = i; j < 6; j++)
                omega[(j*(j+1)/2)+i] = 0.0;
        for (j = 3; j < 6; j++) {
            jj = j*(j+1)/2;
            for (i = 0; i < 3; i++) {
                ii = jj + i;
                omega[ii] = U[i][j-3];
            }
        }

    #ifdef __SYM_DBG
        fprintf(stdout, "omega matrix:\n");
        for (i = 0; i < 21; i++)
            fprintf(stdout, "%7.2f ", omega[i]);
        fprintf(stdout, "\n");
    #endif

        i = 6;  /* dimension of omega matrix */
        j = 0;  /* both, eigenvalues and eigenvectors are calculated */
        eigen(omega, eve_omega, &i, &j);

        for (i = 0; i < 6; i++)
            eva_omega[i] = omega[i*(i+1)/2+i];

    #ifdef __SYM_DBG
        fprintf(stdout, "Eigenvalues:\n");
        for (i = 0; i < 6; i++)
            fprintf(stdout, "%7.2f ", eva_omega[i]);
        fprintf(stdout, "\n");

        ii = 0;
        fprintf(stdout, "Eigenvectors:\n");
        for (j = 0; j < 6; j++) {     /* ----- elements of eigenvector i ------ */
            for (i = 0; i < 6; i++)   /* ----  loop for eigenvectors !!! ----- */
                fprintf(stdout, "%7.2f ", eve_omega[ii++]);
            fprintf(stdout, "\n");
        }

    #endif


        if (det_U < 0.0){
            if (ABS(eva_omega[1] - eva_omega[2]) < SMALL) {
                fprintf(stderr, "determinant of U < 0 && degenerated eigenvalues");
                return false;
            }
        }

        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++) {
                H[i][j] = M_SQRT2 * eve_omega[j*6+i];
                K[i][j] = M_SQRT2 * eve_omega[j*6+i+3];
            }

        sigma = (H[1][0]*H[2][1] - H[2][0]*H[1][1]) * H[0][2] +
                (H[2][0]*H[0][1] - H[0][0]*H[2][1]) * H[1][2] +
                (H[0][0]*H[1][1] - H[1][0]*H[0][1]) * H[2][2];

        if (sigma <= 0.0) {
            for (i = 0; i < 3; i++) {
                H[i][2] = -H[i][2];
                K[i][2] = -K[i][2];
            }
        }

        /* --------- DETERMINE R AND ROTATE X ----------- */
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                R[j][i] = K[j][0]*H[i][0] + K[j][1]*H[i][1] +
                        sign_detU * K[j][2]*H[i][2];


    #ifdef __SYM_DBG
        fprintf(stdout, "Rotation matrix:\n");
        for (j = 0; j < 3; j++) {
            for (i = 0; i < 3; i++)
                fprintf(stdout, "%10.4f", R[i][j]);
            fprintf(stdout, "\n");
        }
    #endif


        for (n = 0; n < nr_atoms; n++) if (weigth[n]!=0) 
        {
            dr.x = R[0][0]*r1[n].x + R[0][1]*r1[n].y + R[0][2]*r1[n].z;
            dr.y = R[1][0]*r1[n].x + R[1][1]*r1[n].y + R[1][2]*r1[n].z;
            dr.z = R[2][0]*r1[n].x + R[2][1]*r1[n].y + R[2][2]*r1[n].z;
            r1[n].x = dr.x;
            r1[n].y = dr.y;
            r1[n].z = dr.z;
        }


        /* ----- calculate RMSD when required -------- */

        rmsd = 0.0;
        for (n = 0; n < nr_atoms; n++) if (weigth[n]!=0) {
        /*printf("%10.4f%10.4f", r1[n].x , r2[n].x);
        printf("%10.4f%10.4f", r1[n].y , r2[n].y);
        printf("%10.4f%10.4f\n", r1[n].z , r2[n].z);*/
            dr.x = r1[n].x - r2[n].x;
            dr.y = r1[n].y - r2[n].y;
            dr.z = r1[n].z - r2[n].z;
            dr_sqrlength = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;

            rmsd +=  dr_sqrlength;
        }
        rmsd /= nr_atoms;
        rmsd = sqrt(rmsd);


        Rot->el[0][0]=R[0][0];Rot->el[0][1]=R[0][1];Rot->el[0][2]=R[0][2];
        Rot->el[1][0]=R[1][0];Rot->el[1][1]=R[1][1];Rot->el[1][2]=R[1][2];
        Rot->el[2][0]=R[2][0];Rot->el[2][1]=R[2][1];Rot->el[2][2]=R[2][2];
        Rot->el[3][0]=0;Rot->el[3][1]=0;Rot->el[3][2]=0;Rot->el[3][3]=1.0;
        Rot->el[0][3]=0;Rot->el[1][3]=0;Rot->el[2][3]=0;

        return true;
    }

    void eigen(double *a, double *r__, int *n, int *mv) {
        int    i__1, i__2, i__3;
        double  d__1;

        double  cosx, sinx, cosx2, sinx2;
        int    i__, j, k, l, m;
        double  x, y, range, anorm, sincs, anrmx;
        int    ia, ij, il, im, ll, lm, iq, mm, jq, lq, mq, ind, ilq, imq, ilr, imr;
        double  thr;


    /* CCCCC W.F. VAN GUNSTEREN, CAMBRIDGE, JUNE 1979 CCCCCCCCCCCCCCCCCCCCCCCC
                                                                        C 
        SUBROUTINE EIGEN (A,R,N,MV)                                      C 
                                                                        C 
            EIGEN COMPUTES EIGENVALUES AND EIGENVECTORS OF THE REAL      C
        SYMMETRIC N*N MATRIX A, USING THE DIAGONALIZATION METHOD         C 
        DESCRIBED IN "MATHEMATICAL METHODS FOR DIGITAL COMPUTERS", EDS.  C 
        A.RALSTON AND H.S.WILF, WILEY, NEW YORK, 1962, CHAPTER 7.        C 
        IT HAS BEEN COPIED FROM THE IBM SCIENTIFIC SUBROUTINE PACKAGE.   C 
                                                                        C 
        A(1..N*(N+1)/2) = MATRIX TO BE DIAGONALIZED, STORED IN SYMMETRIC C 
                        STORAGE MODE, VIZ. THE I,J-TH ELEMENT (I.GE.J) C 
                        IS STORED AT THE LOCATION K=I*(I-1)/2+J IN A;  C 
                        THE EIGENVALUES ARE DELIVERED IN DESCENDING    C 
                        ORDER ON THE DIAGONAL, VIZ. AT THE LOCATIONS   C 
                        K=I*(I+1)/2                                    C 
        R(1..N,1..N) = DELIVERED WITH THE CORRESPONDING EIGENVECTORS     C 
                        STORED COLUMNWISE                                 C 
        N = ORDER OF MATRICES A AND R                                    C 
        MV = 0 : EIGENVALUES AND EIGENVECTORS ARE COMPUTED               C 
            = 1 : ONLY EIGENVALUES ARE COMPUTED                           C 
                                                                        C 
    CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    */

        /* Parameter adjustments */
        --r__;
        --a;



        range = 1e-12;
        if (*mv != 1) {
            iq = -(*n);
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                iq += *n;
                i__2 = *n;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    ij = iq + i__;
                    r__[ij] = 0.;
                    if (i__ - j == 0)
                        r__[ij] = 1.;
                }
            }
        }

    /* *****COMPUTE INITIAL AND FINAL NORMS (ANORM AND ANRMX) */
        anorm = 0.0;
        i__2 = *n;
        for (i__ = 1; i__ <= i__2; ++i__) {
        i__1 = *n;
        for (j = i__; j <= i__1; ++j) {
            if (i__ - j != 0) {
                    ia = i__ + (j * j - j) / 2;
                    anorm += a[ia] * a[ia];
                }
        }
        }
        if (anorm > 0.0) {
            anorm = sqrt(anorm) * 1.414;
            anrmx = anorm * range / (double) (*n);

    /* *****INITIALIZE INDICATORS AND COMPUTE THRESHOLD, THR */
            thr = anorm;

            do {  /* ----- while (thr - anrmx > 0.0) ------ */
                thr /= (double) (*n);

                ind = true;
                while (ind) {
                    ind = false;

                    l = 1;
                    for (l = 1; l <= *n-1; l++) {
                        for (m = l+1; m <= *n; m++) {

    /* *****COMPUT SIN AND COS */

                            mq = (m * m - m) / 2;
                            lq = (l * l - l) / 2;
                            lm = l + mq;
                            d__1 = a[lm];
                            if (ABS(d__1) - thr >= 0.0) {
                                ind = true;
                                ll = l + lq;
                                mm = m + mq;
                                x = (a[ll] - a[mm]) * .5;
                                y = -a[lm] / sqrt(a[lm] * a[lm] + x * x);
                                if (x < 0.0)
                                    y = -y;
                                sinx = y / sqrt((sqrt(1.0 - y * y) + 1.0) * 2.0);
                                sinx2 = sinx * sinx;
                                cosx = sqrt(1.0 - sinx2);
                                cosx2 = cosx * cosx;
                                sincs = sinx * cosx;

    /* *****ROTATE L AND M COLUMNS */
                                ilq = *n * (l - 1);
                                imq = *n * (m - 1);
                                i__1 = *n;
                                for (i__ = 1; i__ <= i__1; ++i__) {
                                    iq = (i__ * i__ - i__) / 2;
                                    if (i__ - l != 0) {

                                        i__2 = i__ - m;
                                        if (i__2 != 0) {
                                            if (i__2 < 0)
                                                im = i__ + mq;
                                            else
                                                im = m + iq;

                                            if (i__ - l >= 0)
                                                il = l + iq;
                                            else
                                                il = i__ + lq;
                                            x = a[il] * cosx - a[im] * sinx;
                                            a[im] = a[il] * sinx + a[im] * cosx;
                                            a[il] = x;
                                        } /* ---- (i__2 != 0) ---- */
                                    } /* ------ if (i__ - l != 0) ---- */

                                    if (*mv != 1) {
                                        ilr = ilq + i__;
                                        imr = imq + i__;
                                        x = r__[ilr] * cosx - r__[imr] * sinx;
                                        r__[imr] = r__[ilr] * sinx + r__[imr] * cosx;
                                        r__[ilr] = x;
                                    }
                                }
                                x = a[lm] * 2. * sincs;
                                y = a[ll] * cosx2 + a[mm] * sinx2 - x;
                                x = a[ll] * sinx2 + a[mm] * cosx2 + x;
                                a[lm] = (a[ll] - a[mm]) * sincs + a[lm] * (cosx2 - sinx2);
                                a[ll] = y;
                                a[mm] = x;

                            } /* --- if ((d__1 = a[lm], ABS(d__1)) - thr >= 0.0) -----
    */
    /* *****TESTS FOR COMPLETION */

    /* *****TEST FOR M = LAST COLUMN */
                        } /* ---- for (m = l+1; m <= *n; m++)  ---- */

    /* *****TEST FOR L = SECOND FROM LAST COLUMN */
                    } /* ----- for (l = 1; l < *n-1; l++) ------ */

                } /* --- while (ind) --- */
    /* *****COMPARE THRESHOLD WITH FINAL NORM */

            } while (thr > anrmx);


        } /* ---- if (anorm > 0) ------ */

    /* *****SORT EIGENVALUES AND EIGENVECTORS */

        iq = -(*n);
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
        iq += *n;
        ll = i__ + (i__ * i__ - i__) / 2;
        jq = *n * (i__ - 2);
        i__2 = *n;
        for (j = i__; j <= i__2; ++j) {
            jq += *n;
            mm = j + (j * j - j) / 2;
            if (a[ll] - a[mm] < 0.0) {
                    x = a[ll];
                    a[ll] = a[mm];
                    a[mm] = x;
                    if (*mv != 1) {
                        i__3 = *n;
                        for (k = 1; k <= i__3; ++k) {
                            ilr = iq + k;
                            imr = jq + k;
                            x = r__[ilr];
                            r__[ilr] = r__[imr];
                            r__[imr] = x;
                        }
                    }
                } /* ---- if (a[ll] - a[mm] < 0.0) ---- */
        }
        }

    } /* eigen */

    bool Lsq_fit_Eigen_version(int nr_atoms, std::vector<Eigen::Vector3d> &r1, std::vector<Eigen::Vector3d> &r2, std::vector<double> &weight, double &rmsd, matrix4 *Rot) {
        Eigen::Matrix3d U = Eigen::Matrix3d::Zero();
        double total_weight = 0.0;

        // Calculate matrix U
        for (int n = 0; n < nr_atoms; n++) {
            if (weight[n] > 0) {
                U += weight[n] * r1[n] * r2[n].transpose();
                total_weight += weight[n];
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

        // Diagonalize omega
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
                rmsd += weight[n] * dr.squaredNorm();
            }
        }
        rmsd /= total_weight;
        rmsd = std::sqrt(rmsd);

        // Set rotation matrix
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                Rot->el[i][j] = R(i, j);
            }
        }
        Rot->el[3][3] = 1.0;
        for (int i = 0; i < 3; i++) {
            Rot->el[i][3] = 0.0;
            Rot->el[3][i] = 0.0;
        }

        return true;
    }

}
