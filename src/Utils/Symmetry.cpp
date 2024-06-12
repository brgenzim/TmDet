#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <Types/Residue.hpp>
#include <ValueObjects/TmdetStruct.hpp>
#include <Utils/Symmetry.hpp>

#include <iostream>
#include <cstring>
#include <cstdio>

using namespace std;

namespace Tmdet::Utils {

    static int Match_Sequence(char *seq1,char *seq2,int *pos1,int *pos2)
    {
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

    std::vector<std::vector<_symmetryData>> Symmetry::CheckSymmetry(Tmdet::ValueObjects::TmdetStruct &tmdetVO) {
        int i,j,k;
        char *seq1, *seq2;
        gemmi::Vec3 t1, t2;
        // double rmsd, dist, *w2, *weight;
        // matrix4 R;
        std::vector<std::vector<_symmetryData>> sim;


        int nc = tmdetVO.chains.size();
        sim.resize(nc);
        for (int i=0;i<nc;i++) {
            sim[i].resize(nc);
            for (j=0;j<nc;j++) {
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
            //     nall=MIN(nca1-pos1,nca2-pos2);
            //     weight=calloc(nall,sizeof(double));
            //     for (k=0;k<nall;k++) {
            //         if ((w1[k+pos1]!=0) && (w2[k+pos2]!=0))
            //             weight[k]=1;
            //         else
            //             weight[k]=0;
            //     }

            //     CM_Translate(nall,&koord1[pos1], weight, &t1);
            //     CM_Translate(nall,&koord2[pos2], weight, &t2);
            //     Lsq_fit(nall,&koord1[pos1],&koord2[pos2],weight, &rmsd,&R);
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

}
