#include "tmall.h"
#include "tmdet_lib.h"
#include "optim.h"

#define CH_CONT 100
#define DIST_CUTOFF  12.0
#define AAN 20

extern int FORCE_NO_DEL_CHAIN;

double Chain_Sort(pdbProtein p, float width, POINT *axis, POINT *memp, double limit)
{
int **ch_clust, cluster, max_clust;
int i, *nr_stAA;
pdbChain ch1,ch2;
pdbRes res;
SIM_OP **sim;
double q;

    sim = gt_Check_Symmetry(p);
	ch_clust = gt_Cluster_Chain(p,sim);
	nr_stAA=calloc(p->nc,sizeof(int));
	for (i=0; i<p->nc; i++) nr_stAA[i]=0;
	max_clust=-1;
	cluster=-1;

	for (ch1=p->fch,i=0; ch1!=NULL; ch1=ch1->next,i++) if (ch1->sel)
	{	for (res=ch1->fres; res!=NULL; res=res->next) if (res->sel)
			if (pdb_isStdAminoAcid(p,ARES(res))) nr_stAA[i]++;
	}

	for (i=0; i<p->nc; i++)
	{	if ((ch_clust[i][i]>max_clust) && (nr_stAA[i]>0))
		{	max_clust=ch_clust[i][i];
			cluster=i;
		}
	}

	if (cluster!=-1)
	if (!FORCE_NO_DEL_CHAIN)
	for (ch1=p->fch,i=0; ch1!=NULL; ch1=ch2,i++)
	{	ch2=ch1->next;
		if (ch1->sel)
		{	if ((ch1->nr>0)&&(i!=cluster)&&(ch_clust[cluster][i]!=cluster))
			{	fprintf(stderr,"Deleting chain %5d %5d %5d\n",i,ch1->nr,ch1->nm);
				//pdb_deleteChain(p,ch1);
				ch1->sel=0;
				if (p->bioMatrix==NULL)
					p->bioMatrix=pdb_newBioMatrix();
				pdb_unselBioMatrix(p->bioMatrix,ch1->id);
			}
			else
			{	//fprintf(stderr,"Selected chain %5d %5d %5d\n",i,ch1->nr,ch1->nm);
			}
		}
	}

	for (i=0; i<p->nc; i++) free(ch_clust[i]);
	free(ch_clust);

	q=Test_Sym_Axis(p,sim,width,axis,memp,limit);

	for (i=0; i<p->nc; i++) free(sim[i]);
	free(sim);
	return q;
}

int **gt_Chain_Contact(pdbProtein p)
{
int i,j;
pdbAtom a1,a2;
pdbChain ch1,ch2;
int **contact;

	contact = IMatrix(p->nc,p->nc);
	for (i=0; i<p->nc; i++)
	{	for (j=0; j<p->nc; j++)
			contact[i][j]=0;
		contact[i][i]=1;
	}

	for (ch1=p->fch,i=0; ch1!=NULL; ch1=ch1->next,i++) if (ch1->sel)
	for (ch2=p->fch,j=0; j<i; ch2=ch2->next,j++) if (ch2->sel)
	{	for (a1=ch1->fa; a1!=NULL; a1=a1->next) if (a1->sel)
		if (!strcmp(ANAME(a1),"CA"))
		for (a2=ch2->fa; a2!=NULL; a2=a2->next) if (a2->sel)
		if (!strcmp(ANAME(a2),"CA"))
		{	if (vec_diff_2(a1->coord,a2->coord) < DIST_CUTOFF*DIST_CUTOFF)
			{	contact[i][j]++;
				contact[j][i]++;
			}
		}
	}
	return contact;
}


int **gt_Cluster_Chain(pdbProtein p, SIM_OP **sim)
{
int i,j;
int nc;
int **ch_cont, **id_clust, *member, *head;
int max, p1, p2;
int stop ;

    nc = p->nc;

    id_clust = IMatrix(nc,nc);
    member = malloc(nc*sizeof(int));
    head = malloc(nc*sizeof(int));
    for (i=0;i<nc;i++) head[i]=member[i]=0;



    ch_cont=gt_Chain_Contact(p);

    for (i=0;i<nc;i++) {
		ch_cont[i][i] = 1;
		head[i] = 1;
		for (j=0;j<i;j++) {
			if (sim[i][j].id == -1) {
				head[i]=0;
			}
			id_clust[i][j] = id_clust[j][i] = sim[i][j].id;
			sim[i][j].cont = sim[j][i].cont = ch_cont[i][j];
		}
		member[i]=i;
    }



    stop=0;
    while (stop==0)
    {


	 max=0;p1=p2=0;

	 for (i=0;i<nc;i++) if (head[i]>0)
	 for (j=0;j<nc;j++) if ((head[j]>0) && (i!=j))
	 {
	      if ((max<ch_cont[i][j]) && (sim[i][j].id!=-1))
	      {
	          max=ch_cont[i][j];
		  p1=i;p2=j;
	      }

	 }

	 /*printf("max %5d i%5d j%5d %5d\n",max,p1,p2,id_clust[p1][p2]);*/

     	 if (((max<CH_CONT)&&((id_clust[p1][p2])==1)) || (max<10))
	     {stop=1; break;}

	 for (i=0;i<nc;i++) if (head[i]>0)
	 {
	      ch_cont[p1][i]+=ch_cont[p2][i];
	      ch_cont[p2][i]=0;
	      ch_cont[i][p1]+=ch_cont[i][p2];
	      ch_cont[i][p2]=0;
	 }
    	 ch_cont[p1][p2]=ch_cont[p2][p1]=0;
	 head[p1]+=head[p2];
	 head[p2]=0;
	 for (i=0;i<nc;i++) if (member[i]==p2) member[i]=p1;


    }

    for (i=0;i<nc;i++)
    {
        for (j=0;j<i;j++)
        {
             if (member[i]==member[j]) ch_cont[i][j]=ch_cont[j][i]=member[j];
	     else ch_cont[i][j]=ch_cont[j][i]=-1;

        }
    }
    for (i=0;i<nc;i++)
    {
        for (j=0;j<nc;j++)  if (member[j]==i) sim[j][i].cl=sim[i][j].cl=i;
	sim[i][i].cl=head[i];
    }

    return (ch_cont);

}

// Some comments ChatGPT-generated
int **gt_Cluster_Chain(pdbProtein p, SIM_OP **sim)
{
    int i, j;
    int nc;
    int **ch_cont, **id_clust, *member, *head;
    int max, p1, p2;
    int stop;

    nc = p->nc; // A láncok száma

    id_clust = IMatrix(nc, nc); // Integer mátrix létrehozása
    member = malloc(nc * sizeof(int)); // Memória foglalás
    head = malloc(nc * sizeof(int)); // Memória foglalás
    for (i = 0; i < nc; i++) head[i] = member[i] = 0; // Inicializálás

    ch_cont = gt_Chain_Contact(p); // Kontakt mátrix létrehozása

    for (i = 0; i < nc; i++) {
        ch_cont[i][i] = 1; // Átló elemek beállítása 1-re
        head[i] = 1; // Fej beállítása - na jó, de mi az a fej?
        for (j = 0; j < i; j++) {
            // ha a sim.id -1, akkor nincs elfogadható transzformáció
            // a két láncra
            if (sim[i][j].id == -1) {
                head[i] = 0;
            }
            id_clust[i][j] = id_clust[j][i] = sim[i][j].id;
            // kontaktusok száma két lánc között
            sim[i][j].cont = sim[j][i].cont = ch_cont[i][j];
        }
        member[i] = i; // Láncok inicializálása
    }

    stop = 0;
    while (stop == 0) {
        max = 0; p1 = p2 = 0;

        // Maximum kontaktok keresése, azokra a lánc-párokra, ahol
        // van egymásba transzformálhatóság symop által.
        for (i = 0; i < nc; i++)
            if (head[i] > 0)
                for (j = 0; j < nc; j++)
                    if ((head[j] > 0) && (i != j)) {
                        if ((max < ch_cont[i][j]) && (sim[i][j].id != -1)) {
                            max = ch_cont[i][j];
                            p1 = i; p2 = j;
                        }
                    }

        // Stop feltétel ellenőrzése
        if (((max < CH_CONT) && ((id_clust[p1][p2]) == 1)) || (max < 10)) {
            stop = 1;
            break;
        }

        // Kontaktusok frissítése
        for (i = 0; i < nc; i++)
            if (head[i] > 0) {
                ch_cont[p1][i] += ch_cont[p2][i];
                ch_cont[p2][i] = 0;
                ch_cont[i][p1] += ch_cont[i][p2];
                ch_cont[i][p2] = 0;
            }
        ch_cont[p1][p2] = ch_cont[p2][p1] = 0;
        head[p1] += head[p2];
        head[p2] = 0;
        for (i = 0; i < nc; i++) {
            if (member[i] == p2) {
                member[i] = p1;
            }
        }
    }

    // Klaszterek beállítása a kontakt mátrixban
    for (i = 0; i < nc; i++) {
        for (j = 0; j < i; j++) {
            if (member[i] == member[j]) {
                ch_cont[i][j] = ch_cont[j][i] = member[j];
            } else {
                ch_cont[i][j] = ch_cont[j][i] = -1;
            }
        }
    }
    for (i = 0; i < nc; i++) {
        for (j = 0; j < nc; j++) {
            // ha tagja az i. lánchoz tartozó cluster-nek
            if (member[j] == i)
                // akkor beállítjuk, melyik lánchoz tartozik
                sim[j][i].cl = sim[i][j].cl = i;
        }
        sim[i][i].cl = head[i];
    }

    return (ch_cont);
}




int Match_Sequence(char *seq1,char *seq2,int *pos1,int *pos2)
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


void Get_Sim_Op(matrix4 R, POINT t1, POINT t2, SIM_OP *simij)
{
matrix4  LR, Rot;
POINT  av, ori;
double tg;



	/*Show_Matrix(R,"Rotori");*/

     //   Rinv= inverseMat(R);


	if ((R.el[0][0]<0.999)&&(R.el[1][1]<0.999) /*&&(R.el[0][0]>-0.999)*/)
	{
	    /* Find vector which is invariant of the rotation
	       lets fix its size by setting z to 1.0
	    */
	     simij->norm.z=1.0;
	     simij->norm.y=
		((1.0-R.el[0][0])*R.el[1][2]+R.el[1][0]*R.el[0][2])/
		((1.0-R.el[1][1])*(1.0-R.el[0][0])-R.el[0][1]*R.el[1][0]);
	     if (R.el[0][0]<1.0)
	         simij->norm.x=(R.el[0][1]*simij->norm.y+R.el[0][2])/
		      (1.0-R.el[0][0]);
	     else simij->norm.x=0;
	 }
 	 else
	 {
	     if (R.el[0][0]>0.999)
	        simij->norm.x=1;simij->norm.y=0;simij->norm.z=0;
	     if (R.el[1][1]>0.999)
	        simij->norm.x=0;simij->norm.y=1;simij->norm.z=0;

	 }
	 simij->norm=vec_norm(simij->norm,1.0);
	 LR=Rotate_Z(simij->norm);

	 /*Show_Matrix(LR,"LR"); */

	 Rot=matrix_multiply(inverseMat(LR),matrix_multiply(R,(LR)));

	 /*Show_Matrix(Rot,"Rot"); */


	if (Rot.el[0][0]>1) Rot.el[0][0]=1.0;
	if (Rot.el[0][0]<-1) Rot.el[0][0]=-1.0;

	if (Rot.el[0][0]!=1.0)
	    tg=Rot.el[0][1]/(1-Rot.el[0][0]);
        else tg=1.0;

	av= vec_sub(midpoint(t1,t2),t1);
	ori=vec_prd(av,simij->norm);
	if (tg!=0) ori=vec_scl(ori,tg);

	simij->ori=vec_add(midpoint(t1,t2),ori);

	simij->angle = RAD2DEG*acos(Rot.el[0][0]);



	/*write_point(simij->ori,"ori");
	printf("%6.2f %6.2f %6.2f    ",
	    simij->norm.x,simij->norm.y,simij->norm.z);
	printf("by %10.f %10.f%10.f%10.f    degree\n",
	   RAD2DEG*acos(Rot.el[0][0]),RAD2DEG*acos(Rot.el[1][1]),
	   RAD2DEG*asin(Rot.el[0][1]),RAD2DEG*asin(-Rot.el[1][0]));
	*/

}

double Test_Sym_Axis(pdbProtein p, SIM_OP **sim, float width, POINT *memaxis, POINT *memmidp, double limit) {
int i, j, k, nc;
POINT mp, axis;
float HPPRO, HPBEST;
pdbChain ch,chn;
int ndel, *delch, *nr_stAA;
double shift;
int sym_tested;
int max_clust, cluster, cont, stop, max;
pdbChain ch1;
pdbRes res;

	ndel=0;
	nc=p->nc;
	delch=calloc(nc,sizeof(int));

	nr_stAA=calloc(nc,sizeof(int));
	for (i=0; i<nc; i++) nr_stAA[i]=0;
	for (ch1=p->fch,i=0;ch1!=NULL;ch1=ch1->next,i++) if (ch1->sel)
	for (res=ch1->fres;res!=NULL;res=res->next)
		if (pdb_isStdAminoAcid(p,ARES(res))!=0) nr_stAA[i]++;

	HPBEST=0;
	sym_tested=0;

	optim_init(p);
	for (k=0;k<nc;k++) if (sim[k][k].cl!=0)
	{	for (i=0;i<nc;i++) if ((sim[k][i].cl==k)||(i==k))
		for (j=0;j<i;j++)  if ((sim[k][j].cl==k)||(j==k))
		{	if (sim[i][j].id==1)
			{	sym_tested++;
				optim_setDistances(p,sim[i][j].norm,sim[i][j].ori);
				shift=0.0;
				HPPRO=100*optim_findBestSlice(p,width,&shift);
				if (HPPRO>HPBEST)
				{	axis.x=sim[i][j].norm.x;
					axis.y=sim[i][j].norm.y;
					axis.z=sim[i][j].norm.z;
					axis=vec_norm(axis,1.0);
					shift+=vec_dot(axis,sim[i][j].ori);
					mp.x=axis.x*shift;
					mp.y=axis.y*shift;
					mp.z=axis.z*shift;
					HPBEST=HPPRO;
				}
			}
		}
	}
	optim_end(p);

	if (sym_tested==0)
		fprintf(stderr,"No internal symmetry\n");
	else
	{	if (HPBEST>limit)
		{	fprintf(stderr,"Symetry axis is likely membrane axis\n");
			memaxis->x=axis.x;
			memaxis->y=axis.y;
			memaxis->z=axis.z;
			memmidp->x=mp.x;
			memmidp->y=mp.y;
			memmidp->z=mp.z;
		}
		else
		{	fprintf(stderr,"Symetry axis is NOT membrane axis\n");
	//		fprintf(stderr,"Keeping a single copy of each chain \n");
			max_clust=0;
			cluster=-1;
			for (i=0;i<nc;i++)
			{	cont=0;
				for (j=0;j<nc;j++) if (sim[i][j].id!=1)
				{	if ((sim[i][j].cont>0) || (i==j))
						cont+=nr_stAA[j];
				}
				if (cont>max_clust)
				{	max_clust=cont;
					cluster=i;
				}
			}
			for (i=0;i<nc;i++) delch[i]=-1;
			delch[cluster]=0;
			stop=0;
			while (!stop)
			{	max_clust=-1;
				max=-1;
				for (i=0;i<nc;i++) if (delch[i]==-1)
				{	for (j=0;j<nc;j++) if (delch[j]==0)
					{	if (sim[i][j].id==1)
						delch[i]=1;
						break;
					}
					if (delch[i]==1) continue;
					for (j=0,cont=0;j<nc;j++) if (delch[j]==0)
						cont+=sim[i][j].cont;
					if (cont>max_clust)
					{	max_clust=cont;
						max=i;
					}
				}
				if (max==-1) {stop=1;break;}
				delch[max]=0;
			}

/*			for (j=0;j<nc;j++)
			{	if (((sim[cluster][j].cont>0)&&(sim[cluster][j].id!=1))||(j==cluster))
					delch[j]=0;
				else
					delch[j]=1;
			}	*/
		}
	}

	if (!FORCE_NO_DEL_CHAIN)
	for (ch=p->fch,i=0; ch!=NULL; ch=chn,i++)
	{	chn=ch->next;
		if  (delch[i]==1)
		{	fprintf(stderr,"Deleting chain %5d %5d %5d\n",i,ch->nr,ch->nm);
			//	pdb_deleteChain(p,ch);
			ch->sel=0;
			ndel++;
			if (p->bioMatrix==NULL)
				p->bioMatrix=pdb_newBioMatrix();
			pdb_unselBioMatrix(p->bioMatrix,ch->id);
		}
	}
	if (ndel>0)
	{	pdb_surface(p);
		pdb_setOutside(p);
	}
	return HPBEST;
}

SIM_OP **gt_Check_Symmetry(pdbProtein p)
{
int i,j,k;
char *seq1, *seq2;
POINT *koord1, *koord2,  t1, t2;
int nc, nall, nca1, nca2, nall1, nall2;
double rmsd, dist, *w1, *w2, *weight;
int pos1, pos2;
pdbChain ch1,ch2;
pdbAtom a1,a2;
pdbRes r1,r2;
matrix4 R;
SIM_OP **sim;


	nc=p->nc;
	sim=calloc(nc,sizeof(SIM_OP *));
	for (i=0;i<nc;i++)
	{	sim[i]=calloc(nc,sizeof(SIM_OP));
		for (j=0;j<nc;j++) sim[i][j].id=0;
	}


     for (ch1=p->fch,i=0; ch1!=NULL; ch1=ch1->next,i++) if (ch1->sel) {

		printf("Working on chain %c\n",ch1->id);
		nall1=ch1->nm + ch1->nr;
		w1=calloc((nall1),sizeof(double));
		koord1=malloc((nall1)*sizeof(POINT));
		seq1=calloc((nall1+1),sizeof(char));
		nca1=0;

		for (r1=ch1->fres; r1!=NULL; r1=r1->next) {
			a1=pdb_findAtomInRes(r1,"CA");
			if (a1!=NULL) {
				koord1[nca1].x=x(a1->coord);
				koord1[nca1].y=y(a1->coord);
				koord1[nca1].z=z(a1->coord);
				w1[nca1]=1;
			}
			else
				w1[nca1]=0;
			seq1[nca1]=RA1(r1);
			nca1++;
		}
		seq1[nca1]='\0';
		if (nca1<4)
			continue;

		for (ch2=p->fch,j=0; j<i; ch2=ch2->next,j++) if (ch2->sel) {
			printf("Working on chains %c %c\n",ch1->id,ch2->id);
			nall2=ch2->nm + ch2->nr;
			koord2=malloc((nall2)*sizeof(POINT));
			seq2=malloc(sizeof(char)*(nall2+1));
			w2=calloc((nall2),sizeof(double));
			nca2=0;

			for (r2=ch2->fres; r2!=NULL; r2=r2->next) {
				a2=pdb_findAtomInRes(r2,"CA");
				seq2[nca2]=RA1(r2);
				if (a2!=NULL) {
					koord2[nca2].x=x(a2->coord);
					koord2[nca2].y=y(a2->coord);
					koord2[nca2].z=z(a2->coord);
					w2[nca2]=1;
				}
				else
					w2[nca2]=0;
				nca2++;
			}
			seq2[nca2]='\0';
			if (nca2<4)
				continue;

		nca1=0;
		for (r1=ch1->fres; r1!=NULL; r1=r1->next) {
			a1=pdb_findAtomInRes(r1,"CA");
			if (a1!=NULL) {
				koord1[nca1].x=x(a1->coord);
				koord1[nca1].y=y(a1->coord);
				koord1[nca1].z=z(a1->coord);
			}
			nca1++;
		}

			if (Match_Sequence(seq1,seq2,&pos1,&pos2)==0) {
				free(seq2);
				free(koord2);
				continue;
			}
			nall=MIN(nca1-pos1,nca2-pos2);
			weight=calloc(nall,sizeof(double));
			for (k=0;k<nall;k++) {
				if ((w1[k+pos1]!=0) && (w2[k+pos2]!=0))
					weight[k]=1;
				else
					weight[k]=0;
			}

			CM_Translate(nall,&koord1[pos1], weight, &t1);
			CM_Translate(nall,&koord2[pos2], weight, &t2);
			Lsq_fit(nall,&koord1[pos1],&koord2[pos2],weight, &rmsd,&R);
			Get_Sim_Op( R, t1, t2, &sim[i][j]);
			dist = vec_diff_2(t1,t2);

			if ((dist>2.0) && (rmsd <3))
				sim[i][j].id=sim[j][i].id=1;

			if ((dist<2.0)&&((R.el[0][0]>0.98)&&(R.el[1][1]>0.98)&&(R.el[2][2]>0.98))) {
				fprintf(stderr,"Chain %5d and %5d are on top of each other \n",i,j);
				sim[i][j].id=sim[j][i].id=-1;
			}
			else if (fabs(vec_cos(vec_sub(t1,t2),sim[i][j].norm))>0.15) {
				fprintf(stdout,"Chain %5d and %5d are not simply rotated\n",i,j);
				sim[i][j].id=sim[j][i].id=-1;
			}

			free(seq2);
			free(koord2);
		}
		free(seq1);
		free(koord1);
	}
	fprintf(stdout,"End of check symmetry\n");
	return sim;
}


