/*****************************************************************************/
/*                                                                           */
/*  optim.c                                                                  */
/*  Copyright (c) Gabor E. Tusnady (tusi), 2003                              */
/*                                                                           */
/*  You may copy, modify, redistribute the source code of optim.c            */
/*  as while as this copyright notice is unchanged.                          */
/*                                                                           */
/*****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "pdb.h"
#include "tmdet_lib.h"


#define OPTIM_BALL_DIST 0.15
#define DIST(a) (*(float *)(a->any))

void optim_init(pdbProtein p);

struct _optimTmp
{	int win;
	int turn;
	double straight;
};
typedef struct _optimTmp *optimTmp;
#define OPTIMTMP_SIZE sizeof(struct _optimTmp)

#define TURN(r) ((optimTmp)(r->any))->turn
#define WIN(r) ((optimTmp)(r->any))->win
#define STRAIGHT(r) ((optimTmp)(r->any))->straight

int OPTIM_RUN=0;
float APOL_RES[]={0,0,0,0,1,1,0,1,0,1,1,0,0,0,0,0,0,1,1,1,0.5};
float HSC[]={1.6,2.0,-9.2,-8.2,3.7,1.0,-3.0,3.1,-8.8,2.8,3.4,
	-4.8,-0.2,-4.1,-12.3,0.6,1.2,2.6,1.9,-0.7,-4.3};

double optim_rotating(pdbProtein p, float width, POINT *axis, POINT *memp);
int optim_numCross(pdbProtein p, float beg, float end, tmdetTypes type,
	int *loop, int *lv);

void optim_setDistances(pdbProtein p, POINT axis, POINT memp)
{
pdbChain c;
pdbAtom a;
float d;

	if (!OPTIM_RUN)
		optim_init(p);
	for_all_atom(p,c,a)
	{	d=axis.x*(a->coord.x-memp.x);
		d+=axis.y*(a->coord.y-memp.y);
		d+=axis.z*(a->coord.z-memp.z);
		DIST(a)=d;
	}
}

double optim_findBestSlice(pdbProtein p, float width, double *shift)
{
pdbChain c;
pdbRes res,pr,ppr,pppr,nr,nnr,nnnr,er,br,sr;
pdbAtom a,pa,ppa,pppa,na,nna,nnna,ba,ea;
double min,max;
int i,j,sdb,r;
double *apol,*surf,*turn,*straight,*lv,*tsurf,*h,q,qs,best;
int *db,*cadb,*hdb,bestpos,empty,len;
double aapol,asurf;
char aa20[]="ACDEFGHIKLMNPQRSTVWY";
double we,swe;

	*shift=0; bestpos=0; best=-1e30;
	min=10000; max=-10000;
	for_all_atom(p,c,a)
 	{	if (min>DIST(a)) min=DIST(a);
 		if (max<DIST(a)) max=DIST(a);
 	}
 	if (max-min<2*width+2)
 		return 0.0;

 	sdb=(int)(max-min+1);
 	apol=calloc(sdb,sizeof(double));
 	surf=calloc(sdb,sizeof(double));
 	turn=calloc(sdb,sizeof(double));
 	tsurf=calloc(sdb,sizeof(double));
 	straight=calloc(sdb,sizeof(double));
 	h=calloc(sdb,sizeof(double));
 	lv=calloc(sdb,sizeof(double));

 	db=calloc(sdb,sizeof(int));
	cadb=calloc(sdb,sizeof(int));
	hdb=calloc(sdb,sizeof(int));

 	for (i=0; i<sdb; i++)
 	{	apol[i]=0;
 		surf[i]=0;
 		turn[i]=0;
 		tsurf[i]=0;
 		straight[i]=0;
 		h[i]=0;
 		lv[i]=0;
 		db[i]=0;
 		cadb[i]=0;
 		hdb[i]=0;
 	}

 	for (c=p->fch; c!=NULL; c=c->next) if (c->sel)
 	for (res=c->fres; res!=NULL; res=res->next) if (res->sel)
 	{	TURN(res)=0;
 		STRAIGHT(res)=0;
 		if ((a=pdb_findAtomInRes(res,"CA"))!=NULL)
		if (a->any!=NULL)
 		if ((pr=res->prev)!=NULL)
 		if ((pa=pdb_findAtomInRes(pr,"CA"))!=NULL)
		if (pa->any!=NULL)
 		if ((ppr=pr->prev)!=NULL)
 		if ((ppa=pdb_findAtomInRes(ppr,"CA"))!=NULL)
		if (ppa->any!=NULL)
 		if ((pppr=ppr->prev)!=NULL)
 		if ((pppa=pdb_findAtomInRes(pppr,"CA"))!=NULL)
		if (pppa->any!=NULL)
 		if ((nr=res->next)!=NULL)
 		if ((na=pdb_findAtomInRes(nr,"CA"))!=NULL)
		if (na->any!=NULL)
 		if ((nnr=nr->next)!=NULL)
 		if ((nna=pdb_findAtomInRes(nnr,"CA"))!=NULL)
		if (nna->any!=NULL)
 		if ((nnnr=nnr->next)!=NULL)
 		if ((nnna=pdb_findAtomInRes(nnnr,"CA"))!=NULL)
		if (nnna->any!=NULL)
 		{	if ((q=pdb_dist(nnna->coord,pppa->coord))>1)
 			{	if (DIST(pppa)<DIST(a)&&DIST(a)<DIST(nnna))
 					STRAIGHT(res)=(DIST(nnna)-DIST(pppa))/q;
 				if (DIST(pppa)>DIST(a)&&DIST(a)>DIST(nnna))
 					STRAIGHT(res)=(DIST(pppa)-DIST(nnna))/q;
 			}
			if (DIST(pppa)<DIST(a)&&DIST(a)>DIST(nnna))
 				TURN(res)=1;
 			if (DIST(pppa)>DIST(a)&&DIST(a)<DIST(nnna))
 				TURN(res)=1;
 		}
 	}


 	for (c=p->fch; c!=NULL; c=c->next) if (c->sel)
	if (c->type&(PDB_CHAIN_CA|PDB_CHAIN_BB))
	for (res=c->fres; res!=NULL; res=res->next) if (res->sel)
	{	er=res;
		if (res!=NULL)
		if (res->any != NULL)
 		for (br=res; (br!=NULL && (br->any != NULL) && STRAIGHT(br)==0); br=br->next);
 		if (br!=NULL)
 		{	for (er=br,len=0; (er!=NULL && (er->any != NULL) && STRAIGHT(er)>0); er=er->next,len++);
 			if (er!=NULL)
 			{	q=0;
 				if (len>=15) {
					ba=pdb_findAtomInRes(br,"CA");
 					ea=pdb_findAtomInRes(er->prev,"CA");
 					if (ba!=NULL && ea!=NULL &&ba->sel && ea->sel) {
 					 	q=fabs(DIST(ba)-DIST(ea))/pdb_dist(ba->coord,ea->coord);
					}
 				}
 				for (sr=br; sr!=er; sr=sr->next) STRAIGHT(sr)=q;
 				er=er->prev;
 			}
 			else er=c->lres;
 		}
 		res=er;
 	}

 	for (c=p->fch; c!=NULL; c=c->next) if (c->sel)
 	{	if (c->type&(PDB_CHAIN_CA|PDB_CHAIN_BB))
 		{	for (res=c->fres; res!=NULL; res=res->next)
 			if ((a=pdb_findAtomInRes(res,"CA"))!=NULL)
			if (a->any != NULL)
 			{	r=(int)(DIST(a)-min);
 				db[r]++;
 				cadb[r]++;
 				turn[r]+=TURN(res);
 				tsurf[r]++;
 				straight[r]+=STRAIGHT(res);
 	 			if (res==c->fres||res==c->lres)
 					lv[r]++;
 				if (pdb_isStdAminoAcid(p,RA3(res)))
 				{	h[r]+=(HSC[(int)(strchr(aa20,RA1(a->res))-aa20)]+12.3)/16.0;
 					hdb[r]++;
 				}
 			}
 		}
 		else
 		{	for (res=c->fres; res!=NULL; res=res->next)
 			for (a=res->fa, i=0; (a!=NULL&&i<res->na); i++, a=a->next) if (a->sel)
 			{	r=(int)(DIST(a)-min);
 				db[r]++;
 				if (!strcmp(ANAME(a),"CA"))
 				{	cadb[r]++;
 					if (p->type&(PDB_CHAIN_CA|PDB_CHAIN_BB))
 					{	turn[r]+=TURN(res);
 						tsurf[r]++;
 					}
 					else
 					{	turn[r]+=TURN(res)*res->surf;
 						tsurf[r]+=res->surf;
 					}
 					straight[r]+=STRAIGHT(res);
 	 				if (res==c->fres||res==c->lres)
 						lv[r]++;
 				}
 				if (strchr(aa20,RA1(res))==NULL)
					j=20;
				else
					j=(int)(strchr(aa20,RA1(res))-aa20);
				apol[r]+=APOL_RES[j]*a->surf;
				surf[r]+=a->surf;
			}
 		}
 	}


 	max=0;
 	for (i=0; i<sdb; i++)
 	{	if (cadb[i]>0)
 		{	straight[i]/=cadb[i];
	 		lv[i]/=cadb[i];
	 	}
		else
		{	straight[i]=0;
			lv[i]=0;
		}
		if (tsurf[i]>1e-5)
			turn[i]/=tsurf[i];
		else
			turn[i]=0;
		if (hdb[i]>0)
			h[i]/=hdb[i];
		if (max<h[i])
			max=h[i];
 	}
 	if (max<1e-5)
 		for (i=0; i<sdb; i++) h[i]=1;

	for (i=width+1; i<sdb-width; i++)
	{	empty=0;
		aapol=asurf=0;
		for (j=-width, empty=0; j<=width; j++)
		{	aapol+=apol[i+j];
			asurf+=surf[i+j];
			if (db[i+j]==0) empty++;
		}
		if (empty<3)
		{	if (p->type&(PDB_CHAIN_CA|PDB_CHAIN_BB))
					q=1;
			else
			{	if (asurf>1e-5)
					q=aapol/asurf;
				else
					q=0;
			}
			qs=0;
			for (j=-width, swe=0; j<=width; j++)
			{	we=abs(j); we/=width; we=1-we; we=1;
				qs+=straight[i+j]*(1-turn[i+j])*(1-lv[i+j])*h[i+j]*we;
				swe+=we;
			}
			qs/=swe;
			qs*=q;
			if (best<qs)
			{	best=qs;
				bestpos=i;
			}
		}
	}
	*shift=min+bestpos;
	return best;
}

double optim_rotating(pdbProtein p, float width, POINT *axis, POINT *memp)
{
POINT sp,ax,mp;
double best,a_dist,q,qq,a1,a2,score,step;
double shift=0;

	sp.x=memp->x;
	sp.y=memp->y;
	sp.z=memp->z;
	best=-1e30;
	a_dist=2*M_PI/OPTIM_BALL_DIST;
	for (a1=0; a1<M_PI/2; a1+=M_PI/a_dist)
	{	q=sin(a1);
		qq=cos(a1);
		if (q>1e-10) step=OPTIM_BALL_DIST/q; else
		{	step=2*M_PI; q=0;}
		for (a2=0; a2<2*M_PI; a2+=step)
		{	ax.x=cos(a2)*q;
			ax.y=sin(a2)*q;
			ax.z=qq;
			mp.x=sp.x;
			mp.y=sp.y;
			mp.z=sp.z;
			optim_setDistances(p,ax,mp);
			score=optim_findBestSlice(p,width,&shift);
			if (score>best)
			{	axis->x=ax.x;
				axis->y=ax.y;
				axis->z=ax.z;
				memp->x=mp.x+ax.x*shift;
				memp->y=mp.y+ax.y*shift;
				memp->z=mp.z+ax.z*shift;
				best=score;
			}
		}
	}
	return best;
}

void optim_init(pdbProtein p)
{
pdbChain c;
pdbRes r;
pdbAtom a;

	for_all_atom(p,c,a)
 		a->any=malloc(sizeof(float));
	for (c=p->fch; c!=NULL; c=c->next) if (c->sel)
	for (r=c->fres; r!=NULL; r=r->next) {
		r->any=malloc(OPTIMTMP_SIZE);
		TURN(r)=0;
 		STRAIGHT(r)=0;
		WIN(r) = 0;
	}
	OPTIM_RUN=1;
}

void optim_end(pdbProtein p)
{
pdbChain c;
pdbRes r;
pdbAtom a;

	for_all_atom(p,c,a)
 	{	free(a->any);
 		a->any=NULL;
 	}
	for (c=p->fch; c!=NULL; c=c->next) if (c->sel)
	for (r=c->fres; r!=NULL; r=r->next)
	{	free(r->any);
		r->any=NULL;
	}
	OPTIM_RUN=0;
}

double optim_findTM(pdbProtein p, float width, POINT *axis, POINT *memp)
{
pdbChain c;
pdbAtom a;
float q;
int db=0;

	if (!OPTIM_RUN)
		optim_init(p);
	memp->x=memp->y=memp->z=0;
	for_all_atom(p,c,a)
 	{	memp->x+=x(a->coord);
 		memp->y+=y(a->coord);
 		memp->z+=z(a->coord);
 		db++;
 	}
 	memp->x/=db;
 	memp->y/=db;
 	memp->z/=db;
 	q=optim_rotating(p,width,axis,memp);
	optim_end(p);
 	return 100*q;

}

double optim_findTMSlice(pdbProtein p, float width, POINT axis, POINT *memp)
{
double score,shift=0;
pdbChain c;
pdbAtom a;
int db=0;
POINT mp;

	if (!OPTIM_RUN)
		optim_init(p);
	mp.x=mp.y=mp.z=0;
	for_all_atom(p,c,a)
 	{	mp.x+=x(a->coord);
 		mp.y+=y(a->coord);
 		mp.z+=z(a->coord);
 		db++;
 	}
 	mp.x/=db;
 	mp.y/=db;
 	mp.z/=db;
	optim_setDistances(p,axis,mp);
	score=optim_findBestSlice(p,width,&shift);
	memp->x=mp.x+axis.x*shift;
	memp->y=mp.y+axis.y*shift;
	memp->z=mp.z+axis.z*shift;
	return score*100;
}

tmdetTypes optim_setType(pdbProtein p, float beg, float end, POINT axis, POINT memp)
{
pdbChain c;
pdbRes res;
pdbAtom a;
int adb,bdb,cdb;

	adb=bdb=cdb=0;
	for (c=p->fch; c!=NULL; c=c->next) if (c->sel)
 	for (res=c->fres; res!=NULL; res=res->next)
 	if ((a=pdb_findAtomInRes(res,"CA"))!=NULL)
 	if (a->sel)
 	if (DIST(a)>beg&&DIST(a)<end)
 	{	if ((res->ss=='H')||(res->ss=='G')||(res->ss=='I'))
 		{	adb++;
 		}
 		else
 		{	if (res->ss=='E') bdb++; else cdb++;}
 	}
 	if (cdb>adb*2&&cdb>bdb*2) return TMDET_COIL;
 	if (bdb>adb) return TMDET_BETA;
	return TMDET_ALPHA;
}


int optim_setWidth(pdbProtein p, POINT *axis, POINT *memp, tmdetTypes type, float *w)
{
int dbe,db,loop,lv,noteq;
float beg,end,move,width;

#define WST 0.5

	beg=-(*w); end=(*w); loop=-1; lv=-1;

	db=dbe=optim_numCross(p,beg,end,type,&loop,&lv);
	while (db<0&&beg<-7.0)
	{	beg++; end--; loop=-1; lv=-1;
		db=dbe=optim_numCross(p,beg,end,type,&loop,&lv);
	}
	if (dbe>0)
	{	for (beg-=WST, end+=WST, db=dbe,noteq=0; (noteq<6&&beg>-25&&end<25&&end-beg<40); beg-=WST, end+=WST) {
			db=optim_numCross(p,beg,end,type,&loop,&lv);
			if (db!=dbe)
				noteq++;
			else
				noteq=0;
			if (db<11&&type==TMDET_BETA) {
				if (end-beg>22)
					noteq++;
			}
		}
		beg+=8*WST;
		end-=8*WST;
		for (end+=WST, db=dbe, noteq=0; (noteq<4&&end<25&&end-beg<40); end+=WST) {
			db=optim_numCross(p,beg,end,type,&loop,&lv);
			if (db!=dbe)
				noteq++;
			else
				noteq=0;
			if (db<11&&type==TMDET_BETA) {
				if (end-beg>22)
					noteq++;
			}
		}
		end-=6*WST;
		for (beg-=WST, db=dbe,noteq=0; (noteq<4&&beg>-25&&end-beg<40); beg-=WST) {
			db=optim_numCross(p,beg,end,type,&loop,&lv);
			if (db!=dbe)
				noteq++;
			else
				noteq=0;
			if (db<11&&type==TMDET_BETA) {
				if (end-beg>22)
					noteq++;
			}
		}
		beg+=6*WST;
	}
	end-=WST;
	beg+=WST;
	move=(beg+end)/2;
	memp->x+=move*axis->x;
	memp->y+=move*axis->y;
	memp->z+=move*axis->z;
	width=(end-beg)/2;
	axis->x*=width;
	axis->y*=width;
	axis->z*=width;
	*w=width;
	return dbe;
}

int optim_numCross(pdbProtein p, float beg, float end,
	tmdetTypes type, int *dbloop, int *dblv)
{
pdbChain c;
pdbRes res,br,er,pr;
pdbRes first,last;
pdbAtom a;
int db,dbp,ok,len,i;
int in,mem,out;
int loop,lv;
int wp,wn,nn;

	in=mem=out=loop=lv=0;
	for (c=p->fch; c!=NULL; c=c->next) if (c->sel)
	{	for (res=c->fres,first=last=NULL; res!=NULL; res=res->next) if (res->sel)
		{	if (first==NULL) first=res;
			last=res;
			WIN(res)=0;
			for (a=res->fa, i=0; (i<res->na&&a!=NULL); a=a->next,i++) if (a->sel)
				{	if (!strcmp(ANAME(a),"CA")||!strcmp(ANAME(a),"C")||
						!strcmp(ANAME(a),"N")||!strcmp(ANAME(a),"O"))
						if (DIST(a)>beg&&DIST(a)<end)
							WIN(res)=1;
				}
			if (!WIN(res))
			{	if (res->fa!=NULL)
				if (res->fa->any!=NULL) {
					if (DIST(res->fa)<0)
					{	WIN(res)=-2;
						in++;
					}
					else
					{	WIN(res)=2;
						out++;
					}
				}
			}
			else mem++;

			if (WIN(res)==1&&res->surf<2.0) WIN(res)=0;
//			printf("%c %d %.2f %d\n",res->chid,res->rsn,res->surf,WIN(res));
		}
		if (first!=NULL)
			if (WIN(first)==1) lv++;
		if (last!=NULL)
			if (WIN(last)==1) lv++;

	}
//	if (in<2||out<2||mem<5)
//		return -1;
	db=0;
	for (c=p->fch; c!=NULL; c=c->next) if (c->sel)
	for (res=c->fres; res!=NULL; res=res->next) if (res->sel)
	if (WIN(res)==0) {
		for (br=res; (br!=NULL&&WIN(br)==0); br=br->next);
		wp=wn=0;
		if (res->prev!=NULL) wp=WIN(res->prev);
		if (br!=NULL) wn=WIN(br);
		if (wp==1||wn==1)
			nn=1;
		else {
			if (wp==2) nn=2;
			else
				nn=-2;
		}
		if (wp!=0)
			for (er=res; er!=br; er=er->next) WIN(er)=nn;
		if (wn!=0)
			for (er=res; er!=br; er=er->next) WIN(er)=nn;
		if (br!=NULL)
			res=br->prev;
	}
	for (c=p->fch; c!=NULL; c=c->next) if (c->sel)
	{	for (res=c->fres; res!=NULL; res=res->next) if (res->sel)
		{	for (br=res; (br!=NULL&&WIN(br)!=1); br=br->next);
			if (br==NULL)
			{	res=c->lres; }
			else
			{	res=br;
				for (er=br->next,len=0; (er!=NULL&&WIN(er)==1); er=er->next)
				if (er->sel)
					len++;
				if (er==NULL) res=c->lres;
				else
				{	res=er;
					pr=br->prev;
					ok=0;
					dbp=1;
					if (pr!=NULL)
					{	if (WIN(pr)*WIN(er)<0) ok=1;
						else
						{	if (type==TMDET_BETA) {
								if (len>15) {
									ok=1; dbp=2; loop++;
								}
							}
							else if (type==TMDET_ALPHA) {
								if (len>5) {
									ok=1; dbp=2; loop++;
								}
							}
							else
								loop++;
						}
					}
					if (ok)
					{	db+=dbp;}
				}
			}
		}
 	}
	if (*dbloop!=-1) {
 		if (loop!=*dbloop) {
//			printf("Loop not ok: %d %d\n",loop,*dbloop);
			db=-1;
		}
	}
 	else {
		*dbloop=loop;
//		printf("Loop a vegen: %d\n",loop);
	}
 	if (*dblv!=-1) {
 		if (lv!=*dblv) {
//			printf("Lv not ok: %d %d\n",lv,*dblv);
			db=-1;
		}
	}
 	else {
		*dblv=lv;
//		printf("Lv a vegen: %d\n",lv);
	}
 	return db;
}

void tm_rotateToZ(pdbProtein p, POINT axis, POINT memp, pdbMatrix rot)
{
pdbMatrix m,rotn;
double a,b,c,ca,sa,d;
int i,j,k;

	m=pdb_newMatrix();
	a=axis.x;
	b=axis.y;
	c=axis.z;
	d=sqrt(b*b+c*c);
	if (d>1e-5)
	{	sa=c/d;
		ca=b/d;
		m->m[0][0]=d; m->m[0][1]=-a*ca; m->m[0][2]=-a*sa; m->m[0][3]=-memp.x;
		m->m[1][0]=0; m->m[1][1]=sa;    m->m[1][2]=-ca;  m->m[1][3]=-memp.y;
		m->m[2][0]=a; m->m[2][1]=d*ca;  m->m[2][2]=d*sa; m->m[2][3]=-memp.z;
	}
	else
	{	m->m[0][0]=0; m->m[0][1]=0; m->m[0][2]=1; m->m[0][3]=-memp.x;
		m->m[1][0]=0; m->m[1][1]=1; m->m[1][2]=0; m->m[1][3]=-memp.y;
		m->m[2][0]=1; m->m[2][1]=0; m->m[2][2]=0; m->m[2][3]=-memp.z;
	}

	pdb_transformProtein(p,m);
	rotn=pdb_newMatrix();
	for (i=0; i<4; i++)
	{	for (j=0; j<4; j++)
		for (k=0, rotn->m[i][j]=0; k<4; k++)
			rotn->m[i][j]+=m->m[i][k]*rot->m[k][j];
 	}
	for (i=0; i<4; i++) for (j=0; j<4; j++)
		rot->m[i][j]=rotn->m[i][j];

	pdb_freeMatrix(rotn);
	pdb_freeMatrix(m);
}
