/*****************************************************************************/
/*                                                                           */
/*  util.c                                                                   */
/*  Copyright (c) Gabor E. Tusnady (tusi), 2003                              */
/*                                                                           */
/*  You may copy, modify, redistribute the source code of util.c             */
/*  as while as this copyright notice is unchanged.                          */
/*                                                                           */
/*****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "../../PdbLib/include/pdb.h"
#include "pdb.h"
#include "tmdet_lib.h"
#include "uniprot.h"
#include "util.h"
#include "xml.h"

struct _tmType
{	tmdetTypes t;
	int clust;
};
typedef struct _tmType *tmType;
#define TMTYPE_SIZE sizeof(struct _tmType)

#define TMT(a) ((tmType)(a->any))->t
#define TMB(a) ((tmType)(a->any))->clust


float tm_norm(POINT *v);
void tm_pymolSet(FILE *f,char c, pdbRes b, pdbRes e, char type, char *color);
void tm_oavSet(FILE *f,char c, pdbRes b, pdbRes e, char type, char *color);
void tm_pvSet(FILE *f, char id, pdbRes b, pdbRes e, char type, char *color);
void tm_countBeta(pdbProtein p);
void tm_checkLoop(pdbRes re, pdbRes rk, pdbRes rn, double h);
void tm_splitHelix(pdbRes re, pdbRes rk, pdbRes rn);

char tmpart_chains[100];
double IFH_MEMBRDIST=9.0;
double IFH_ANGLE_LIMIT=30.0;


void tm_setTypes(pdbProtein p, tmDet tm) {
pdbChain c;
pdbRes r,re,rs,rn,rsp,rspp,rsppp,rspppp,rss;
pdbAtom a;
int i,count,slen,ok,hi,sign,psign,turn,cdb;
float h,d;
tmdetTypes t,tp,tn;
int len;
POINT v;
tmdetChain ch,lch;
tmdetRegion reg,lreg;
char s[100];
char pattern[100];
int tmpart;

	v.x=tm->axis.x;
	v.y=tm->axis.y;
	v.z=tm->axis.z;

	h=tm_norm(&v);
	for (c=p->fch; c!=NULL; c=c->next) if (c->sel)
	{	//printf("\n<---------- %c ---------->\n", c->id);
		if (strchr(tmpart_chains,c->id)==NULL)
			tmpart=0;
		else
			tmpart=1;
		for (r=c->fres; r!=NULL; r=r->next)
		{	r->any=malloc(TMTYPE_SIZE);
			TMT(r)=TMDET_UNKNOWN;
			TMB(r)=-1;
			r->hbe[0]=-1000;
			if ((a=pdb_findAtomInRes(r,"CA"))!=NULL) {
				d=v.x*(a->coord.x-tm->memp.x);
				d+=v.y*(a->coord.y-tm->memp.y);
				d+=v.z*(a->coord.z-tm->memp.z);
				r->hbe[0]=d;
				if (d<-h) TMT(r)=TMDET_SIDE1;
				else {
					if (d<=h) {
						TMT(r)=tm->type;
						if (tm->type==TMDET_CA_TM) TMT(r)=TMDET_ALPHA;
					}
					else
						TMT(r)=TMDET_SIDE2;
				}
				if (tmpart)
					TMT(r)=TMDET_SIDE1;
			}
			//fprintf(stderr,"%c%5d%5d\n",r->chid, r->rsn, TMT(r));
		}
		tp=TMDET_CBEG;
		for (r=c->fres; r!=NULL;)
		{	t=TMT(r);
			for (re=r,len=0; (re!=NULL&&TMT(re)==t); re=re->next,len++);
			if (re==NULL)
				tn=TMDET_CEND;
			else
				tn=TMT(re);
			if ((t==TMDET_ALPHA||t==TMDET_BETA||t==TMDET_COIL)&&
				(((tp==tn)&&tp!=TMDET_UNKNOWN)||tp==TMDET_CBEG||tn==TMDET_CEND
				|| tp==TMDET_UNKNOWN || tn==TMDET_UNKNOWN))
			{
				if (len<=10)
				{	if (tp==tn)
					{	for (rs=r; rs!=re; rs=rs->next) TMT(rs)=tn;
					}
					else
					{	if (tp==TMDET_UNKNOWN) {
							for (rs=r; rs!=re; rs=rs->next) TMT(rs)=tn;
						}
						else if (tn==TMDET_UNKNOWN) {
							for (rs=r; rs!=re; rs=rs->next) TMT(rs)=tp;
						}
						else {
							for (rs=r; rs!=re; rs=rs->next) TMT(rs)=TMDET_INSIDE;
						}
					}
				}
				else
				{	if (tp==tn) {
						for (rs=r; rs!=re; rs=rs->next) TMT(rs)=TMDET_LOOP;
					}

				}
			}
			r=re;
			tp=t;

		}

		/*Loop correction*/
		/* TODO */
		/* somehow this section changes the the type of the residue on G in 1be3 on the longest helix */

		if (tm->type==TMDET_ALPHA)
		for (r=c->fres; r!=NULL;) {
			for (re=r; (re!=NULL && !tmdet_AlphaType[TMT(re)]); re=re->next);
			if (re!=NULL) {
				for (rn=re,len=0; (rn!=NULL && tmdet_AlphaType[TMT(rn)]); rn=rn->next,len++);
				if (rn==NULL) rn=c->lres;
				if (len>10) {
					rsp=rspp=rsppp=rspppp=NULL;
					psign=0;
					sign=0;
					turn=0;
					for (rs=re; rs!=rn; rs=rs->next)
						rs->hb[0]=0;
					for (rss=rs=re; rs!=rn; rs=rs->next) {
						if (rspppp!=NULL) {
							if (rspppp->hbe[0]<rs->hbe[0])
								sign=1;
							else
								sign=-1;
							if (psign*sign<0) {
								if (fabs(rss->hbe[0]-rsp->hbe[0])>5) {
									rsp->hb[0]=1;
									turn++;
									rss=rsp;
								}
							}
						}
						rspppp=rsppp;
						rsppp=rspp;
						rspp=rsp;
						rsp=rs;
						psign=sign;
					}
					if (re->hbe[0]*rn->hbe[0]<0)
						sign=-1;
					else
						sign=1;
					if (sign==1&&turn==0)
						t=TMDET_LOOP;
					else
						t=TMDET_ALPHA;
					for (rs=re; rs!=rn; rs=rs->next) {
						TMT(rs)=t;
						if (rs->hb[0]==1 || rs==re || rs==rn) {
							if (rs->hbe[0]<0)
								TMT(rs)=TMDET_SIDE1;
							else
								TMT(rs)=TMDET_SIDE2;
						}
					}
					//("Turns: %d (%.2f), ",re->rsn,re->hbe[0]);
					//for (rs=re; rs!=rn; rs=rs->next)
					//	if (rs->hb[0]==1)
					//		printf("%d (%.2f), ",rs->rsn,rs->hbe[0]);
					//printf("%d (%.2f)\n",rn->rsn,rn->hbe[0]);
					for (rsp=rs=re; rs!=rn; rs=rs->next)
						if (rs->hb[0]==1) {
							for (rss=rs,rs=rs->next; (rs!=rn && rs->hb[0]!=1); rs=rs->next);
							//printf("%d %d %d\n",rsp->rsn,rss->rsn,rs->rsn);
							if (rss->rsn-rsp->rsn<15&&rs->rsn-rss->rsn<15)
								tm_checkLoop(rsp,rss,rs,h);
							else {
								if (rs->rsn-rss->rsn<=10) {
									for (rspp=rss;rspp!=rs;rspp=rspp->next)
										TMT(rspp)=TMT(rss);
									//printf("Type changed to %s for %d-%d\n",
									//	tmdet_TypeNames[TMT(rss)],rss->rsn,rs->rsn);
								}
								if (rss->rsn-rsp->rsn<=10) {
									for (rspp=rsp;rspp!=rss;rspp=rspp->next)
										TMT(rspp)=TMT(rsp);
									//printf("Type changed to %s for %d-%d\n",
									//	tmdet_TypeNames[TMT(rsp)],rsp->rsn,rss->rsn);
								}
							}
							rsp=rss;
							if (rs==rn) rs=rn->prev; else rs=rs->prev;
						}
				}
			}
			else rn=c->lres;
			r=rn->next;
		}

		/*End of loop correction*/


		for (re=c->fres; (re!=NULL&&TMT(re)==TMDET_UNKNOWN); re=re->next);
		for (r=re,len=hi=0,ok=1; (r!=NULL && ok); r=r->next,len++) {
			if (TMT(r)>=TMDET_GLOB)
				ok=0;
			else
				hi=1;
		}
		if (len<10&&hi&&re!=NULL&&r!=NULL)
			for (rs=re; rs!=r; rs=rs->next) TMT(rs)=TMT(r);
		for (re=c->lres; (re!=NULL&&TMT(re)==TMDET_UNKNOWN); re=re->prev);
		for (r=re,len=hi=0,ok=1; (r!=NULL && ok); r=r->prev,len++) {
			if (TMT(r)>=TMDET_GLOB)
				ok=0;
			else
				hi=1;
		}
		if (len<10&&hi&&re!=NULL&&r!=NULL)
			for (rs=re; rs!=r; rs=rs->prev) TMT(rs)=TMT(r);
	}


	/* Interfacial helix detection */
	//if (tm->type==TMDET_ALPHA)
	//	tm_detectIFH(p,tm);

	if (tm->type==TMDET_BETA)
		tm_countBeta(p);

	lch=NULL;
	for (c=p->fch; c!=NULL; c=c->next) if (c->sel)
	{	ch=tmdet_newChain(lch);
		if (lch==NULL)
			tm->c=ch;
		lch=ch;
		ch->id=c->id;
		count=0;
		for (r=c->fres; r!=NULL;)
		{	for (re=r->next; (re!=NULL&&TMT(re)==TMT(r)); re=re->next);
			if (TMT(r)==TMDET_BETA||TMT(r)==TMDET_ALPHA||TMT(r)==TMDET_COIL)
			{	count++;
			}
			r=re;
		}
		strcpy(s,"error");
		if (tm->type==TMDET_ALPHA||tm->type==TMDET_CA_TM)
			strcpy(s,"alpha");
		else if (tm->type==TMDET_BETA)
			strcpy(s,"beta");
		else if (tm->type==TMDET_COIL)
			strcpy(s,"coil");

		if (c->sel) ch->sel=True; else ch->sel=False;
		ch->numtm=count;
		if (tm->type==TMDET_BETA&&(count%2)==1) {
			printf("!!!!!WARNING!!!!\n");
			printf("!!Invalide beta barrel structure detected for %s.\n",tm->code);
			printf("!!Please check it\n");
		}

		if (count==0) {
			if (tm->type==TMDET_TMPART) {
				if (strchr(tmpart_chains,c->id)==NULL)
					ch->sel=False;
				else
					strcpy(s,"tm_part");
			}
			else
				strcpy(s,"non_tm");
		}

		if (c->sel)
		{	ch->type=calloc(strlen(s)+1,sizeof(char));
			strcpy(ch->type,s);
			ch->rseq=calloc(c->nr+1,sizeof(char));
			for (i=0; i<c->nr; i++)
				ch->rseq[i]=RA1(c->res[i]);
			ch->rseq[i]='\0';
			lreg=(tmdetRegion)NULL;
			for (r=c->fres,cdb=0; r!=NULL;)
			{	for (rs=r,re=r->next,slen=0; (re!=NULL&&TMT(re)==TMT(r)&&re->icode==r->icode&&(abs(re->rsn-rs->rsn)<=1)); rs=re,re=re->next) slen++;
				if (re==NULL) rs=c->lres; else rs=re->prev;
				if (r!=rs) {
					//printf("Region: %d-%d: %p\n",r->rsn,rs->rsn,lreg);
					reg=tmdet_newRegion(lreg);
					if (lreg==NULL)
						ch->reg=reg;
					lreg=reg;
					reg->beg=r->rsn;
					reg->begi=r->icode;
					reg->end=rs->rsn;
					reg->endi=rs->icode;
					reg->type=TMT(r);
					reg->rbeg=r->relsn+1;
					reg->rend=rs->relsn+1;	
					pattern[cdb++]=tmdet_TypeCodes[reg->type];
				}
				r=re;
        	}
			pattern[cdb]='\0';
			if (strstr(pattern,"1H1")!=NULL ||
				strstr(pattern,"2H2")!=NULL ||
				strstr(pattern,"12")!=NULL ||
				strstr(pattern,"21")!=NULL ||
				strstr(pattern,"1B1")!=NULL ||
				strstr(pattern,"2B2")!=NULL) {
					printf("!!Invalide structure has been detected: %s >%c< (%s)\n",
						tm->code,ch->id,pattern);
				}
      	}
   }
   for (c=p->fch; c!=NULL; c=c->next) if (c->sel)
	{	for (r=c->fres; r!=NULL; r=r->next)
		{	free(r->any);
			r->any=NULL;
		}
	}
	//tm_detectAndStoreIFH(p, tm);
}


void tm_checkStructure(tmDet tm) {
char pattern[1000];
unsigned int cdb;
tmdetChain c;
tmdetRegion reg;

	for (c=tm->c; c!=NULL; c=c->next) if (c->sel) {
		for (reg=c->reg,cdb=0; reg!=NULL; reg=reg->next)
			pattern[cdb++]=tmdet_TypeCodes[reg->type];
		pattern[cdb]='\0';
		if (strstr(pattern,"1H1")!=NULL ||
			strstr(pattern,"2H2")!=NULL ||
			strstr(pattern,"12")!=NULL ||
			strstr(pattern,"21")!=NULL ||
			strstr(pattern,"1B1")!=NULL ||
			strstr(pattern,"2B2")!=NULL) {
				printf("!!Invalide structure has been detected: %s >%c< (%s)\n",
					tm->code,c->id,pattern);
		}
	}
}

void tm_helixDirection(pdbAtom A, pdbAtom B, pdbAtom C, double res[3]) {
double a[3],b[3],x[3],y[3],z[3];
double ah,bh,xh,yh,zh;
double alpha = M_PI * 100.0 / 180.0;
int i;

	ah=0;
	a[0]=B->coord.x-A->coord.x; ah+=a[0]*a[0];
	a[1]=B->coord.y-A->coord.y; ah+=a[1]*a[1];
	a[2]=B->coord.z-A->coord.z; ah+=a[2]*a[2];
	ah=sqrt(ah);

	bh=0;
	b[0]=C->coord.x-B->coord.x; bh+=b[0]*b[0];
	b[1]=C->coord.y-B->coord.y; bh+=b[1]*b[1];
	b[2]=C->coord.z-B->coord.z; bh+=b[2]*b[2];
	bh=sqrt(bh);

	xh=yh=0;
	for (i=0; i<3; i++) {
		a[i]/=ah;
		b[i]/=bh;
		x[i]=b[i]-a[i];
		y[i]=b[i]+a[i];
		xh+=x[i]*x[i];
		yh+=y[i]*y[i];
	}
	xh=sqrt(xh);
	yh=sqrt(yh);

	for (i=0; i<3; i++) {
		x[i]/=xh;
		y[i]/=yh;
	}

	zh=0;
	z[0]=x[1]*y[2]-x[2]*y[1]; zh+=z[0]*z[0];
	z[1]=x[2]*y[0]-x[0]*y[2]; zh+=z[1]*z[1];
	z[2]=x[0]*y[1]-x[1]*y[0]; zh+=z[2]*z[2];
	zh=sqrt(zh);

	for (i=0; i<3; i++) {
		z[i]/=zh;
		res[i]= y[i]*cos(alpha) + z[i]*sin(alpha);
//		printf("%10.2f",res[i]);
	}
//	printf("\n");
}

void tm_detectIFH(pdbProtein p, tmDet tm) {
pdbRes r,rr;
pdbRes ifh_beg, ifh_end, beg, end;
pdbAtom cylinder[6];
int i,j,hl;
float h, d, dd;
POINT u,v,w;
tmdetChain ch;
tmdetRegion reg;
double dir[100][3];
int dird;
pdbAtom A, B, C;

	v.x=tm->axis.x;
	v.y=tm->axis.y;
	v.z=tm->axis.z;
	h=tm_norm(&v);

	for (ch=tm->c; ch!=NULL; ch=ch->next) if (ch->sel==True) {
		for (reg=ch->reg; reg!=NULL; reg=reg->next) {
			beg=pdb_findResidue(p,ch->id,reg->beg);
			end=pdb_findResidue(p,ch->id,reg->end);
			if (beg!=NULL&&end!=NULL) {

				for (r=beg; r!=end->next; r=r->next) {
					r->any=malloc(TMTYPE_SIZE);
					TMT(r)=reg->type;
					if (r->ss=='G')
						r->ss='H';
					if (r->ss=='I')
						r->ss='H';
				}

				if ( reg->type==TMDET_SIDE1 || reg->type==TMDET_SIDE2 ) {
					for (r=beg; ( r!=NULL && r!=end->next);) {
						ifh_beg = ifh_end = NULL;
						// looking for helices
						for( ; ( r!=NULL && r->ss!='H' && r!=end ); r=r->next);
						if (r==end) break;
						ifh_beg = r;
						//printf("Begin: %d\n",r->rsn);
						for ( ; (r!=NULL && r!=end->next); r=r->next ) {
							if ( r->ss!='H' ) break;
							ifh_end = r;
						}
						hl = ifh_end->rsn - ifh_beg->rsn + 1;
						//printf("Helix: %c %d %d length: %d \n ", ifh_beg->chid, ifh_beg->rsn, ifh_end->rsn, hl);
						if ( hl < 6 ) continue;

						int err=0;
						//Calculating the helix directions for each ca in the helix
						for (dird=1, rr=ifh_beg->next; rr!=ifh_end; rr=rr->next, dird++ ) {
							A = pdb_findAtomInRes( rr->prev, "CA" );
							B = pdb_findAtomInRes( rr, "CA" );
							C = pdb_findAtomInRes( rr->next, "CA" );
							if ( A == NULL || B == NULL || C == NULL )
								err=1;
							else {
								tm_helixDirection(A, B, C, dir[dird] );
							}
						}
						//err means missing ca atom
						if ( err )
							continue;

						//calculating straightness of the helix
							for (i=2, rr=ifh_beg->next->next; i<dird; rr=rr->next, i++) {
								double d = 0;
								double ah=0,bh=0;
								for (j=0; j<3;j++) {
									d+=dir[i][j]*dir[i-1][j];
									ah+=dir[i][j]*dir[i][j];
									bh+=dir[i-1][j]*dir[i-1][j];
								}
								ah=sqrt(ah);
								bh=sqrt(bh);
								//printf("%10.2f%10.2f%10.2f\n",dir[i-1][0],dir[i-1][1],dir[i-1][2]);
								//printf("%10.2f%10.2f%10.2f\n",dir[i][0],dir[i][1],dir[i][2]);
								//printf("Helix alkotok: %2d %.2f\n",i, 180*acos(d/(ah*bh))/M_PI);
								if ( 180*acos(d/(ah*bh))/M_PI > 80.0 )
									rr->ss = ' ';
							}
					}

					//printf("Second round, dealing straight helices\n");
					for (r=beg; r!=end->next;) {
						ifh_beg = ifh_end = NULL;
						// looking for helices
						while( r->ss!='H' && r!=end ) r=r->next;
						if (r==end) break;
						ifh_beg = r;
						for ( ; r!=end->next; r=r->next ) {
							if ( r->ss!='H' ) break;
							ifh_end = r;
						}

						// check helix length
						hl = ifh_end->rsn - ifh_beg->rsn + 1;
						//printf("Helix: %c %d %d length: %d \n ", ifh_beg->chid, ifh_beg->rsn, ifh_end->rsn, hl);
						if ( hl < 6 ) continue;
						// else turn the averages of first and last 3 CA to vector
						else {
							// TODO should it be continuous ????
							for (i=0; i<6; i++) cylinder[i] = (pdbAtom)NULL;
							for (i=0, rr=ifh_beg; rr!=ifh_end && i<3; rr=rr->next, i++)
								cylinder[i] = pdb_findAtomInRes( rr, "CA" );
							if ( i<3 ) {printf("First three are missing\n"); continue;}
							for (i=0, rr=ifh_end; rr!=ifh_beg && i<3; rr=rr->prev, i++)
								cylinder[5-i] = pdb_findAtomInRes( rr, "CA" );
							if ( i<3 ) {printf("Last three are missing\n"); continue;}

							// check if all CAs are found
							for (i=0; i<6; i++) if ( cylinder[i]==NULL ) break;
							if ( i<6 ) {printf("First and last overlap\n"); continue;}

							w.x = ( cylinder[0]->coord.x + cylinder[1]->coord.x + cylinder[2]->coord.x - cylinder[3]->coord.x - cylinder[4]->coord.x - cylinder[5]->coord.x )/3.;
							w.y = ( cylinder[0]->coord.y + cylinder[1]->coord.y + cylinder[2]->coord.y - cylinder[3]->coord.y - cylinder[4]->coord.y - cylinder[5]->coord.y )/3.;
							w.z = ( cylinder[0]->coord.z + cylinder[1]->coord.z + cylinder[2]->coord.z - cylinder[3]->coord.z - cylinder[4]->coord.z - cylinder[5]->coord.z )/3.;
							u.x = ( cylinder[0]->coord.x + cylinder[1]->coord.x + cylinder[2]->coord.x + cylinder[3]->coord.x + cylinder[4]->coord.x + cylinder[5]->coord.x )/6.;
							u.y = ( cylinder[0]->coord.y + cylinder[1]->coord.y + cylinder[2]->coord.y + cylinder[3]->coord.y + cylinder[4]->coord.y + cylinder[5]->coord.y )/6.;
							u.z = ( cylinder[0]->coord.z + cylinder[1]->coord.z + cylinder[2]->coord.z + cylinder[3]->coord.z + cylinder[4]->coord.z + cylinder[5]->coord.z )/6.;

						}

						// criteria
						dd = fabs(v.x*(u.x-tm->memp.x) + v.y*(u.y-tm->memp.y) + v.z*(u.z-tm->memp.z));
						//printf("\tDistance: %lf ", fabs(h-dd));
						if ( fabs(h-dd)>IFH_MEMBRDIST ) continue;
						d = (v.x*w.x + v.y*w.y + v.z*w.z) / ( sqrt( v.x*v.x + v.y*v.y + v.z*v.z) * sqrt( w.x*w.x + w.y*w.y + w.z*w.z));
						d = fabs(acos(d)/M_PI*180.-90.);
						//printf("\tAngle: %lf\n", d);
						if ( d > IFH_ANGLE_LIMIT ) continue;

						// set type
						for (rr=ifh_beg; rr!=ifh_end->next; rr=rr->next)
							TMT(rr) = TMDET_IFH;
						//printf("IntFacHel: %d - %d (%f)\n", ifh_beg->rsn, ifh_end->rsn, dd);
					}
				}
			}
		}
	}
}

void tm_detectAndStoreIFH(pdbProtein p, tmDet tm) {
pdbChain c;
pdbRes r,re, rs;
int slen,cdb;
tmdetChain ch;
tmdetRegion reg,lreg;
tmdetModification mod;
char date[12];
char pattern[100];
time_t ct;
struct tm *ttm;

	if (p==(pdbProtein)NULL) return;
	if (tm==(tmDet)NULL) return;
	if (tm->type!=TMDET_ALPHA) return;

	//tm_detectIFH(p, tm);

	for (ch=tm->c; ch!=NULL; ch=ch->next) if (ch->sel==True) {
		if ((c = pdb_findChain(p->fch, ch->id)) == NULL ) {
			printf("Could not find chain: >>%c<<\n",ch->id);
			exit(EXIT_FAILURE);
		}

		if (c->sel) {
			for (reg=ch->reg; reg!=NULL; reg=lreg) {
				lreg=reg->next;
				tmdet_freeRegion(reg);
				reg = (tmdetRegion)NULL;
			}
			lreg=(tmdetRegion)NULL;
			for (r=c->fres,cdb=0; r!=NULL;) {
				if (r->sel&&r->any!=NULL)
				{	for (re=r->next,slen=0; (re!=NULL&&re->any!=NULL&&r->any!=NULL&&TMT(re)==TMT(r)); re=re->next) if (re->sel) slen++;
					if (re==NULL) rs=c->lres; else rs=re->prev;
					reg=tmdet_newRegion(lreg);
					if (lreg==NULL)
						ch->reg=reg;
					lreg=reg;
					reg->beg=r->rsn;
					reg->end=rs->rsn;
					reg->type=TMT(r);
					reg->rbeg=r->relsn+1;
					reg->rend=rs->relsn+1;
					r=re;
					pattern[cdb++]=tmdet_TypeCodes[reg->type];
				}
				else
					r=r->next;
			}
			pattern[cdb]='\0';
			if (strstr(pattern,"1H1")!=NULL ||
				strstr(pattern,"2H2")!=NULL ||
				strstr(pattern,"12")!=NULL ||
				strstr(pattern,"21")!=NULL ||
				strstr(pattern,"1B1")!=NULL ||
				strstr(pattern,"2B2")!=NULL) {
					printf("!!Invalide structure has been detected: %s >%c< (%s)\n",
						tm->code,ch->id,pattern);
				}
		}
   }

	for (c=p->fch; c!=NULL; c=c->next) if (c->sel) {
		for (r=c->fres; r!=NULL; r=r->next) {
			if (r->any!=NULL) free(r->any);
			r->any=NULL;
		}
	}


	for (ch=tm->c; ch!=NULL; ch=ch->next) if (ch->sel==True)
		for (reg=ch->reg; reg!=NULL; reg=reg->next) {
			if (reg->type==TMDET_IFH) {
				time(&ct);
				ttm=gmtime(&ct);
				sprintf(date,"%4d-%02d-%02d",1900+ttm->tm_year,1+ttm->tm_mon,ttm->tm_mday);

				for (mod=tm->mod; mod!=NULL; mod=mod->next) {
					if ( !strcmp(mod->descr, "Interfacial helix detected.") ) {
						strcpy(mod->date, date);
						return;
					}
				}
				tmdet_addModification(tm,date,"Interfacial helix detected.");
				return;
			}
		}
}

void tm_checkLoop(pdbRes re, pdbRes rk, pdbRes rn, double h) {
tmdetTypes t;
pdbRes r;
	if (fabs(re->hbe[0]-rk->hbe[0])>4.0 && fabs(rk->hbe[0]-rn->hbe[0])>4.0) {
		t=TMDET_LOOP;
	}
	else {
		if (rk->hbe[0]<0)
			t=TMDET_SIDE1;
		else
			t=TMDET_SIDE2;
	}
	for (r=re; r!=rn; r=r->next)
		TMT(r)=t;
	if (re->hbe[0]<0)
		TMT(re)=TMT(rn)=TMDET_SIDE1;
	else
		TMT(re)=TMT(rn)=TMDET_SIDE2;
}

void tm_splitHelix(pdbRes re, pdbRes rk, pdbRes rn) {
pdbRes r;

	for (r=re; r!=rn; r=r->next)
		TMT(r)=TMDET_ALPHA;
	if (rk->hbe[0]<0)
		TMT(rk)=TMDET_SIDE1;
	else
		TMT(rk)=TMDET_SIDE2;
}

void tm_countBeta(pdbProtein p)
{
pdbChain c;
pdbRes r,re,sr,to;
int clust,cnum,count[100],sheetcount,i,j,max,nclust;
tmdetTypes t,t1;


	for (c=p->fch; c!=NULL; c=c->next) if (c->sel)
	{	for (r=c->fres; r!=NULL;r=r->next)
			if (TMT(r)==TMDET_LOOP)
				TMT(r)=TMDET_INSIDE;
		cnum=0;
		for (r=c->fres; r!=NULL; r=r->next)
		{	if (TMT(r)==TMDET_BETA)
	      {	if (r->ss!='E')
	      		TMT(r)=TMDET_INSIDE;
	      }
	      if (TMT(r)==TMDET_LOOP||TMT(r)==TMDET_INSIDE)
	      {	if (r->ss=='E')
	      		TMT(r)=TMDET_BETA;
	      	else
	      		TMT(r)=TMDET_INSIDE;
	      }
		}
		for (r=c->fres; r!=NULL;)
		{	for (re=r->next; (re!=NULL&&TMT(re)==TMT(r)); re=re->next);
			if (TMT(r)==TMDET_BETA)
	      {	clust=-1;
	      	for (sr=r; sr!=re; sr=sr->next)
	      		if (TMB(sr)!=-1)
	      				clust=TMB(sr);
	      	if (clust==-1)
	      		clust=cnum++;
	      	for (sr=r; sr!=re; sr=sr->next)
	      		TMB(sr)=clust;
	      	nclust=-1;
	      	for (sr=r; sr!=re; sr=sr->next)
	      	{	if (sr->hb[0]!=-1)
	      		{	to=c->res[sr->hb[0]];
	      			if (TMB(to)!=-1) nclust=TMB(to);
	      		}
	      		if (sr->hb[1]!=-1)
	      		{	to=c->res[sr->hb[1]];
	      			if (TMB(to)!=-1) nclust=TMB(to);
	      		}
	      	}
	      	if (nclust!=-1)
	      	if (nclust!=clust)
	      	{	for (sr=c->fres; sr!=NULL; sr=sr->next)
	      		if (TMB(sr)==clust) TMB(sr)=nclust;
	      		clust=nclust;
	      	}
	      	for (sr=r; sr!=re; sr=sr->next)
	      	{	if (sr->hb[0]!=-1)
	      		{	to=c->res[sr->hb[0]];
	      			TMB(to)=clust;
	      		}
	      		if (sr->hb[1]!=-1)
	      		{	to=c->res[sr->hb[1]];
	      			TMB(to)=clust;
	      		}
	      	}
	      }
        	r=re;
		}
		if (cnum>0)
		{	for (i=0; i<cnum; i++) count[i]=0;
			for (r=c->fres; r!=NULL; r=r->next)
				if (TMB(r)!=-1) count[TMB(r)]++;
			max=0; j=-1;
			for (i=0; i<cnum; i++)
			{	if (max<count[i])
				{	max=count[i]; j=i;}
			}
			sheetcount=0;
			for (r=c->fres; r!=NULL;)
			{	for (re=r->next; (re!=NULL&&TMT(re)==TMT(r)); re=re->next);
				if (TMT(r)==TMDET_BETA)
				{	clust=-1;
					for (sr=r; sr!=re; sr=sr->next)
						if (TMB(sr)!=-1)
							clust=TMB(sr);
					if (clust!=j)
					{	for (sr=r; sr!=re; sr=sr->next) TMT(sr)=TMDET_INSIDE;}
					else
						sheetcount++;
				}
				r=re;
			}
			if (sheetcount>=8)
			{	printf("Chain \'%c\' has a %d-stranded beta barrel\n",
					c->id,sheetcount);
			}
		}
		for (r=c->fres; r!=NULL;)
		{	for (re=r->next,i=1; (re!=NULL&&TMT(re)==TMT(r)); re=re->next, i++);
			if (TMT(r)==TMDET_INSIDE)
			{	t=t1=TMDET_UNKNOWN;
				if (re!=NULL)
				if (i<=8)
				{	if (r->prev!=NULL)
						if (TMT(r->prev)==TMDET_BETA)
							t=TMT(r->prev);
					if (re!=NULL)
						if (TMT(re)==TMDET_BETA)
							t1=TMT(re);
					if (t==TMDET_BETA||t1==TMDET_BETA)
					{	if (i>2&&(t==t1||re==NULL))
							t=TMDET_INSIDE;
						else
							t=TMDET_BETA;
					}
				}
				if (t!=TMDET_UNKNOWN) {
					for (sr=r; sr!=re; sr=sr->next)
						TMT(sr)=t;
				}
			}
			r=re;
		}
		for (r=c->fres; r!=NULL;)
		{	for (re=r->next,i=1; (re!=NULL&&TMT(re)==TMT(r)); re=re->next, i++);
			if (TMT(r)==TMDET_BETA)
			{	if (i<3)
				{	t=TMDET_UNKNOWN;
					if (r->prev!=NULL)
						if (TMT(r->prev)==TMDET_SIDE1||TMT(r->prev)==TMDET_SIDE2)
							t=TMT(r->prev);
					if (re!=NULL)
						if (TMT(re)==TMDET_SIDE1||TMT(re)==TMDET_SIDE2)
							t=TMT(re);
					if (t!=TMDET_UNKNOWN)
						for (sr=r; sr!=re; sr=sr->next)
							TMT(sr)=t;
				}
				if (re!=NULL)
				if (re->prev!=NULL) {
					if (fabs(r->hbe[0]-re->prev->hbe[0])<3.0) {
						printf("Invalid beta loop: %d %d\n",r->rsn,re->prev->rsn);
						for (sr=r; sr!=re; sr=sr->next)
							TMT(sr)=TMDET_INSIDE;
					}
				}
			}
			r=re;
		}
	}
}
