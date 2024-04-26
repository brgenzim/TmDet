#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "pdb.h"
#include "protein.h"
#include "tmall.h"



matrix4 inverseMat(matrix4 m)
{
matrix4 i;

   i.el[0][0]=m.el[0][0];i.el[0][1]=m.el[1][0];i.el[0][2]=m.el[2][0];i.el[0][3]=-m.el[0][3];
   i.el[1][0]=m.el[0][1];i.el[1][1]=m.el[1][1];i.el[1][2]=m.el[2][1];i.el[1][3]=-m.el[1][3];
   i.el[2][0]=m.el[0][2];i.el[2][1]=m.el[1][2];i.el[2][2]=m.el[2][2];i.el[2][3]=-m.el[2][3];
   i.el[3][0]=0;         i.el[3][1]=0;         i.el[3][2]=0;         i.el[3][3]=0;

   return i;
}

matrix4 matrix_multiply(matrix4 a, matrix4 b)
{
matrix4 c;

    int i,j,k;
    
    for (i=0;i<4;i++)
    {
        for (j=0;j<4;j++)
	{
	     c.el[i][j] = 0.0;
	     for (k=0;k<4;k++)
	     {
	         c.el[i][j] += a.el[i][k]*b.el[k][j];
	     }
	}
    }
    
    return c;
}



/* Rotates the Molecule that Z axis will be T */
matrix4 Rotate_Z( POINT T)
{
float  a, b, c, d, cosa, sina;
matrix4 R1, R2, R;
POINT new;
    
    new = vec_norm(T, 1.0);  

    a=new.x; b=new.y; c=new.z;
      
    d=sqrt( b*b + c*c);
    if (d!=0)     
    {
        sina=c/d; cosa=b/d;
    /*printf("abcd %10.4f %10.4f %10.4f %10.4f \n",a,b,c,d);*/
    
    R1.el[0][0]=1;   R1.el[0][1]=0;    R1.el[0][2]=0;	  R1.el[0][3]=0;
    R1.el[1][0]=0;   R1.el[1][1]=sina; R1.el[1][2]=cosa;  R1.el[1][3]=0;
    R1.el[2][0]=0;   R1.el[2][1]=-cosa; R1.el[2][2]=sina; R1.el[2][3]=0;
    R1.el[3][0]=0;   R1.el[3][1]=0;    R1.el[3][2]=0;	  R1.el[3][3]=1;
    
             
    R2.el[0][0]=d;   R2.el[0][1]=0;    R2.el[0][2]=a;     R2.el[0][3]=0;
    R2.el[1][0]=0;   R2.el[1][1]=1;    R2.el[1][2]=0;	  R2.el[1][3]=0;
    R2.el[2][0]=-a;  R2.el[2][1]=0;    R2.el[2][2]=d;	  R2.el[2][3]=0;
    R2.el[3][0]=0;   R2.el[3][1]=0;    R2.el[3][2]=0;	  R2.el[3][3]=1;

    R=matrix_multiply(R1,R2);
    }
    else
    {
         R.el[0][0]=0;R.el[0][1]=1;R.el[0][2]=0;R.el[0][3]=0;
         R.el[1][0]=0;R.el[1][1]=0;R.el[1][2]=1;R.el[1][3]=0;
         R.el[2][0]=1;R.el[2][1]=0;R.el[2][2]=0;R.el[2][3]=0;
         R.el[3][0]=0;R.el[3][1]=0;R.el[3][2]=0;R.el[3][3]=1;
    
    }
            
    return R;
}

POINT midpoint(POINT a, POINT b)
{
POINT mp;

    mp.x = 0.5*(a.x + b.x );
    mp.y = 0.5*(a.y + b.y );
    mp.z = 0.5*(a.z + b.z );
    
    return mp;
}

/* ---------------- vec_diff_2 --------------------------------
 * Return the squared difference between two POINTS'
 */
float vec_diff_2 ( POINT a, POINT b)
{
    float x, y, z;
    x = a.x - b.x;
    x *= x;
    y = a.y - b.y;
    y *= y;
    z = a.z - b.z;
    z *= z;
    return ( x + y + z);
}

/* ---------------- vec_len   ---------------------------------
 */
float vec_len (POINT v)
{
    return sqrt (v.x * v.x + v.y * v.y + v.z * v.z);
}

/* ---------------- vec_scl   ---------------------------------
 */
POINT vec_scl (POINT v, float a)
{
    POINT tmp;
    tmp.x = v.x * a;
    tmp.y = v.y * a;
    tmp.z = v.z * a;
    return (tmp);
}

/* ---------------- vec_sub   ---------------------------------
 */
POINT vec_sub (POINT a, POINT b)
{
    POINT tmp;
    tmp.x = a.x - b.x;
    tmp.y = a.y - b.y;
    tmp.z = a.z - b.z;
    return (tmp);
}

/* ---------------- vec_add   ---------------------------------
 */
POINT vec_add (POINT a, POINT b)
{
    POINT tmp;
    tmp.x = a.x + b.x;
    tmp.y = a.y + b.y;
    tmp.z = a.z + b.z;
    return (tmp);
}

/* ---------------- vec_norm   ---------------------------------
 */
POINT vec_norm (POINT v, float r)
{
    float scl = r / vec_len (v);
    return ( vec_scl (v, scl));
}

/* ---------------- vec_prd   ---------------------------------
 */
POINT vec_prd (POINT u, POINT v)
{
    POINT tmp;
    tmp.x = u.y * v.z - u.z * v.y;
    tmp.y = u.z * v.x - u.x * v.z;
    tmp.z = u.x * v.y - u.y * v.x;
    return (tmp);
}

/* ---------------- vec_dot   ---------------------------------
 */
float vec_dot (POINT a, POINT b)
{
     float vd;
     vd = a.x*b.x + a.y*b.y + a.z*b.z;
     /*printf("vd %10.4f %10.4f %10.4f %10.4f\n",a.x,a.y,a.z,vd);*/
     return vd;

}

/* ---------------- cosinus angle between two vector   --------
 */
float vec_cos (POINT a, POINT b)
{
     float vd, la,lb, cos;
     vd = a.x*b.x + a.y*b.y + a.z*b.z;
     la=vec_len(a);lb=vec_len(b);
     cos = vd/(la*lb);
     /*printf("cos %10.4f \n",cos);*/
     return cos;

}


