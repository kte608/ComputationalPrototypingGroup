/*!\page LICENSE LICENSE
 
Copyright (C) 2003 by the Board of Trustees of Massachusetts Institute of Technology, hereafter designated as the Copyright Owners.
 
License to use, copy, modify, sell and/or distribute this software and
its documentation for any purpose is hereby granted without royalty,
subject to the following terms and conditions:
 
1.  The above copyright notice and this permission notice must
appear in all copies of the software and related documentation.
 
2.  The names of the Copyright Owners may not be used in advertising or
publicity pertaining to distribution of the software without the specific,
prior written permission of the Copyright Owners.
 
3.  THE SOFTWARE IS PROVIDED "AS-IS" AND THE COPYRIGHT OWNERS MAKE NO
REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED, BY WAY OF EXAMPLE, BUT NOT
LIMITATION.  THE COPYRIGHT OWNERS MAKE NO REPRESENTATIONS OR WARRANTIES OF
MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT THE USE OF THE
SOFTWARE WILL NOT INFRINGE ANY PATENTS, COPYRIGHTS TRADEMARKS OR OTHER
RIGHTS. THE COPYRIGHT OWNERS SHALL NOT BE LIABLE FOR ANY LIABILITY OR DAMAGES
WITH RESPECT TO ANY CLAIM BY LICENSEE OR ANY THIRD PARTY ON ACCOUNT OF, OR
ARISING FROM THE LICENSE, OR ANY SUBLICENSE OR USE OF THE SOFTWARE OR ANY
SERVICE OR SUPPORT.
 
LICENSEE shall indemnify, hold harmless and defend the Copyright Owners and
their trustees, officers, employees, students and agents against any and all
claims arising out of the exercise of any rights under this Agreement,
including, without limiting the generality of the foregoing, against any
damages, losses or liabilities whatsoever with respect to death or injury to
person or damage to property arising from or out of the possession, use, or
operation of Software or Licensed Program(s) by LICENSEE or its customers.
 
*/


#include "mulGlobal.h"

int num_calcp_grid; 
#define NGRID3 ((gridsize+1)*(gridsize+1)*(gridsize+1))
#define NGRID (gridsize+1)
#define DGRID (2.0 / (double) gridsize)

double llength;
double *Qwts; 

double **Q2G(qchgs, numqchgs, xc, yc, zc, Colloc, csize, invG2C, calc, gridsize)
charge **qchgs;
int numqchgs, csize;
double xc, yc, zc;
double **Colloc, **invG2C;
int calc;
int gridsize; 
{
  static double *p = NULL;
  double **q2g;
  int i, j, k;
  double x, y, z;
  double calcp(); 
  double dtmp; 
  double tilelength(); 

  if(p == NULL) {
    CALLOC(p,NGRID3,double,ON,AQ2M); 
  }

  CALLOC(q2g,NGRID3,double *, ON, AQ2M); 
  for(i=0; i < NGRID3; i++) 
    CALLOC(q2g[i],numqchgs,double, ON, AQ2M); 
  for(i = 0; i < numqchgs; i++) {
      for(j=0; j < csize; j++) {
	  x = Colloc[j][0] + xc;
	  y = Colloc[j][1] + yc;
	  z = Colloc[j][2] + zc;
	  p[j] = calcp(qchgs[i], x, y, z, NULL);
          num_calcp_grid++; 
      }
      for (j=csize; j < NGRID3; j++) {
        p[j] = 0.0; 
      }
      for(j=0; j < NGRID3; j++) {
	  for(k=0; k < csize; k++) {
            q2g[j][i] += invG2C[j][k] * p[k]; 
	  }
      }
    dtmp = 0.0; 
    for (j=0; j<NGRID3; j++) {
      dtmp += q2g[j][i] * q2g[j][i]; 
    }
    dtmp = sqrt(dtmp); 

  }

  return(q2g);
}

precorrect(q2p, dq2g, sq2g, numd, nums, g2g, gridsize, length)
double **q2p, **dq2g, **sq2g, **g2g, length;
int numd, nums, gridsize;
{
  int i, j, k, l;
  double gp;
  FILE *fp, *fopen(); 

  /* Compute q2p - dq2p^t g2g sq2g. */

  for(j=0; j < nums; j++) {
      for(l=0; l < NGRID3; l++) { 
	  for(gp = 0, k=0; k < NGRID3; k++) {
	      gp += g2g[l][k] * sq2g[k][j]; 
	  }
	  for(i=0; i < numd; i++) {
	      q2p[i][j] -= dq2g[l][i] * gp;
	  }
      }
  }
}

static double radstor[6] = {11, 11, 3.5, 3.5, 3.0, 2.6}; 

double **initColloc(pcsize, length, gridsize)
int *pcsize, gridsize;
double length;
{
  int index, i, j, k;
  double **col; 
  int on_surface(); 
  double t;
  double r, s, u, v, x0,x1,x2; 
  double A, B, C; 
  double rad; 

  CALLOC(col,NGRID3,double *, ON, AQ2M); 

/* get rad for colloc sphere */ 
  if (gridsize < 7)
    rad = 0.5 * length * radstor[gridsize-1]; 
  else
    rad = 0.5 * length * 2.3; 

  if (gridsize == 1 ) { 
    index = 0;
    for(i= 0; i < NGRID; i++) {
      for(j= 0; j < NGRID; j++) {
        for(k= 0; k < NGRID; k++) {
          if(on_surface(i,j,k,NGRID)) {
            CALLOC(col[index],3,double,ON,AQ2M); 
            col[index][0] = rad * (-1.0 + (double)i * DGRID); 
            col[index][1] = rad * (-1.0 + (double)j * DGRID); 
            col[index][2] = rad * (-1.0 + (double)k * DGRID); 
            index++;
          }
        }
      }
    }
  }
  else if (gridsize == 2) {
    CALLOC(Qwts, 26, double, ON, AMSC); 
    for (i=0; i < 26; i++) {
      CALLOC(col[i],3,double,ON,AMSC); 
    }
    r = 1.0; 
    s = 0.5; 
    t = 1.0/3.0; 

    s = sqrt(s); t = sqrt(t); 
    index = 0;
    for (i=0,x0=t; i<2; i++,x0*=-1.0) {
      for (j=0,x1=t; j<2; j++,x1*=-1.0) {
        for (k=0,x2=t; k<2; k++,x2*=-1.0) {
          col[index][0] = x0; 
          col[index][1] = x1; 
          col[index][2] = x2;
          index++; 
        }
      }
    }
    
    for (i=0,x0=s; i<2; i++,x0*=-1.0) {
      for (j=0,x1=s; j<2; j++,x1*=-1.0) {
          col[index][0] = x0; 
          col[index][1] = x1; 
          col[index][2] = 0.0;
          index++; 
      }
    }
    for (i=0,x0=s; i<2; i++,x0*=-1.0) {
      for (j=0,x1=s; j<2; j++,x1*=-1.0) {
          col[index][0] = x0; 
          col[index][1] = 0.0; 
          col[index][2] = x1;
          index++; 
      }
    }
    for (i=0,x0=s; i<2; i++,x0*=-1.0) {
      for (j=0,x1=s; j<2; j++,x1*=-1.0) {
          col[index][0] = 0.0; 
          col[index][1] = x0; 
          col[index][2] = x1;
          index++; 
      }
    }


    col[index][0] = r; col[index][1] = 0.0; col[index][2] = 0.0; index++; 
    col[index][0] = -r; col[index][1] = 0.0; col[index][2] = 0.0; index++; 
    col[index][0] = 0.0; col[index][1] = r; col[index][2] = 0.0; index++; 
    col[index][0] = 0.0; col[index][1] = -r; col[index][2] = 0.0; index++; 
    col[index][0] = 0.0; col[index][1] = 0.0; col[index][2] = r; index++; 
    col[index][0] = 0.0; col[index][1] = 0.0; col[index][2] = -r; index++;  

    for (i=0; i < 26; i++) {
      for (j=0; j<3; j++) {
        col[i][j] *= rad; 
      }
    }
  }
  else if (gridsize == 3) {
    CALLOC(Qwts, 56, double, ON, AMSC); 
    for (i=0; i < 56; i++) {
      CALLOC(col[i],3,double,ON,AMSC); 
    }
    r = (9. - 4. * sqrt(3.))/33.0; 
    s = (15.0 + 8 * sqrt(3.0))/33.0; 
    u = (15.0 - 8 * sqrt(3.0))/33.0; 
    v = (9 + 4*sqrt(3.0))/33.0; 
    t = sqrt(1.0/3.0); 
    A = 9.0/560.0; 
    B = (122+9.0*sqrt(3.0))/6720.0; 
    C = (122-9.0*sqrt(3.0))/6720.0; 

    r = sqrt(r); s = sqrt(s); u = sqrt(u); v = sqrt(v); 
    index = 0;
    for (i=0,x0=t; i<2; i++,x0*=-1.0) {
      for (j=0,x1=t; j<2; j++,x1*=-1.0) {
        for (k=0,x2=t; k<2; k++,x2*=-1.0) {
          col[index][0] = x0; 
          col[index][1] = x1; 
          col[index][2] = x2;
          Qwts[index] = A; 
          index++; 
        }
      }
    }
    x0 = u; x1 = v; x2 = v; 
    for (i=0; i<2; i++,x0*=-1.0) {
      for (j=0; j<2; j++,x1*=-1.0) {
        for (k=0; k<2; k++,x2*=-1.0) {
          col[index][0] = x0; 
          col[index][1] = x1; 
          col[index][2] = x2;
          Qwts[index] = C; 
          index++; 
        }
      }
    }

    x0 = v; x1 = u; x2 = v; 
    for (i=0; i<2; i++,x0*=-1.0) {
      for (j=0; j<2; j++,x1*=-1.0) {
        for (k=0; k<2; k++,x2*=-1.0) {
          col[index][0] = x0; 
          col[index][1] = x1; 
          col[index][2] = x2;
          Qwts[index] = C; 
          index++; 
        }
      }
    }
    x0 = v; x1 = v; x2 = u; 
    for (i=0; i<2; i++,x0*=-1.0) {
      for (j=0; j<2; j++,x1*=-1.0) {
        for (k=0; k<2; k++,x2*=-1.0) {
          col[index][0] = x0; 
          col[index][1] = x1; 
          col[index][2] = x2;
          Qwts[index] = C; 
          index++; 
        }
      }
    }
    x0 = s; x1 = r; x2 = r; 
    for (i=0; i<2; i++,x0*=-1.0) {
      for (j=0; j<2; j++,x1*=-1.0) {
        for (k=0; k<2; k++,x2*=-1.0) {
          col[index][0] = x0; 
          col[index][1] = x1; 
          col[index][2] = x2;
          Qwts[index] = B; 
          index++; 
        }
      }
    }
    x0 = r; x1 = s; x2 = r; 
    for (i=0; i<2; i++,x0*=-1.0) {
      for (j=0; j<2; j++,x1*=-1.0) {
        for (k=0; k<2; k++,x2*=-1.0) {
          col[index][0] = x0; 
          col[index][1] = x1; 
          col[index][2] = x2;
          Qwts[index] = B; 
          index++; 
        }
      }
    }
    x0 = r; x1 = r; x2 = s; 
    for (i=0; i<2; i++,x0*=-1.0) {
      for (j=0; j<2; j++,x1*=-1.0) {
        for (k=0; k<2; k++,x2*=-1.0) {
          col[index][0] = x0; 
          col[index][1] = x1; 
          col[index][2] = x2;
          Qwts[index] = B; 
          index++; 
        }
      }
    }
    for (i=0; i < 56; i++) {
      for (j=0; j<3; j++) {
        col[i][j] *= rad; 
      }
    }
  }
  else if (gridsize == 4) { 
    {
      double z[7], u[7], v[7], w[7]; 
      CALLOC(Qwts, 72, double, ON, AMSC); 
      for (i=0; i < 72; i++) {
        CALLOC(col[i],3,double,ON,AMSC); 
      }
      r = ((5.0-sqrt(5.0))/10.0); r = sqrt(r); 
      s = ((5.0+sqrt(5.0))/10.0); s = sqrt(s); 
      z[1] = 0.91206379026291; 
      z[2] = 0.74883416366819; 
      z[3] = 0.64178606967939; 
      z[4] = 0.38470531769520; 
      z[5] = 0.21149786319039; 
      z[6] = 0.05261322061079; 

      u[1] = (-z[3]+z[4])/(2.0*s); 
      u[2] = (-z[5]+z[2])/(2.0*s); 
      u[3] = (-z[2]+z[6])/(2.0*s); 
      u[4] = (-z[6]+z[3])/(2.0*s); 
      u[5] = (-z[4]+z[5])/(2.0*s); 

      v[1] = (z[5]+z[6])/(2.0*s); 
      v[2] = (z[6]+z[4])/(2.0*s); 
      v[3] = (z[3]+z[5])/(2.0*s); 
      v[4] = (z[4]+z[2])/(2.0*s); 
      v[5] = (z[2]+z[3])/(2.0*s); 

      w[1] = (z[1]+z[2])/(2.0*s); 
      w[2] = (z[1]+z[3])/(2.0*s); 
      w[3] = (z[1]+z[4])/(2.0*s); 
      w[4] = (z[1]+z[5])/(2.0*s); 
      w[5] = (z[1]+z[6])/(2.0*s); 

      index = 0;
      col[index][0] = r; col[index][1] = s; col[index][2] = 0.0; index++; 
      col[index][0] = r; col[index][1] = -s; col[index][2] = 0.0; index++; 
      col[index][0] = -r; col[index][1] = s; col[index][2] = 0.0; index++; 
      col[index][0] = -r; col[index][1] = -s; col[index][2] = 0.0; index++; 

      col[index][0] = 0.0; col[index][1] = r; col[index][2] = s; index++; 
      col[index][0] = 0.0; col[index][1] = r; col[index][2] = -s; index++; 
      col[index][0] = 0.0; col[index][1] = -r; col[index][2] = s; index++; 
      col[index][0] = 0.0; col[index][1] = -r; col[index][2] = -s; index++; 

      col[index][0] = s; col[index][1] = 0.0; col[index][2] = r; index++; 
      col[index][0] = s; col[index][1] = 0.0; col[index][2] = -r; index++; 
      col[index][0] = -s; col[index][1] = 0.0; col[index][2] = r; index++; 
      col[index][0] = -s; col[index][1] = 0.0; col[index][2] = -r; index++; 

      for (i=1; i<=5; i++) {
        col[index][0] = u[i]; col[index][1] = v[i]; col[index][2] = w[i]; index++; 
        col[index][0] = u[i]; col[index][1] = -v[i]; col[index][2] = -w[i]; index++; 
        col[index][0] = -u[i]; col[index][1] = -v[i]; col[index][2] = w[i]; index++; 
        col[index][0] = -u[i]; col[index][1] = v[i]; col[index][2] = -w[i]; index++; 

        col[index][0] = v[i]; col[index][1] = w[i]; col[index][2] = u[i]; index++; 
        col[index][0] = v[i]; col[index][1] = -w[i]; col[index][2] = -u[i]; index++; 
        col[index][0] = -v[i]; col[index][1] = -w[i]; col[index][2] = u[i]; index++; 
        col[index][0] = -v[i]; col[index][1] = w[i]; col[index][2] = -u[i]; index++; 

        col[index][0] = w[i]; col[index][1] = u[i]; col[index][2] = v[i]; index++; 
        col[index][0] = w[i]; col[index][1] = -u[i]; col[index][2] = -v[i]; index++; 
        col[index][0] = -w[i]; col[index][1] = -u[i]; col[index][2] = v[i]; index++; 
        col[index][0] = -w[i]; col[index][1] = u[i]; col[index][2] = -v[i]; index++; 

      }
    }
    for (i=0; i < 72; i++) {
      for (j=0; j<3; j++) {
        col[i][j] *= rad; 
      }
    }
  }
  else {
    fprintf(stderr, "Grid order > 5 not implemented.\n"); 
    abort(); 
  }

  *pcsize = index;

  return(col);
}


int on_surface(i,j,k,ngrid)
int i,j,k,ngrid; 
{
  int flag; 
  if (i > 0 && i < ngrid-1 && j > 0 && j < ngrid-1 && k > 0 && k < ngrid-1) { 
    return(0); 
  }
  else {
    return(1); 
  }
}


Real **getGrid(centerj, centerk, centerl, data, gridsize,size1,size2,size3)
int centerj, centerk, centerl, gridsize;
Real ***data;
int size1,size2,size3;
{
  Real **grid;
  int j, k, l;
  
  CALLOC(grid,NGRID3,Real *, ON, AQ2M); 

  for(j= 0; j < NGRID; j++) {
    for(k= 0; k < NGRID; k++) {
      for(l= 0; l < NGRID; l++) {
        grid[gindex(j,k,l,gridsize)] = 
          &(data[centerj+j][centerk+k][centerl+l]);
      }
    }
  }
  return(grid);
}

double **invG2C(col, length, csize, gridsize, zoff)
int csize, gridsize;
double length, **col;
double zoff; 
{
  int i,j,k,l;
  double **invg2c, dx, dy, dz;
  static int *reorder=NULL;
  FILE *fp, *fopen(); 

  llength = length; 

  if(reorder == NULL) {
    CALLOC(reorder,NGRID3,int,ON,AMSC); 
    
  }
  CALLOC(invg2c,NGRID3,double *, ON, AQ2M); 
  
  for(i=0; i < NGRID3; i++) {
    CALLOC(invg2c[i],NGRID3,double,ON,AQ2M);
  }

  for(l=0; l < csize; l++) {
    for(i= 0; i < NGRID; i++) {
      for(j= 0; j <  NGRID; j++) {
	for(k= 0; k < NGRID; k++) {
	  dx = (col[l][0] - 0.5 * length * (-1.0 + i * DGRID));
          dy = (col[l][1] - 0.5 * length * (-1.0 + j * DGRID));
	  dz = (col[l][2] - 0.5 * length * (-1.0 + k * DGRID));
	  invg2c[l][gindex(i,j,k,gridsize)] = 
	    1.0/sqrt(dx * dx + dy * dy + dz * dz); 
	}
      }
    }
  }

  pinv(invg2c,csize,NGRID3); 

  return(invg2c);
}



double ***grid2grid(gridsize, length, zoff)
int gridsize;
double length;
double zoff; 
{
  double ***g2g;
  double r, delta = length/gridsize, dx, dy, dz;
  int i, j, k, id, jd, kd, is, js, ks, cubeindex, dindex;
  int cubes = (2 * NNBRS + 1) * (2 * NNBRS + 1)  * (2 * NNBRS + 1);

  CALLOC(g2g,cubes,double **, ON, AM2M); 
  for(i= -NNBRS; i <= NNBRS; i++) {
    for(j= -NNBRS; j <= NNBRS; j++) {
      for(k= -NNBRS; k <= NNBRS; k++) {
        cubeindex = cindex(i,j,k);
	CALLOC(g2g[cubeindex],NGRID3,double *, ON, AM2M); 
	for(id= 0; id < NGRID; id++) {
	  for(jd= 0; jd < NGRID; jd++) {
	    for(kd= 0; kd < NGRID; kd++) {
	      dindex = gindex(id,jd,kd,gridsize);
	      CALLOC(g2g[cubeindex][dindex],NGRID3,double,ON,AM2M); 
	      for(is= 0; is < NGRID; is++) {
		for(js= 0; js < NGRID; js++) {
		  for(ks= 0; ks < NGRID; ks++) {
		    dx = (i * gridsize + (id - is)) * delta;
		    dy = (j * gridsize + (jd - js)) * delta;
		    dz = (k * gridsize + (kd - ks)) * delta;
		    r = sqrt(dx * dx + dy * dy + dz * dz); 
		    if(r != 0.0) r = 1.0/r; 
                    

		    g2g[cubeindex][dindex][gindex(is,js,ks,gridsize)] = r;


#if GNDPLN == ON 
                    dz = k * gridsize + kd + ks; 
                    dz = delta * dz + 2.0 * zoff; 
		    r = sqrt(dx * dx + dy * dy + dz * dz); 
                    g2g[cubeindex][dindex][gindex(is,js,ks,gridsize)] += Rcoef/r; 

#endif 
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  return(g2g);
}


gindex(i,j,k,gridsize)
int i,j,k,gridsize;
{
  int limit;

  ASSERT(i >= 0 && i <= gridsize+1); 
  ASSERT(j >= 0 && j <= gridsize+1); 
  ASSERT(k >= 0 && k <= gridsize+1); 

  limit = i*(NGRID)*(NGRID) 
    + j*(NGRID) + k;
  return limit; 

}

cindex(i,j,k)
int i,j,k;
{
  int limit = NNBRS, size = 2 * NNBRS + 1;

  ASSERT((-limit <= i) && (i <= limit));
  ASSERT((-limit <= j) && (j <= limit));
  ASSERT((-limit <= k) && (k <= limit));
  return((i+limit)*(size*size) 
	 + (j+limit)*(size) + k+limit);
}

createFFT(pdata, pspeq, s1, s2, s3)
Real ****pdata, ***pspeq;
int s1, s2, s3;
{
    Real *datas, ***data, **speq;
    int i, j;

    CALLOC(data,s1+1,Real **, ON, AM2M); 
    CALLOC(datas, (s1+1)*(s2+1)*(s3+1), Real, ON, AM2M); 
    for(i=1; i <= s1; i++) {
	CALLOC(data[i], s2+1, Real *, ON, AMSC); 
	for(j=1; j <= s2; j++) 
	  data[i][j] = &(datas[1 + (i-1)*(s2 * s3) + (j-1)*s3]);
    }
    CALLOC(speq, s1+1, Real *, ON, AM2M); 
    CALLOC(datas, (s1+1)*(2*s2+2), Real, ON, AM2M); 
    for(i=1; i <= s1; i++) 
      speq[i] = &(datas[1+(i-1) * (2 * s2)]);
    *pspeq = speq;
    *pdata = data;
}



ft_reverse(data,speq,size1,size2,size3)
Real ***data, **speq;
int size1, size2, size3; 
{
/* compute the FFT of a sequence reversed about the Z-direction 
   from the transform of the original sequence 
   this does 
   (1) an implicit shift 
   (2) multiplication by the phase factor W 
   (3) implicit reversal; for real-data
       means take conjugate 
*/ 
  int i,j,k; 
  static double *W = NULL; 
  static int lastsize = 0; 
  double re, im, s,re2, im2; 
  if (size3 > lastsize) {
    if (W) free(W); 
    CALLOC(W, 2*size3, double, ON, AM2M);   
  }
  if (size3 != lastsize) {
    for (i=0; i< size3; i++) {
      k = size3 - i; 
      W[2*k-1] = cos(2.0*M_PI*(double)(i+1)/size3); 
      W[2*k] = sin(2.0*M_PI*(double)(i+1)/size3); 
    }
  }
  lastsize = size3;  

  /* the f3=fc component gets mult by -1 */ 
  for (i=1; i<=size1; i++) {
    for (j=1; j<=2*size2; j++) {
      speq[i][j] *= -1.0; 
    }
  }

  /* the f3=0 (DC) component is unchanged */ 
  s = -1.0; 
  for (k=2; k<= size3/2; k++) {

    for (i=1; i<=size1/2+1; i += size1/2) {
      for (j=1; j <=size2/2+1; j += size2/2) {
        re = W[2*k-1] * data[i][j][2*k-1] + W[2*k] * data[i][j][2*k]; 
        im = W[2*k] * data[i][j][2*k-1] - W[2*k-1] * data[i][j][2*k]; 

        data[i][j][2*k-1] = s * re; 
        data[i][j][2*k] = s * im; 
      }
    }


    for (i=1; i<=size1/2+1; i += size1/2) {
      for (j=2; j<=size2/2; j++) {
        re = W[2*k-1] * data[i][j][2*k-1] + W[2*k] * data[i][j][2*k]; 
        im = W[2*k] * data[i][j][2*k-1] - W[2*k-1] * data[i][j][2*k]; 

        re2 = W[2*k-1] * data[i][size2-j+2][2*k-1] + W[2*k] * data[i][size2-j+2][2*k]; 
        im2 = W[2*k] * data[i][size2-j+2][2*k-1] - W[2*k-1] * data[i][size2-j+2][2*k]; 
        data[i][j][2*k-1] = s * re2; 
        data[i][j][2*k] = s * im2; 
      
        data[i][size2-j+2][2*k-1] = s * re; 
        data[i][size2-j+2][2*k] = s * im; 

      }
    }
    for (j=1; j<=size2/2+1; j+= size2/2) {
      for (i=2; i<=size1/2; i++) {
        re = W[2*k-1] * data[i][j][2*k-1] + W[2*k] * data[i][j][2*k]; 
        im = W[2*k] * data[i][j][2*k-1] - W[2*k-1] * data[i][j][2*k]; 

        re2 = W[2*k-1] * data[size1-i+2][j][2*k-1] + W[2*k] * data[size1-i+2][j][2*k]; 
        im2 = W[2*k] * data[size1-i+2][j][2*k-1] - W[2*k-1] * data[size1-i+2][j][2*k]; 
        data[i][j][2*k-1] = s * re2; 
        data[i][j][2*k] = s * im2; 
      
        data[size1-i+2][j][2*k-1] = s * re; 
        data[size1-i+2][j][2*k] = s * im; 

      }
    }

    for (i=2; i<=size1; i++) {
      for (j=2; j<=size2/2; j++) {
        if (i != size1/2+1) {
        re = W[2*k-1] * data[i][j][2*k-1] + W[2*k] * data[i][j][2*k]; 
        im = W[2*k] * data[i][j][2*k-1] - W[2*k-1] * data[i][j][2*k]; 
        
        re2 = W[2*k-1] * data[size1-i+2][size2-j+2][2*k-1] + W[2*k] * data[size1-i+2][size2-j+2][2*k]; 
        im2 = W[2*k] * data[size1-i+2][size2-j+2][2*k-1] - W[2*k-1] * data[size1-i+2][size2-j+2][2*k]; 

        data[i][j][2*k-1] = s * re2; 
        data[i][j][2*k] = s * im2; 

        data[size1-i+2][size2-j+2][2*k-1] = s * re; 
        data[size1-i+2][size2-j+2][2*k] = s * im; 
        }
      }
    }
    /* now shift it by N/2 spaces to account for zero padding */ 
    s *= -1.0;  
  }

}

double **G2P_poly(qchgs, numqchgs, xc, yc, zc, length, csize, gridsize)
charge **qchgs;
int numqchgs, csize;
double xc, yc, zc;
double length;
{
  static double *p = NULL;
  double **g2p;
  int i, j;
  double x, y, z;
  double calcp(); 
  double xi, yi, zi, xk, tmp; 
  int i1,j1,k1,i2; 
  
  if(p == NULL) {
    CALLOC(p,NGRID3,double,ON,AQ2M); 
  }
  CALLOC(g2p,NGRID3,double *, ON, AQ2M); 
  for(i=0; i < NGRID3; i++) 
    CALLOC(g2p[i],numqchgs,double, ON, AQ2M); 

  for(i = 0; i < numqchgs; i++) {
      /* loop over grid points in cube */ 

      x = qchgs[i]->x - xc; y = qchgs[i]->y - yc; z = qchgs[i]->z - zc; 
      x = 2.0 * x / length; y = 2.0 * y / length;  z = 2.0 * z / length; 

      for (i1 = 0; i1 < NGRID; i1++) { 
        for (j1 = 0; j1 < NGRID; j1++) { 
          for (k1 = 0; k1 < NGRID; k1++) { 
            j = gindex(i1,j1,k1, gridsize); 
            /* j is the grid index, i the charge index */ 
            /* put in normalized cube coords */ 

            xi = (double) i1; yi = (double) j1; zi = (double) k1; 
            xi = (-1.0 + i1 * DGRID);
            yi = (-1.0 + j1 * DGRID);
            zi = (-1.0 + k1 * DGRID);
            tmp = 1.0; 
            for (i2 = 0; i2 < NGRID; i2 ++) { 
              xk = (-1 + i2 * DGRID); 
              if (i1 != i2) {
                tmp *= (x-xk) / (xi-xk); 
              }
            }
            for (i2 = 0; i2 <  NGRID; i2 ++) { 
              xk = (-1 + i2 * DGRID); 
              if (j1 != i2) {
                tmp *= (y-xk) / (yi-xk); 
              }
            }
            for (i2 = 0; i2 < NGRID; i2 ++) { 
              xk = (-1 + i2 * DGRID); 
              if (k1 != i2) {
                tmp *= (z-xk) / (zi-xk); 
              }
            }
            g2p[j][i] = tmp; 
          }
        }
      }

    }
  return(g2p); 
}

double **G2P(qchgs, numqchgs, xc, yc, zc, Colloc, csize, invG2C, calc, gridsize)
charge **qchgs;
int numqchgs, csize;
double xc, yc, zc;
double **Colloc;
double **invG2C;
int calc;
int gridsize; 
{
  static double *p = NULL;
  double **g2p;
  int i, j, k;
  double x, y, z;
  double calcp(); 
  double xn,yn,zn; 
  double R; 
  if(p == NULL) {
      CALLOC(p,NGRID3,double,ON,AMSC); 
  }
  CALLOC(g2p,NGRID3,double *, ON, AMSC); 
  for(i=0; i <= NGRID3; i++) 
    CALLOC(g2p[i],numqchgs,double, ON, AMSC); 
  for(i = 0; i < numqchgs; i++) {
    for (j=0; j <NGRID3; j++) { 
      g2p[j][i] = 0.0; 
    }
    for(j=0; j < csize; j++) {
      x = Colloc[j][0] + xc;
      y = Colloc[j][1] + yc;
      z = Colloc[j][2] + zc;
      xn = qchgs[i]->x - x; 
      yn = qchgs[i]->y - y; 
      zn = qchgs[i]->z - z; 
      R = sqrt(xn*xn+yn*yn+zn*zn); 
      if (R > 0.0)
        p[j] = 1.0 / R; 
      else 
        p[j] = 0.0; 
    }
    for (j=csize; j < NGRID3; j++) {
      p[j] = 0.0; 
    }
    for(j=0; j < NGRID3; j++) {
      for(k=0; k < csize; k++) {
        g2p[j][i] += invG2C[j][k] * p[k]; 
      }
    }
  }
  return(g2p);
}
