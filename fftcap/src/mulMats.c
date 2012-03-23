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

double **Q2P(), **Q2PDiag();
double **mulMulti2P(), **mulQ2Multi(), **mulMulti2Multi();
double **mulLocal2Local(), **mulLocal2P(), **mulQ2Local(), **mulMulti2Local();

int *localcnt, *multicnt, *evalcnt;	/* counts of builds done by level */

int **Q2Mcnt, **Q2Lcnt, **Q2Pcnt, **L2Lcnt; /* counts of xformation mats */
int **M2Mcnt, **M2Lcnt, **M2Pcnt, **L2Pcnt, **Q2PDcnt;

/*
MulMatDirect creates the matrices for the piece of the problem that is done
directly exactly.
*/
mulMatDirect(sys)
ssystem *sys;
{
  cube *nextc, *nextnbr;
  int i, nummats, **temp;
  extern double lutime, dirtime;

  /* First count the number of matrices to be done directly. */
  for(nextc=sys->directlist; nextc != NULL; nextc = nextc->dnext) {
    for(nummats=1, i=0; i < nextc->numnbrs; i++) {
      nextnbr = nextc->nbrs[i];
      ASSERT(nextnbr->upnumvects > 0);
      nummats++;
    }

  /* Allocate space for the vects and mats. */
    nextc->directnumvects = nummats;
    if(nummats > 0) {
      CALLOC(nextc->directq, nummats, double*, ON, AQ2P);
      CALLOC(temp, nummats, int*, ON, AQ2P);
      CALLOC(nextc->directnumeles, nummats, int, ON, AQ2P);
      CALLOC(nextc->directmats, nummats, double**, ON, AQ2P);
    }
  }

/* Now place in the matrices. */
  for(nextc=sys->directlist; nextc != NULL; nextc = nextc->dnext) {
    nextc->directq[0] = nextc->upvects[0];
    nextc->directnumeles[0] = nextc->upnumeles[0];

    starttimer;
    nextc->directmats[0] 
	= Q2PDiag(nextc->chgs, nextc->upnumeles[0], TRUE);

        stoptimer;
    dirtime += dtime;

#if DSQ2PD == ON
    dumpQ2PDiag(nextc);
#endif

#if DMTCNT == ON
    Q2PDcnt[nextc->level][nextc->level]++;
#endif

    starttimer;
    for(nummats=1, i=0; i < nextc->numnbrs; nummats++, i++) {
      nextnbr = nextc->nbrs[i];
      ASSERT(nextnbr->upnumvects > 0);
      nextc->directq[nummats] = nextnbr->upvects[0];
      nextc->directnumeles[nummats] = nextnbr->upnumeles[0];
      nextc->directmats[nummats] = Q2P(nextnbr->chgs, 
				       nextnbr->upnumeles[0], 
				       nextc->chgs, nextc->upnumeles[0],
				       TRUE);

      
#if DMTCNT == ON
      Q2Pcnt[nextc->level][nextnbr->level]++;
#endif
    }
    stoptimer;
    dirtime += dtime;
  }
}

/* In this routine, the 3x3x3 grid per cube is built in. */
mulMatGrid(sys)
ssystem *sys;
{
  int i, j, k, l, dj, dk, dl, gridj, gridk, gridl, size1, size2, size3;
  int maxj=0, minj=sys->side, maxk=0, mink=sys->side, maxl=0, minl=sys->side;
  int centerj, centerk, centerl;
  int gridsize, nummats;
  cube *nextc, *nextnbr;
  Real ***kdata, ***data, **kspeq, **speq, *datas, *speqs;
  double dx, dy, dz, length;
  double **Q2G(), **invG2C(), **initColloc(), ***grid2grid(), **G2P_poly(), **G2P();
  Real **getGrid();
  FILE *fp, *fopen(); 
  double fac, *sp1, *sp2; 
  double r,im; 
  double zoff; 
  int xindex, ncub; 
  extern int num_cond_panels; 

  gridsize = sys->gridsize; 

  /* Find smallest power of 2 grid that covers the domain. */
  for(nextc=sys->directlist; nextc != NULL; nextc = nextc->dnext) {
    maxj = MAX(maxj, nextc->j);
    minj = MIN(minj, nextc->j);
    maxk = MAX(maxk, nextc->k);
    mink = MIN(mink, nextc->k);
    maxl = MAX(maxl, nextc->l);
    minl = MIN(minl, nextc->l);
  }
  sys->minj = minj; sys->mink = mink; sys->minl = minl; 

  /* count # of cubes with charge */ 
  xindex = 0; 
  for (j=minj; j<=maxj; j++) {
    for (k=mink; k<= maxk; k++) {
      for (l=minl; l<=maxl; l++) {
        if (nextc = sys->cubes[sys->depth][j][k][l]) {
          if (nextc->upnumvects) {
            xindex++; 
          }
        }
      }
    }
  }
  /* # of cubes in the domain */ 
  ncub = (maxj-minj+1) * (maxk-mink+1)*(maxl-minl+1); 

  dj = (maxj - minj +1) * gridsize + 1;
  dk = (maxk - mink +1) * gridsize + 1;
  dl = (maxl - minl +1) * gridsize + 1;
  for(gridj=1; gridj < dj; gridj *= 2);
  for(gridk=1; gridk < dk; gridk *= 2);
  for(gridl=1; gridl < dl; gridl *= 2);

  size1 = 2 * gridj; size2 = 2 * gridk; size3 = 2 * gridl;

  /* Do the FFT on a grid twice as long (padding to avoid aliasing). */
  sys->sizej = size1;
  sys->sizek = size2;
  sys->sizel = size3;

  /* Create the data for the fft of the kernel. */
  createFFT(&(sys->kdata), &(sys->kspeq), size1, size2, size3);

  /* Create the data for the fft of the grid potentials. */
  createFFT(&(sys->data), &(sys->speq), size1, size2, size3/2);

#if GNDPLN == ON 
/* data for Hankel-like part of kernel */ 

  createFFT(&(sys->k2data), &(sys->k2speq), size1, size2, size3); 
  createFFT(&(sys->data2), &(sys->speq2), size1, size2, size3); 
#endif 

  /* Get collocation points. */
  sys->colloc = initColloc(&(sys->csize), sys->length, sys->gridsize);

  /* Create and Invert the grid to collocation map. */
  CALLOC(sys->invg2c, maxl+1, double **, ON, AQ2M); 
  CALLOC(sys->g2g, maxl+1, double ***, ON, AM2M); 
  sys->invg2c[0] = invG2C(sys->colloc, sys->length, sys->csize, sys->gridsize,0.0);

#if GNDPLN == ON 
  for (l = minl; l <=maxl; l++) {
    zoff = GroundPlaneHeight*sys->length + ((double) l) * sys->length; 
    sys->g2g[l] = grid2grid(gridsize, sys->length, zoff);
  }
#else 
  /* Create the grid to grid mapping matrices. */
  sys->g2g[0] = grid2grid(gridsize, sys->length,0.0);
  for (l = minl; l <=maxl; l++) {
    sys->g2g[l] = sys->g2g[0]; 
  }
#endif 

  /* Calculate the Fourier transform of the kernel. */
  dx = sys->length/sys->gridsize;
  dy = sys->length/sys->gridsize;
  dz = sys->length/sys->gridsize;

  kdata = sys->kdata;
  centerj = size1/2;
  centerk = size2/2;
  centerl = size3/2;

  kdata[centerj][centerk][centerl] = 0.0;
  for(j = 1-centerj; j <= centerj; j++) {
    for(k = 1-centerk; k <= centerk; k++) {
      for(l = 1-centerl; l <= centerl; l++) {    
	if((j != 0) || (k != 0) || (l !=0 ))  {
	  length = 1.0/sqrt((j*j*dx*dx) + (k*k*dy*dy) + (l*l*dz*dz));
	}
	else {
	  length = 0.0;
	}

#define CIRC(xi, size) ( ( (xi) >= 0) ? (1+(xi)) : (1+(size)+(xi)) ) 
	kdata[CIRC(j,size1)][CIRC(k,size2)][CIRC(l,size3)] = length;  
      }
    }
  }

  rfft3(sys->kdata,sys->kspeq,size1,size2,size3);     

#if GNDPLN == ON
  kdata = sys->k2data;
  centerj = size1/2;
  centerk = size2/2;
  centerl = size3/2;
  kdata[centerj][centerk][centerl] = 0.0;
  for(j = 1-centerj; j <= centerj; j++) {
    for(k = 1-centerk; k <= centerk; k++) {
      for(l = 1-centerl; l <= centerl; l++) {     
        if (l >= 0) 
          length = sqrt((j*j*dx*dx) + (k*k*dy*dy) + dz * dz * 
                        (l + size3/2 +1.0 -2.0 + 2.0 * (sys->gridsize)*GroundPlaneHeight) *
                        (l + size3/2 +1.0 -2.0 + 2.0 * (sys->gridsize)*GroundPlaneHeight)) ; 
        else 
          length = sqrt((j*j*dx*dx) + (k*k*dy*dy) + dz * dz * 
                        (l +  size3/2 +1.0  -2.0 + 2.0 * (sys->gridsize)*GroundPlaneHeight) * 
                        (l +  size3/2 +1.0  -2.0 + 2.0 * (sys->gridsize)*GroundPlaneHeight)); 
#define CIRC(xi, size) ( ( (xi) >= 0) ? (1+(xi)) : (1+(size)+(xi)) ) 
	kdata[CIRC(j,size1)][CIRC(k,size2)][CIRC(l,size3)] = Rcoef/length;  
      }
    }
  }
  rfft3(sys->k2data,sys->k2speq,size1,size2,size3); 

#endif 

  /* Compute the panel to grid matrices. */
  for(nextc=sys->directlist; nextc != NULL; nextc = nextc->dnext) {
    nextc->directq[0] = nextc->upvects[0];
    nextc->q2g = Q2G(nextc->chgs, nextc->upnumeles[0], nextc->x, 
		     nextc->y, nextc->z, sys->colloc, sys->csize, 
		     sys->invg2c[0], TRUE, gridsize);

    nextc->g2p = nextc->q2g; 
    nextc->g2p = G2P(nextc->chgs, nextc->upnumeles[0], nextc->x, 
		     nextc->y, nextc->z, sys->colloc, sys->csize, 
		     sys->invg2c[0], TRUE, gridsize);
  }

  /* Precorrect the direct part. */ 
  for(nextc=sys->directlist; nextc != NULL; nextc = nextc->dnext) { 
    for(nummats=0, i= -1; i < nextc->numnbrs; nummats++, i++) { 
      if(i >= 0) nextnbr = nextc->nbrs[i]; 
      else nextnbr = nextc; 
      ASSERT(nextnbr->upnumvects > 0);
      j = nextc->j - nextnbr->j; 
      k = nextc->k - nextnbr->k; 
      l = nextc->l - nextnbr->l; 
 
      precorrect(nextc->directmats[nummats], nextc->g2p, nextnbr->q2g, 
		 nextc->upnumeles[0], nextnbr->upnumeles[0],
		 sys->g2g[nextnbr->l][cindex(j, k, l)], sys->gridsize,
		 sys->length);

 
    }
  }
}


