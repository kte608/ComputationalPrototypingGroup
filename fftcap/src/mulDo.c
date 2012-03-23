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

#if OPCNT == ON
static directops=0, upops=0, downops=0, evalops=0;
#endif

/* 
Compute the direct piece. 
*/
mulDirect(sys)
ssystem *sys;
{
int i, j, k, dsize;
double pc, *p, *q, *qn, *pn, **mat;
cube *nextc;
double *p1, *p2; 

/* Assumes the potential vector has been zero'd!!!! */
  for(nextc=sys->directlist; nextc != NULL; nextc = nextc->dnext) {
    dsize = nextc->directnumeles[0];  /* Equals number of charges. */

  /* Through all nearest nbrs. */
    for(i=nextc->directnumvects - 1; i >= 0; i--) {
      for(j = dsize - 1, p = &(nextc->eval[dsize-1]); j >= 0; j--, p--) {  
        for(k = nextc->directnumeles[i] - 1, p1 = nextc->directq[i], p2 = nextc->directmats[i][j]; 
            k >= 0; k--, p1++, p2++) {   
          (*p) += (*p2) * (*p1); 
	}
      }
    }
  }
}

  
mulGrid(sys)
ssystem *sys;
{
  int size1 = sys->sizej, size2 = sys->sizek, size3 = sys->sizel;
/*  int fac = size1 * size2 * size3;  */ 
  int i, j, k, ip, jp, kp;
  double *q, *p, **q2g, r;
  Real *sp1, *sp2, **g;
  cube *nextc;
  Real fac = 2.0 / (Real)(size1*size2*size3);  
  Real ***data = sys->data, ***kdata = sys->kdata;
  double delta;
  Real ***kern;
  FILE *fp, *fopen(); 
  double im; 
  int cj,ck,cl, tj, tk, tl; 
  double *p2, *p1, **pp, *px, *py; 
  int gridsize = sys->gridsize; 

  /* First zero out the data. */
  sp1 = &(sys->data[1][1][1]);
  for(j=1; j <= (size1 * size2 * size3/2); j++, sp1 += 1) {
    *sp1 = 0.0;
  }

  /* Extrapolate local charges onto the grid. */

  for(nextc=sys->directlist; nextc != NULL; nextc = nextc->dnext) {
    tj = gridsize * (nextc->j - sys->minj) +1 + gridsize; 
    tk = gridsize * (nextc->k - sys->mink) +1 + gridsize; 
    tl = gridsize * (nextc->l - sys->minl) +1; 
    pp = nextc->q2g; 
    
    for (cj = gridsize; cj >=0; cj--) { 
      for (ck = gridsize; ck >=0; ck--) { 
        p2 = &(sys->data[tj-cj][tk-ck][tl]);  
        for (cl = gridsize; cl >=0; cl--, p2++) { 
          p1 = *pp++; 
          p = nextc->directq[0];  
          for(i = nextc->directnumeles[0] - 1; i >= 0; i--, p++) {
            *(p2) += (*p1++) * (*p); 
          }
        }
      }
    }

  }

#if GNDPLN == ON 
  convolve(&(sys->data[1][1][1]), sys->kdata, sys->kspeq, -1, size1, size2, size3/2, sys->k2data, sys->k2speq); 
#else 
  convolve(&(sys->data[1][1][1]), sys->kdata, sys->kspeq, -1, size1, size2, size3/2, NULL, NULL); 
#endif 

  /* Interpolate potentials from the grid. */
  for(nextc=sys->directlist; nextc != NULL; nextc = nextc->dnext) {
    tj = gridsize * (nextc->j - sys->minj) +1 + gridsize; 
    tk = gridsize * (nextc->k - sys->mink) +1 + gridsize; 
    tl = gridsize * (nextc->l - sys->minl) +1; 
    pp = nextc->g2p; 
    
    for (cj = gridsize; cj >=0; cj--) { 
      for (ck = gridsize; ck >=0; ck--) { 
        p2 = &(sys->data[tj-cj][tk-ck][tl]);  
        for (cl = gridsize; cl >=0; cl--, p2++) { 
          p1 = *pp++; 
          p = nextc->eval;  
          for(i = nextc->directnumeles[0] - 1; i >= 0; i--, p++) {
            *p += (*p1++) * (*p2); 
          }
        }
      }
    }

  }
}

