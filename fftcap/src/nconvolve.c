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

struct cplx {
  double r, i;
};
typedef struct cplx cplx;

/*
* Bitrev takes a integer vector of length N and fills it so the i^th entry
* is i bit reversed.
*/
atic void bitrev(x, N)
int x[];
unsigned long N;
{
  int i, j, k, m, s, q;

  for(k = 0; k < N; k++) x[k] = k;
  for(k = 0; k < N; k++) {
    for(j = 0, m = k, q = 1; q < N; q *= 2) {
      s = m >> 1;
      j = 2 * j + (m - 2 * s);
      m = s;
    }
    if (j > k) {
      m = x[j];
      x[j] = x[k];
      x[k] = m;
    }
  }
}


static cplx *calcwf(N)
unsigned long N;
{
  int i, j, Ls, index;
  cplx *w;
  double angle, pi;

  pi = 2 * acos(0.0);

  w = (cplx *) calloc(N, sizeof(cplx));
  for(Ls = 1, index = 0; Ls < N; Ls <<= 1) {
    angle = pi / Ls;
    for(j=0; j < Ls; j++) {
      w[index].r = cos(j * angle);
      w[index].i = -sin(j * angle);
      w[index].i = sin(j * angle);
      index++;
    }
  }
  return(w);
}

/* wr(n) = exp(- 2 * pi * sqrt(-1) * n /(2 * N)); */

static cplx *calcwr(N)
unsigned long N;
{
  int i;
  cplx *w;
  double angle;

  angle = 2 * acos(0.0) / N; /* 2 * acos(0) = pi */ 

  w = (cplx *) calloc(N, sizeof(cplx));
  for(i=0; i < N; i++) {
    w[i].r = cos(i * angle);
    w[i].i = - sin(i * angle);
    w[i].i = sin(i * angle);
  }
  return(w);
}

/*
* xm[N1][N2][N3/2] is a complex array of data.  The convolution will
* be returned in this array.
*
* kernel[N1][N2][N3/2] is a complex array of fourier transformed Kernel data.
*
* datatype < 0 indicates that xm and kernel were packed real data.
* 
* N1, N2, and N3 are the sizes of the padded data.
*/
void convolve(xm, kernel, kspec, datatype, N1, N2, N3, kernel2, kspec2)
cplx *xm, ***kernel;
double **kspec; 
int datatype;
long N1, N2, N3;
cplx ***kernel2; 
double **kspec2; 
{
  register cplx *x, *xa, *xb, *xplane, *wp;
  register double wr, wi, tempr, tempi;
  double temppr, temppi;
  int j, k, Ls, LsN, L, kL, plane, row, jN;
  static int allocN1=0, allocN2=0, allocN3=0;
  static cplx *w1, *w2, *w3, *wreal, *workvector, *workvector2;
  register cplx *y; 
  int rowx, rowy; 
  static int *bitr1, *bitr2, *bitr3; 
  int wrapx, wrapy, ct, wrap; 
  cplx *z; 
  double ***dat, **speq; 
  cplx *zd; 
  int i; 
  double fac; 
  cplx *zs; 
  cplx * zp, *zs2; 
  cplx s0, s3; 
  static double *xr, *xi; 

#if GNDPLN == ON 
  static cplx *rrev, *wrev; 
  static double *W = NULL; 
  static int lastsize = 0; 
  { 
    int size3 = N3 * 2; 
    if (size3 > lastsize) {
      if (W) free(W); 
      W = (double *) calloc(2*size3+ 2, sizeof(double)); 
    }
    if (size3 != lastsize) {
      for (i=0; i< size3; i++) {
        k = size3 - i; 
        W[2*k-1] = cos(2.0*M_PI*(double)(i+1)/size3); 
        W[2*k] = sin(2.0*M_PI*(double)(i+1)/size3); 
      }
    }
    lastsize = size3;  
  }
#endif 

  if((N1 != allocN1) || (N2 != allocN2) || (N3 !=allocN3)) {
    /* Allocate a bit reversed vector for N3. */
    bitr1 = (int *) calloc(N1, sizeof(int));
    bitrev(bitr1,N1);

    bitr2 = (int *) calloc(N2, sizeof(int));
    bitrev(bitr2,N2);

    bitr3 = (int *) calloc(N3, sizeof(int));
    bitrev(bitr3,N3);
    /* Allocate the w vector for each dimension. */
    w1 = calcwf(N1);
    w2 = calcwf(N2);
    w3 = calcwf(N3);
    
    /* Allocate the wreal vector for unpacking the real data. */
    wreal = calcwr(N3);

    /* Allocate a temporary work vector of length max(N1,N2,N3). */
    j = N1; if(j < N2) j = N2; if(j < N3) j = N3;
    workvector = (cplx *) calloc(j+1, sizeof(cplx));
    workvector2 = (cplx *) calloc(j+1, sizeof(cplx));

    allocN1 = N1;
    allocN2 = N2;
    allocN3 = N3;
#if GNDPLN == ON 
    rrev = (cplx *) calloc(j+1, sizeof(cplx));
    wrev = (cplx *) calloc(j+1, sizeof(cplx));
#endif 

  }

  /* Note, data[planes][rows][cols], N1 planes, N2 rows, N3/2 cols.
   * The steps are organized so that the rows are done last, this is
   * done so that if the data was real, the conversion can be done easily.
   */ 

  /* In each of N1/2 planes, do all the N3/2 lines of length N2 in parallel. 
   * Effectively do all the columns. 
   * Maybe someday exploit the fact that half the entries are zero. */
  for(plane = 0; plane < N1/2; plane++) {  
    xplane = xm + plane * (N2 * N3 / 2);
    for(L = N2, Ls = N2 >> 1; Ls > 0; L = Ls, Ls >>= 1) { /* FFT in a Col. */
      LsN = Ls * N3/2;
      for (j=0; j < Ls; j++) {
	wr = w2[Ls-1 + j].r;
	wi = w2[Ls-1 + j].i;
	jN = j * N3/2;
	for (kL=0; kL < N2; kL += L) { /* Rest of FFT in a Col. */
	  xa = xplane + (kL * N3/2) + jN;
	  x = xa + N3/2;
	  for(xb = xa + LsN; xa < x; xa++, xb++) { /* Zip down N3/2 columns. */
	    tempr = xa->r - xb->r;
	    tempi = xa->i - xb->i;
	    xa->r += xb->r;
	    xa->i += xb->i;
	    xb->r = wr * tempr - wi * tempi;
	    xb->i = wr * tempi + wi * tempr;
	  }
	}
      }
    }
  }
  /* Do (N2 * N3/2) lines of length N1. Here, copy length N1/2 vector into
   * a vector of length N1 with zero padding.  Then transform the line.
   * Maybe someday exploit the fact that half the entries are zero. 
   * This vector copying seems to be faster than trying to do planes
   * simultaneously. I think the issue is that doing planes simultaneously
   * touches ALL the data at every one of log n steps.  It might
   * be faster to copy a N3/2 lines of length N1 into a plane, and then
   * follow the step 1 routine.
   */
  x = workvector;

  for(row = 0; row < (N2 * N3/2); row++) {
    /* Copy the vector and zero pad. */
    for(j = 0; j < N1/2; j++) x[j]  = xm[row + j * (N2 * N3/2)];
    for(; j < N1; j++) { x[j].r  = 0.0; x[j].i = 0.0; }
    /* Fourier transform the vector. */
    for(L = N1, Ls = N1 >> 1; Ls > 0; L = Ls, Ls >>= 1) {
      for (j=0; j < Ls; j++) {
	wr = w1[Ls-1 + j].r;
	wi = w1[Ls-1 + j].i;
	for (kL=0; kL < N1; kL += L) {
	  xa = x + kL + j;
	  xb = xa + Ls;
	  tempr = xa->r - xb->r;
	  tempi = xa->i - xb->i;
	  xa->r += xb->r;
	  xa->i += xb->i;
	  xb->r = wr * tempr - wi * tempi;
	  xb->i = wr * tempi + wi * tempr;
	}
      }
    }
    /* Put the vector back in the array. */
    for(j = 0; j < N1; j++) xm[row + j * (N2 * N3/2)] = x[j];
  }

  /* Here is the combined last step transform, multiply, 
   * first step inverse transform.  It is done this way to
   * save a factor of two in  memory.  Also, the rows are done
   * in this last step to make the conversion from packed real data easier.
   */


  for (rowx = 0; rowx < N1; rowx++) {

    wrapx = (rowx == 0 ?  0 : N1 - rowx); 

    for (rowy = 0; rowy < N2; rowy++) {
      
      wrapy = (rowy == 0 ?  0 : N2 - rowy); 
      
      /* the rows we are working on are 
         (rowx, rowy) and (wrapx, wrapy) 
         but they are stored by bit-reversed addresses, so decode first 
       */

      row = bitr2[rowy] + bitr1[rowx] * N2; 
      wrap = bitr2[wrapy] + bitr1[wrapx] * N2; 
      if (row > wrap) 
        continue;


      for (ct = 0, x=workvector; 
           ct<2; ct++, x=workvector2, row = wrap)  { 
        /* Copy the vector and zero pad. */
        for(j = 0; j < N3/2; j++) x[j]  = xm[j + row * N3/2];
        for(; j < N3; j++) { x[j].r  = 0.0; x[j].i = 0.0; }

        /* Fourier transform the vector. */
        for(L = N3, Ls = N3 >> 1; Ls > 0; L = Ls, Ls >>= 1) {
          for (j=0; j < Ls; j++) {
            wr = w3[Ls-1 + j].r;
            wi = w3[Ls-1 + j].i;
            for (kL=0; kL < N3; kL += L) {
              xa = x + kL + j;
              xb = xa + Ls;
              tempr = xa->r - xb->r;
              tempi = xa->i - xb->i;
              xa->r += xb->r;
              xa->i += xb->i;
              xb->r = wr * tempr - wi * tempi;
              xb->i = wr * tempi + wi * tempr;
            }
          }
        }
        if (row == wrap)
          break; 
      }

      x = workvector; 
      /* Correct for transform of real data packed into a complex vector. */
      if(datatype < 0) {
        row = bitr2[rowy] + bitr1[rowx] * N2; 
        wrap = bitr2[wrapy] + bitr1[wrapx] * N2; 
        x = workvector;
        if (row == wrap) 
          y = x; 
        else
          y = workvector2; 

        tempr = 0.5 * (x[0].r + y[0].r); 
        tempi = 0.5 * (x[0].i - y[0].i); 
        temppr = 0.5 * (x[0].r - y[0].r); 
        temppi = 0.5 * (x[0].i + y[0].i); 

        x[0].r = tempr + temppi; 
        x[0].i = tempi - temppr; 

        x[N3].r = tempr - temppi; 
        x[N3].i = tempi + temppr; 

        if (row != wrap) {
          tempr = tempr; 
          tempi = -tempi; 
          temppr =  -temppr; 
          temppi = temppi; 

          y[0].r = tempr + temppi; 
          y[0].i = tempi - temppr; 

          y[N3].r = tempr - temppi; 
          y[N3].i = tempi + temppr; 
        }

        for(L = 1; L < N3/2 ; L++) {
          k = bitr3[L];
          j = bitr3[N3 - L];
          tempr = 0.5 * (x[k].r + y[j].r);
          tempi = 0.5 * (x[k].i - y[j].i);
          temppr = 0.5 * (x[k].r - y[j].r); 
          temppi = 0.5 * (x[k].i + y[j].i);
          x[k].r = tempr + wreal[L].i * temppr + wreal[L].r * temppi;
          x[k].i = tempi + wreal[L].i * temppi - wreal[L].r * temppr;
          y[j].r = tempr - wreal[N3-L].i * temppr + wreal[N3-L].r * temppi;
          y[j].i = -tempi + wreal[N3-L].i * temppi + wreal[N3-L].r * temppr; 
          
          if (row != wrap) {
            tempr = 0.5 * (x[j].r + y[k].r);
            tempi = 0.5 * (x[j].i - y[k].i);
            temppr = 0.5 * (x[j].r - y[k].r); 
            temppi = 0.5 * (x[j].i + y[k].i);

            x[j].r = tempr + wreal[N3-L].i * temppr + wreal[N3-L].r * temppi;
            x[j].i = tempi + wreal[N3-L].i * temppi - wreal[N3-L].r * temppr;
            y[k].r = tempr - wreal[L].i * temppr + wreal[L].r * temppi;
            y[k].i = -tempi + wreal[L].i * temppi + wreal[L].r * temppr; 
          }
        }
      }
      x = workvector; y= workvector2; 

      fac = 1.0 / (double) (N1*N2*N3); 

#if GNDPLN == ON 
/* copy data, compute reversed transform in preparation for 
   Hankel-like part of kernel  */ 
      /* fc = 0 (DC) component unchanged */ 
      { 
        double s= -1.0; 
        int k; 
        double re, im, re2, im2; 
        int t; 

        if (row != wrap) y = workvector2; else y = workvector; 
        rrev[0] = x[0]; wrev[0] = y[0]; 

        for (k = 1; k<N3; k++) { 
          /* t = where is this component */ 
          t = bitr3[k]; 
          re = W[2*k-1+2] * x[t].r + W[2*k+2] * x[t].i; 
          im = W[2*k+2] * x[t].r - W[2*k-1+2] * x[t].i; 
          
          re2 = W[2*k-1+2] * y[t].r + W[2*k+2] * y[t].i; 
          im2 = W[2*k+2] * y[t].r - W[2*k-1+2] * y[t].i; 
          
          rrev[t].r = s * re2; rrev[t].i = s * im2; 
          wrev[t].r = s * re;  wrev[t].i = s * im; 
          s *= -1.0; 
        }
        /* N3 comp mult by -1 */ 
        rrev[N3].r = -1 * x[N3].r;  wrev[N3].r = -1 * y[N3].r; 
        rrev[N3].i = -1 * x[N3].i;  wrev[N3].i = -1 * y[N3].i; 
      }

#endif 

      /* Multiply by the kernel. */
      z = (cplx *) &(kernel[rowx+1][rowy+1][0].i); 
      for(j = 0; j < N3; j++)  { 
        tempr = x[j].r * z[bitr3[j]].r - x[j].i * z[bitr3[j]].i; 
        tempi = x[j].r * z[bitr3[j]].i + x[j].i * z[bitr3[j]].r; 
        x[j].r  = tempr * fac; 
        x[j].i  = tempi * fac; 
      }
        tempr = kspec[rowx+1][2*(rowy+1)-1]; tempi = kspec[rowx+1][2*(rowy+1)]; 
        temppr = x[N3].r * tempr - x[N3].i * tempi; 
        temppi = x[N3].r * tempi + x[N3].i * tempr; 
        x[N3].r = fac * temppr; x[N3].i = fac * temppi; 

      if (row != wrap) {
        z = (cplx *) &(kernel[wrapx+1][wrapy+1][0].i); 
        x = workvector2; 
        for(j = 0; j < N3; j++)  { 
          tempr = x[j].r * z[bitr3[j]].r - x[j].i * z[bitr3[j]].i; 
          tempi = x[j].r * z[bitr3[j]].i + x[j].i * z[bitr3[j]].r; 
          x[j].r  = tempr * fac; 
          x[j].i  = tempi * fac; 
        }
        tempr = kspec[wrapx+1][2*(wrapy+1)-1]; tempi = kspec[wrapx+1][2*(wrapy+1)]; 
        temppr = x[N3].r * tempr - x[N3].i * tempi; 
        temppi = x[N3].r * tempi + x[N3].i * tempr; 
        x[N3].r = fac * temppr; x[N3].i = fac * temppi; 
      }


#if GNDPLN == ON 
      /* Multiply by the Hankel-kernel  and add to Toeplitz part */
      z = (cplx *) &(kernel2[rowx+1][rowy+1][0].i); 
      x = rrev; 
      for(j = 0; j < N3; j++)  { 
        tempr = x[j].r * z[bitr3[j]].r - x[j].i * z[bitr3[j]].i; 
        tempi = x[j].r * z[bitr3[j]].i + x[j].i * z[bitr3[j]].r; 
        x[j].r  = tempr * fac; 
        x[j].i  = tempi * fac; 
      }
      tempr = kspec2[rowx+1][2*(rowy+1)-1]; tempi = kspec2[rowx+1][2*(rowy+1)]; 
      temppr = x[N3].r * tempr - x[N3].i * tempi; 
      temppi = x[N3].r * tempi + x[N3].i * tempr; 
      x[N3].r = fac * temppr; x[N3].i = fac * temppi; 
      for(j = 0; j <= N3; j++)  {   
        workvector[j].r += x[j].r; 
        workvector[j].i += x[j].i; 
      }

      if (row != wrap) {
        z = (cplx *) &(kernel2[wrapx+1][wrapy+1][0].i); 
        x = wrev; 
        for(j = 0; j < N3; j++)  { 
          tempr = x[j].r * z[bitr3[j]].r - x[j].i * z[bitr3[j]].i; 
          tempi = x[j].r * z[bitr3[j]].i + x[j].i * z[bitr3[j]].r; 
          x[j].r  = tempr * fac; 
          x[j].i  = tempi * fac; 
        }
        
        tempr = kspec2[wrapx+1][2*(wrapy+1)-1]; tempi = kspec2[wrapx+1][2*(wrapy+1)]; 
        temppr = x[N3].r * tempr - x[N3].i * tempi; 
        temppi = x[N3].r * tempi + x[N3].i * tempr; 
        x[N3].r = fac * temppr; x[N3].i = fac * temppi; 
        for(j = 0; j <= N3; j++)  {   
          workvector2[j].r += x[j].r; 
          workvector2[j].i += x[j].i; 
        }
      }
        
#endif /*GNDPLN */ 
      
      /* Unscramble the packed real data. */

      if(datatype < 0) {
          x = workvector; 
        if (row == wrap) 
          y = x; 
        else
          y = workvector2; 

        tempr = 0.5 * (x[0].r + y[N3].r); 
        tempi = 0.5 * (x[0].i - y[N3].i); 
        temppr = 0.5 * (x[0].r - y[N3].r); 
        temppi = 0.5 * (x[0].i + y[N3].i); 


        x[0].r = tempr - temppi; 
        x[0].i = tempi + temppr; 

        y[N3].r = tempr + temppi; 
        y[N3].i = -tempi + temppr; 

        if (row != wrap) {
        
          tempr = 0.5 * (x[N3].r + y[0].r); 
          tempi = 0.5 * (x[N3].i - y[0].i); 
          temppr = 0.5 * (x[N3].r - y[0].r); 
          temppi = 0.5 * (x[N3].i + y[0].i); 
            
          x[N3].r = tempr + temppi; 
          x[N3].i = tempi - temppr; 

          y[0].r = tempr - temppi; 
          y[0].i = -tempi - temppr; 
        }
        
        for(L = 1; L < N3/2 ; L++) {
          k = bitr3[L];
          j = bitr3[N3 - L];
          tempr = 0.5 * (x[k].r + y[j].r);
          tempi = 0.5 * (x[k].i - y[j].i);
          temppr = -0.5 * (x[k].r - y[j].r); 
          temppi = -0.5 * (x[k].i + y[j].i);
          x[k].r = tempr - wreal[L].i * temppr + wreal[L].r * temppi;
          x[k].i = tempi - wreal[L].i * temppi - wreal[L].r * temppr;
          y[j].r = tempr + wreal[N3-L].i * temppr + wreal[N3-L].r * temppi;
          y[j].i = -tempi - wreal[N3-L].i * temppi + wreal[N3-L].r * temppr; 

          if (row != wrap) {
            tempr = 0.5 * (x[j].r + y[k].r);
            tempi = 0.5 * (x[j].i - y[k].i);
            temppr = -0.5 * (x[j].r - y[k].r); 
            temppi = -0.5 * (x[j].i + y[k].i);
            
            x[j].r = tempr - wreal[N3-L].i * temppr + wreal[N3-L].r * temppi;
            x[j].i = tempi - wreal[N3-L].i * temppi - wreal[N3-L].r * temppr;
            y[k].r = tempr + wreal[L].i * temppr + wreal[L].r * temppi;
            y[k].i = -tempi - wreal[L].i * temppi + wreal[L].r * temppr; 
          }
        }
      }

    /* Inverse transform the product of x with the kernel. Note w*. */
          
      row = bitr2[rowy] + bitr1[rowx] * N2; 
      wrap = bitr2[wrapy] + bitr1[wrapx] * N2; 
      for (ct = 0, x=workvector; 
           ct<2; ct++, x=workvector2, row = wrap)  { 


        for(L = 2, Ls = 1; Ls < N3; Ls = L, L <<= 1) {
          wp = w3 + Ls - 1;
          for (kL=0; kL < N3; kL += L) {
            xa = x + kL;
            xb = xa + Ls;
            for (j=0; j < Ls; j++) {
              tempr = wp[j].r * xb[j].r + wp[j].i * xb[j].i;
              tempi = wp[j].r * xb[j].i - wp[j].i * xb[j].r; 
              xb[j].r = xa[j].r - tempr;
              xb[j].i = xa[j].i - tempi;
              xa[j].r += tempr;
              xa[j].i += tempi;
            }
          }
        }
        /* Copy back the result (only half the data is needed). */

        for(j = 0; j < N3/2; j++) xm[j + row * N3/2] = x[j];
        if (row == wrap) 
          break; 
      }
    }
  }
  /* Do (N2 * N3/2) lines of length N1. Here, copy length N1 length vector.
       * Maybe someday exploit the fact that half of the resulting N1 entries
       * are unneeded. Note w*.
       */
  x = workvector;

  for(row = 0; row < (N2 * N3/2); row++) {

    /* Copy the vector. */
    for(j = 0; j < N1; j++) x[j]  = xm[row + j * (N2 * N3/2)];
    /* Inverse transform */
    for(L = 2, Ls = 1; Ls < N1; Ls = L, L <<= 1) {
      wp = w1 + Ls - 1;
      for (kL=0; kL < N1; kL += L) {
        xa = x + kL;
        xb = xa + Ls;
        for (j=0; j < Ls; j++) { 
          tempr = wp[j].r * xb[j].r + wp[j].i * xb[j].i;
          tempi = wp[j].r * xb[j].i - wp[j].i * xb[j].r; 
          xb[j].r = xa[j].r - tempr;
          xb[j].i = xa[j].i - tempi;
          xa[j].r += tempr;
          xa[j].i += tempi;
        }
      }
    }
    /* Copy back the result (only half the data is really needed). */
    for(j = 0; j < N1; j++) xm[row + j * (N2 * N3 / 2)] = x[j];
  }

  /* In N1/2 planes, inverse trans N3/2 lines of length N2 in parallel. 
       * Maybe one day exploit that fact that half the N2 entries are unneeded.
       * Note w*.
       */
  for(plane = 0; plane < N1/2; plane++) {  
    xplane = xm + plane * (N2 * N3/2);
    for(L = 2, Ls = 1; Ls < N2; Ls = L, L <<= 1) {
      LsN = Ls * N3/2;
      for (j=0; j < Ls; j++) {
        wr = w2[Ls-1 + j].r;
        wi = w2[Ls-1 + j].i;
        jN = j * N3/2;
        for (kL=0; kL < N2; kL += L) {
          xa = xplane + (kL * N3/2) + jN;
          x = xa + N3/2;
          for(xb = xa + LsN; xa < x; xa++, xb++) { 
            tempr = wr * xb->r + wi * xb->i;
            tempi =  wr * xb->i - wi * xb->r;
            xb->r = xa->r - tempr;
            xb->i = xa->i - tempi;
            xa->r += tempr;
            xa->i += tempi;
          }
        }
      }
    }
  }
}
