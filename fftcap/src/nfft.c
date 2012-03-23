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

#include <stdio.h>
#include <math.h>

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
  if (N < 1) return 0; 
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
* date[N1][N2][N3] is a real array of data. 
*
* it is psychologically complex [N1][N2][N3/2] in much of this routine 
* data must be contiguous storage (!) 
*/
void rfft3(data, speq, N1, N2, N3) 
double ***data;
double **speq; 
int N1, N2, N3;
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
  int i; 
  double fac; 
  cplx ctmp; 
  int N32 = N3/2; 

  cplx *xm; cplx*xs; 
  xm = (cplx *) &(data[1][1][1]); 
  xs = (cplx *) &(speq[1][1]); 
  
  /* Allocate stuff if the sizes have changed. */

  if((N1 != allocN1) || (N2 != allocN2) || (N3 !=allocN3)) {
    /* Allocate a bit reversed vector for N3. */
    bitr1 = (int *) calloc(N1, sizeof(int));
    bitrev(bitr1,N1);

    bitr2 = (int *) calloc(N2, sizeof(int));
    bitrev(bitr2,N2);

    bitr3 = (int *) calloc(N3, sizeof(int));
    bitrev(bitr3,N32); 
    /* Allocate the w vector for each dimension. */
    w1 = calcwf(N1);
    w2 = calcwf(N2);
    w3 = calcwf(N3/2);
    
    /* Allocate the wreal vector for unpacking the real data. */
    wreal = calcwr(N32);

    /* Allocate a temporary work vector of length max(N1,N2,N3). */
    j = N1; if(j < N2) j = N2; if(j < N3) j = N3;
    workvector = (cplx *) calloc(j+1, sizeof(cplx));
    workvector2 = (cplx *) calloc(j+1, sizeof(cplx));

    
    allocN1 = N1;
    allocN2 = N2;
    allocN3 = N3;

  }

  /* Note, data[planes][rows][cols], N1 planes, N2 rows, N3/2 cols.
   * The steps are organized so that the rows are done last, this is
   * done so that if the data was real, the conversion can be done easily.
   */ 

  /* In each of N1 planes, do all the N3/2 lines of length N2 in parallel. 
   * Effectively do all the columns. 
   * Maybe someday exploit the fact that half the entries are zero. */

  for(plane = 0; plane < N1; plane++) {  
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
  

  /* Do (N2 * N3/2) lines of length N1. Here, copy length N1 vector into
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
    for(j = 0; j < N1; j++) x[j]  = xm[row + j * (N2 * N3/2)];

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

  /* figure out the bit-reverse in x-y */ 

  for (rowx = 0; rowx < N1; rowx++)  {
    for (rowy = 0; rowy < N2; rowy++) {
      row = bitr2[rowy] + bitr1[rowx] * N2; 
      wrap = rowy + rowx * N2; 
      if (row > wrap ) { 
        for(L = 0; L < N32 ; L++) {
          ctmp = xm[row * N32 + L]; xm[row* N32+L]=xm[wrap*N32+L]; xm[wrap*N32+L] = ctmp; 
        }   
        ctmp = xs[row]; xs[row] = xs[wrap]; xs[wrap] = ctmp; 
      }
    }
  }

  /* the z-direction transform */ 

  for (rowx = 0; rowx < N1; rowx++) {

    wrapx = (rowx == 0 ?  0 : N1 - rowx); 

    for (rowy = 0; rowy < N2; rowy++) {
      
      wrapy = (rowy == 0 ?  0 : N2 - rowy); 
      
      /* the rows we are working on are 
         (rowx, rowy) and (wrapx, wrapy) 
       */


      row = rowy + rowx * N2; 
      wrap = wrapy + wrapx * N2; 
      if (row > wrap) continue; 

      for (ct = 0, x=workvector; 
           ct<2; ct++, x=workvector2, row = wrap)  { 

        /* Copy the vector and zero pad. */
        for(j = 0; j < N3/2; j++) x[j]  = xm[j + row * N3/2];

        /* Fourier transform the vector. */
        for(L = N32, Ls = N32 >> 1; Ls > 0; L = Ls, Ls >>= 1) {
          for (j=0; j < Ls; j++) {
            wr = w3[Ls-1 + j].r;
            wi = w3[Ls-1 + j].i;
            for (kL=0; kL < N3/2; kL += L) {
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

      /* figure out z-bit-reverse */ 

      /* bit-reverse in z */ 
      x = workvector; y  = workvector2; 

      row = rowy + rowx * N2; 
      wrap = wrapy + wrapx * N2; 

      for (i=0; i<N32; i++) { 
        if (bitr3[i] > i) { 
          ctmp = x[i]; x[i] = x[bitr3[i]]; x[bitr3[i]] = ctmp; 
          if (row != wrap) { 
            ctmp = y[i]; y[i] = y[bitr3[i]]; y[bitr3[i]] = ctmp; 
          }
        }
      }

      /* Correct for transform of real data packed into a complex vector. */

      row = rowy + rowx * N2; 
      wrap = wrapy + wrapx * N2; 

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

      x[N32].r = tempr - temppi; 
      x[N32].i = tempi + temppr; 

      if (row != wrap) {
        tempr = tempr; 
        tempi = -tempi; 
        temppr =  -temppr; 
        temppi = temppi; 

        y[0].r = tempr + temppi; 
        y[0].i = tempi - temppr; 

        y[N32].r = tempr - temppi; 
        y[N32].i = tempi + temppr; 
      }

      for(L = 1; L < N32/2; L++) {
        k = L;
        j = N32 - L;
        tempr = 0.5 * (x[k].r + y[j].r);
        tempi = 0.5 * (x[k].i - y[j].i);
        temppr = 0.5 * (x[k].r - y[j].r); 
        temppi = 0.5 * (x[k].i + y[j].i);
        x[k].r = tempr + wreal[L].i * temppr + wreal[L].r * temppi;
        x[k].i = tempi + wreal[L].i * temppi - wreal[L].r * temppr;
        y[j].r = tempr - wreal[N32-L].i * temppr + wreal[N32-L].r * temppi;
        y[j].i = -tempi + wreal[N32-L].i * temppi + wreal[N32-L].r * temppr; 
          
        if (row != wrap) {
          tempr = 0.5 * (x[j].r + y[k].r);
          tempi = 0.5 * (x[j].i - y[k].i);
          temppr = 0.5 * (x[j].r - y[k].r); 
          temppi = 0.5 * (x[j].i + y[k].i);

          x[j].r = tempr + wreal[N32-L].i * temppr + wreal[N32-L].r * temppi;
          x[j].i = tempi + wreal[N32-L].i * temppi - wreal[N32-L].r * temppr;
          y[k].r = tempr - wreal[L].i * temppr + wreal[L].r * temppi;
          y[k].i = -tempi + wreal[L].i * temppi + wreal[L].r * temppr; 
        }
      }

      /* copy into place */ 

      row = rowy + rowx * N2; 
      wrap = wrapy + wrapx * N2; 

      xs[row] = x[N32]; 
      if (row != wrap) { 
        xs[wrap] = y[N32]; 
      }
      for (ct = 0, x=workvector; 
           ct<2; ct++, x=workvector2, row = wrap)  { 
        for (j=0; j<N32; j++)  
          xm[j+row*N32] = x[j]; 
        if (row == wrap) 
          break; 
      }

    }
  }

}
