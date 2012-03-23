/*
Copyright (c) 2000 Massachusetts Institute of Technology, Cambridge, MA.
All rights reserved.

This Agreement gives you, the LICENSEE, certain rights and obligations.
By using the software, you indicate that you have read, understood, and
will comply with the terms.

Permission to use, copy and modify for internal, noncommercial purposes
is hereby granted.  Any distribution of this program or any part thereof
is strictly prohibited without prior written consent of M.I.T.

Title to copyright to this software and to any associated documentation
shall at all times remain with M.I.T. and LICENSEE agrees to preserve
same.  LICENSEE agrees not to make any copies except for LICENSEE'S
internal noncommercial use, or to use separately any portion of this
software without prior written consent of M.I.T.  LICENSEE agrees to
place the appropriate copyright notice on any such copies.

Nothing in this Agreement shall be construed as conferring rights to use
in advertising, publicity or otherwise any trademark or the name of
"Massachusetts Institute of Technology" or "M.I.T."

M.I.T. MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.  By
way of example, but not limitation, M.I.T. MAKES NO REPRESENTATIONS OR
WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR
THAT THE USE OF THE LICENSED SOFTWARE COMPONENTS OR DOCUMENTATION WILL
NOT INFRINGE ANY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS.
M.I.T. shall not be held liable for any liability nor for any direct,
indirect or consequential damages with respect to any claim by LICENSEE
or any third party on account of or arising from this Agreement or use
of this software.
*/
/*
* FILE: nfft.c
*
*  This file contains routines that operate on the data used in 
*  the FFT part of the code. 
*
* JKW: Okay, so the code in this file indexes arrays starting 
* from one, and everywhere else the indices start from zero.  
* So sue me.
*
* BB: OK. It's your hash ;) 
* Seriously: DO TAKE NOTE that the arrays in this part of the code 
* are made from 1 to N (i.e. allocating N+1 entries, but not using 
* the first one).
* Historically, the choice was made to reuse the FFT part of the 
* old code (FFTCap). The extra amount of memory needed is probably 
* negligible - even though it is not a constant (probably the 
* largest part of extra memory is O((s1+1)*(s2+1)*(s3+1)-s1*s2*s3)
* = O(s1*s2+s1*s3+s2*s3) where "si" is (twice) the number of grid 
* points in the i'th direction. Anyway, the primary down side is 
* the "inconsistency" when compared with the rest of the code.
* An exception to this rule is, of course(?), the grid point 
* indices, which are numbered from zero to N-1 as in the rest 
* of the code. (There may be more cases like this one.)
* Also, rfft3 seems to use the convention 0 to N-1, at least in 
* some cases.
*/

# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <assert.h>
# include <float.h>

# include "Global.h"
# include "Calcp.h"
# include "Nfft.h"

#define CIRC(xi, size) ( ( (xi) >= 0) ? (1+(xi)) : (1+(size)+(xi)) ) 

/*#define INDEX(i,j,k,s1,s2,s3) ((i-1)*s2*s3 + (j-1)*s3 + k)*/
/* Given (i,j,k) and maximum size (s1,s2,s3) calculate index n into 
 * the grid-data array                                                 */
#define INDEX(i,j,k,s1,s2,s3) ((i-1)*s2*s3 + (j-1)*s3 + k-1)
/* Given the index n into the grid-data array and the maximum 
 * size (s1,s2,s3), calculate (i,j,k). This should be the inverse 
 * mapping of INDEX defined above. Note, that it is exploited that 
 * integer division (by the ANSI standard) truncates any fractional 
 * part.                                                               */
#define INDEX_IJK_NEW(n, i,j,k, s1,s2,s3){              \
	  (i)      = 1 + (n)/((s2)*(s3));              \
	  (j)      = 1 + ((n)%((s2)*(s3))) /(s3);      \
	  (k)      = 1 + ((n)%(s3));                   \
}
/* This "old" definition may be slighly faster, but it is also harder 
 * to see what is going on. It should be absolutely equivalent to 
 * INDEX_IJK_NEW above                                                 */
#define INDEX_IJK_OLD(n, i,j,k, s1,s2,s3){int remain;\
	  i      = 1 + n/(s2*s3);                    \
	  remain = n - (i-1)*s2*s3;                  \
	  j      = 1+  remain/s3;                    \
	  k      = remain - (j-1)*s3 +1;             \
}
/* Feel free to swap between "old" and "new" versions                  */
#define INDEX_IJK(n, i,j,k, s1,s2,s3)             \
        INDEX_IJK_NEW(n, i,j,k, s1,s2,s3)



/* Size settings (increase) for the data/kernel FFT grid (s1,s2,s3) 
 * compared to the "physical" model problem grid (nx,ny,nz).
 * This was moved to a macro since it is used several places in 
 * the code (both when estimating memory size and when actually 
 * allocating). */
#define S1S2S3INCREASE(s1,s2,s3,Type,fctName) {  s1 *= 2;s2 *= 2;\
   if(strcmp(Type, "kernel") == 0) {s3 *= 2;}\
   else if(strcmp(Type, "data") != 0) {\
    fprintf(stderr,"%s: ERROR. Unknown type %s\n"\
	    "    Please use \"kernel\" or \"data\"\n"\
	    "    Exiting\n",\
	    fctName,Type);\
    exit(1);\
  }\
}

struct cplx {
  double r, i;
};
typedef struct cplx cplx;

struct gridloc {
  int i, j, k;
};
typedef struct gridloc gridloc;

struct fftGridData {
  double ***data;     /* The grid data. 3D array [i][j][k]             */
  double **speq;      /* Extra work vector used by the fft.            */
  int s1, s2, s3;     /* Number of grid points in each direction.      */
  int c1, c2, c3;     /* Center point in each direction. 
		       * At the moment these are s1/2, s2/2, s3/2, but 
		       * that could change, so they will be stored.     */
};
typedef struct fftGridData fftGridData;

/* STATIC routines for this file: */
static void bitrev(int *x, unsigned long N);
/*static void bitrev(int *x, int N);*/

/* 
* ==================================================================== 
* Creates the FFT data areas for the real data / real kernel FFT's 
* and  convolutions.
* ==================================================================== */
void *createFFT(int s1, int s2, int s3, char *Type)
{
  char fctName[] = "createFFT";
  fftGridData *griddata;
  double *datas, ***data, **speq;
  int i, j;

  griddata = (fftGridData *) calloc(1, sizeof(fftGridData));

  /* Make grid 2x larger in every dimension if kernel. 
   * If data, then make larger only in two of the three dimensions.    
  s1 *= 2;
  s2 *= 2;
  if(strcmp(Type, "kernel") == 0) {
    s3 *= 2;
  }
  else if(strcmp(Type, "data") != 0) {
    fprintf(stderr,"%s: ERROR. Unknown type %s\n"
	    "    Please use \"kernel\" or \"data\"\n"
	    "    Exiting\n",
	    fctName,Type);
    exit(1);
  }*/
  S1S2S3INCREASE(s1,s2,s3,Type,fctName);

  /* Store sizes (maximum entry value) of grid data in each dimension  */
  griddata->s1 = s1;
  griddata->s2 = s2;
  griddata->s3 = s3;

  /* Allocate the space as one long vector, then set up pointers.      */
  data = (double ***) calloc(s1+1, sizeof(double **));
  datas = (double *) calloc(s1*s2*s3+1, sizeof(double));
  /* Keep the starting point of "datas" around, so that it can be freed 
   * if necessary. Store it in data[0], which is not used anyway, since 
   * JKW has chosen to use a 1 to N convention, rather than 0 to N-1.  */
  data[0] = (double **) datas;
  /* Add pointers into the long vector of data, such that data[i][j][k] 
   * will be in the array allocated as "datas" above. Make sure that 
   * data[i][j][s3] is followed immediately by data[i][j+1][1] (j<s2).
   * This implies that data[i][j+1] = data[i][j+1][0] = data[i][j][s3] */
  for(i=1; i <= s1; i++) {
    data[i] = (double **) calloc(s2+1, sizeof(double *));
    for(j=1; j <= s2; j++) 
      data[i][j] = &(datas[INDEX(i,j,0,s1,s2,s3)+1]);
  }

  /* Extra work vector used during the fft. */
  speq = (double **) calloc(s1+1, sizeof(double *));
  datas = (double *) calloc((s1+1)*(2*s2+2), sizeof(double));
  /* Store beginning of "datas" to enable freeing the memory later     */
  speq[0] = datas;
  for(i=1; i <= s1; i++) 
    speq[i] = &(datas[1+(i-1) * (2 * s2)]);

  griddata->speq = speq;
  griddata->data = data;

  return (void *) griddata;
} /* End of routine createFFT */

/* 
* ==================================================================== 
* Free the FFT data areas for the real data / real kernel FFT's 
* and  convolutions which is set up by createFFT.
* Note:
*   This routine has been tested in an infinite loop with createFFT. 
*   After 10^6 successive allocations and freeing, no increase in the 
*   memory was observed (1kb resolution). (Allocating but not freeing an 
*   equal sized problem used >100Mb in less than 4000 allocations.) 
* ==================================================================== */
void freeFFT(void *dataIn)
{
  /*char fctName[] = "freeFFT";*/
  fftGridData *griddata;
  int i;

  /* Cast data to correct type: */
  griddata = (fftGridData*) dataIn;

  /* Free main data areas: */
  FREE(griddata->data[0]);
  FREE(griddata->speq[0]);

  /* Free data and speq pointers */
  for(i=1; i <= griddata->s1; i++) {
    FREE(griddata->data[i]);
  }
  FREE(griddata->speq);
  FREE(griddata->data);

  /* Free the griddata structure itself */
  FREE(griddata);

  return;
} /* End of routine freeFFT */

/* 
* ==================================================================== 
* Converts external grid indices [i][j][k] to vector offsets. 
* Note that griddata->s1, griddata->s2, griddata->s3 are the 
* NUMBER OF grid points in each dimension. 
* Thus, e.g. 0 < i < griddata->s1 + 1. 
* ==================================================================== */
int dataIndex(void *gridDataIn, int i, int j, int k)
{
  fftGridData *griddata;
  /* Cast input data to correct format                                 */
  griddata = (fftGridData *) gridDataIn;

  /* Bump i, j and k by one, as in here grids go from 1->n, not 0->n-1.*/
  i += 1;
  j += 1;
  k += 1;

  /* Make sure that we are inside the grid structure, i.e. that an
   * existing point is queried for (the alternative is a point outside 
   * the grid, which is generally a bad idea.                          */
  assert(i>0);
  assert(j>0);
  assert(k>0);
  assert(i <= griddata->s1);
  assert(j <= griddata->s2);
  assert(k <= griddata->s3);
  
  return(INDEX(i,j,k,griddata->s1, griddata->s2, griddata->s3));
}

/* 
* ==================================================================== 
* fillKernel fills in the kernel data with the kernel values. 
* Must be careful about the definitions of everything so that 
* the fft comes out right.  Not yet checked. 
* ==================================================================== */
void fillKernel(void *griddataIn, double dx, double dy, double dz)
{
  fftGridData *g;
  int centerj, centerk, centerl, j, k, l;
  double val;

  g = (fftGridData *) griddataIn;

  /* Centre of the Greens function; r=(0,0,0) */
  g->c1 = g->s1/2;
  g->c2 = g->s2/2;
  g->c3 = g->s3/2;

  /* Local shorthands */
  centerj = g->c1;
  centerk = g->c2;
  centerl = g->c3;
  /* BB: Why are the values needed for j=centerj, k=centerk and l=centerl? */
  for(j = 1-centerj; j <= centerj; j++) {
    for(k = 1-centerk; k <= centerk; k++) {
      for(l = 1-centerl; l <= centerl; l++) {
	val = kernel(j * dx, k * dy, l * dz);
	g->data[CIRC(j,g->s1)][CIRC(k,g->s2)][CIRC(l,g->s3)] = val;
      }
    }
  }
 return;
} /* End of routine fillKernel */

/* 
* ==================================================================== 
* Calculate correction factor for a single entry in the direct matrix.
* 
* NOTE: This routine is called once per entry in the direct matrix.
* For relatively high accuracy (maybe also for medium accuracy) 
* this is a lot of times. An efficient implementation of this routine 
* seems to be very important. At the time of writing, the present 
* routine in some cases uses more than 60% of the total CPU time 
* (not including the solution process).
*
* This routine is now obsolete and kept only to check timings against 
* this "old" method (preconditioning element-by-element). Also, this 
* routine may come in handy in later changes in implementation.
* ==================================================================== */
double precorrectVal(void *kernelIn, void *griddataIn, 
	  double *projectVals, int *projectGpoints, int projectPoints,
	  double *interpVals, int *interpGpoints, int interpPoints)
{
  /*char fctName[]="precorrectVal";*/  /* Name of this routine             */
  /* Input parameters: 
   *
   * kernelIn=g :      Kernel (Greens function) data on grid.
   * projectVals :     Values of the projection matrix. One value for 
   *                   each grid point.
   * projectGpoints:   Grid point indices. These are row numbers of a 
   *                   particular column of the projection matrix.
   * projectPoints:    Number of grid points in above.
   * interpVals :      Values of the interpolation matrix. One value for 
   *                   each grid point.
   * interpGpoints:    Grid point indices. These are row numbers of a 
   *                   particular column of the interpolation matrix.
   * interpPoints:     Number of grid points in above.
   */
  fftGridData *kernel, *griddata;
  double ***kdata, coeff, ival, g2g;
  int i, j;
#define TMPSIZE_precorrectVal  1000
  static int 
    iP[TMPSIZE_precorrectVal],
    jP[TMPSIZE_precorrectVal],
    kP[TMPSIZE_precorrectVal];
  /*register*/ 
  int 
    iInterp, jInterp, kInterp, 
    kernelS1,kernelS2,kernelS3, 
    dataS1,  dataS2,  dataS3;

  assert(projectPoints<=TMPSIZE_precorrectVal);

  /* Cast input data kernel to correct type */
  kernel   = (fftGridData *) kernelIn;
  griddata = (fftGridData *) griddataIn;

  /* Local shorthands: */
  kdata = kernel->data; /* kdata[i][j][k] contain the kernel. Note that at 
		   * this point the kernel still hasn't been convolved */

  kernelS1 = kernel->s1;
  kernelS2 = kernel->s2;
  kernelS3 = kernel->s3;
  dataS1   = griddata->s1;
  dataS2   = griddata->s2;
  dataS3   = griddata->s3;

  /* Setup of index_ijk for inner loop to avoid doing the job 
   * multiple times                                            */
  for(j = 0; j < projectPoints; j++) {
    INDEX_IJK(projectGpoints[j], iP[j],jP[j],kP[j],
	      dataS1,dataS2,dataS3);
  }
  /* Now do multiplication, I^t * G2G * P.  */
  for(coeff = 0.0, i = 0; i < interpPoints; i++) {
    INDEX_IJK(interpGpoints[i], iInterp,jInterp, kInterp,
	      dataS1,dataS2,dataS3);
    ival = interpVals[i];
    for(j = 0; j < projectPoints; j++) {
      /* This is known to work: */
      g2g = kdata[CIRC((iInterp-iP[j]),kernelS1)]
	         [CIRC((jInterp-jP[j]),kernelS2)]
	         [CIRC((kInterp-kP[j]),kernelS3)];
      /* Interpolation * kernel * projection */
      coeff += ival * g2g * projectVals[j];
    }
  }
  return coeff;
} /* End of routine precorrectVal */

/* 
* ==================================================================== 
* Given lists of projection and interpolation points 
* return kernel values for all combinations of point pairs.
* This assumes that the kernel has not been convolved.
* ==================================================================== */
void supplyKernelDataValues(void *kernelIn, void *griddataIn,
			    int nProjectPoints, int *projectPoints,
			    int nInterpPoints, int *interpPoints,
			    double **kernelValuesArray){
  /*char fctName[]="supplyKernelDataValues";*/  /* Name of this routine    */
  /* Input parameters: 
   *
   * kernelIn=g :      Kernel (Greens function) data on grid.
   * projectVals :     Values of the projection matrix. One value for 
   *                   each grid point.
   * projectGpoints:   Grid point indices. These are row numbers of a 
   *                   particular column of the projection matrix.
   * projectPoints:    Number of grid points in above.
   * interpVals :      Values of the interpolation matrix. One value for 
   *                   each grid point.
   * interpGpoints:    Grid point indices. These are row numbers of a 
   *                   particular column of the interpolation matrix.
   * interpPoints:     Number of grid points in above.
   */
  fftGridData *kernel, *griddata;
  int i,j;
  double ***kdata;
#define TMPSIZE_supplyKernelDataValues  2000
  static int 
    iI[TMPSIZE_supplyKernelDataValues],
    jI[TMPSIZE_supplyKernelDataValues],
    kI[TMPSIZE_supplyKernelDataValues];
  /*The following could be declared "register variables" */ 
  int 
    iProject, jProject, kProject,
    kernelS1,kernelS2,kernelS3, 
    dataS1,  dataS2,  dataS3;
  /* Cast input data kernel to correct type */
  kernel   = (fftGridData *) kernelIn;
  griddata = (fftGridData *) griddataIn;

  assert(nInterpPoints<=TMPSIZE_supplyKernelDataValues);
 
  /* Local shorthand: */
  kdata = kernel->data; /* kdata[i][j][k] contain the kernel. 
			 * Note that at this point the kernel still 
			 * hasn't been convolved (at least it should 
			 * not be convolved yet).                      */

  /* Set up shorthands */
  kernelS1 = kernel->s1;
  kernelS2 = kernel->s2;
  kernelS3 = kernel->s3;
  dataS1   = griddata->s1;
  dataS2   = griddata->s2;
  dataS3   = griddata->s3;

  /* Setup of index_ijk for inner loop to avoid doing the job 
   * multiple times                                            */
  for(j = 0; j < nInterpPoints; j++) {
    INDEX_IJK(interpPoints[j], iI[j],jI[j],kI[j],
		dataS1,dataS2,dataS3);  
  }

  /* Now fill in the grid-to-grid interactions:                         */
  for(i = 0; i < nProjectPoints; i++) {
    INDEX_IJK(projectPoints[i], iProject,jProject,kProject,
	      dataS1,dataS2,dataS3);
    for(j = 0; j < nInterpPoints; j++) {
      /*INDEX_IJK(interpPoints[j], iInterp,jInterp, kInterp,
	dataS1,dataS2,dataS3);*/
      /* This is "g2g":                                                */
      kernelValuesArray[i][j] = 
	kdata[CIRC((iI[j]-iProject),kernelS1)]
	     [CIRC((jI[j]-jProject),kernelS2)]
	     [CIRC((kI[j]-kProject),kernelS3)];
    }
  }
  return;
} /* End of routine supplyKernelDataValues */

/*
* Bitrev takes a integer vector of length N and fills it so the i^th entry
* is i bit reversed.
*/
static void bitrev(x, N)
int x[];
unsigned long N;
{
  int j, k, m, s, q;

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


/* BB: N must be 2^n - 1 where i is an integer (N=2^n will work also) */
static cplx *calcwf(N)
unsigned long N;
{
  int j, Ls, index;
  cplx *w;
  double angle, pi;

  pi = 2 * acos(0.0);

  w = (cplx *) calloc(N, sizeof(cplx));
  /* Ls = 1,2,4,8,... "Ls <<= 1;" is equivalent to "Ls *= 2;". */
  for(Ls = 1, index = 0; Ls < N; Ls <<= 1) {
    angle = pi / Ls;
    /* First pass fills cos(0*pi/1) and sin(0*pi/1) into w[0].
     * Second pass fills cos(0*pi/2) and sin(0*pi/1) into w[1] and
     *                   cos(1*pi/2) and sin(1*pi/1) into w[2].
     * And so forth.
     * If N=2^i (i being a natural number) then the last entry 
     * w[N-1] will not be given a new value.
     */
    for(j=0; j < Ls; j++) {
      w[index].r = cos(j * angle);
      /*w[index].i = -sin(j * angle);*/
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
    /*w[i].i = - sin(i * angle);*/
    w[i].i = sin(i * angle);
  }
  return(w);
}

/*
* FFT transform kernel data. 
* This is used for calculating the discrete convolution product of 
* the kernel data and the "source" data using the FFT
* (See e.g. Oppenheim and Schafer, 1989, 
*  "Discrete-Time Signal Processing", p. 548):
*
* data[N1][N2][N3] is a real array of data. 
*
* it is psychologically complex [N1][N2][N3/2] in much of this routine 
* data must be contiguous storage (!) 
*/
void rfft3(void *fftdata)
{
  fftGridData *g;
  double ***data, **speq;
  int N1, N2, N3, N32;
  register cplx *x, *xa, *xb, *xplane;
  register double wr, wi, tempr, tempi;
  double temppr, temppi;
  int j, k, Ls, LsN, L, kL, plane, row, jN;
  static int allocN1=0, allocN2=0, allocN3=0;
  static cplx *w1=NULL, *w2=NULL, *w3=NULL, *wreal=NULL, 
    *workvector=NULL, *workvector2=NULL;
  register cplx *y; 
  int rowx, rowy; 
  static int *bitr1=NULL, *bitr2=NULL, *bitr3=NULL; 
  int wrapx, wrapy, ct, wrap; 
  int i; 
  cplx ctmp; 
  cplx *xm, *xs; 

  /* Cast input to correct type */
  g = (fftGridData *) fftdata;

  /* Local shorthands */
  N1 = g->s1;
  N2 = g->s2;
  N3 = g->s3;
  N32 = N3/2; 

  data = g->data;
  speq = g->speq;

  xm = (cplx *) &(data[1][1][1]); 
  xs = (cplx *) &(speq[1][1]); 
  
  /* Allocate stuff if the sizes have changed. */

  if((N1 != allocN1) || (N2 != allocN2) || (N3 !=allocN3)) {
    /* Allocate a bit reversed vector for N3. */
    if (bitr1 != NULL) FREE(bitr1);
    if (bitr2 != NULL) FREE(bitr2);
    if (bitr3 != NULL) FREE(bitr3);
    bitr1 = (int *) calloc(N1, sizeof(int));
    bitrev(bitr1,N1);
    bitr2 = (int *) calloc(N2, sizeof(int));
    bitrev(bitr2,N2);
    bitr3 = (int *) calloc(N3, sizeof(int));
    bitrev(bitr3,N32); 
    /* Allocate the w vector for each dimension. */
    if (w1 != NULL) FREE(w1);
    if (w2 != NULL) FREE(w2);
    if (w3 != NULL) FREE(w3);
    w1 = calcwf(N1);
    w2 = calcwf(N2);
    w3 = calcwf(N3/2);
    
    /* Allocate the wreal vector for unpacking the real data. */
    if (wreal != NULL) FREE(wreal);
    wreal = calcwr(N32);

    /* Allocate a temporary work vector of length max(N1,N2,N3). */
    j = N1; if(j < N2) j = N2; if(j < N3) j = N3;
    if (workvector  != NULL) FREE(workvector );
    if (workvector2 != NULL) FREE(workvector2);
    workvector  = (cplx *) calloc(j+1, sizeof(cplx));
    workvector2 = (cplx *) calloc(j+1, sizeof(cplx));


    allocN1 = N1;
    allocN2 = N2;
    allocN3 = N3;

  } /* The work vectors now have the right size */

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

} /* End of routine rfft3 */

/*
* xm[N1][N2][N3/2] is a complex array of data.  The convolution will
* be returned in this array.
*
* kernel[N1][N2][N3/2] is a complex array of fourier transformed 
* Kernel data.
*
* datatype < 0 indicates that xm and kernel were packed real data.
* 
* N1, N2, and N3 are the sizes of the padded data.
*/
void convolve(void *xIn, void *kernelIn, int datatype)
{
  fftGridData *xdata, *kerneldata;
  register cplx *x, *xa, *xb, *xplane, *wp;
  register double wr, wi, tempr, tempi;
  double temppr, temppi;
  int j, k, Ls, LsN, L, kL, plane, row, jN;
  static int allocN1=0, allocN2=0, allocN3=0;
  static cplx *w1=NULL, *w2=NULL, *w3=NULL, *wreal=NULL, 
    *workvector=NULL, *workvector2=NULL;
  register cplx *y; 
  int rowx, rowy; 
  static int *bitr1=NULL, *bitr2=NULL, *bitr3=NULL; 
  int wrapx, wrapy, ct, wrap; 
  cplx *z; 
  double fac; 
  cplx *xm, ***kernelVals;
  double **kspec;
  int N1, N2, N3;
  
  xdata = (fftGridData *) xIn;
  kerneldata = (fftGridData *) kernelIn;
  xm = (cplx *) (&(xdata->data[1][1][1]));
  kernelVals = (cplx ***) kerneldata->data;
  kspec = kerneldata->speq;
  N1 = xdata->s1;
  N2 = xdata->s2;
  N3 = xdata->s3;

  if((N1 != allocN1) || (N2 != allocN2) || (N3 !=allocN3)) {
    /* Allocate a bit reversed vector for N3. */
    if (bitr1 != NULL) FREE(bitr1);
    if (bitr2 != NULL) FREE(bitr2);
    if (bitr3 != NULL) FREE(bitr3);
    bitr1 = (int *) calloc(N1, sizeof(int));
    bitrev(bitr1,N1);
    bitr2 = (int *) calloc(N2, sizeof(int));
    bitrev(bitr2,N2);
    bitr3 = (int *) calloc(N3, sizeof(int));
    bitrev(bitr3,N3);
    /* Allocate the w vector for each dimension. */
    if (w1 != NULL) FREE(w1);
    if (w2 != NULL) FREE(w2);
    if (w3 != NULL) FREE(w3);
    w1 = calcwf(N1);
    w2 = calcwf(N2);
    w3 = calcwf(N3);
    
    /* Allocate the wreal vector for unpacking the real data. */
    if (wreal != NULL) FREE(wreal);
    wreal = calcwr(N3);

    /* Allocate a temporary work vector of length max(N1,N2,N3). */
    j = N1; if(j < N2) j = N2; if(j < N3) j = N3;
    if (workvector  != NULL) FREE(workvector );
    if (workvector2 != NULL) FREE(workvector2);
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

      /* Multiply by the kernel. */
      z = (cplx *) &(kernelVals[rowx+1][rowy+1][0].i); 
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
        z = (cplx *) &(kernelVals[wrapx+1][wrapy+1][0].i); 
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
} /* End of routine convolve */


/*
* ==================================================================== 
* Project a vector (of source-element values) onto the grid. 
* This routine extracts the necessary data from the griddata structure 
* and passes them on (to an sparse.c routine) for further processing, 
* since the grid doesn't know what a sparse matrix is.
* ==================================================================== */ 
void projectOntoGridData(void *gridDataIn,double *src,void *projectMat)
{
  fftGridData *gridData;
  double *dest;
  int matNum=0, i,imax;

  /* Cast input data to correct type */
  gridData = (fftGridData *) gridDataIn;

  /* Extract information on the starting point of the array of 
   * grid data. The first element in dest should be 
   * dest[0]=data[1][1][1]                                            */
  dest = &(gridData->data[1][1][1]);

  /* Make sure that the values of the grid data are zero before making
   * the sparse matrix-vector product. Recall that the matrix-vector 
   * multiply routine will ADD the product to the values on the 
   * destination vector. 
   * Exploit that the gridData are stored contiguous in memory. Thus 
   * a SINGLE pass through dest[i] can be used, rather than a tripple 
   * loop through gridData->data[i][j][k]:                            */
  imax = gridData->s1*gridData->s2*gridData->s3;
  for (i=0;i<imax;i++){
    dest[i] = 0.0;
  }

  /* Make a sparse matrix-vector multiply (projctMat times src)
   * Store result in dest, i.e. in gridData->data                     */
  spMatMultiply(projectMat,src,dest,matNum);
  /* The data should now be projected onto the grid                   */
  return;
} /* End of routine projectOntoGridData */

/*
* ==================================================================== 
* Interpolate the griddata values to a vector (of evaluation-element 
* values).
* This routine extracts the necessary data from the griddata structure 
* and passes them on (to an sparse.c routine) for further processing, 
* since the grid doesn't know what a sparse matrix is.
* ==================================================================== */ 
void interpolateFromGridData(void *gridDataIn,double *dest,void *interpMat)
{
  fftGridData *gridData;
  double *src;
  int matNum=0;

  /* Cast input data to correct type */
  gridData = (fftGridData *) gridDataIn;
  /* Extract information on the starting point of the array of 
   * grid data. The first element in dest should be 
   * dest[0]=data[1][1][1]                                            */
  src = &(gridData->data[1][1][1]);
  /* Make a sparse matrix-vector multiply (projctMat times src)
   * Store result in dest, i.e. in gridData->data                     */
  spTransposeMatMultiply(interpMat,src,dest,matNum);
  /* The data should now be projected onto the grid                   */
  return;
}

/*
* ==================================================================== 
* Estimate the amount of memory that will be required to realize the 
* grid data and kernel data for a grid of a specified size.
* ==================================================================== */
double griddataMemoryEstimate(int nx, int ny, int nz, int verbose)
{
  char fctName[] = "griddataMemoryEstimate";
  int nnx,nny,nnz;
  double thisMem, totMem;
  /* Initialize memory counters.                                       */
  thisMem = totMem = 0.0;

  /* Estimate memory needs for grid data                               */
  nnx=nx;
  nny=ny;
  nnz=nz;
  S1S2S3INCREASE(nnx,nny,nnz,"data",fctName);
  thisMem = ( ((double)nnx)*((double)nny)*((double)nnz)
	      *((double)sizeof(double)) );
  if (verbose) printf("%s: Estimated memory for griddata   %10.3e Mb\n",
		      fctName, thisMem/pow(2.0,20.0));
  totMem+=thisMem;

  /* Estimate memory needs for kernel data                             */
  nnx=nx;
  nny=ny;
  nnz=nz;
  S1S2S3INCREASE(nnx,nny,nnz,"kernel",fctName);
  thisMem = ( ((double)nnx)*((double)nny)*((double)nnz)
	      *((double)sizeof(double)) );
  if (verbose) printf("%s: Estimated memory for gridkernel %10.3e Mb\n",
		      fctName, thisMem/pow(2.0,20.0));
  totMem+=thisMem;

  return totMem;
}




/* Test functions and stuff that are not really needed (but
 * which need access to the variables in this file) can go
 * in the following file: */
#include "nfft.testfunctions.c"



