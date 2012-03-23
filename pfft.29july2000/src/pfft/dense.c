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
/* FILE: dense.c
* 
* Routines for definitions of and operations on/with dense
* matrices, e.g. LU-factorization and back-solves for direct
* solutions of problems.
*
*/

#include <stdlib.h>
#include <stdio.h>
#include <float.h>

#include "Global.h"
#include "Dense.h"

/* For interaction with cLapack routine: */
#include "dgesvd.h"

/* The routines from fftcap have been used to some extend during the 
 * implementation of the present routines. However, fftcap uses the 
 * disk for storage, and defines special storage for triangular 
 * matrices (L or U parts). In the present code, an approach like 
 * the one in LAPACK is more convenient. However, the storage convention
 * for a dense matrix (SQDEX) is adapted from fftcap
 *
 * At the present LU-factorization WITHOUT PIVOTING has been implmented.
 * It is recommended to implement a version with partial pivoting.
 */

/* 
* ==================================================================== 
* Makes an LU-factorization witout pivoting. 
* A minimum of testing is performed to see only if a pivot is exactly 
* zero. It is generally a bad idea to use this routine, since it is 
* unstable in nature (see e.g. Trefethen and Bau), and the results 
* can be unreliable even when no pivots are exactly zero.
* This routine was taking in large from "blkLudecomp" in fftcap. In my 
* version of fftcap "blkLudecomp" is located in the file "blkDirect.c".
* ==================================================================== */
void LUdecompNoPivoting(double *mat,int n)
{
  char fctName[]="LUdecompNoPivoting";
  double factor;
  int i, j, k;
#if defined(MYDEBUGFLAG)
  int verbose = 0;
#endif

#if defined(MYDEBUGFLAG)
  if (verbose)
      fprintf(stderr, "%s: progress (%d rows total):\n",fctName,n);
#endif
  /* For all rows */
  for(k = 0; k < n-1; k++) {
    /* Check if the diagonal element is nonzero - it will beused to 
     * produce zeros in the entries below the diagonal in this column  */
    if(mat[SQDEX(k, k, n)] == 0.0) {
      fprintf(stderr, "%s: Zero diagonal element in row %d\n",fctName,k+1);
      exit(0);
    }
#if defined(MYDEBUGFLAG)
    if (verbose)
      fprintf(stderr,"%d ", k); /* Print progress */
#endif
    for(i = k+1; i < n; i++) { /* loop on remaining rows */
    
      factor = (mat[SQDEX(i, k, n)] /= mat[SQDEX(k, k, n)]);
      for(j = k+1; j < n; j++) { /* loop on remaining columns */
	mat[SQDEX(i, j, n)] -= (factor*mat[SQDEX(k, j, n)]);
      }
    }
  }
#if defined(MYDEBUGFLAG)
    if (verbose)
      fprintf(stderr, "\n");
#endif
  return;
}

/* 
* ==================================================================== 
* Uses an existing LU-factorization witout pivoting to solve the 
* system A*x = L*U*x = b
* On return, b is overwritten with x.
* ==================================================================== */
void LUbacksolveNoPivoting(double *mat,int n,double* b)
{
#if defined(MYDEBUGFLAG)
  char fctName[]="LUbacksolveNoPivoting";
  int verbose = 0;
#endif
  int i, j;

#if defined(MYDEBUGFLAG)
  if (verbose)
    fprintf(stderr, "%s: progress (%d rows total):\n",fctName,n);
#endif
  /* Solve L y = b using forward elimination.
   * Overwrite input b with y.                                         */
  /* For all rows */
  for(i = 0; i < n; i++) {
#if defined(MYDEBUGFLAG)
    if (verbose)
      fprintf(stderr,"%d ", i); /* Print progress */
#endif
    /* Find y[i] using forward elimination with already obtained 
     * unknowns (j<i): */
    for (j=0; j<i; j++){
      b[i] -= mat[SQDEX(i,j,n)] * b[j];
    }
  }
  /* Solve U x = y using backward elimination
   * Overwrite input b with x.                                         */
  for(i = n-1; i>-1 ; i--) {
    /* Find x[i] using backward elimination with already obtained 
     * unknowns (j>i): */
    for (j=i+1; j<n; j++){
      b[i] -= mat[SQDEX(i,j,n)] * b[j];
    }
    /* Normalize by diagonal entry */
    b[i] = b[i]/mat[SQDEX(i,i,n)];
  }
#if defined(MYDEBUGFLAG)
    if (verbose)
      fprintf(stderr, "\n");
#endif
  return;
}

/*
* ==================================================================== 
* Make a (pseudo)inverse of a dense matrix using CLAPACK SVD.
*
* Note the different definitions of a matrix here (double **) and 
* in other routines (double *).
*
* The matrix A must be of size at least [max(rows,cols)][max(rows,cols)] 
* since the pseudoinverse (and not its transpose) is returned in A.
*
* Singular value decomposition is used to calculate the psuedo-inverse.
* If #rows=#cols then (an approximation of) the inverse is found. 
* Otherwise an approximation of the pseudoinverse is obtained.
*
* The return value of this routine is the estimated rank of the system.
* ==================================================================== */ 
int pinv(double **A, int rows, int cols) 
{
  char fctName[] = "pinv_new";
  int i,j,k;
  int nsv, rank;
  double tol;
  /* Variables needed for interaction with CLAPACK routines 
   * (FORTRAN style)                                                   */
  double *amat, *svals, *U,  *V, *work;
  long int lda,        ldu, ldvt, lwork, m,n,info;

  /* Short-hands for this routine only [undef'ed at end of routine]    */
#define AA(i,j)  amat[i + j * (int) lda]
#define UU(i,j)  U[i + j * (int) ldu]
#define VV(i,j)  V[i + j * (int) ldvt]

  /* Set up arrays for call to CLAPACK routine                         */
  /* Use leading dimensions large enough to calculate both under- 
   * and over-determined systems of equations. The way this is 
   * calculated, we need to store U*(Z^-1^T) in U. The array U must 
   * have at least cols columns to perform this operation, and at 
   * least rows columns to store the initial U. However, for an under-
   * determined system the last cols-rows columns of U*(Z^-1^T) will be 
   * zero. Exploiting this, the size of the array U is just mxm        */
  m = (long int) rows;
  n = (long int) cols;
  lda  = m;
  ldu  = m;
  ldvt = n;

  nsv = MIN(rows,cols);
  /* Allocate memory for:
   * amat  : The matrix to factor
   * svals : Vector of singular values 
   * U     : mxm unitary matrix
   * V     : nxn unitary matrix                                        */
  amat  = (double *) calloc((int) m*n+1, sizeof(double)); 
  svals = (double *) calloc((int) m+n,   sizeof(double)); 
  U     = (double *) calloc((int) ldu*m+1,sizeof(double));
  V     = (double *) calloc((int) ldvt*n+1, sizeof(double));
  /* In matrix form amat = U * diag(svals) * V^T                       */
  /* Leading dimensions [number of rows] of the above arrays           */
  /* Work storage */
  lwork = 100*(m+n);
  work  = (double *) calloc(lwork, sizeof(double));

  /* Copy A to temporary storage                                       */
  for (j=0; j<cols; j++) {
    for (i=0; i<rows; i++) {
      AA(i,j) = A[i][j];
    }
  }

  {
    char jobu,jobvt;
    jobu  = 'A'; /* Calculate all columns of matrix U                  */
    jobvt = 'A'; /* Calculate all columns of matrix V^T                */

    /* Note that dgesvd_ returns V^T rather than V                     */
    dgesvd_(&jobu, &jobvt, &m, &n, 
	    amat, &lda, svals, U, &ldu, V, &ldvt, 
	    work, &lwork, &info);

    /* Test info from dgesvd */
    if (info) {
      if (info<0) {
	printf("%s: ERROR: clapack routine 'dgesvd_' complained about\n"
	       "    illegal value of argument # %d\n",fctName,(int) -info);
	_EXIT_;
      }
      else{
	printf("%s: ERROR: clapack routine 'dgesvd_' complained that\n"
	       "    %d superdiagonals didn't converge.\n",
	       fctName,(int) info);
	_EXIT_;
      }
    }
  }

  /* Test the the singular values are returned in correct ordering 
   * (This is done because there were problems with this with a former 
   * implementation [using f2c'ed linpack] when high optimization 
   * was used)                                                         */
  if (svals[0]<0.0) {
    printf("%s: ERROR: First singular value returned by clapack \n"
	   "   is negative: %16.6e.\n",fctName,svals[0]);
    _EXIT_;
  }
  for (i=1; i<nsv; i++){
    if ( svals[i] > svals[i-1] ) {
      printf("%s: ERROR: Singular values returned by clapack \n"
	     "  are not approprately ordered!\n"
	     "  svals[%d] = %16.6e > svals[%d] = %16.6e",
	     fctName,i,svals[i],i-1,svals[i-1]);
      _EXIT_;
    }
  }

  /* Test rank of matrix by examining the singular values.             */
  /* The singular values of the matrix is sorted in decending order, 
   * so that svals[i] >= svals[i+1]. The first that is zero (to some 
   * precision) yields information on the rank (to some precision)     */
  rank = nsv;
  tol  = DBL_EPSILON * svals[0] * MAX(rows,cols); 
  for (i=0; i<nsv; i++){
    if ( svals[i] <= tol ){
      rank = i;
      break;
    }
  }

  /* Compute (pseudo-) inverse matrix using the computed SVD:
   *   A   = U*S*V'
   *   A^+ = V * (Zi) U'
   *       = V * (U * (Zi)')'
   *       = { (U * Zi') * V' }'
   * Here Zi is the "inverse" of the diagonal matrix S, i.e. Zi is 
   * diagonal and has the same size as S'. The non-zero entries in 
   * Zi is calculated from the non-zero entries in S as Zi_ii=1/S_ii
   * Note that Zi' is of the same size as S.        
   * The last line here is used in the present computation. This 
   * notation avoids any need to transpose the output from the CLAPACK
   * rotines (which deliver U, S and V).                               */

  /* Inverse of [non-zero part of] diagonal matrix                     */
  tol = 1.0e-10 * svals[0]; 
  for (i=0; i< nsv; i++) {
    if (svals[i] < tol)
      svals[i] = 0.0; 
    else 
      svals[i] = 1.0 / svals[i]; 
  }

  /* Calculate  UZ = U * Zi', ie. scale COLUMN j in U by 
   * the j'th singular value [nsv columns only - since the diagonal 
   * matrix in general is not square]. If rows>cols then the last 
   * columns in UZ wil be zero (no need to compute).                   */
  for (j=0; j<nsv; j++){
    for (i=0; i<rows; i++){
      UU(i,j) *= svals[j];
    }
  }
  /* U*Zi' is stored in array U. It has size (rows x nsv). 
   * If cols>rows, then it should be though of as the larger matrix 
   * (rows x cols) with zero columns added on the right.               */

  /* Zero out the full array A to avoid confusion upon return.
   * This is not abosolutely necessary, zeroin out A[j][i] could be
   * part of the next loop.                                            */
  for (i=0; i< MAX(cols,rows); i++) {
    for (j=0; j< MAX(cols,rows); j++) {
      A[i][j] = 0.0;
    }
  }

  /* Matrix-matrix multiply  (U*Zi') * V'. 
   * Only the first nsv columns in U*Zi' are non-zero, so the inner 
   * most loop will go to k=nsv (-1). 
   * The result will be the transpose of the (psuedo-) inverse of A, 
   * so store directly in A'=A[j][i]. 
   *  A[j][i] = sum_k (U*Zi)[i][k] * V'[k][j] :                        */
  for (i=0; i< rows; i++) {
    for (j=0; j< cols; j++) {
      /*A[j][i] = 0.0; */ /* Include if A is not zeroed out above      */
      for (k=0;k<nsv;k++){
	A[j][i] += UU(i,k)*VV(k,j);
      }
    }
  }

  /* Free memory allocated in this routine */
  FREE(amat);
  FREE(svals);
  FREE(U);
  FREE(V);
  FREE(work);

  return rank;

  /* Undefine macros for this routine */
#undef AA
#undef UU
#undef VV
} /* End of routine pinv */



/* 
* ==================================================================== 
* Dump a matrix to stdout.
* ==================================================================== */
void dumpDenseMatrix(double *mat,int n)
{
  int i,j;
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      fprintf(stdout," %g",mat[SQDEX(i,j,n)]);
    }
    fprintf(stdout,"\n");
  }
  return;
}

/* ==== TEST FUNCTIONS (not really needed) ==== */ 
/* Test functions and stuff that are not really needed (but
 * which need access to the variables in this file) can go
 * in the following file: */
#include "dense.testfunctions.c"

