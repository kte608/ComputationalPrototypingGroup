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
/* FILE: gmres.c
*
* This file contains routines necessary to perform an iterative 
* solution of a system of linear algebraic equations using a 
* Generalized Minimum Residual method. 
*
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
# include <time.h> /* timing routines may eventually not be needed ? */

#include "Global.h"
#include "Gmres.h"

/*
* ====================================================================
* Solve a system of linear algebraic equations A*x=b using a 
* Generalized Minimum Residual method.
* Note that this routine does not deal with preconditioning. 
* Preconditioning is left for the calling routine and the 
* routine in charge of the matrix-vector products.
*
* This routine originates from fftcap. Multiple people have been 
* involved in its making.
*
* Note that the vector indices in this routine go from 0 to N-1 
* ("C-style"), rather than from 1 to N ("Fortran-style" or 
* "black-board style"), the latter being the convenvetion in the 
* gmres routine in fftcap.
*
* The calling routine will supply 
*  solutionVector: an initial guess vector - Solution vector on return,
*  rightHandSide:  a right-hand-side, 
*  size:           the system size,
*  tol:            a maximum tolerance on the solution,
*  maxiter:        a maximum iteration count, and
*  matmulFunction: a routine name for performing matrix-vector products.
*
* However, it must be the responsibility of the GMRES routine to 
* calculate all the needed residuals - including the first one(!) - 
* by performing the necessary matrix-vector products (i.e. calling 
* the named routine to do so).
*
* The matrix-vector-multiply routine should perform a product 
* dest := dest + A * src, dest and src being vectors, src being 
* unrelated to the so-called "source elements".
*
* By default this routine will setup necessary work space when it is 
* called, and free the memory upon return. However, this behaviour 
* may be changed later.
*
* The following references are of value if you try to figure out 
* what this routine is doing:
*
*    Yousef Saad (1996) 
*    Iterative Methods for Sparse Linear Systems, 
*    The PWS Series in Computer Science,
*    PWS Publishing Company, Boston 
* (on GMRES and practical implementation issues, pp 158-165)
*
*    Trefethen and Bau (1997), 
*    Numerical Linear Algebra, 
*    SIAM, 
* (on Arnoldi Iteration p250-256, and GMRES p266-275)
*
* ==================================================================== */
/*Original call (sys, q, p, r, ap, bv, bh, size, maxiter, tol, chglist)*/
int gmres(double *solutionVector, double *rightHandSide,
	  int size, double tol, int maxiter,
	  void (*matmulFunction)(double *src,
				 double *dest,
				 int size,
				 void *mData),
	  void *matmulData)
     /*
ssystem *sys;
double *q, *p, *ap, tol;
double **bv, **bh;
int size, maxiter;
charge *chglist;
     */
{
  char fctName[] = "gmres";
  int iter, i, j;
  double rnorm, norm;
  double hi, hip1, length;

  /* These could be moved to static variables outside this routine, 
   * such that they can be easily deallocated?                         */
  double 
    *c=NULL,          /* Diagonal non-unity entry in #iter rotation 
		       * matrix "Omega"_iter, i.e. "c_i" in Saad (1996)
		       * p. 161-165.                                   */
    *s=NULL,          /* Off-diagonal non-zero entry in #iter rotation 
		       * matrix "Omega"_iter, i.e. "s_i" in Saad (1996)
		       * p. 161-165.                                   */
    *y=NULL,
    /* New names */
    *residual,        /* Residual vector                               */
    *projection,      /* Projection vector - (normalized)              */
    *rnorms,          /* Norm of (previous) residual vectors(?)        */
    *av,              /* Work vector - When doing matrix-vector 
		       * products this vector is often the output.
		       * (ap is a shorthand for "A * vector")          */
    **backVectors,    /* Set of orthogonal projection vectors used     */
    **arnoldiHessenberg; /* Hessenberg matrix growing in the Arnoldi 
			  * iteration.                                 */
  /* Temporary - for timing of matrix-vector products: */
  double timeA, timeB, totalTime;
  int numProducts, countTime=1;

  totalTime=0.0;
  numProducts=0;

  /* Allocation first time through. */
  if (c == NULL) {
    c = (double*) calloc(size,sizeof(double)); /* -> ? */
    s = (double*) calloc(size,sizeof(double)); /* -> ? */
    /*g = (double*) calloc(size,sizeof(double));  -> "rnorms"   */
    y = (double*) calloc(size,sizeof(double)); /* -> ? */
    /*r = (double*) calloc(size,sizeof(double));  -> "residual" */
    /*p = (double*) calloc(size,sizeof(double));  -> "projection" (?)*/
    /* New names: */
    residual    = (double*) calloc(size,sizeof(double)); /* "r"  */
    projection  = (double*) calloc(size,sizeof(double)); /* "p"  */
    rnorms      = (double*) calloc(size+1,sizeof(double)); /* "g"  */
    av          = (double*) calloc(size,sizeof(double)); /* "ap" */

    backVectors 
      = (double**) calloc(maxiter,sizeof(double*));      /* bv */
    arnoldiHessenberg 
      = (double**) calloc(maxiter,sizeof(double*));      /* bh */
  }
  /* Calculate initial residual                                        */
  if (countTime) timeA = ((double) clock());
  matmulFunction(solutionVector,av,size,matmulData);
  if (countTime){ 
    timeB = ((double) clock()); 
    totalTime+= timeB-timeA;
    numProducts++;
  }
  VMINUS(residual,rightHandSide,av,size);
  INNER(rnorm, residual, residual, size);
  rnorm = sqrt(rnorm);

  if(rnorm<=tol){
    /* In this case the initial guess is really good!
     * This will (probably) never, happen, but we should be able to
     * capture this special case! In any case, we definitely need an
     * if-statement here, just to capture rnorm==0!                    */
    /* Just return the initial guess (solutionVector) unmodified. */
    return 0;
  }
  /* Set up v^1 and g^0. */
  for(i=0; i<size; i++) {
    projection[i] = residual[i] / rnorm;
    rnorms[i] = 0.0;
  }
  rnorms[0] = rnorm;

#if ITRDAT == ON
  fprintf(stdout, "||res|| = %g\n", rnorm); /* initial guess residual norm */
#endif


  for(iter = 0; 
      (iter < maxiter) && (rnorm > tol) && (iter < size); 
      iter++) {
    
    /* allocate the back vectors if they haven't been already          */
    if(backVectors[iter] == NULL) {
      backVectors[iter]       =  (double*) calloc(size,sizeof(double));
      arnoldiHessenberg[iter] =  (double*) calloc(iter+2,sizeof(double));
    }

    /* Save projection as the v{iter} (i.e. backVector number iter).
     * This is the projection vector calculated in the previous loop.  */
    for(i=0; i < size; i++) backVectors[iter][i] = projection[i];

    /* The following is an Arnoldi iteration step (modified 
     * Gram-Schmidt iteration), see e.g.: 
     *  Trefethen and Bau (1997), p252 
     *  "Algorithm 33.1. Arnoldi Iteration".
     * It should be noted that "arnoldiHessenberg" is stored 
     * column-by-column, i.e. arnoldiHessenberg should be read 
     * as "arnoldiHessenberg[column][row]".                            */
    /* Form Av{iter}, i.e. the matrix times the previous projection 
     * vector                                                          */
    VZERO(av,size);
    if (countTime) timeA = ((double) clock());
    matmulFunction(projection,av,size,matmulData);
    if (countTime){ 
      timeB = ((double) clock()); 
      totalTime+= timeB-timeA;
      numProducts++;
    }
    /* Initialize the projection vector to Av^{iter}, 
     * i.e. initiialize the new projection vector for this loop        */
    for(i=0; i < size; i++) projection[i] = av[i];
    /* Make the projection vector orthogonal to all the previously 
     * obvained projection vectors (v^{j}, 0 <= j <= iter), 
     * one at a time.                                                  */
    for(j=0; j <= iter; j++) {
      INNER(hi, av, backVectors[j], size);
      for(i=0; i < size; i++) projection[i] -= hi * backVectors[j][i];
      arnoldiHessenberg[iter][j] = hi;
    }
    /* Normalize the new projection.                                   */
    INNER(norm, projection, projection, size);    
    norm = sqrt(norm);
    if (norm > 0.0) 
      for(i=0; i<size; i++) projection[i] /= norm;
    /* Last entry in column #iter of the Hessenberg matrix             */
    arnoldiHessenberg[iter][iter+1] = norm;

    /* Apply earlier rotations "Omega_i" (i=0 ... iter-1) 
     * to the new column in the Hessenberg matrix.                     */
    /*for(i=0; i<iter-1; i++) { This is wrong! */
    for(i=0; i<iter; i++) { /* This is correct! */
      hi   = arnoldiHessenberg[iter][i];
      hip1 = arnoldiHessenberg[iter][i+1];
      arnoldiHessenberg[iter][i]   = c[i] * hi   - s[i] * hip1;
      arnoldiHessenberg[iter][i+1] = c[i] * hip1 + s[i] * hi;
    }

    /* Compute new rotations. 
     * For the rotation matrix "Omega", these are for the only entries 
     * which differ from the entries in the identity matrix.
     * See e.g. Saad (1996) pp 161-165                                 */
    hi = arnoldiHessenberg[iter][iter];
    hip1 = arnoldiHessenberg[iter][iter+1];
    length = sqrt(hi * hi + hip1 * hip1);
    c[iter] = hi/length;
    s[iter] = -hip1/length;

    /* Apply new rotations. */
    arnoldiHessenberg[iter][iter]   = c[iter] * hi - s[iter] * hip1;
    arnoldiHessenberg[iter][iter+1] = c[iter] * hip1 + s[iter] * hi;
    hi = rnorms[iter];
    rnorms[iter] = c[iter] * hi;
    rnorms[iter+1] = s[iter] * hi;

    rnorm = ABS(rnorms[iter+1]);

#if ITRDAT == ON
    fprintf(stdout, "||res|| = %g\n", rnorm);
#else
    fprintf(stdout, "%d ", iter);
    if((iter) % 15 == 0 && iter != 0) fprintf(stdout, "\n");
#endif
    fflush(stdout);


  }

  /* Decrement from the last increment. */
  iter--;
printf(" Used iter = %d\n",iter);
  /* Compute solution, note, the Hessenberg matrix 
   * is arnoldiHessenberg[col][row].                                   */
  for(i=0; i <= iter; i++) y[i] = rnorms[i];
  for(i = iter; i >= 0; i--) {
    y[i] /=  arnoldiHessenberg[i][i];
    for(j = i-1; j >= 0; j--) {
      y[j] -= arnoldiHessenberg[i][j]*y[i];
    }
  }
  for(i=0; i < size; i++) {
    solutionVector[i] = 0.0;
    for(j=0; j <= iter; j++) {
      solutionVector[i] += y[j] * backVectors[j][i];
    }
  }

  /* Free memory */
  for (i=0; i<maxiter; i++){
    if(backVectors[i]!=NULL) FREE(backVectors[i]);
    if(arnoldiHessenberg[i]!=NULL) FREE(arnoldiHessenberg[i]);
  }
  FREE(backVectors);
  FREE(arnoldiHessenberg);
  FREE(c);
  FREE(s);
  FREE(y);
  FREE(residual);
  FREE(projection);
  FREE(rnorms);
  FREE(av);

  if(rnorm > tol) {
    fprintf(stdout, "\n%s: WARNING exiting without converging\n",fctName);
  }

  if(countTime) printf(" %d matrix-vector products in %g secs.\n"
		       " Time per product estimated to %g secs.\n",
		       numProducts, totalTime/((double) CLOCKS_PER_SEC),
		       totalTime/(numProducts * ((double) CLOCKS_PER_SEC)));

printf("%s: Returning\n",fctName);
  return iter+1;
} /* End of routine gmres */
