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
/* FILE: solve.c
*
*  This file contains routines that deals with solving 
*  the system of equations for a given right hand side.
*  Calls will be made to routines in other files when
*  specific information is needed.
*  In principle, very little information on specifics 
*  is needed here.
*
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
# include <time.h> /* timing routines may eventually not be needed ? */

#include "Global.h"

/* Structure for bundling data for passing through the gmres solver */
struct fullSolveBundle {
  /* Pointers to matrices: */
  void *directMat, *interpMat, *projectMat, *precondMat; 
  /* Preconditioner type ("none"/"right"/"left") */
  char *precondType;
  /* Pointer to the grid */
  void *grid;
  /* Number of columns and rows in the direct part. Also the 
   * length of the values of sources and evals, respectively. */
  int numColsDirect, numRowsDirect;
  /* Matrix number to use (direct + precond) */
  int matNum;
  /* Work vectors */
  double *workSrc, *workEval;
};

/* 
* ====================================================================
* Solve the full system of equations including direct and grid parts.
*
* Note: This routine does not "need to know" about definitions of 
* grid, sparse matrices, elements etc. 
*
* This is still in TEST PHASE. I need to implement routines for dealing 
* with boundary conditions (in "elements.c"?) and setting up a 
* right-hand-side. These data should be passed on by the calling routine!
*
* Maybe later implement restarted GMRES to decrease memory consumption?
* ==================================================================== */
void solveFull(void *directMat, void *interpMat, void *projectMat,
	       void *precondMat, char *precondType, void *grid,
	       double *sourceVector, double *destVector,
	       int matNum, 
	       int numColsDirect, int numRowsDirect)
{
  char fctName[] = "solveFull";
  int testSolution = 0; /* Perform a final matrix-vector product to 
			 * confirm the soltion accuracy.               */
  int i, iter;
  double norm1, norm2, tol2;
  struct fullSolveBundle bundleData[1];
  double err, *testVector=NULL, *solutionVectorTmp, *rightHandSideTmp;

  /* These should really be passed from the calling routine            */
  double tol = 1.0E-4, *solutionVector, *rightHandSide;
  int modifyTol=1; /* 1 if tol may be modified based on the RHS norm */
  int maxiter  = 100; /* Maximum number of iteration in the gmres      */
  /* Never do more iterations than the number of unknowns              */
  maxiter = MIN(maxiter, numColsDirect);

  /*
  solutionVector   = (double*) calloc(numColsDirect,sizeof(double));
  rightHandSide    = (double*) calloc(numColsDirect,sizeof(double));
  */
  solutionVector   = destVector;
  rightHandSide    = sourceVector;
  


/* JKW DUMMY = unit vector right-hand-side: */
/*
  fprintf(stderr,"%s: WARNING: Using dummy right-hand-side!\n",fctName);
  for (i=0; i<numColsDirect; i++)
    rightHandSide[i] = ((double) rand())/((double) RAND_MAX);
  for (i=0; i<numColsDirect; i++)
    rightHandSide[i] = 0.0;
  rightHandSide[0] = 1.0;
*/

  printf("%s: Solving system with preconditioner: %s\n",
	 fctName,precondType);

  /* Bundle data for the "matrix multiply routine". Use a structure 
   * for the bundle, so that it cannot be messed up easily.            */
  bundleData->directMat     = directMat;
  bundleData->interpMat     = interpMat;
  bundleData->projectMat    = projectMat;
  bundleData->precondMat    = precondMat;
  bundleData->precondType   = precondType;
  bundleData->numColsDirect = numColsDirect;
  bundleData->numRowsDirect = numRowsDirect;
  bundleData->grid          = grid;

  /* Presently, only work on the first matrix. 
   * This *could* be changed at a later point. */
  bundleData->matNum = matNum;
  /* Allocate working vector for the "matrix multiply routine" */
  /*  Work vectors for eval values for direct and grid parts           */
  if (strcmp(precondType,"right")==0){
    /* Make an extra work vector to store source values                */
    bundleData->workSrc  = (double *) calloc(numColsDirect,sizeof(double));
  }
  else if (strcmp(precondType,"left")==0){
    /* Make an extra work vector to store eval values                    */
    bundleData->workEval  = (double *) calloc(numRowsDirect,sizeof(double));
  }
  else if (strcmp(precondType,"none")!=0){
    /* This is everything else, which is unknown = bad! 
     * Print error and exit */
    fprintf(stderr,"%s: ERROR. Unknown preconditioner type\n"
	    "    Please use \"none\", \"left\" or \"right\"\n"
	    "    Exiting\n",
	    fctName);
    exit(1);
  }

  /* === Call iterative solver === */
  /* Possibly repeat the call to the iterative solver.
   * Note that the memory consumption is O(maxiter^2*N) and the 
   * CPU consumption is O(maxgmres*maxiter^2*N), so it may pay out 
   * to use restarted gmres.                                           */
  if (strcmp(precondType,"none")==0){
    /*  NO precond, solve:     Ax  = b */
      iter = gmres(solutionVector, rightHandSide,
		   numColsDirect, tol, maxiter,
		   fullMatrixMultiply,
		   (void*) bundleData);
  }
  else if (strcmp(precondType,"right")==0){
    /*  RIGHT precond, solve:     APy  = b */
    /* Modify stop criterion based on the 2-norm of the 
     * right-hand-side vector                                         */
    for (i=0,norm1=0.0;i<numColsDirect;i++)
      norm1 += rightHandSide[i]*rightHandSide[i];
    norm1 = sqrt(norm1);
    tol2 = tol;
    if (norm1>0.0 && modifyTol){
      fprintf(stdout,"%s: Scaling tolerance with ||b|| = %g\n",
	      fctName,norm1);
      tol2 *= norm1;
      fprintf(stdout,
	      "%s: Aiming for tolerance %g, rather than the specified %g\n",
	      fctName,tol2,tol);
    }
    solutionVectorTmp = (double*) calloc(numColsDirect,sizeof(double));
    iter = gmres(solutionVectorTmp, rightHandSide,
		 numColsDirect, tol2, maxiter,
		 fullMatrixMultiply,
		 (void*) bundleData);
    /*  Then calculate x = Py. */
    spMatMultiply(precondMat, solutionVectorTmp, 
		  solutionVector, matNum);
    FREE(solutionVectorTmp); 
  }
  else if (strcmp(precondType,"left")==0){
    /* LEFT preconditioner - calculate y := Pb. */
    rightHandSideTmp = (double*) calloc(numColsDirect,sizeof(double));
    spMatMultiply(precondMat,rightHandSide,rightHandSideTmp,matNum);
    /* Note that rightHandSideTmp and rightHandSide do not have the 
     * same norm. In fact they may differ considerably, so make sure 
     * this is taken into account when selecting the accuracy of the 
     * solution. At the moment I'm slightly at a loss on how to do 
     * this accurately. For now, the tolerance will be scaled by 
     * ||Pb||/||b||, where P is the preconditioner and b is the right 
     * hand side (Pb is the new right hand side). 
     * Issue a warning statement about this                        */
    fprintf(stderr,"%s: WARNING: Scaling tolerance with ||Pb||/||b||!\n",
	    fctName);
    for (i=0,norm1=0.0,norm2=0.0;i<numColsDirect;i++){
      norm1 += rightHandSide[i]*rightHandSide[i];
      norm2 += rightHandSideTmp[i]*rightHandSideTmp[i];
    }
    norm1 = sqrt(norm1);
    norm2 = sqrt(norm2);
    if(norm1>0.0 && modifyTol){
      tol2 = tol*norm2/norm1; 
      fprintf(stderr,
	      "%s: Aiming for tolerance %g, rather than the specified %g!\n",
	      fctName,tol2,tol);
    }
    else tol2 = tol;
    /* Now solve PAx  = y */
    iter = gmres(solutionVector, rightHandSideTmp,
		 numColsDirect, tol2, maxiter,
		 fullMatrixMultiply,
		 (void*) bundleData);
    FREE(rightHandSideTmp);
  }
  printf("%s: GMRES used %d iterations\n",fctName,iter);

  /* Test the solution making a matrix-vector product 
   * (no preconditioning) A*x and comparing with the right-hand-side   */
  if (testSolution) {
    bundleData->precondType   = "none";
    testVector = (double*) calloc(numColsDirect,sizeof(double));
    fullMatrixMultiply(solutionVector,testVector,numColsDirect,
		       (void*)bundleData);
    for (i=0, err=0.0; i<numColsDirect; i++)
      err += (testVector[i]-rightHandSide[i])
	*(testVector[i]-rightHandSide[i]);
    err=sqrt(err);
    printf(" Estimated ||err||:  %g\n",err);
  }

  /* TEST SECTION */
  if (0){
    if (0){
      for (i=0; i<numColsDirect; i++){
	printf("%g %g\n",solutionVector[i],rightHandSide[i]);
      }
    }
    /*dumpSparseMatrixTranspose(precondMat,0);*/
    if (0){
      int j;
      for (j=0;j<numColsDirect;j++){
	for (i=0; i<numColsDirect; i++) testVector[i]=solutionVector[i]=0.0;
	testVector[j]=1;
	fullMatrixMultiply(testVector,solutionVector,numColsDirect,
			   (void*)bundleData);
	printf("Column %d: \n",j);
	for (i=0; i<numColsDirect; i++)
	  printf(" %12.6e",solutionVector[i]);
	printf("\n");
      }
    }
    if (0){/* Dump transpose matrix to file: */
      int j;
      FILE *dumpFile;
      dumpFile = fopen("mat.out","w");
      printf("Writing to mat.out\n");

      for (j=0;j<numColsDirect;j++){
	for (i=0; i<numColsDirect; i++) testVector[i]=solutionVector[i]=0.0;
	testVector[j]=1;
	fullMatrixMultiply(testVector,solutionVector,numColsDirect,
			   (void*)bundleData);
	for (i=0; i<numColsDirect; i++)
	  fprintf(dumpFile," %15.8e",solutionVector[i]);
	fprintf(dumpFile,"\n");
      }
      fclose(dumpFile);
   }
    if (0){ /* Calculate and write A*e_1: */
      int j;
      FILE *dumpFile;
      dumpFile = fopen("a_times_e1.out","w");

      j=0;
      for (i=0; i<numColsDirect; i++) testVector[i]=solutionVector[i]=0.0;
      testVector[j]=1;
      fullMatrixMultiply(testVector,solutionVector,numColsDirect,
			 (void*)bundleData);
      printf("Writing to %s\n","a_times_e1.out");
      for (i=0; i<numColsDirect; i++)
	fprintf(dumpFile," %12.6e\n",solutionVector[i]);
    }
  }
  /* Test a preconditioner. For a small system P*A=I, so PA*b=b 
  for (i=0; i<numColsDirect; i++) testVector[i]=solutionVector[i]=0.0;
  bundleData->precondType   = "none";
  testVector = (double*) calloc(numColsDirect,sizeof(double));
  fullMatrixMultiply(rightHandSide,testVector,numColsDirect,
		     (void*)bundleData);
  spMatMultiply(precondMat,testVector,solutionVector,matNum);
  for (i=0; i<numColsDirect; i++)
    printf(" %16.8e %16.8e %16.8e \n",rightHandSide[i],solutionVector[i],testVector[i]);
  */



  /* Free working vector and other stuff used. */
  IFFREE(testVector);
  /*FREE(solutionVector);
    FREE(rightHandSide);*/
  if (strcmp(precondType,"right")==0){
    FREE(bundleData->workSrc);
  }
  else if (strcmp(precondType,"left")==0){
    FREE(bundleData->workEval);
  }
      
  return; /* Should probably return (a pointer to) the solution vector? */
} /* End of routine solveFull */

/* 
* ====================================================================
* This is simply a wrapper for "fullMatrixMultiply".
* TODO: This routine has too many similarities with "solveFull". 
* Put the similar stuff in static routines callable by both the present 
* routine and by solveFull.
* ==================================================================== */
void systemMatrixMultiply(void *directMat, void *interpMat, 
			  void *projectMat, void *precondMat, 
			  char *precondType, void *grid,
			  double *sourceVector, double  *destVector,
			  int matNum, 
			  int numColsDirect, int numRowsDirect)
{
  char fctName[] = "systemMatrixMultiply";
  struct fullSolveBundle bundleData[1];
  /* Bundle data for the "matrix multiply routine". Use a structure 
   * for the bundle, so that it cannot be messed up easily.            */
  bundleData->directMat     = directMat;
  bundleData->interpMat     = interpMat;
  bundleData->projectMat    = projectMat;
  bundleData->precondMat    = precondMat;
  bundleData->precondType   = precondType;
  bundleData->numColsDirect = numColsDirect;
  bundleData->numRowsDirect = numRowsDirect;
  bundleData->grid          = grid;
  bundleData->matNum        = matNum;
    /* Allocate working vector for the "matrix multiply routine" */
  /*  Work vectors for eval values for direct and grid parts           */
  if (strcmp(precondType,"right")==0){
    /* Make an extra work vector to store source values                */
    bundleData->workSrc  = (double *) calloc(numColsDirect,sizeof(double));
  }
  else if (strcmp(precondType,"left")==0){
    /* Make an extra work vector to store eval values                    */
    bundleData->workEval  = (double *) calloc(numRowsDirect,sizeof(double));
  }
  else if (strcmp(precondType,"none")!=0){
    /* This is everything else, which is unknown = bad! 
     * Print error and exit */
    fprintf(stderr,"%s: ERROR. Unknown preconditioner type\n"
	    "    Please use \"none\", \"left\" or \"right\"\n"
	    "    Exiting\n",
	    fctName);
    exit(1);
  }

  /* Call routine to do the dirty work */
  fullMatrixMultiply(sourceVector,destVector,
		     (int) NULL, (void*) bundleData);

  /* Free allocated memory */
  if (strcmp(precondType,"right")==0){
    FREE(bundleData->workSrc);
  }
  else if (strcmp(precondType,"left")==0){
    FREE(bundleData->workEval);
  }
  return;
}
/* 
* ====================================================================
* "Matrix-vector multiply" routine for iterative solution.
* For a test vector "src" return the matrix-vector product 
* dest = A * src calculated as the sum of grid and direct parts. 
* Note that "size" is not actually used. It is included for 
* possible future use.
* ==================================================================== */
void fullMatrixMultiply(double *src,double *dest,
			int size,
			void *inputData)
{
  char fctName[] = "fullMatrixMultiply";
  struct fullSolveBundle *bundleData;
  void *directMat, *interpMat, *projectMat, *precondMat, *grid;
  char *precondType;
  int matNum, numRowsDirect, numColsDirect, i;
  double  *workSrc, *workEval;
  double timeA, timeB;
  static double totalTime=0.0;
  int countConvolveTime = 0; /* Set to 1 to see timings for convolution part */

  /* Cast input data to correct type */
  bundleData = (struct fullSolveBundle*) inputData;
  /* Local shorthands. Get information from the bundled data */
  directMat     = bundleData->directMat;
  interpMat     = bundleData->interpMat;
  projectMat    = bundleData->projectMat;
  precondMat    = bundleData->precondMat;
  precondType   = bundleData->precondType;
  numRowsDirect = bundleData->numRowsDirect;
  numColsDirect = bundleData->numColsDirect;
  grid          = bundleData->grid;
  matNum        = bundleData->matNum;

  /* Make sure that the destination vector contains only zeros. 
   * This is needed since the sparse matrix-vector product routines 
   * ADD a matrix-vector product to a vector.                          */
  if (strcmp(precondType,"none")==0){
    for (i=0;i<numRowsDirect;i++){
      dest[i] = (double) 0.0;
    }
    /* No working vectors needed in this case:                         */
    workEval = dest;
    workSrc  = src;
  }  
  else if (strcmp(precondType,"right")==0){
    workSrc    = bundleData->workSrc;
    for (i=0;i<numRowsDirect;i++){
      dest[i] = (double) 0.0;
    }
    for (i=0;i<numColsDirect;i++){
      workSrc[i] = (double) 0.0;
    }
    /* No EVAL working vector needed in this case: */
    workEval = dest;
  }
  else if (strcmp(precondType,"left")==0){
    workEval = bundleData->workEval;
    for (i=0;i<numRowsDirect;i++){
      dest[i] = (double) 0.0;
      workEval[i] = (double) 0.0;
    }
    /* No SOURCE working vector needed in this case: */
    workSrc  = src;
  }
  else {
    /* This is everything else, which is unknown = bad! 
     * Print error and exit */
    fprintf(stderr,"%s: ERROR. Unknown preconditioner type\n"
	    "    Please use \"none\", \"left\" or \"right\"\n"
	    "    Exiting\n",
	    fctName);
    exit(1);
  }

  /* --- PRECONDITIONING --- */
  if (strcmp(precondType,"right")==0){
    /* If RIGHT preconditioning is used, then start out by 
     * multiplying the trial vector (src) by the preconditioner        */
    /* Calculate P*src, store in workSrc:                              */
    workSrc    = bundleData->workSrc;
    spMatMultiply(precondMat,src,workSrc,matNum);
  }

  /* --- GRID PART ---  */
  /* Projection step. Project the source element values given 
   * by src to the grid. Stores result in griddata, which is an 
   * array internal to the grid. 
   * Maybe this routine should have been called "src2gridCalculate"?   */
  projectOntoGrid(grid, workSrc, projectMat);

 
  if (countConvolveTime) timeA = ((double) clock());
  /* Grid-to-grid part. This is (at the moment of writing) a 
   * convolution of the grid data with the kernel data using the FFT.  */
  grid2gridCalculate(grid);
  if (countConvolveTime) { 
    timeB = ((double) clock());
    totalTime += timeB-timeA;
    printf("%s: Total time for grid2grid is now: %g  secs. (%g)\n",
	   fctName,
	   totalTime/((double) CLOCKS_PER_SEC),
	   (timeB-timeA)/((double) CLOCKS_PER_SEC));
  }
  /* Interpolation step. Interpolate grid values to values at the 
   * evaluation elements. Add the interpolated values to workEval.
   * Maybe this routine should have been called "grid2evalCalculate"?  */
  interpolateFromGrid(grid,workEval,interpMat);

  /* --- DIRECT PART --- */
  /*  ADD to the grid part the direct matrix-vector product:           */
  spMatMultiply(directMat,workSrc,workEval,matNum);

  /* --- PRECONDITIONING --- */
  if (strcmp(precondType,"left")==0){
  /* If LEFT preconditioning is used, then end by multiplying 
   * the solution vector (workDest) by the preconditioner.
   * Store in the dest destination (output) array.                     */
    spMatMultiply(precondMat,workEval,dest,matNum);
  }

  return;
} /* End of routine fullMatrixMultiply */

/* 
* ====================================================================
* Allocate vector either for e.g. right-hand-side or solution.
* Note that other pieces of the code may use information of the 
* structure (i.e. be assuming that it is a double precision array).
* ==================================================================== */
void *allocateSystemVector(int nEntries){
  return (void *) calloc(nEntries, sizeof(double));
}
/* 
* ====================================================================
* Deallocate system vector.
* ==================================================================== */
void freeSystemVector(void *vector){
  double *v;
  v = (double *) vector;
  FREE(v);
  return;
}
