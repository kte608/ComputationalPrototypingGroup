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
/* FILE: gridData.c
*
*  NOTE: The routines for the grid interaction.
*
* To ease the implementation, these routines will "know" 
* the sparse matrix structures. Note however, that these 
* routines do not "know" the grid structure itself. Only 
* the structure of the DATA used on the grid.
*/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "Global.h"
#include "Sparse.h"
#include "Dense.h"

/* Prototypes for static routines. */
static void calcNumNonzeroElementGrid(void *grid,
				      void *element,
				      char *Type,
				      int i, int j, int k, 
				      void *colsData);

static void calcElementGrid(void *grid,
			    void *element,
			    char *Type,
			    int i, int j, int k, 
			    void *matIn);

/* ==== Special structures needed in this file only ==== */
/* Structure for calling calcElementGrid */
struct cegBundle{
  void *mat;
  char *systemType;
};

/*
* ==================================================================== 
* Initialize the nonsquare sparse matrix which projects sources to the 
* grid or interpolates the grid onto evals.
* Figure out how much storage is needed, and subsequently allocate
* the needed memory. Then return a pointer to the matrix. 
* NOTE: the interpolation and projection matrices are not necessarily
* structurally symmetric, and can be allocated by separate calls to
* this function.
*
* The matrices are implemented with one column per element (source or 
* eval). The number of nonzeros needed in each column is then the 
* number of unknowns (degrees of freedom) related to the element times 
* the number of grid points in the interpolation/projection stencil.
*
* grid    - Grid pointer, passed through
* Type    - Typically indicates source or eval (project or interp). 
*            Passed through
* n       - Number of columns in the matrix.
* numMats - Number of matrices to allocate. Passed through to spMatAlloc
* ==================================================================== */ 
void *callocElementGridSparse(void *grid, char *Type, int n, int numMats)
{
  int *numNonzeroInCols;          /* # nonzeros in sparse matrix cols. */
  void *mat;
  int nCols;

  numNonzeroInCols  = (int *) calloc(nCols=n,sizeof(int));
  
  /* Count the number of nonzeros in each column by querying the grid. */
  elementAndGridOperate(grid, Type, calcNumNonzeroElementGrid,
			(void*) numNonzeroInCols);

  /* Allocate the matrix. */
  mat  = spMatAlloc(nCols, numNonzeroInCols, numMats);

  /* Free the generated vector. */
  FREE(numNonzeroInCols);

  /* Return generic pointer to the matrix. */
  return mat;
}

/*
* ==================================================================== 
* Calculates the number of nonzeros in a set of columns in the sparse
* matrix for mapping between elements and the grid.
* Works by asking the element how many unknowns it has, and multiplies
* by the size of the interpolation/projection stencil.
*
* grid    - Grid pointer.
* Type    - Typically indicates  source or eval (project or interp). 
*            Passed through.
* element - Element pointer.
* i,j,k   - Grid point index.
* ==================================================================== */ 
static void calcNumNonzeroElementGrid(void *grid,
			       void *element,
			       char *Type,
			       int ip, int jp, int kp, 
			       void *colsData)
{
  char fctName[]="calcNumNonzeroElementGrid";
  int numUnknowns, numStencilPoints;
  int *indices;
  int *numNonzeroInCols;
  int ii;

  if(strcmp(Type,"source") == 0) {
    numUnknowns = numSourceUnknowns(element);
    indices = sourceIndices(element);
    /* Get the number of points in the projection stencil              */
    /* BB: Note that later the interpolation and projection stencils may 
     * differ. Thus, different routines should be called for "source" 
     * and "eval". */
    numStencilPoints = numInterpPoints(grid);
  }
  else if(strcmp(Type,"eval") == 0) {
    numUnknowns = numEvalUnknowns(element);
    indices = evalIndices(element);
    /* Get the number of points in the interpolation stencil           */
    numStencilPoints = numInterpPoints(grid);
  }
  else {
    /* This is everything else - which is unknown=bad. 
     * Report error and exit                                           */
    fprintf(stderr,"%s: ERROR. Unknown type for element-grid stencils %s\n"
	    "    Please use \"source\" or \"eval\"\n"
	    "    Exiting\n",
	    fctName,Type);
    exit(1);
  }

  numNonzeroInCols = (int*) colsData;
  for(ii = 0; ii < numUnknowns; ii++) {
    numNonzeroInCols[indices[ii]] = numStencilPoints;
  }
  return;
} /* End of routine calcNumNonzeroElementGrid */

/*
* ==================================================================== 
* Calculate the values for the sparse matrix for
* interpolation/projection, and put those values in the matrix 
* (interpolation or projection respectively).
* This routine just calls a grid-controlled routine which makes 
* repeated calls to calcElementGrid for each element.
* grid - pointer to grid (passed through)
* Type - type of element (char *) passed through
* mat - pointer to the sparse matrix passed through
* ==================================================================== */ 
void calcElementGridSparse(void *grid, char *elementType, 
			   char *systemType, void *mat)
{
  struct cegBundle bundle[1];
  /* Bundle data for calculations */
  bundle->mat        = mat;
  bundle->systemType = systemType;

  /* Calculate the matrix elements. */
  elementAndGridOperate(grid, elementType, calcElementGrid, bundle);
  return;
}

/*
* ==================================================================== 
* Actually calculate the values for the sparse matrix for
* interpolation/projection, and puts those values in the matrix
* (interpolation or projection respectively).
* Works by first getting the interpolation weights matrix, which
* relates grid point values to polynomial coefficients.  
* Then, integrals over elements of polys are computed.  
* The matrix entries are then the product of the weights matrix and the 
* integral (or moments) matrix.
*
* grid       - generic pointer to grid (passed through)
* element    - generic pointer to element
* Type       - type of element (char *) passed through
* ip, jp, kp - grid indices for the grid points associated with element.
* matIn      - pointer to the sparse matrix 
* ==================================================================== */ 
static void calcElementGrid(void *grid,
			    void *element,
			    char *elementType,
			    int ip, int jp, int kp, 
			    void *bundleIn)
{
  char fctName[]="calcElementGrid";
  struct cegBundle *bundle;
  int i,j;
  int *is, *js, *ks; /* Stencil point indices. */
  int *imonoOrder,*jmonoOrder,*kmonoOrder;
  int nGridPoints, unknown;
  double **wts;
  sparseMat *mat;
  char *systemType;
  int *rowIndex, iColStart;
  double **matVals, temp;
  int numUnknowns, *indices;
  double x, y, z, **moments, polyScale;
  int nterms, polyOrder, maxOrder1D;

  /* Get data from bundle */
  bundle = (struct cegBundle *) bundleIn;
  systemType = bundle->systemType;

  /* Yeah, we know the sparse matrix data structure here. */
  mat = (sparseMat *) bundle->mat;

  /* Find the number of points (nGridPoints) in the interpolation or 
   * projection stencil, their relative position (is, js, ks), the 
   * number of polynomial terms (nterms) used and the contribution of 
   * each point to each polynomial term (wts) (wts = "weights" is a 
   * # gridpts times # poly terms matrix). Note that pointers to the 
   * data in grid are returned, so please do not corrupt these data.
   * No allocation or deallocation is needed. In fact, deallocating 
   * this storage would be a very bad idea!                            */
  nGridPoints = getGridInterp(grid, elementType, &is, &js, &ks, &wts, 
			      &imonoOrder, &jmonoOrder, &kmonoOrder,
			      &(nterms), &(polyScale),
			      &(polyOrder), &(maxOrder1D));

  /* Get the number of terms and their indices in the 
   * direct computations.                                              */
  if(strcmp(elementType,"eval") == 0) {
    numUnknowns = numEvalUnknowns(element);
    indices = evalIndices(element);
    /* Do not change the values of these indices[] !                   */
  }
  else if(strcmp(elementType,"source") == 0) {
    numUnknowns = numSourceUnknowns(element); 
    indices = sourceIndices(element);
    /* Do not change the values of these indices[] !                   */
  }
  else {
    /* This is everything else - which is unknown=bad. 
     * Report error and exit                                           */
    fprintf(stderr,"%s: ERROR. Unknown type for element-grid stencils %s\n"
	    "    Please use \"source\" or \"eval\"\n"
	    "    Exiting\n",
	    fctName,elementType);
    exit(1);
  } 

  /* Allocate temporary storage for the element moments. 
   * Array is # unknowns on element times # terms in basis (of 
   * interpolation or projection)                                      */
  moments = (double **) calloc(numUnknowns, sizeof(double *));
  for(i = 0; i < numUnknowns; i++) {
    moments[i] = (double *) calloc(nterms, sizeof(double));
  }

  /* Get grid absolute location (center of stencil). */
  getGridLocation(grid, ip, jp, kp, &x, &y, &z);

  /* Calculate element moments about given center. */
  calcElementMoments(element, elementType, systemType,
		     x, y, z, polyScale, moments, 
		     nterms, polyOrder, maxOrder1D,
		     imonoOrder, jmonoOrder, kmonoOrder);

  /* Short hands */
  rowIndex  = mat->rowIndex;
  matVals   = mat->matVals;
  /* We exploit knowing the sparse matrix data structure here.         */
  for(unknown = 0; unknown < numUnknowns; unknown++) {
    iColStart = mat->colIndex[indices[unknown]];

    for (i=0;i<nGridPoints;i++){
      /* Write the row index of this entry:                            */
      rowIndex[iColStart+i]  
	= gridDataIndex(grid, ip+is[i], jp+js[i], kp+ks[i]);
      for(temp = 0.0, j=0; j < nterms; j++) {
	temp += wts[j][i] * moments[unknown][j];
      }
      matVals[0][iColStart+i] = temp;
    }
  }

  /* Free temporary storage for the element moments. */
  for(i = 0; i < numUnknowns; i++) {
    FREE(moments[i]);
  }
  FREE(moments);
} /* End of routine calcElementGrid */


