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
/* FILE: direct.c
*
*  This file contains routines that deals with setting up 
*  the direct interaction. These routines may need 
*  information regarding the elements (for integration) and 
*  the sparse matrix structure (for storing the information), 
*  but hopefully information with regards to the grid can 
*  be avoided.
*
*  NOTE: The routines for the direct interaction, which need 
*  only information on grid or on sparse go into grid.c and 
*  sparse.c, respectively.
*
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "Global.h"
#include "Calcp.h"
#include "Direct.h"
#include "Elements.h"
#include "Sparse.h"
#include "Dense.h"

#include "GridBundles.h"

/* ======= Prototypes for static routines ======= */
static
void calcpWrapper(struct constantElement *source,
		  struct constantElement **evals, 
		  int nEvals,
		  double **infCoefs);
static 
void calcDirectMatCalloc(int nEvals,int nMats);
static
void calcDirectMatRealloc(int nEvals);
static
void calcDirectMatFree(void);

static 
void calcPrecondMatCalloc(int nEvals,int nMats);
static
void calcPrecondMatRealloc(int nEvals);
static
void calcPrecondMatFree(void);

static
void calcpWrapperCalloc(int nEvals);
static
void calcpWrapperRealloc(int nEvals);
static
void calcpWrapperFree(void);

static
void addToDiagonal(void *mat, 
		   int matNum, 
		   double addConst, 
		   char *addType);

/* ======= MACROS for this file only ======= */
/* 
* It is important for efective use of the preconditioner that
* the preconditioner use the same setup as the direct part, 
* such that a true "local part" of the larger problem is considered.
* To minimize the risc of accidently calling two different routines 
* for doing the integration, the call is implemented as a macro
* set here, and used both when setting up the direct part, and when 
* (if) calculating remaining values during the preconditioning 
* step. At the moment, calcp is used for the actual calculations, 
* and a "wrapper" has been implemented to handle to interaction 
* with calcp. All we need to do here, is to make sure that the 
* direct part and the preconditioner calls the same wrapper:           */
#define CALCPWRAPPER(ARG1,ARG2,ARG3,ARG4) \
        calcpWrapper(ARG1,ARG2,ARG3,ARG4)

/* ==== Special structures needed in this file only ==== */
/* Structure for calling calcDirectMat_1 */
struct cd1Bundle{
  void *mat;
  int nSources;
  char *type;
};
/* Structure for calling calcPrecondSparse_1 */
struct cpm1Bundle{
  void *dMat;
  void *pMat;
  char *type;
  int nDirectCols;
};
/* Structure for calling precorrectDirectMat_1 */
struct pd1Bundle{
  void *dMat;
  void *iMat;
  void *pMat;
  void *grid;
  int nSources;
  char *type;
};



/*
* ==================================================================== 
* Initialize the (sparse) direct matrix based on the current grid 
* and element distributions.  
* Figure out how much storage is needed, and subsequently allocate
* the needed memory. Then return.
* ==================================================================== */ 
void *callocDirectSparse(void *grid, int nCols, int nRows, 
			 int numDirectMats)
{
  char fctName[] = "callocDirectSparse";
  int *numNonzeroInColsDirect;        /* Number of nonzero entries in 
				       * each column of the (sparse)
				       * direct matrix                 */
  void *directMat;

  /* It the moment the numbers of sources and evals must match to build
   * a square system. Later that may not be so (one could consider e.g. 
   * building an overdetermined system...)                             */
  assert(nCols==nRows);

  /* Allocate memory for counting the number of nonzero entries 
   * in each column for both the direct matrix and the preconditioner
   * matrix                                                            */
  numNonzeroInColsDirect  = (int *) calloc(nCols,sizeof(int));

  /* Calculate the number of nonzero entries in each row of the 
   * sparse direct matrix. 
   * This call modifies numNonzeroInColsDirect, filling in the number 
   * of nonzero entries needed in each column.                         */
  printf("%s: Finding number of nonzero entries \n"
	 "       in each column of the direct part\n",fctName);
  /*sourceAndNeighboursOperate(grid,
			     "direct",
			     calcNumNonzeroDirectEntries,
			     (void*) numNonzeroInColsDirect);*/
  sourcesAndNeighboursOperate(grid,
			      "direct", 
			      0, /* No, don't use multiple sources! */
			      0, /* No limit */
			      0, /* No thanks, keep your weird data */
			      calcNumNonzeroDirectEntries,
			      (void*) numNonzeroInColsDirect);
  /* Initialize the sparse matrix for the direct part.
   * This does not actually calculate the direct influences.
   * For now just store room for one direct matrix, 
   * later the number of matrices should be problem-dependent.         */
  printf("%s: Allocating memory for the direct part\n",fctName);
  directMat  = spMatAlloc(nCols,
			  numNonzeroInColsDirect,
			  numDirectMats);
  /* Free memory allocated in the present routine                      */
  FREE(numNonzeroInColsDirect);
  /* Return a generic pointer to the direct matrix                     */
  return directMat;
} /* End of routine setupDirectSparse */

/*
* ==================================================================== 
* Calculate the (sparse) direct matrix. 
* The procedure will be the following:
* - Set up data to inform the integration function what is to be done.
* - Call the grid routine which operates on sources and neighbours.
* - The grid routine in turn calls the integration routine and 
*    supplies list of neighbours and the additional data set up here.
* - The intgration routine calculates the needed influence 
*    coefficients and stores them in the direct matrix.
*
* The first version will be to set up to integrate either the monopole 
* or the dipole coeficients, storing in the first entries of the matrix.
* Using calcp it could be efficient in respect to cpu time (if not in 
* memory) to calculate both monopole and dipole coefficients in one go, 
* storing in the first two entries in the direct matrix.
*
* As an alternative to this behaviour a linear combination of the 
* monopole and dipole coefficients will be considered as well. This 
* latter approach will be somewhat more expensive with respect to 
* CPU time, but with savings in memory. The linear combination can 
* then yield first the monopoles on the entire boundary and secondly 
* the dipoles on the entire boundary. Or, a linear combination 
* corresponding to first the right-hand-side of the mixed problem, 
* and secondly the left-hand-side could also be considered. In the 
* latter case, the linear combination will depend on what boundary 
* the evaluation element is. The easiest way to get this information 
* to the integration routine is probably to let the information be 
* "part of" the (evaluation) element itself.
*
* It should be fairly easy to swap between these different methods by 
* providing different "target function" and "target data" to the grid 
* routine. The two approaches will then be nearly de-coupled, but will 
* look extremely similar in implementation.
* 
* Note: This routine does not "need to know" about definitions of 
* grid, sparse matrices, elements etc.
*
* The parameter 'type' indicates wether the matrix should be used for 
* the system matrix (left hand side) [type="lhs"] or for constructing 
* the right-hand-side [type="rhs"].
* The value type="both" may be used to set up left-hand-side and 
* right-hand-side simultaneouosly, thus saving CPU time, but at the 
* cost of increased memory consumption. This option (type) has not 
* yet been implemented.
* ==================================================================== */
void calcDirectSparse(void *grid,
		      void *directMat,
		      int nSources,
		      char *type)
{
  char fctName[] = "calcDirectSparse";
  struct cd1Bundle bundle[1];

  /* Bundle data for calculations */
  bundle->mat      = directMat;
  bundle->nSources = nSources;
  bundle->type     = type;

  /* Calculate the direct interactions and store in the 
   * sparse direct matrix. 
   * At the moment I guess the direct part might involve:
   * 1) 1/r singular kernel parts (G) 
   * 2) d(1/r)/dn singular kernel parts (Gn) 
   * 3) Frequency dependent or other regular parts                     */
  printf("%s: Calculating direct part entries\n",fctName);
  /* This part is for setting up the sources(?) over the 
   * entire boundary.                                                  */
  sourcesAndNeighboursOperate(grid,
			      "direct", 
			      0, /* No, don't use multiple sources! */
			      0, /* No limit */
			      0, /* No thanks, keep your weird data */
			      calcDirectMat_1,
			      (void*) bundle);
  /* Free memory used in the calcp wrapper (it may be used again in a 
   * precondition step, but then it will just be allocated once more). */
  calcpWrapperFree();
  /* Free memory used in the calculation of the direct matrix          */
  calcDirectMatFree();
  /* For a Neumann-part add 2*pi to the diagonal of the matrix
   * NOTE: This should depend both on the context of the direct matrix 
   * (for solve r setting up right-hand-side) AND on the boundary 
   * condition applied at each panel. This should be taken care of in
   * calcp(?). No need to do it here also...                           */
  return;
} /* End of routine calcDirectSparse */

/*
* ==================================================================== 
* Precorrect the (sparse) direct matrix. 
* The procedure will be the following:
*  - For each source and eval configuration use precomputed 
*    grid-to-grid (G2G) matrices to compute I*(G2G*P) as 
*    the precorrection value. 
*  - (G2G*P) can be recycled for all the evals that are associated 
*    with a single grid point. 
*  - For neighbouring evals, which cannot be reached "through the grid" 
*    use compression of the interpolation points and compute the 
*    (larger) G2G matrix for all these evals together. Then sort 
*
* Type should either be "many small" or "one large".
* ==================================================================== */
void calcPrecorrectSparse(void *directMat,
			  void *projectMat,
			  void *interpMat,
			  void *grid,
			  int nSources,
			  char *type)
{
  char fctName[] = "calcPrecorrectSparse";
  struct pd1Bundle bundle[1];

  /* Bundle data for calculations */
  bundle->dMat     = directMat;
  bundle->iMat     = interpMat;
  bundle->pMat     = projectMat;
  bundle->grid     = grid;
  bundle->nSources = nSources;
  bundle->type     = type;

  /* Calculate the direct interactions and store in the 
   * sparse direct matrix. 
   * At the moment I guess the direct part might involve:
   * 1) 1/r singular kernel parts (G) 
   * 2) d(1/r)/dn singular kernel parts (Gn) 
   * 3) Frequency dependent or other regular parts                     */
  printf("%s: Calculating precorrection entries\n",fctName);
  /* This part is for setting up the sources(?) over the 
   * entire boundary.                                                  */
  if (strcmp(type,"many small")==0)
    sourcesAndNeighboursOperate(grid,
				"direct", 
				0, /* No, don't use multiple sources! */
				0, /* No limit */
				1, /* Yes please, I need data set 1   */
				precorrectDirectMat_1,
				(void*) bundle);
  else if (strcmp(type,"one large")==0)
    sourcesAndNeighboursOperate(grid,
				"direct", 
				0, /* No, don't use multiple sources! */
				0, /* No limit */
				1, /* Yes please, I need data set 1   */
				precorrectDirectMat_2, 
				(void*) bundle);
  else {
    /* This is unknown = bad */
    fprintf(stderr,"%s: ERROR: Unknown type specifier for precorrection.\n"
	    "     Please use \"many small\" or \"one large\"\n",
	    fctName);
    _EXIT_;
  }
  
  /* Free memory used in the calcp wrapper (it may be used again in a 
   * precondition step, but then it will just be allocated once more). */
  calcpWrapperFree();
  /* Free memory used in the calculation of the direct matrix          */
  calcDirectMatFree();
  return;
} /* End of routine calcPrecorrectSparse */

/*
* ==================================================================== 
* Initialize the (sparse) preconditioner for the direct matrix based 
* on the current grid and element distributions.  
* Figure out how much storage is needed, and subsequently allocate
* the needed memory. Then return.
* ==================================================================== */
void *callocPrecondSparse(void *grid, void *directMat, 
			 int maxNeighbours,
			 int nCols, int nRows, 
			 int numPrecondMats)
{
  char fctName[] = "callocPrecondSparse";
  int *numNonzeroInColsPrecond;       /* Number of nonzero entries in 
				       * each column of the (sparse)
				       * direct matrix                 */
  void *precondMat;

  /* It the moment the numbers of sources and evals must match to build
   * a square system. Later that may not be so (one could consider e.g. 
   * building an overdetermined system...)                             */
  assert(nCols==nRows);

  /* Allocate memory for counting the number of nonzero entries 
   * in each column for both the direct matrix and the preconditioner
   * matrix                                                            */
  numNonzeroInColsPrecond  = (int *) calloc(nCols,sizeof(int));

  /* Calculate the number of nonzero entries in each row of the 
   * sparse direct matrix. 
   * This call modifies numNonzeroInColsDirect, filling in the number 
   * of nonzero entries needed in each column.                         */
  printf("%s: Finding number of nonzero entries \n"
	 "       in each column of the preconditioning part\n",fctName);
  /*sourcesAndNeighboursOperateOld(grid, 
			      "precond",
			      maxNeighbours,
			      calcNumNonzeroPrecondEntries,
			      (void*) numNonzeroInColsPrecond);*/
  sourcesAndNeighboursOperate(grid, 
			      "precond",
			      1,   /* Yes, use multiple sources */
			      maxNeighbours,
			      0,   /* No, don't add your weird data   */
			      calcNumNonzeroPrecondEntries, 
			      (void*) numNonzeroInColsPrecond);
  /* Initialize the sparse matrix for the preconditioning part.
   * For now just store room for one direct matrix, 
   * later the number of matrices should be problem-dependent.         */
  printf("%s: Allocating memory for the preconditioning part\n",fctName);
  precondMat  = spMatAlloc(nCols,
			   numNonzeroInColsPrecond,
			   numPrecondMats);

  /* Free memory allocated in the present routine                      */
  FREE(numNonzeroInColsPrecond);

  /* Return a generic pointer to the preconditioning matrix            */
  return precondMat;
} /* End routine callocPrecondSparse */

/*
* ==================================================================== 
* Setup the preconditioning matrix.
*
* Note: This routine does not "need to know" about definitions of 
* grid, sparse matrices, elements etc. 
* ==================================================================== */
void calcPrecondSparse(void *grid, void *directMat, void *precondMat,
			int maxNeighbours)
{
  /*char fctName[] = "calcPrecondSparse";*/
  struct cpm1Bundle bundle[1];

  /* Bundle the direct matrix and the preconditioning matrix, 
   * so that the entries in directMat can be recycled when 
   * constructing precondMatfirst.                                     */
  bundle->dMat = directMat;
  bundle->pMat = precondMat;
  /*sourcesAndNeighboursOperateOld(grid, 
			      "precond",
			      maxNeighbours,
			      calcPrecondMat_1,
			      (void*) bundle);*/
  sourcesAndNeighboursOperate(grid, 
			      "precond",
			      1, /* Yes, use multiple sources */
			      maxNeighbours,
			      0, /* No, don't add your weird data   */
			      calcPrecondMat_1, /* target function  */
			      (void*) bundle);  /* target data      */
  /* If the wrapper to calcp has been used then free the allocated 
   * memory (it may be used again later?, but then it will just be 
   * allocated once more).                                             */
  calcpWrapperFree();
  /* Free memory used in the precondition calculation:                 */
  calcPrecondMatFree();

#ifdef NORECYCLE
  /* If everything is recomputed for the preconditioner, then we need
   * to modify the diagonal. Otherwise the diagonal should not need 
   * modification!                                                     */
  /* For a Neumann-part add 2*pi to the diagonal of the matrix
   * NOTE: This should depend both on the context of the direct matrix 
   * (for solve r setting up right-hand-side) AND on the boundary 
   * condition applied at each panel.
   * calcp should take care of this. No need to do this stuff I guess! 
   addToDiagonal(directMat, 0, -2.0*PI, "all");
*/
#endif

  return;
} /* End routine calcPrecondSparse */



/* 
* ====================================================================
* The following routines and static variables go together as parts of 
* building the direct matrix.
* At the present calcDirectMat_1, the direct calculation routine itself, 
* checks if the storage has been allocated - and takes care of the 
* allocation. If the storage needs to be released for other use, 
* then calcDirectMatFree needs to be called. calcDirectMat wont 
* do this, since it doesn't know if it will be called again.
* Use:
*  1)  Allocate needed storage (once):  calcDirectMatCalloc
*  2)  Calls to calcp (multiple times): calcDirectMat_1
*  2b) Reallocate storage if needed:    calcDirectMatRealloc
*  3)  Free the used memory (once):     calcDirectMatFree
*
* Different preconditioners may be implemented - all using the same 
* memory for working space.
* Pointers are initially explicitly set to NULL, so that anything else 
* will indicate that memory has been allocated.
* ==================================================================== */
static int 
  isAllocated_calcDirectMat=0, /* Flag for allocation status.
				* 0: Not allocated
				* 1: Allocated                         */
  nEvalsMax_calcDirectMat=50,  /* Size of arrays being allocated, set 
				* to a minimum default value           */
  nMats_calcDirectMat,         /* Used for keeping track of memory 
				* allocated for the influence 
				* coefficients - needed in the 
				* deallocation step.                   */
  *isortEvals_calcDirectMat=NULL,/* Arrays for sorting after global    */
  *iorigEvals_calcDirectMat=NULL;/* index the passed list of 
				  * neighbouring evaluation elements.  */


static struct constantElement 
  **sortEvals_calcDirectMat=NULL;/* Pointers to neighbouring evaluation 
				  * elements, sorted after global index  */
static double 
  **infCoefs_calcDirectMat=NULL; /* Storage for those influence 
				  * coefficients that are calculated 
				  * in the preconditioning step.         */
double *g2gp_calcDirectMat=NULL, /* These arrays are needed in a newer   */
  *precVals_calcDirectMat=NULL,  /* version of the precorrection.        */
  *precVals_sort_calcDirectMat=NULL;/* See routine "precorrectDirectMat_1"*/
int *dRows_calcDirectMat=NULL;

/*
* --------------------------------------------------------------------
* Given a source and a list of neighbours, store the influence 
* coefficients in the direct matrix (also supplied).
*
* This version deals only with triangular and quadrilateral 
* "constant panels", as used in e.g. wamit. Thus "calcp" can be
* used for the actual integration.
*
* Note that this function must conform to the layout of a "target 
* function" of "sourceAndNeighboursOperate" (see e.g. grid.c).
* This function is a companion to "calcPrecondMat_1". These two routines
* should be modified together for consistency.
* -------------------------------------------------------------------- */
void calcDirectMat_1(void **sourcesIn, 
		     void **evalsIn, 
		     int nSources,
		     int nTotNeighbours,
		     void *bundleIn)
{
  char fctName[] = "calcDirectMat_1";
  int verbose =1;
  struct constantElement *source, **evals;
  struct cd1Bundle *bundle;
  char *type;
  void *dmat;

  int i, iCol, matNum;
  double *values;
  char mark = '>';
  static int nSourcesDone;

  /* Only ONE source at a time here! */
  assert(nSources==1);

  /* Get data from bundle */
  bundle = (struct cd1Bundle *) bundleIn;
  dmat   = bundle->mat;
  type   = bundle->type;

  /* Check that type can be recognized                                 */
  if (  strcmp(type,"rhs") &&
        strcmp(type,"lhs")  ){
    printf("%s: Unrecognized value of parameter type: \"%s\"\n",
	   fctName,type);
    printf("     Please use \"rhs\", or \"lhs\".\n");
    _EXIT_;
  }

  /* Convert generic pointers given in input to "true" types:          */
  source = (struct constantElement*)  sourcesIn[0];
  evals  = (struct constantElement**) evalsIn; /* Is this dangerous(?) */

  /*{ 
    int *indices;
    indices=sourceIndices(source);
    printf("%s called, source %d, %d evals\n",
	   fctName,indices[0],nTotNeighbours);
	   }*/
  /* Allocate memory needed for this routine.                          */
  if(!isAllocated_calcDirectMat) {
    calcDirectMatCalloc(nTotNeighbours,2); /* Allocate memory for 
					    * storing both G and Gn    */
    /* Initiate plotting progress marks */
    if (verbose) {
      /* Nothing done yet. */
      plotProgressMark(nSourcesDone=0, 1, 0, bundle->nSources,1,mark);
    }
  }
  /* Check that enough memory has been allocated                       */
  else if(nTotNeighbours>nEvalsMax_calcDirectMat)
    calcDirectMatRealloc(nTotNeighbours);

  /* The evaluation elements are delivered unsorted (except for a small 
   * part), so it may be an advantage to sort the entries. This would 
   * mean that the column entries are not sorted, i.e. the sparse matrix 
   * may "say" something like:
   * Column 5 have nonzeros in rows 5, 8, 2 and 1, 
   * rather than:
   * Column 5 have nonzeros in rows 1, 2, 5 and 8.
   * For performing a matrix-vector product this does not matter (much),
   * but if someone wanted to see if two elements (eval=i,source=j) are
   * neighbours by checking if entry [i,j] is nonzero, then we would 
   * have to check *all* nonzero entries of column j to see if one of 
   * them is i, rather than performing a binary search. Also, when 
   * setting up the *preconditioner* it may be advantageous to
   * lookup for entries in the existing matrix before calculation.
   * This could mean that a lot of searches within columns may be 
   * needed, and, consequently, that sorting could be an advantage.
   * 
   * For now, it has been decided to sort the entries                  */
  /* The present calling routine recycles at least part of the 
   * list of neighbours. Thus, a second list is allocated here to 
   * avoid corruption of the initial list (sortEvals).                 */

  /* Sort the entries                                                  */
  for (i=0;i<nTotNeighbours;i++){
    isortEvals_calcDirectMat[i] = evals[i]->directIndex;
    iorigEvals_calcDirectMat[i] = i;
  }
  intSort2(isortEvals_calcDirectMat,
	   iorigEvals_calcDirectMat,
	   nTotNeighbours);
  for (i=0;i<nTotNeighbours;i++){
    sortEvals_calcDirectMat[i]  = evals[iorigEvals_calcDirectMat[i]];
  }


  /* TEST block: 
  {
    int multEvals;
    multEvals=0;
    for (i=1;i<nTotNeighbours;i++){
      if (sortEvals_calcDirectMat[i]==sortEvals_calcDirectMat[i-1]){
	printf("%s:  %d and %d are actually the same eval!\n",fctName,i-1,i);
	multEvals=1;
	break;
      }
    }
    if (!multEvals)
      printf("%s: No multiple evals found\n",fctName);
  }
  */

  /* Calculate influence coefficients. Use a (macro for) a wrapper 
   * for calcp                                                         */
  CALCPWRAPPER(source,sortEvals_calcDirectMat,nTotNeighbours,
	       infCoefs_calcDirectMat);

  /* Store information in direct matrix                                */
  iCol     = source->directIndex;

  /* Test only to see if the sorting is OK at this point.
   * This can safely be removed when the test is done a few times      
  for (i=0;i<nTotNeighbours;i++) 
    assert(sortEvals_calcDirectMat[i]->directIndex
	   ==isortEvals_calcDirectMat[i]);*/

  /* Figure out which values that should be stored                     */
  if ((int)source->bcType==0){
    /* This is a Neumann condition boundary, i.e. PHI_n is known.
     * Thus G goes on the right-hand-side (with PHI_n) and Gn
     * goes on the system (left-hand) side (with PHI).                 */
    if (strcmp(type,"rhs")==0) 
      values = infCoefs_calcDirectMat[0]; /* This is G */
    else if (strcmp(type,"lhs")==0) 
      values = infCoefs_calcDirectMat[1]; /* This is Gn */
  }
  else if ((int)source->bcType==1){
    /* This is a Diriclet condition boundary, i.e. PHI is known.
     * Thus Gn goes on the right-hand-side (with PHI) and G
     * goes on the system (left-hand) side (with PHI_n).               */
    if (strcmp(type,"rhs")==0) 
      values = infCoefs_calcDirectMat[1]; /* This is Gn */
    else if (strcmp(type,"lhs")==0) 
      values = infCoefs_calcDirectMat[0]; /* This is G */
  }
  else {
    printf("%s: What kind of boundary condition type is this: %d\n",
	   fctName,(int)source->bcType);
    _EXIT_;
  }

  /* Store entries as a column in the direct matrix                    */
  spMatColumnStore(dmat, 
		   matNum=0,                  /* Matrix values number  */
		   iCol,                      /* Column number         */
		   nTotNeighbours,            /* # rows                */
		   isortEvals_calcDirectMat,  /* row indices           */
		   values,                    /* row values            */
		   "overwrite",               /* New row structure     */
		   "discard");                /* New row values        */

  if (verbose) { /* Plot progress */
    if (++nSourcesDone == bundle->nSources) /* All done */ 
      plotProgressMark(nSourcesDone, -1, 0, bundle->nSources,0,mark);
    else  
      plotProgressMark(nSourcesDone,  0, 0, bundle->nSources,0,mark);
  }

  return;
} /* End of routine calcDirectMat_1                                    */

/*
* --------------------------------------------------------------------
* Given a source and a list of neighbours as well as some other 
* information, calculate the precorrection values for the direct
* matrix. 
*
* Note that information on the actual element type is used here!!
*
* It is assumed that the grid-routine hand over a (pointer to) an
* array of precomputed small G2G matrices and a vector with an index
* for each evaluation element telling which matrix to use in that 
* particular case.
*
* NOTE: It is assumed in this routine that all elements are projected 
* onto the grid using ONE projection scheme ("on level 1") and 
* interpolated using ONE interpolation scheme (also "on level one").
* If projection and interpolation of "large" elements is implemented, 
* then the fact that some elements are "special" MUST BE PASSED ON
* TO THIS ROUTINE!! For now this flexibility is not implemented 
* (but it would certainly not be impossibly to implement it!). 
*
* Note that this function must conform to the layout of a "target 
* function" of "sourceAndNeighboursOperate" (see e.g. grid.c).
* This function should be VERY similar to "calcDirectMat_1" and 
* "calcPrecondMat_1". These three routines should be modified 
* together for consistency.
*
* Procedure: 
* Get the column number for the source.
* Get the projection column for the source.
*  Get the next (or first) g2g index
*  Calculate p=G2G*P
*  For each eval with this direct index (they should 
*   be sorted correctly)
*   Get interpolation column and store directr entry index.
*   Calculate c = p*Interp and store this value as precorrection 
*    value. 
*  Next direct index
* For index <0 calculate precorrection the expensive way, by
*  calling a grid routine to do the dirty work (projection and 
*  interpolation columns can be supplied). Store direct entry 
*  index and correction value.
* Sort the correction values by row entry number as done in the 
*  direct-part computations.
* Call a sparse-matrix routine to add (subtract) the values to
*  the existing entries. Make a check on the indexing.
* -------------------------------------------------------------------- */
/* These inclusions are for test purposes only. */
#include "Grid.h"
#include "Sparse.h"
void precorrectDirectMat_1(void **sourcesIn, 
			   void **evalsIn, 
			   int nSources,
			   int nTotNeighbours,
			   void *bundleIn)
{
  char fctName[] = "precorrectDirectMat_1";
  int verbose = 1;
  char mark = '>'; /* Progression mark to use (not influencing computations) */
  struct constantElement *source, **evals;
  struct pd1Bundle *bundle;
  struct gridBundle_SANO_1 *gridBundle;
  char *type;
  void *dmat, *pmat, *imat, *grid;
  double ***g2g; /* small G2G matrices */
  int *g2gIndex; /* Index to g2g. One for each eval */
  int  ng2g, nInterp, nProj, ievl, iSource, iEval,
    startEvl, doPartialColumns;
  int *projectRowIndices, *interpRowIndices;
  double *projectValues, *interpValues, **thisG2G, 
    *thisG2Grow;
  int thisG2Gindex, lastG2Gindex;
  static int nSourcesDone;
  register int mg2g, i, j;
  register double addIt, addIt2;
  double  *precValsPointer[1];

  /* Only ONE source at a time here! */
  assert(nSources==1);

  /* Get data from grid bundle */
  gridBundle = (struct gridBundle_SANO_1 *) bundleIn;
  g2g        = gridBundle->nearFieldG2G;
  g2gIndex   = gridBundle->g2gIndex;
  ng2g       = gridBundle->G2Gn;
  mg2g       = gridBundle->G2Gm;

  /* Get data from bundle */
  bundle = (struct pd1Bundle *) gridBundle->origData;
  dmat   = bundle->dMat;
  imat   = bundle->iMat;
  pmat   = bundle->pMat;
  grid   = bundle->grid;
  type   = bundle->type;

  /* Convert generic pointers given in input to "true" types:          */
  source = (struct constantElement*)  sourcesIn[0];
  evals  = (struct constantElement**) evalsIn; /* Is this dangerous(?) */

  /* Allocate memory needed for this routine.                          */
  if(!isAllocated_calcDirectMat) {
    printf(">>%s: REMOVE HARDCODED ARRAY SIZES!!\n", fctName);
    calcDirectMatCalloc(nTotNeighbours,2); /* Allocate memory          */
    /* Initiate plotting progress marks */
    if (verbose) {
      /* Nothing done yet. */
      plotProgressMark(nSourcesDone=0, 1, 0, bundle->nSources,1,mark);
    }
  }
  /* Check that enough memory has been allocated                       */
  else if(nTotNeighbours>nEvalsMax_calcDirectMat)
    calcDirectMatRealloc(nTotNeighbours);

  /* Column number in the direct matrix (and projetion matrix)
   *  for this source element                                          */
  iSource = source->directIndex;
  /* Get this column from the projection matrix                        */
  spMatColumnGet(pmat,  /* Projection matrix */
		 0,     /* First matrix (sets of values) */
		 iSource, /* Column # */
		 &nProj, /* Number of rows (#proj. points) */
		 &projectRowIndices, /* Row numbers */
		 &projectValues,  /* Entry values */
		 NULL); /* Last argument not used */
  /* Make sure that the right number of grid points are used 
   * in the projection                                                */
  assert(nProj==mg2g);

  /* Don't do partial stuff unless we have to (see later) */
  doPartialColumns=0; 

  for (ievl=0, lastG2Gindex=-1; ievl<nTotNeighbours; ievl++){
    iEval = evals[ievl]->directIndex;
    dRows_calcDirectMat[ievl] = iEval; /* Store row number */
    thisG2Gindex = g2gIndex[ievl];
    if (thisG2Gindex!= lastG2Gindex){
      /* Make new matrix-vector product G2G*P                          */      
      if (thisG2Gindex<0){
	/* The rest evals are not found through the grid               */
	startEvl = ievl;
	/* So in this case there ARE some evals that we need to handle 
	 * separately. Set a flag accordingly */
	doPartialColumns=1; 
	break;
      }
      assert(thisG2Gindex>=0 && thisG2Gindex<gridBundle->numG2G);
      /* Shorthand to this G2G matrix:                                 */
      thisG2G = g2g[thisG2Gindex];
      for (i=0; i<ng2g; i++){
	thisG2Grow = thisG2G[i];
	for (j=0, addIt=0.0; j<mg2g; j++){
	  addIt += projectValues[j]*thisG2Grow[j];
	}
	g2gp_calcDirectMat[i] = addIt; /* Store in data array */
      }
      lastG2Gindex = thisG2Gindex; /* New last index */
    }
    /* Get column from interpolation matrix */
    spMatColumnGet(imat,  /* Interpolation matrix */
		   0,     /* First matrix (sets of values) */
		   iEval, /* Column # */
		   &nInterp, /* Number of rows (#proj. points) */
		   &interpRowIndices, /* Row numbers */
		   &interpValues,  /* Entry values */
		   NULL); /* Last argument not used */
    /* Calculate I*G2G*P as vector dot product: I*(G2G*P)              */
    for (i=0, addIt2 = 0.0; i<ng2g; i++){
      addIt2 += interpValues[i]*g2gp_calcDirectMat[i];
    }
    precVals_calcDirectMat[ievl] = addIt2;
  }
  if(doPartialColumns){
    /* We still need to precorrect the remaining entries 
     * corresponding to evals that cannot be reached "through the grid",
     * (but then must be in a source-local list).                      */
    /* Set direct indices of the remaining evals                       */
    for (/*just continue counting*/;ievl<nTotNeighbours; ievl++){
      dRows_calcDirectMat[ievl] = evals[ievl]->directIndex;
    }

    /* Get precorrection values for the remaining evals                  */
    /*precPointer = &(precVals[startEvl]);*/
    precValsPointer[0] = &(precVals_calcDirectMat[startEvl]);
    precorrectPartialColumns(1,          /* One source only */
			     &(iSource), /* Source index as "array" */
			     nTotNeighbours-startEvl, /* #remaining evals */
			     &(dRows_calcDirectMat[startEvl]),   /* indices of remaining evals */
			     precValsPointer,/* Storage for return values */
			     dmat, pmat, imat, grid);
  }
  /* Some little thing seems to be wrong. 
   * All the computed values have checked out OK, but even so there is
   * a discrepancy between earlier obtained results! 
   * I'll try to dump the direct matrix for comparison? */

  /* Sort the correction values by row entry number as done in the 
   *  direct-part computations.                                        */
  for (i=0;i<nTotNeighbours;i++){
    isortEvals_calcDirectMat[i] = evals[i]->directIndex;
    iorigEvals_calcDirectMat[i] = i;
  }
  intSort2(isortEvals_calcDirectMat,
	   iorigEvals_calcDirectMat,
	   nTotNeighbours);
  for (i=0;i<nTotNeighbours;i++){
    precVals_sort_calcDirectMat[i]  
      = precVals_calcDirectMat[iorigEvals_calcDirectMat[i]];
  }

  /* Double check entries */
  for (ievl=0; ievl<nTotNeighbours; ievl++){
    iEval = isortEvals_calcDirectMat[ievl];
    /*assert(iEval == evals[ievl]->directIndex);*/
    /* Put test here */
    if (0){/* Test vs. existing routine */
      struct grid *g;
      struct sparseMat *PMAT,*IMAT;
      double cval, cval2, ee;
      int pStart,pStop,iStart,iStop;
      static double maxee=0.0;

      cval=precVals_sort_calcDirectMat[ievl];
      g = (struct grid*) grid;
      PMAT = (struct sparseMat *) pmat;
      IMAT = (struct sparseMat *) imat;

      pStart = PMAT->colIndex[iSource];
      pStop  = PMAT->colIndex[iSource+1];
      iStart = IMAT->colIndex[iEval];
      iStop  = IMAT->colIndex[iEval+1];
      
      cval2 = precorrectVal(g->kernel,g->griddata,
			    &(PMAT->matVals[0][pStart]), 
			    &(PMAT->rowIndex[pStart]), 
			    pStop-pStart,
			    &(IMAT->matVals[0][iStart]), 
			    &(IMAT->rowIndex[iStart]), 
			    iStop-iStart);

      ee = ABS(cval-cval2)/ABS(cval2);

      if (ee>1e-6){
	printf("DISCREPANCY ");
	printf("(%d,%d,%d-%d): %.12e %.12e %.6e\n",
	       iSource,iEval,iStart,iStop-1,cval,cval2,ee);
	assert(ee<1e-14);
      }
      else if (0){
	printf("OK> ");
	printf("(%d,%d,%d-%d): %.16e %.16e %.6e\n",
	       iSource,iEval,iStart,iStop-1,cval,cval2,ee);
      }
      if (ee>maxee){
	maxee = ee;
	printf(" Maximum error encountered is %.8g\n",maxee);
      }
    }
  }

  /* Subtract precorrection values from the correct column of the 
   * direct matrix                                                     */
  spMatColumnStore(dmat, 0, 
		   iSource, 
		   nTotNeighbours,      /* Total number of rows */
		   isortEvals_calcDirectMat, /* Row numbers (sorted) */
		   precVals_sort_calcDirectMat, /* Entry values */
		   "check",             /* Verify that row numbering is consistent! */
		   "subtractFrom");     /* Add to existing values */


  if (verbose) { /* Plot progress */
    if (++nSourcesDone == bundle->nSources) /* All done */ 
      plotProgressMark(nSourcesDone, -1, 0, bundle->nSources,0,mark);
    else  
      plotProgressMark(nSourcesDone,  0, 0, bundle->nSources,0,mark);
  }
  return;
} /* End of routine precorrectDirectMat_1 */

/* Precorrect using one long skinny grid-to-grid matrix */
void precorrectDirectMat_2(void **sourcesIn, 
			   void **evalsIn, 
			   int nSources,
			   int nTotNeighbours,
			   void *bundleIn)
{
  char fctName[] = "precorrectDirectMat_2";
  int verbose = 1;
  char mark = '>'; /* Progression mark to use (not influencing computations) */
  struct constantElement *source, **evals;
  struct pd1Bundle *bundle;
  struct gridBundle_SANO_1 *gridBundle;
  char *type;
  void *dmat, *pmat, *imat, *grid;
  double **g2gunion; /* Large skinny G2G near-field matrix */
  double *gpWork;    /* Work storage for G2G*PROJ */
  int *g2gIndex; /* Index to g2gIndex. One for each eval */
  int  ng2gCols, ng2gRows, ng2g, nInterp, nProj, ievl, iSource, iEval,
    startEvl, doPartialColumns;
  int *projectRowIndices, *interpRowIndices, *compressedIndices, **unionG2Gindices;
  double *projectValues, *interpValues, **thisG2G, 
    *thisG2Grow;
  int thisG2Gindex, lastG2Gindex;
  static int nSourcesDone;
  register int i, j, nProjReg;
  register double addIt, addIt2;
  double *precValsPointer[1];

  /* Only ONE source at a time here! */
  assert(nSources==1);

  /* Get data from grid bundle */
  gridBundle = (struct gridBundle_SANO_1 *) bundleIn;
  g2gunion   = gridBundle->nearFieldUnionG2G;
  g2gIndex   = gridBundle->g2gIndex;
  unionG2Gindices = gridBundle->nearFieldUnionG2Gindices;
  ng2gRows   = gridBundle->nRowsNearFieldUnion;
  ng2gCols   = gridBundle->G2Gm;
  gpWork     = gridBundle->unionWork;

  /* Get data from bundle */
  bundle = (struct pd1Bundle *) gridBundle->origData;
  dmat   = bundle->dMat;
  imat   = bundle->iMat;
  pmat   = bundle->pMat;
  grid   = bundle->grid;
  type   = bundle->type;

  /* Convert generic pointers given in input to "true" types:          */
  source = (struct constantElement*)  sourcesIn[0];
  evals  = (struct constantElement**) evalsIn; /* Is this dangerous(?) */

  /* Allocate memory needed for this routine.                          */
  if(!isAllocated_calcDirectMat) {
    calcDirectMatCalloc(nTotNeighbours,2); /* Allocate memory for 
					    * storing both G and Gn    */
    /* Initiate plotting progress marks */
    if (verbose) {
      /* Nothing done yet. */
      plotProgressMark(nSourcesDone=0, 1, 0, bundle->nSources,1,mark);
    }
  }
  /* Check that enough memory has been allocated                       */
  else if(nTotNeighbours>nEvalsMax_calcDirectMat)
    calcDirectMatRealloc(nTotNeighbours);

  /* TEST: Just do the matrix-vector products to see if 
   * this method is worth pursuing! */
  iSource = source->directIndex;
  /* Get this column from the projection matrix */
  spMatColumnGet(pmat,  /* Projection matrix */
		 0,     /* First matrix (sets of values) */
		 iSource, /* Column # */
		 &nProj, /* Number of rows (#proj. points) */
		 &projectRowIndices, /* Row numbers */
		 &projectValues,  /* Entry values */
		 NULL); /* Last argument not used */
  assert(nProj==ng2gCols);

  /* Compute G2G*PROJ for a large skinny near-field G2G matrix. 
   * This may result in too much work done, but at least no 
   * interpolation point will be accessed more than once.              */
  nProjReg = nProj;
  for (i=0; i<ng2gRows; i++){
    thisG2Grow = g2gunion[i];
    for (j=0, addIt = 0.0; j<nProjReg; j++){
      addIt += projectValues[j]*thisG2Grow[j];
    } 
    /* Store in data array */
    gpWork[i] = addIt;
  }

  /* Don't do partial stuff unless we have to (see later) */
  doPartialColumns=0; 

  /* For each eval, calculate a linear combination of the pre-computed
   * G2G*PROJ to obtain the precorrection value INTERP*G2G*PROJ        */
  for (ievl=0, lastG2Gindex=-1; ievl<nTotNeighbours; ievl++){
    iEval = evals[ievl]->directIndex;
    dRows_calcDirectMat[ievl] = iEval; /* Store row number */
    thisG2Gindex = g2gIndex[ievl];
    if (thisG2Gindex!= lastG2Gindex){
      /* Make vector inner (dot) product with selected values of the 
       * G2G*P vector computed above                                   */      
      if (thisG2Gindex<0){
	/* The rest evals are not found through the grid               */
	startEvl = ievl;
	/* So in this case there ARE some evals that we need to handle 
	 * separately. Set a flag accordingly */
	doPartialColumns=1; 
	break;
      }
      assert(thisG2Gindex>=0 && thisG2Gindex<gridBundle->numG2G);
      /* Shorthand to this G2G indices:                                */
      compressedIndices = unionG2Gindices[thisG2Gindex];
      lastG2Gindex = thisG2Gindex; /* New last index */
    }
    /* Get column from interpolation matrix */
    spMatColumnGet(imat,  /* Interpolation matrix */
		   0,     /* First matrix (sets of values) */
		   iEval, /* Column # */
		   &nInterp, /* Number of rows (#proj. points) */
		   &interpRowIndices, /* Row numbers */
		   &interpValues,  /* Entry values */
		   NULL); /* Last argument not used */
    /* Calculate I*G2G*P as vector dot product: I*(G2G*P)              */
    for (i=0, addIt2 = 0.0; i<nInterp; i++){
      addIt2 += interpValues[i]*gpWork[compressedIndices[i]];
    }
    precVals_calcDirectMat[ievl] = addIt2;
  }
  if(doPartialColumns){
    /* We still need to precorrect the remaining entries 
     * corresponding to evals that cannot be reached "through the grid",
     * (but then must be in a source-local list).                      */
    /* Set direct indices of the remaining evals                       */
    for (/*just continue counting*/;ievl<nTotNeighbours; ievl++){
      dRows_calcDirectMat[ievl] = evals[ievl]->directIndex;
    }

    /* Get precorrection values for the remaining evals                  */
    /*precPointer = &(precVals[startEvl]);*/
    precValsPointer[0] = &(precVals_calcDirectMat[startEvl]);
    precorrectPartialColumns(1,          /* One source only */
			     &(iSource), /* Source index as "array" */
			     nTotNeighbours-startEvl, /* #remaining evals */
			     &(dRows_calcDirectMat[startEvl]),   /* indices of remaining evals */
			     precValsPointer,/* Storage for return values */
			     dmat, pmat, imat, grid);
  }
  /* Some little thing seems to be wrong. 
   * All the computed values have checked out OK, but even so there is
   * a discrepancy between earlier obtained results! 
   * I'll try to dump the direct matrix for comparison? */

  /* Sort the correction values by row entry number as done in the 
   *  direct-part computations.                                        */
  for (i=0;i<nTotNeighbours;i++){
    isortEvals_calcDirectMat[i] = evals[i]->directIndex;
    iorigEvals_calcDirectMat[i] = i;
  }
  intSort2(isortEvals_calcDirectMat,
	   iorigEvals_calcDirectMat,
	   nTotNeighbours);
  for (i=0;i<nTotNeighbours;i++){
    precVals_sort_calcDirectMat[i]  = precVals_calcDirectMat[iorigEvals_calcDirectMat[i]];
  }

  /* Double check entries */
  for (ievl=0; ievl<nTotNeighbours; ievl++){
    iEval = isortEvals_calcDirectMat[ievl];
    /*assert(iEval == evals[ievl]->directIndex);*/
    /* Put test here */
    if (0){/* Test vs. existing routine */
      struct grid *g;
      struct sparseMat *PMAT,*IMAT;
      double cval, cval2, ee;
      int pStart,pStop,iStart,iStop;
      static double maxee=0.0;

      cval=precVals_sort_calcDirectMat[ievl];
      g = (struct grid*) grid;
      PMAT = (struct sparseMat *) pmat;
      IMAT = (struct sparseMat *) imat;

      pStart = PMAT->colIndex[iSource];
      pStop  = PMAT->colIndex[iSource+1];
      iStart = IMAT->colIndex[iEval];
      iStop  = IMAT->colIndex[iEval+1];
      
      cval2 = precorrectVal(g->kernel,g->griddata,
			    &(PMAT->matVals[0][pStart]), 
			    &(PMAT->rowIndex[pStart]), 
			    pStop-pStart,
			    &(IMAT->matVals[0][iStart]), 
			    &(IMAT->rowIndex[iStart]), 
			    iStop-iStart);

      ee = ABS(cval-cval2)/ABS(cval2);

      if (ee>1e-6){
	printf("DISCREPANCY ");
	printf("(%d,%d,%d-%d): %.12e %.12e %.6e\n",
	       iSource,iEval,iStart,iStop-1,cval,cval2,ee);
	assert(ee<1e-14);
      }
      else if (0){
	printf("OK> ");
	printf("(%d,%d,%d-%d): %.16e %.16e %.6e\n",
	       iSource,iEval,iStart,iStop-1,cval,cval2,ee);
      }
      if (ee>maxee){
	maxee = ee;
	printf(" Maximum error encountered is %.8g\n",maxee);
      }
    }
  }

  /* Subtract precorrection values from the correct column of the 
   * direct matrix                                                     */
  spMatColumnStore(dmat, 0, 
		   iSource, 
		   nTotNeighbours,      /* Total number of rows */
		   isortEvals_calcDirectMat, /* Row numbers (sorted) */
		   precVals_sort_calcDirectMat, /* Entry values */
		   "check",             /* Verify that row numbering is consistent! */
		   "subtractFrom");     /* Add to existing values */


  if (verbose) { /* Plot progress */
    if (++nSourcesDone == bundle->nSources) /* All done */ 
      plotProgressMark(nSourcesDone, -1, 0, bundle->nSources,0,mark);
    else  
      plotProgressMark(nSourcesDone,  0, 0, bundle->nSources,0,mark);
  }
  return;
} /* End of routine */

/*
* --------------------------------------------------------------------
* Allocate needed memory for calcDirectMat .
* A maximum is set on the number of expected evals. If the actual 
* number of evals exceeds this number, the memory will be reallocated
* later.
* -------------------------------------------------------------------- */
static 
void calcDirectMatCalloc(int nEvals, int nMats)
{
  char fctName[]="calcDirectMatCalloc";
  int i;
  /*printf("%s called, %d evals\n",fctName,nEvals);*/
  /* Make absolutely sure that this memory has not been allocated      */
  if (isAllocated_calcDirectMat) {
    fprintf(stderr,"%s: Apparently memory is already allocated.\n"
	    "    (this may indicate a bug. Better fix it!)\n"
	    "    In the mean time I'll stop the program\n",
	    fctName);
    _EXIT_;
  }
  if (isortEvals_calcDirectMat!=NULL ||
      iorigEvals_calcDirectMat!=NULL ||
      sortEvals_calcDirectMat!=NULL  ||
      infCoefs_calcDirectMat!=NULL   ||
      g2gp_calcDirectMat!=NULL       ||
      precVals_calcDirectMat!=NULL   ||
      precVals_sort_calcDirectMat!=NULL||
      dRows_calcDirectMat!=NULL){
    fprintf(stderr,"%s: Memory is already allocated. \n"
	    "    Dont use this routine for re-allocation of memory\n"
	    "    This halts the program\n",
	    fctName);
    _EXIT_;
  }
  /* Set number of matrices - used for deallocating memory later       */
  nMats_calcDirectMat = nMats;
  /* Set the allocation size:                                          */
  nEvalsMax_calcDirectMat = MAX(nEvalsMax_calcDirectMat,nEvals);
  /* Shorthand (to increase the readability of this routine)           */
#define NN nEvalsMax_calcDirectMat
  /* Actually allocate storage:                                        */
  /* Memory for the sorting process                                    */
  isortEvals_calcDirectMat = (int*) calloc(NN, sizeof(int));
  iorigEvals_calcDirectMat = (int*) calloc(NN, sizeof(int));
  /* The present calling routine does not recycle the list of neighbours. 
   * However, to increase the ease of implementation as well as the 
   * similarity with the direct computation, allocate new storage 
   * regardless. There should be room for optimization here.           */
  sortEvals_calcDirectMat  = 
    (struct constantElement**) 
    calloc(NN, sizeof(struct constantElement*));
  /* Influence coefficients                                            */
  infCoefs_calcDirectMat = (double**) calloc(nMats,sizeof(double*));
  for (i=0;i<nMats;i++){
    infCoefs_calcDirectMat[i] = (double*) calloc(NN, sizeof(double));
  }
  /* For precdorrection step                                           */
  g2gp_calcDirectMat          = (double*) calloc(NN, sizeof(double));
  precVals_calcDirectMat      = (double*) calloc(NN, sizeof(double));
  precVals_sort_calcDirectMat = (double*) calloc(NN, sizeof(double));
  dRows_calcDirectMat = (int*) calloc(NN, sizeof(int));
#undef NN
  /* Memory is now allocated - set flag:                               */
  isAllocated_calcDirectMat = 1;
  return;
} /* end of routine calcDirectMatCalloc */
/* 
* --------------------------------------------------------------------
* Re-allocate the memory for calcDirectMat. 
* Make sure this routine is not called too often, by increasing 
* the memory by at least a factor two.
* Note that this routine does *not* preserve the data in the arrays
* being reallocated (as does the routine realloc). Rather, the arrays are
* first freed and then allocated, destroying all the data in the process.
* -------------------------------------------------------------------- */
static
void calcDirectMatRealloc(int nEvals)
{
  char fctName[]="calcDirectMatRealloc";  /* Name of this routine     */
  int newMax, nMats;
  /*printf("%s called, %d evals\n",fctName,nEvals);*/
  /* Test if the memory is already allocated (as it should be)         */
  if (!isAllocated_calcDirectMat ||
      isortEvals_calcDirectMat    == NULL ||
      iorigEvals_calcDirectMat    == NULL ||
      sortEvals_calcDirectMat     == NULL ||
      infCoefs_calcDirectMat      == NULL ||
      g2gp_calcDirectMat          == NULL ||
      precVals_calcDirectMat      == NULL ||
      precVals_sort_calcDirectMat == NULL ||
      dRows_calcDirectMat         == NULL){
    fprintf(stderr,"%s: ERROR.\n"
	    "    Memory should be re-allocated, but isn't allocated.\n"
	    "    (this may indicate a bug. Better fix it!)\n"
	    "    In the meantime I'll stop the program\n",
	    fctName);
    exit(1);
  }
  /* Keep value of nMats:                                              */
  nMats = nMats_calcDirectMat;
  /* Set new max size:                                                 */
  newMax = MAX(2*nEvalsMax_calcDirectMat,nEvals);
  /* Free old memory:                                                  */
  calcDirectMatFree();
  /* Allocate new memory:                                              */
  calcDirectMatCalloc(newMax,nMats);
  return;
} /* End of routine calcDirectMatRealloc */
/*
* --------------------------------------------------------------------
* Free the memory used by calcDirectMat.
* -------------------------------------------------------------------- */
static
void calcDirectMatFree()
{
  char fctName[]="calcDirectMatFree";  /* Name of this routine         */
  int writeWarning, i;
  int verbose = 1;
  /* If verbose and the memory is not allocated, then write a warning  */
  writeWarning=!isAllocated_calcDirectMat;
  /* Free memory. If an array is not allocated, then toggle the 
   * warning to be written                                             */
  /* When deallocating - then let the pointer be reset to NULL, so that
   * other routines will know that it is deallocated!                  */
#define LOCALFREE(ARRAY) if (ARRAY != NULL){ \
                           FREE(ARRAY)}      \
                         else writeWarning=1;
  LOCALFREE(isortEvals_calcDirectMat);
  LOCALFREE(iorigEvals_calcDirectMat);
  LOCALFREE(sortEvals_calcDirectMat);
  LOCALFREE(g2gp_calcDirectMat);
  LOCALFREE(precVals_calcDirectMat);
  LOCALFREE(precVals_sort_calcDirectMat);
  LOCALFREE(dRows_calcDirectMat);
#undef LOCALFREE
  if (infCoefs_calcDirectMat!=NULL){
    for (i=0;i<nMats_calcDirectMat;i++){
      FREE(infCoefs_calcDirectMat[i]);
    }
    FREE(infCoefs_calcDirectMat);
  }
  else writeWarning=1;

  if(verbose && writeWarning){
    fprintf(stderr,"%s: WARNING!\n"
	    "    Some (or all) of the memory could not be deallocated"
	    "    since it was not allocated to begin with!\n",
	    fctName);
  }
  /* In debug mode make sure that all pointers now point to NULL
   * ("FREE" should set the pointer to NULL after deallocating memory) */
  assert(isortEvals_calcDirectMat==NULL);
  assert(iorigEvals_calcDirectMat==NULL);
  assert(sortEvals_calcDirectMat==NULL);
  assert(sortEvals_calcDirectMat==NULL);
  assert(g2gp_calcDirectMat==NULL);
  assert(precVals_calcDirectMat==NULL);
  assert(precVals_sort_calcDirectMat==NULL);
  assert(dRows_calcDirectMat==NULL);

  /* Memory is now deallocated. Set flag to show this:                 */
  isAllocated_calcDirectMat = 0;
  return;
} /* End of routine calcDirectMatFree */

/*
* ====================================================================
* The following routines and static variables go together as parts of 
* the preconditioning step. 
* At the present calcPrecondMat_1, the preconditioner routine itself, 
* checks if the storage has been allocated - and takes care of the 
* allocation. If the storage needs to be released for other use, 
* then calcPrecondMatFree needs to be called. calcPrecondMat wont 
* do this, since it doesn't know if it will be called again.
* Use:
*  1)  Allocate needed storage (once):  calcPrecondMatCalloc
*  2)  Calls to calcp (multiple times): calcPrecondMat_1
*  2b) Reallocate storage if needed:    calcPrecondMatRealloc
*  3)  Free the used memory (once):     calcPrecondMatFree
*
* Different preconditioners may be implemented - all using the same 
* memory for working space.
* ==================================================================== */
static int 
  isAllocated_calcPrecondMat=0, /* Flag for allocation status.
				 * 0: Not allocated
				 * 1: Allocated                        */
  nEvalsMax_calcPrecondMat=50,  /* Size of arrays being allocated, set 
				 * to a minimum default value          */
  nMats_calcPrecondMat,         /* Used for keeping track of memory 
				 * allocated for the influence 
				 * coefficients - needed in the 
				 * deallocation step.                  */
  *isortEvals_calcPrecondMat,   /* Arrays for sorting after global     */
  *iorigEvals_calcPrecondMat,   /* index the passed list of 
				 * neighbouring evaluation elements.   */
  *isortSrcs_calcPrecondMat,    /* Arrays for sorting after global     */
  *iorigSrcs_calcPrecondMat,    /* index the passed list of 
				 * source elements.                    */
  *iCalcEvals_calcPrecondMat,   /* Global index of evaluation elements 
				 * for which influence coefficients 
				 * need to be calculated.              */
  *iCalcEvalsLoc_calcPrecondMat;/* Local row index in a small dense 
				 * matrix - to store an influence 
				 * coefficient once calculated.        */
static struct constantElement 
  **sortEvals_calcPrecondMat,   /* Pointers to neighbouring evaluation 
				 * elements, sorted after global index.*/
  **sortSrcs_calcPrecondMat,    /* Pointers to source elements, 
				 * sorted after global index.          */
  **calcEvals_calcPrecondMat;   /* Pointers to those elements for which 
				 * influence coefficients need to be 
				 * calculated.                         */
static double 
  *mat_calcPrecondMat,          /* Small dense matrix for solving the 
				 * local problem. This matrix will be 
				 * LU-factored for solving each local 
				 * problem.                            */
  *column_calcPrecondMat,       /* Vector for extracing columns of the 
				 * inverse matrix for the local problem.
				 * This is both the "right hand side" 
				 * and the solution upon return from a
				 * backsolve.                          */
  **infCoefs_calcPrecondMat;    /* Storage for those influence 
				 * coefficients that are calculated 
				 * in the preconditioning step.        */
/* 
* --------------------------------------------------------------------
* Given a set of sources and a list of neighbours, set up the part 
* of the preconditioning matrix that relates to these sources.
*
* This version deals only with triangular and quadrilateral 
* "constant panels", as used in e.g. wamit. Thus "calcp" can be
* used for the actual integration of the entries that cannot be 
* taken from the direct part.
*
* Note that this function must conform to the layout of a "target 
* function" of "sourcesAndNeighboursOperate" (see e.g. grid.c).
* This function is a companion to "calcDirectMat_1". These two routines
* should be modified together for consistency.
*
* Given time, the present routine should be rewritten so that it 
* does not require information on the structure of sparse matrices. 
* -------------------------------------------------------------------- */
void calcPrecondMat_1(void **sourcesIn,
		      void **evalsIn, 
		      int nSources,
		      int nTotNeighbours,
		      void *dataIn)
{
  char fctName[] = "calcPrecondMat_1";
  int verbose = 1;
  char mark = '>';
  /* For recasting input data:                                         */
  struct constantElement **sources,**evals;
  sparseMat *directMat, *precondMat;
  /* Shorthands to the direct matrix                                   */
  double **matVals;
  int *rowIndex;
  /* Other local variables                                             */
  int i,ievl,isrc, 
    iColDirect,iRowDirect,
    iColPrecond,
    iColDense,iRowDense,
    iLastRow, iSearchRow, nStore, nCompute, nCompLoc, iColStart;
  static int nColsDone; 
  struct constantElement *source;
  struct cpm1Bundle *bundle;
#ifdef MYDEBUGFLAG
  int verbose=0;
#endif

  sources = (struct constantElement**) sourcesIn;
  evals   = (struct constantElement**) evalsIn;
  bundle = (struct cpm1Bundle*) dataIn;
  directMat  = (sparseMat *) bundle->dMat;
  precondMat = (sparseMat *) bundle->pMat;

  /* Allocate memory needed for this routine.                          */
  if(!isAllocated_calcPrecondMat) {
    calcPrecondMatCalloc(nTotNeighbours,precondMat->nMats);
    /* Initiate plotting progress marks */
    if (verbose) {
      /* Nothing done yet. */
      plotProgressMark(nColsDone = 0, 1, 0, directMat->n,1,mark);
    }
  }
  /* Check that enough memory has been allocated                       */
  else if(nTotNeighbours>nEvalsMax_calcPrecondMat)
    calcPrecondMatRealloc(nTotNeighbours);

  /* Sort the neighbours.                                              */
  for (ievl=0;ievl<nTotNeighbours;ievl++){
    isortEvals_calcPrecondMat[ievl] = evals[ievl]->directIndex;
    iorigEvals_calcPrecondMat[ievl] = ievl;
  }
  intSort2(isortEvals_calcPrecondMat,
	   iorigEvals_calcPrecondMat,
	   nTotNeighbours);
  for (ievl=0;ievl<nTotNeighbours;ievl++){
    sortEvals_calcPrecondMat[ievl] 
      = evals[iorigEvals_calcPrecondMat[ievl]];
  }

  /* Sort the sources. 
   * It is doubtful whether this results in a speedup. However, it 
   * should not slow down the process. Sorting the sources reduces 
   * the search time when copying data from the small dense matrix 
   * to the sparse preconditioning matrix.                             */
  for (isrc=0;isrc<nSources;isrc++){
    isortSrcs_calcPrecondMat[isrc] = sources[isrc]->directIndex;
    iorigSrcs_calcPrecondMat[isrc] = isrc;
  }
  intSort2(isortSrcs_calcPrecondMat,
	   iorigSrcs_calcPrecondMat,
	   nSources);
  for (isrc=0;isrc<nSources;isrc++){
    sortSrcs_calcPrecondMat[isrc] 
      = sources[iorigSrcs_calcPrecondMat[isrc]];
  }

  /* A matrix of size n x n needs to be 
   * 1) allocated
   * 2) constructed 
   * 3) LU-decomposed
   *
   * Rows (columns?) of the inverse matrix are extracted, using 
   * (dense) matrix-vector products, and subsequently stored in
   * the (sparse) preconditioning matrix.                              */

  /* Short hands:                                                      */
  rowIndex  = directMat->rowIndex;
  matVals   = directMat->matVals;

  /* nStore and nCompute are only used to get an idea of how large a 
   * part of the small preconditioning matrix that needs to be computed 
   * and how much can be taken from the direct part. nStore and nCompute 
   * is not really needed for the actual computations.                 */
  nStore = nCompute = 0;
  /* Construct column # iColDense of the small dense matrix.           */
  for (iColDense=0;iColDense<nTotNeighbours;iColDense++){
    /* The corresponding column of the direct matrix.
     * IMPORTANT: This assumes that the sources and the evals 
     * "are the same":                                                 */
    iColDirect = sortEvals_calcPrecondMat[iColDense]->directIndex;
    /* If possible copy entry # iRowDense of this column of the 
     * dense matrix from the direct sparse matrix. 
     * Start searching at the first nonzero (row) entry in this column 
     * in the direct matrix                                            */
    iSearchRow=directMat->colIndex[iColDirect];
    /* Find the number of the last nonzero (row) entry in this column
     * of the direct matrix                                            */
    iLastRow = directMat->rowIndex[directMat->colIndex[iColDirect+1]-1];
    /* Counter for how many entries needs to be computed for 
     * this column                                                     */
    nCompLoc = 0;
    for (iRowDense=0;iRowDense<nTotNeighbours;iRowDense++){
      /* The (global) index number for this evaluation element, 
       * corresponding to the number of the related row in the direct
       * matrix, is:                                                   */
      iRowDirect = sortEvals_calcPrecondMat[iRowDense]->directIndex;
      /* If the entry (iRowDirect,iColDirect) is non-zero (exists) in 
       * the direct part sparse matrix, then we can used it for entry 
       * (iRowDense,iColDense) in the local dense matrix. If the entry 
       * is not in the direct matrix, then we need to calculate it.    */
      /* Since the entries are sorted in the direct matrix it is not 
       * necessary to scan the whole column each time. Just continue 
       * from where we left of the last time.                          */
      while (iSearchRow<directMat->colIndex[iColDirect+1]){
	/* See if row # iSearchRow is the one we want                  */
	if (rowIndex[iSearchRow]==iRowDirect){
#ifndef NORECYCLE /* This is the normal procedure.                     */
	  /* Store the corresponding entry in the small dense matrix   */
	  mat_calcPrecondMat[SQDEX(iRowDense,iColDense,nTotNeighbours)] 
	    = matVals[0][iSearchRow];
	  nStore++;
          break; /* Search no longer - we found the needed entry!      */
#else /* This forces the preconditioner to compute everything! 
       * Use this only as part of debugging! Even though we have now 
       * found a usable value, discard it, and compute a new one!      */
	  iCalcEvalsLoc_calcPrecondMat[nCompLoc] = iRowDense;
	  iCalcEvals_calcPrecondMat[nCompLoc]    = iRowDirect;
	  calcEvals_calcPrecondMat[nCompLoc] 
	    = sortEvals_calcPrecondMat[iRowDense];
	  nCompLoc++;
	  iSearchRow++; /* Go to next entry */
	  break;
#endif
	}
	else if(rowIndex[iSearchRow]>iRowDirect || iRowDirect>iLastRow){
	  /* The needed entry is not in the direct matrix. Thus, store 
	   * the entry among the "stuff to compute". 
	   * Register iRowDense and iRowDirect for later reference.    */
	  iCalcEvalsLoc_calcPrecondMat[nCompLoc] = iRowDense;
	  iCalcEvals_calcPrecondMat[nCompLoc]    = iRowDirect;
	  /* Add the proper element to the (list of) "evals" for 
	   * computation of influence coeficients                      */
	  calcEvals_calcPrecondMat[nCompLoc] 
	    = sortEvals_calcPrecondMat[iRowDense];
	  /* Update number of influence coeficients to compute         */
	  nCompLoc++;
	  break;/* Search no longer - The needed entry is not here!    */
	}
	iSearchRow++;
      } /* Next element in this column of the direct sparse matrix 
	 * (only if rowIndex[iSearchRow] < iRowDirect)                 */

    } /* Next entry in this row in the small dense matrix 
       * (next iRowDense)                                              */

    if (nCompLoc){/* If nCompLoc is nonzero then we need to calculate 
		   * some entries to the small matrix. Not all the 
		   * entries had been calculated for the direct matrix */
      /* Add to the total number of computed entries                   */
      nCompute += nCompLoc;
      /* The source element corresponding to the present column
       * NOTE that this assumes that the sources and the evals "are 
       * the same" - just as when finding the global column index.     */
      source = sortEvals_calcPrecondMat[iColDense];
      /* Calculate the missing entries for this column (row?) 
       * in the small dense matrix using a (macro for a) wrapper for 
       * calcp.                                                        */
      CALCPWRAPPER(source,calcEvals_calcPrecondMat,nCompLoc,
		   infCoefs_calcPrecondMat);
      /* Store the calculated entries in the small dense matrix        */
      for (ievl=0;ievl<nCompLoc;ievl++){
	iRowDense = iCalcEvalsLoc_calcPrecondMat[ievl];
	mat_calcPrecondMat[SQDEX(iRowDense,iColDense,nTotNeighbours)] 
	  = infCoefs_calcPrecondMat[0][ievl];
      }
    }
  } /* Next column (row?) in the small dense matrix (next iColDense)   */

  /* This is a simple check to catch missmatches in the numbers of
   * accessed entries (obvious bugs)                                   */
  if (nStore+nCompute!=nTotNeighbours*nTotNeighbours)
    printf("%s: %d+%d=%d, but I expected %d?? FIX THIS!\n",
	   fctName,nStore,nCompute,nStore+nCompute,
	   nTotNeighbours*nTotNeighbours);
#ifdef MYDEBUGFLAG
  if (verbose)
    printf("%s: Would use %d stored values (%g%%) "
	   "and compute %d values (%g%%)!\n"
	   "   for %d evals (and %d other neighbours) i.e. a %dx%d system\n",
	   fctName,
	   nStore,nStore*100.0/(nTotNeighbours*nTotNeighbours),
	   nCompute,nCompute*100.0/(nTotNeighbours*nTotNeighbours),
	   nSources,nTotNeighbours-nSources,nTotNeighbours,nTotNeighbours); 
  /* Test block
  {
    dumpDenseMatrix(mat_calcPrecondMat,nTotNeighbours);
    dumpSparseMatrixTranspose((void*) directMat,0);
  } 
  */
#endif

  /* LU factor the small dense matrix                                  */
  LUdecompNoPivoting(mat_calcPrecondMat,nTotNeighbours);

  /* Short hands (for preconditioning matrix rather than direct matrix)*/
  rowIndex  = precondMat->rowIndex;
  matVals   = precondMat->matVals;
  /* Extract the needed columns (rows?) from the inverse of the small 
   * matrix and store in the sparse preconditioning matrix             */
  for (isrc=0,ievl=0;isrc<nSources;isrc++){
    /* Find the global number of this source, i.e. the column number 
     * in the preconditioning matrix.                                  */
    iColPrecond = isortSrcs_calcPrecondMat[isrc];
    /* Only for tests:                                                 */
    assert(iColPrecond==sortSrcs_calcPrecondMat[isrc]->directIndex);
    /* Find index for first entry in this column of the preconditioning
     * matrix                                                          */
    iColStart = precondMat->colIndex[iColPrecond];
    /* Find the column number in the local dense matrix                */
    for (;;ievl++){
      if (iColPrecond==isortEvals_calcPrecondMat[ievl]){
	iColDense = ievl;
	ievl++; /* Go to next row number (not very important) */
	break;  /* Exit the for-loop. Found the needed column number.  */
      }
      else if (ievl>=nTotNeighbours){
	/* The control flow should never get inside this if statement.
	 * This "else if" was put in here rather than using a limit 
	 * in the "for" statement in order to catch bugs.              */
	fprintf(stderr,
		"%s: ERROR: Index of evaluation element is out of bounds!\n"
		"     This error is terminal\n",
		fctName);
	exit(1);
      }
    }
    /* Setup the needed right-hand-side:                               */
    for (i=0;i<nTotNeighbours;i++){
      column_calcPrecondMat[i]=0.0;
    }
    column_calcPrecondMat[iColDense]=1.0;
    /* Extract the needed column of the inverse dense matrix           */
    LUbacksolveNoPivoting(mat_calcPrecondMat,
			  nTotNeighbours,
			  column_calcPrecondMat);
    /* Copy this column into the sparse precondition matrix            */
    for(i=0;i<nTotNeighbours;i++){
      /* Copy this entry into the sparse precondition matrix           */
      rowIndex[iColStart+i]   = isortEvals_calcPrecondMat[i];
      matVals[0][iColStart+i] = column_calcPrecondMat[i];
    }
  }
  /* We have now done nSources more columns. Add to total 
   * and update ...? */
  nColsDone += nSources;
  if (verbose) { /* Plot progress */
    if (nColsDone == directMat->n) /* All done */ 
      plotProgressMark(nColsDone, -1, 0, directMat->n,0,mark);
    else  
      plotProgressMark(nColsDone,  0, 0, directMat->n,0,mark);
  }

  return;
}/* End of routine calcPrecondMat_1                                    */
/* 
* --------------------------------------------------------------------
* Allocate needed memory for calcPrecondMat .
* A maximum is set on the number of expected evals. If the actual 
* number of evals exceeds this number, the memory will be reallocated
* later.
* -------------------------------------------------------------------- */
static 
void calcPrecondMatCalloc(int nEvals, int nMats)
{
  char fctName[]="calcPrecondMatCalloc";
  int i;
  /* Make absolutely sure that this memory has not been allocated      */
  if (isAllocated_calcPrecondMat) {
    fprintf(stderr,"%s: Apparently memory is already allocated.\n"
	    "    (this may indicate a bug. Better fix it!)\n"
	    "    In the mean time I'll stop the program\n",
	    fctName);
    exit(1);
  }
  if (isortEvals_calcPrecondMat!=NULL ||
      iorigEvals_calcPrecondMat!=NULL ||
      isortSrcs_calcPrecondMat!=NULL ||
      iorigSrcs_calcPrecondMat!=NULL ||
      iCalcEvals_calcPrecondMat!=NULL ||
      iCalcEvalsLoc_calcPrecondMat!=NULL ||
      sortEvals_calcPrecondMat!=NULL ||
      sortSrcs_calcPrecondMat!=NULL ||
      calcEvals_calcPrecondMat!=NULL ||
      mat_calcPrecondMat!=NULL ||
      column_calcPrecondMat!=NULL ||
      infCoefs_calcPrecondMat!=NULL){
    fprintf(stderr,"%s: Memory is already allocated. \n"
	    "    Dont use this routine for re-allocation of memory\n"
	    "    This halts the program\n",
	    fctName);
    exit(1);
  }
  /* Set number of matrices - used for deallocating memory later       */
  nMats_calcPrecondMat = nMats;
  /* Set the allocation size:                                          */
  nEvalsMax_calcPrecondMat = MAX(nEvalsMax_calcPrecondMat,nEvals);
  /* Shorthand (to increase the readability of this routine)           */
#define NN nEvalsMax_calcPrecondMat
  /* Actually allocate storage:                                        */
  /* Memory for the sorting process                                    */
  isortEvals_calcPrecondMat = (int*) calloc(NN, sizeof(int));
  iorigEvals_calcPrecondMat = (int*) calloc(NN, sizeof(int));
  isortSrcs_calcPrecondMat = (int*) calloc(NN, sizeof(int));
  iorigSrcs_calcPrecondMat = (int*) calloc(NN, sizeof(int));
  /* The present calling routine does not recycle the list of neighbours. 
   * However, to increase the ease of implementation as well as the 
   * similarity with the direct computation, allocate new storage 
   * regardless. There should be room for optimization here.           */
  sortEvals_calcPrecondMat  = 
    (struct constantElement**) 
    calloc(NN, sizeof(struct constantElement*));
  sortSrcs_calcPrecondMat  = 
    (struct constantElement**) 
    calloc(NN, sizeof(struct constantElement*));
  /* Memory for storing indices of those elements that needs to be 
   * computed                                                          */
  iCalcEvals_calcPrecondMat = (int*) calloc(NN, sizeof(int));
  iCalcEvalsLoc_calcPrecondMat = (int*) calloc(NN, sizeof(int));
  /* Similar structure for pointers to elements                        */
  calcEvals_calcPrecondMat = 
    (struct constantElement**) 
    calloc(NN, sizeof(struct constantElement*));
  /* Dense matrix and right-hand-side:                                 */
  mat_calcPrecondMat = (double*) calloc(NN*NN,sizeof(double));
  column_calcPrecondMat = (double*) calloc(NN,sizeof(double));
  /* Influence coefficients                                            */
  infCoefs_calcPrecondMat = (double**) calloc(nMats,sizeof(double*));
  for (i=0;i<2;i++){
    infCoefs_calcPrecondMat[i] = (double*) calloc(NN, sizeof(double));
  }
#undef NN
  /* Memory is now allocated. Set Flag:                                */
  isAllocated_calcPrecondMat = 1;
  return;
} /* end of routine calcPrecondMatCalloc */
/* 
* --------------------------------------------------------------------
* Re-allocate the memory for calcPrecondMat. 
* Make sure this routine is not called too often, by increasing 
* the memory by at least a factor two.
* Note that this routine does *not* preserve the data in the arrays
* being reallocated (as does the routine realloc). Rather, the arrays are
* first freed and then allocated, destroying all the data in the process.
* -------------------------------------------------------------------- */
static
void calcPrecondMatRealloc(int nEvals)
{
  char fctName[]="calcPrecondMatRealloc";  /* Name of this routine     */
  int newMax, nMats;
  /* Test if the memory is already allocated (as it should be)         */
  if (!isAllocated_calcPrecondMat ||
      isortEvals_calcPrecondMat    == NULL ||
      iorigEvals_calcPrecondMat    == NULL ||
      isortSrcs_calcPrecondMat     == NULL ||
      iorigSrcs_calcPrecondMat     == NULL ||
      iCalcEvals_calcPrecondMat    == NULL ||
      iCalcEvalsLoc_calcPrecondMat == NULL ||
      sortEvals_calcPrecondMat     == NULL ||
      sortSrcs_calcPrecondMat      == NULL ||
      calcEvals_calcPrecondMat     == NULL ||
      mat_calcPrecondMat           == NULL ||
      column_calcPrecondMat        == NULL ||
      infCoefs_calcPrecondMat      == NULL){
    fprintf(stderr,"%s: \n"
	    "    Memory should be re-allocated, but isn't allocated.\n"
	    "    (this may indicate a bug. Better fix it!)\n"
	    "    In the meantime I'll stop the program\n",
	    fctName);
    exit(1);
  }
  /* Keep value of nMats:                                              */
  nMats = nMats_calcPrecondMat;
  /* Set new max size:                                                 */
  newMax = MAX(2*nEvalsMax_calcPrecondMat,nEvals);
  /* Free old memory:                                                  */
  calcPrecondMatFree();
  /* Allocate new memory:                                              */
  calcPrecondMatCalloc(newMax,nMats);
  return;
} /* End of routine calcPrecondMatRealloc */
/*
* --------------------------------------------------------------------
* Free the memory used by calcPrecondMat.
* -------------------------------------------------------------------- */
static
void calcPrecondMatFree()
{
  char fctName[]="calcPrecondMatFree";  /* Name of this routine        */
  int writeWarning, i;
  int verbose = 1;
  /* If verbose and the memory is not allocated, then write a warning  */
  writeWarning=!isAllocated_calcPrecondMat;
  /* Free memory. If an array is not allocated, then toggle the 
   * warning to be written                                             */
  /* When deallocating - then let the pointer be reset to NULL,
   * so that other routines will know that it is deallocated!          */
#define LOCALFREE(ARRAY) if (ARRAY != NULL){ \
                           FREE(ARRAY)}      \
                         else writeWarning=1;
  LOCALFREE(isortEvals_calcPrecondMat);
  LOCALFREE(iorigEvals_calcPrecondMat);
  LOCALFREE(isortSrcs_calcPrecondMat);
  LOCALFREE(iorigSrcs_calcPrecondMat);
  LOCALFREE(iCalcEvals_calcPrecondMat);
  LOCALFREE(iCalcEvalsLoc_calcPrecondMat);
  LOCALFREE(sortEvals_calcPrecondMat);
  LOCALFREE(sortSrcs_calcPrecondMat);
  LOCALFREE(calcEvals_calcPrecondMat);
  LOCALFREE(mat_calcPrecondMat);
  LOCALFREE(column_calcPrecondMat);
#undef LOCALFREE
  if (infCoefs_calcPrecondMat!=NULL){
    for (i=0;i<2;i++){
      FREE(infCoefs_calcPrecondMat[i]);
    }
    FREE(infCoefs_calcPrecondMat);
  }
  else writeWarning=1;

  if(verbose && writeWarning){
    fprintf(stderr,"%s: WARNING!\n"
	    "    Some (or all) of the memory could not be deallocated"
	    "    since it was not allocated to begin with!\n",
	    fctName);
  }
  /* Memory is now free - set the flag to reflect this:                */
  isAllocated_calcPrecondMat=0;
  return;
} /* End of routine calcPrecondMatFree */

/*
* ====================================================================
* The following routines and static variables go together as a
* wrapper to the calcp routine. At the present calcpWrapper itself 
* checks if the storage has been allocated - and takes care of the 
* allocation. If the storage needs to be released for other use, 
* then calcpWrapperFree needs to be called. calcpWrapper wont do this, 
* since it doesn't know if it will be called again.
* Use:
*  1)  Allocate needed storage (once):  calcpWrapperCalloc
*  2)  Calls to calcp (multiple times): calcpWrapper
*  2b) Reallocate storage if needed:    calcpWrapperRealloc
*  3)  Free the used memory (once):     calcpWrapperFree
*
* Different wrappers may be implemented - all using the same memory
* for working space.
* ==================================================================== */
static int isAllocated_calcpWrapper=0,nEvalsMax_calcpWrapper=1000;
static double 
  **panel_calcpWrapper,             
  **evalpnts_calcpWrapper, 
  **directions_calcpWrapper,
  *fss_calcpWrapper, 
  *fds_calcpWrapper, 
  *fess_calcpWrapper, 
  *feds_calcpWrapper;
/* 
* --------------------------------------------------------------------
* Given a source and a number of evals, return the wanted influence
* coefficients, by calling calcp.
* Basically, this routine acts as an interface (hopefully a fairly 
* clean one) to calcp.
* To do: The present allocation and deallocation should be substituted 
* by allocating needed memory at the beginning of the simulation (before
* the first call) and then keeping that throughout the calculation, 
* deallocating only once it is not needed again. (An if statement in 
* here may then be needed to check if enough memory is present...)
* -------------------------------------------------------------------- */
static
void calcpWrapper(struct constantElement *source,
		  struct constantElement **evals, 
		  int nEvals,
		  double **infCoefs)
{
  int i;
  point X,Y,Z, *corners;
  double znorm;

  /* Allocate memory if it has not been done before:                   */
  if (!isAllocated_calcpWrapper) calcpWrapperCalloc(nEvals);
  /* Check if enough memory have been allocated, if not then 
   * reallocate memory                                                 */
  else if (nEvals>nEvalsMax_calcpWrapper) calcpWrapperRealloc(nEvals);

  /* In debug mode make a few tests on the status of the memory:       */
  assert(panel_calcpWrapper!=NULL);
  assert(evalpnts_calcpWrapper!=NULL);
  assert(directions_calcpWrapper!=NULL);
  assert(fss_calcpWrapper!=NULL);
  assert(fds_calcpWrapper!=NULL);
  assert(fess_calcpWrapper!=NULL);
  assert(feds_calcpWrapper!=NULL);

  /* Copy data for corner coordinates                                  */
  for (i=0;i<source->shape;i++){
    VASSIGN(panel_calcpWrapper[i],source->corners[i]);
  }

  /* Let the evaluation points (collocation points) be the centroids 
   * of the panels:                                                    */
  for (i=0;i<nEvals;i++){
    VASSIGN(evalpnts_calcpWrapper[i],evals[i]->centroid);
  }
  /*Test output: 
    printf("Calling calcp. src=%d, eval=",source->directIndex);
    for(i=0;i<nEvals;i++)printf(" %d",evals[i]->directIndex);
    printf("\n");*/
  /* Calculate Greens function...                                      */
  calcp(panel_calcpWrapper,            /* "panel"    */
	source->shape,                 /* "verts"    */
	evalpnts_calcpWrapper,         /* "evalpnts" */
	NULL,                          /* No directions */
	nEvals,                        /* numeval    */
	fss_calcpWrapper, 
	fds_calcpWrapper, 
	fess_calcpWrapper, 
	feds_calcpWrapper);
  /* Store information in return array     */
  for (i=0;i<nEvals;i++){
    infCoefs[0][i] = fss_calcpWrapper[i]; /* This is Int G             */
    infCoefs[1][i] = fds_calcpWrapper[i]; /* This is Int Gn            */
  }

  /* Note: A linear combination of G and Gn could be used for each 
   * eval to signify the entry to store, e.g. a*G+b*Gn.
   * Maybe one entry for each eval for each matrix to calculate is 
   * a better way... In this case less CPU time may be wasted, i.e. 
   * memory *can* be sacrificed for CPU time if necessary/convenient. 
   * Actually, it may be pretty easy to make "clones" of the present 
   * routine that do either G, Gn, a linear combination, G and Gn, 
   * multiple linear combinations or whatever is needed. 
   * Using sourceAndNeighboursOperate from setupDirectSparse may make 
   * such an approach attractive and easy to manage. Note that the 
   * data for linear combinations and what matrices to set up etc 
   * must be "bundled" with the sparse matrix in a structure, such 
   * that a generic pointer can be passed through the grid routines.
   * Alternatives:
   *   Store "Diriclet/Neumann" status with the panels, and
   *   calculate G or Gn accordingly (an input parameter can then
   *   flag the calculation of Gn or G, i.e. the other order).
   */

  return;

  /* The remaining part of this routine is not used at the present.
   * It shows how to calculate directions as normals on each panel.    */
  for (i=0;i<nEvals;i++){
    VASSIGN(evalpnts_calcpWrapper[i],evals[i]->centroid);

    /* FOR NOW use the boundary normal vectors at the collocation 
     * points (element centroids) for "directions" in calcp.           */
    /* Calculate the panel coordinate system. */ 
    corners = evals[i]->corners;
    VOP(X, corners[2], -, corners[0]);
    if(evals[i]->shape == 3) {
      VOP(Y, corners[1], -, corners[0]); 
    } else {
      VOP(Y, corners[1], -, corners[3]); 
    }
    /* Calculate Z-axis is normal to two diags, i.e. pointing 
     * normal direction of the corners.                                */
    XPROD(Z, X, Y); 
    /* Normalize normal vector and assign this to be the direction to
     * use for the corresponding collocation point.                    */
    znorm = VNORM(Z);
    /* znorm = sqrt(Z[0]*Z[0]+Z[1]*Z[1]+Z[2]*Z[2]);*/
    directions_calcpWrapper[i][0] = Z[0] / znorm;
    directions_calcpWrapper[i][1] = Z[1] / znorm;
    directions_calcpWrapper[i][2] = Z[2] / znorm;
  }
  return;
} /* End of routine calcpWrapper */
/* 
* --------------------------------------------------------------------
* Allocate needed memory for calcpWrapper routines for calls to calcp.
* A maximum is set on the number of expected evals. If the actual 
* number of evals exceeds this number, the memory will be reallocated
* later.
* -------------------------------------------------------------------- */
static
void calcpWrapperCalloc(int nEvals)
{
  char fctName[]="calcpWrapperCalloc";  /* Name of this routine        */
  int i;
  /* Make absolutely sure that this memory has not been allocated      */
  if (isAllocated_calcpWrapper) {
    fprintf(stderr,"%s: Apparently memory is already allocated.\n"
	    "    (this may indicate a bug. Better fix it!)\n"
	    "    In the mean time I'll stop the program\n",
	    fctName);
    exit(1);
  }
  if (panel_calcpWrapper!=NULL ||
      evalpnts_calcpWrapper!=NULL ||
      directions_calcpWrapper!=NULL ||
      fss_calcpWrapper!=NULL ||
      fds_calcpWrapper!=NULL ||
      fess_calcpWrapper!=NULL||
      feds_calcpWrapper!=NULL ) {
    fprintf(stderr,"%s: Memory is already allocated. \n"
	    "    Dont use this routine for re-allocation of memory\n"
	    "    This halts the program\n",
	    fctName);
    exit(1);
  }
  /* See if the default size is big enough                             */
  nEvalsMax_calcpWrapper = MAX(nEvalsMax_calcpWrapper,nEvals);

  panel_calcpWrapper = (double**) calloc(4, sizeof(double*));
  for (i=0;i<4;i++){
    panel_calcpWrapper[i] = (double*) calloc(3, sizeof(double));
  }
  evalpnts_calcpWrapper   
    = (double**) calloc(nEvalsMax_calcpWrapper, sizeof(double*));
  directions_calcpWrapper 
    = (double**) calloc(nEvalsMax_calcpWrapper, sizeof(double*));
  for (i=0;i<nEvalsMax_calcpWrapper;i++){
    evalpnts_calcpWrapper[i]   = (double*) calloc(3, sizeof(double));
    directions_calcpWrapper[i] = (double*) calloc(3, sizeof(double));    
  }
  fss_calcpWrapper  = calloc(nEvalsMax_calcpWrapper, sizeof(double));
  fds_calcpWrapper  = calloc(nEvalsMax_calcpWrapper, sizeof(double));
  fess_calcpWrapper = calloc(nEvalsMax_calcpWrapper, sizeof(double));
  feds_calcpWrapper = calloc(nEvalsMax_calcpWrapper, sizeof(double));
  /* Memory is now allocated - set flag to show this                    */
  isAllocated_calcpWrapper=1;
  return;
} /* End of routine calcpWrapperCalloc */
/* 
* --------------------------------------------------------------------
* Re-allocate the memory for calcpWrapper that depends on the number
* of passed evaluation points. Make sure this routine is not called 
* too often, by increasing the memory by at least a factor two.
* Note that this routine does *not* preserve the data in the arrays
* being reallocated (as does the routine realloc). Rather, the arrays are 
* first freed and then allocated, destroying all the data in the process.
* -------------------------------------------------------------------- */
static
void calcpWrapperRealloc(int nEvals)
{
  char fctName[]="calcpWrapperRealloc";  /* Name of this routine       */
  int newMax;
  /* Test if the memory is already allocated (as it should be)         */
  if (!isAllocated_calcpWrapper ||
      evalpnts_calcpWrapper==NULL ||
      directions_calcpWrapper==NULL ||
      fss_calcpWrapper==NULL ||
      fds_calcpWrapper==NULL ||
      fess_calcpWrapper==NULL ||
      feds_calcpWrapper==NULL){
    fprintf(stderr,"%s: \n"
	    "    Memory should be re-allocated, but isn't allocated.\n"
	    "    (this may indicate a bug. Better fix it!)\n"
	    "    In the meantime I'll stop the program\n",
	    fctName);
    exit(1);
  }
  /* Set new memory size:                                              */
  newMax = MAX(2*nEvalsMax_calcpWrapper,nEvals);
  /* Free old memory:                                                  */
  calcpWrapperFree();
  /* Allocate new memory:                                              */
  calcpWrapperCalloc(newMax);
  return;
} /* End of routine calcpWrapperRealloc */
/* 
* --------------------------------------------------------------------
* Free the memory used by calcpWrapper.
* -------------------------------------------------------------------- */
static
void calcpWrapperFree()
{
  char fctName[]="calcpWrapperFree";  /* Name of this routine          */
  int i,writeWarning;
  int verbose = 0;
  /* If verbose and the memory is not allocated, then write a warning  */
  writeWarning=!isAllocated_calcpWrapper;
  /* Free memory. If an array is not allocated, then toggle the 
   * warning to be written                                             */
  /* When deallocating - then let the pointer be reset to NULL, so that
   * other routines will know that it is deallocated!                  */
  if (panel_calcpWrapper!=NULL){
    for (i=0;i<4;i++){
      FREE(panel_calcpWrapper[i]);
    }
    FREE(panel_calcpWrapper);
  }
  else writeWarning=1;
  /* If both the "long" arrays are allocated together, 
   * then deallocate them in one loop                                  */
  if(evalpnts_calcpWrapper!=NULL &&
     directions_calcpWrapper!=NULL){
    for (i=0;i<nEvalsMax_calcpWrapper;i++){
      FREE(evalpnts_calcpWrapper[i]);
      FREE(directions_calcpWrapper[i]);   
    }
    FREE(evalpnts_calcpWrapper);
    FREE(directions_calcpWrapper);
  }
  else {
    writeWarning=1;
    if (evalpnts_calcpWrapper!=NULL){
      for (i=0;i<nEvalsMax_calcpWrapper;i++){
	FREE(evalpnts_calcpWrapper[i]);
      }
      FREE(evalpnts_calcpWrapper);
    }
    if (directions_calcpWrapper!=NULL){
      for (i=0;i<nEvalsMax_calcpWrapper;i++){
	FREE(directions_calcpWrapper[i]);   
      }
      FREE(directions_calcpWrapper);
    }
  }
#define LOCALFREE(ARRAY) if (ARRAY != NULL){ \
                           FREE(ARRAY)}      \
                         else writeWarning=1;
  LOCALFREE(fss_calcpWrapper);
  LOCALFREE(fds_calcpWrapper);
  LOCALFREE(fess_calcpWrapper);
  LOCALFREE(feds_calcpWrapper);
#undef LOCALFREE

  if(verbose && writeWarning){
    fprintf(stderr,"%s: WARNING!\n"
	    "    Some (or all) of the memory could not be deallocated\n"
	    "    since it was not allocated to begin with!\n"
	    "    This warning may occur after the preconditioning step\n"
	    "    if no additional influence coefficients were calculated.\n"
	    "    (In that case you may disregard this warning.)\n",
	    fctName);
  }
  /* In debug mode make sure that all pointers now point to NULL:      */
  assert(panel_calcpWrapper==NULL);
  assert(evalpnts_calcpWrapper==NULL);
  assert(directions_calcpWrapper==NULL);
  assert(fss_calcpWrapper==NULL);
  assert(fds_calcpWrapper==NULL);
  assert(fess_calcpWrapper==NULL);
  assert(feds_calcpWrapper==NULL);

  /* Memory is now deallocated. Set flag to show this:                 */
  isAllocated_calcpWrapper=0;
  return;
} /* End of routine calcpWrapperFree */


/*
* ==================================================================== 
* Add a constant to the diagonal of a matrix.
* 
* Eventually, the constant should be added only to some of the diagonal 
* elements. In particular 2*pi should be added only to those entries 
* which correspond to an integration of Gn. Whether these entries 
* correspond to the panels with Neumann conditions or to the panels 
* with Diriclet conditions depend on the context, i.e. whether the 
* matrix in question should be used for building a right-hand-side 
* or for solving the problem.
*
* NOTE: Depending on the numbering scheme used, the entries which needs 
* to be modified may or may not lie on the diagonal. Only the "elements" 
* will know this.
* ==================================================================== */
static
void addToDiagonal(void *inMat, int matNum, double addConst, char *addType)
{
  char fctName[] = "addToDiagonal";
  sparseMat *mat;
  int i, j, n, startRow, endRow, *rowIndex, *colIndex;
  double *vals;


printf(" %s: WARNING: Returning without adding to diagonal!\n",fctName);
return;

  mat = (sparseMat*) inMat;

  n = mat->n;
  rowIndex = mat->rowIndex;
  colIndex = mat->colIndex;
  vals = mat->matVals[matNum];

printf(" %s: Adding %g\n",fctName,addConst);
  if (strcmp(addType,"all")==0){
    for (j=0; j<n; j++){
      startRow = colIndex[j];
      endRow   = colIndex[j+1]; 
      for(i = startRow;  i < endRow; i++) {
	if(rowIndex[i]==j){
	  /* This is the diagonal entry.                               */
	  vals[i] += addConst;
	  break; /* Skip to next column                             */
	}
      }
    }
  }
  else {
    /* This is an unknown type. Eventually, "all", "neumann" 
     * and "diriclet" should be implemented.                           */
    fprintf(stderr,"%s: ERROR. Unknown type for diagonal entry: \"%s\"\n"
	    "    Please use \"all\"\n"
	    "    Exiting\n",
	    fctName,addType);
    exit(1);
  }

  return;
} /* End of routine addToDiagonal */


/*
* ==================================================================== 
* Estimate the amount of memory needed to allocate the direct-part
* matrices: the direct part itself or the preconditioning matrix.
* ==================================================================== */ 
double directPartMemoryEstimate(void *grid, 
				int nCols, int nRows, int numMats,
				char *Type)
{
  char fctName[] = "directPartMemoryEstimate";
  int verbose = 0;
  int *numNonzeroInColsDirect;  /* Number of nonzero entries in 
				 * each column of the (sparse) matrix  */
  int i, nEntries, nLinks,returnZeroIfAssociated;
  double memUse;

  assert(strcmp(Type,"direct")==0 || strcmp(Type,"precond")==0);

  /* At the moment the numbers of sources and evals must match to build
   * a square system. Later that may not be so (one could consider e.g. 
   * building an overdetermined system...)                             */
  assert(nCols==nRows);

  /* Allocate temporary memory for counting the number of nonzero entries 
   * in each column for both the direct matrix and the preconditioner
   * matrix                                                            */
  numNonzeroInColsDirect  = (int *) calloc(nCols,sizeof(int));
  
  /* Calculate the number of nonzero entries in each row of the 
   * sparse direct matrix. 
   * This call modifies numNonzeroInColsDirect, filling in the number 
   * of nonzero entries needed in each column.                         */
  if (verbose)
    printf("%s: Finding number of nonzero entries \n"
	   "       in each column of the %s matrix\n",fctName,Type);
  sourcesAndNeighboursOperate(grid,
			      Type, 
			      0, /* No, don't use multiple sources! */
			      0, /* No limit */
			      0, /* No thanks, keep your weird data */
			      calcNumNonzeroDirectEntries,
			      (void*) numNonzeroInColsDirect);

  /* Add the numbers to find the total number of non-zero elements     */
  for (nEntries=0,i=0; i<nCols; i++) nEntries += numNonzeroInColsDirect[i];

  /* Possibly add memory if neighbour lists have not been set up       */
  nLinks = countPossibleNeighbourLinks(grid,
				       returnZeroIfAssociated=1,
				       Type);
  nEntries += nLinks;

  /* Find the needed memory: query a sparse-matrix routine             */
  memUse = spMatMemoryEstimate(nCols, nEntries, numMats);

  /* Free temporary memory                                             */
  FREE(numNonzeroInColsDirect);

  return memUse;
} /* End of routine directPartMemoryEstimate */
