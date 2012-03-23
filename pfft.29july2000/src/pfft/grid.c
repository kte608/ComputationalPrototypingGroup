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
* FILE: grid.c
*
*  This file contains routines that operate on 
*  or through the grid.
*
*  Note: Since all allocations so far have been made using 
*  "calloc", it is possible to check each pointer in a 
*  previously allocated structure against NULL using "assert"
*  before allocation of the pointer in the structure to verify 
*  that the pointer has not yet been allocated.
*/

# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <assert.h>
# include <float.h>

# include "Global.h"
# include "Grid.h" 
# include "Nfft.h" 

# include "GridBundles.h"

/* ======= Prototypes for static routines ======= */
static 
void setupGridElements(void *sources,int nSources,
		       void *evals,int nEvals,
		       struct grid *grid);
static
void setupGridSizeNomN(struct grid *grid,int nNomPow,int iDir,
		       point minCoordinates,point maxCoordinates,
		       double extraLayers, int verbose);
static
void setupGridSizeNi(struct grid *grid,int nNomPow,int iDir,
		     point minCoordinates, point maxCoordinates,
		     double extraLayers, int verbose);
static
int setupGridSizeNNN(struct grid *grid,
		     int nNomPow0,int nNomPow1,int nNomPow2,
		     point minCoordinates, point maxCoordinates,
		     double extraLayers, int verbose);

static
void setupGridElementNeighbours(struct grid *grid, 
				struct gridStencil *stencil);

static 
void allocateGridElements(void *sources,int nSources,
			  void *evals,int nEvals,
			  struct grid *grid);
static
void associateGridElements(void *sources,
			   void *evals,
			   struct grid *grid);
static
void findGridElementBoundingSpheres(struct grid *grid);

static
void setupBasicGridValues(struct grid *grid);

static
void setupGridStencils(struct grid *grid, int refreshGridSizeOnly);

static
void setupGridDimensions(struct grid *grid,
			 int gridMethod,
			 int p2x, int p2y, int p2z,
			 int iDir,
			 int verbose);

static
void setupGridPointNeighbourStencil(struct grid *grid,
				    struct gridStencil *stencil,
				    int stencilType,
				    double stencilSize,
				    double separationDistanceFactor,
				    int refreshGridSizeOnly);
static 
int findStencilIndices(int stencilType,
		       double stencilSize,
		       int storeInArray,
		       int refreshGridSizeOnly,
		       double dx,
		       double dy,
		       double dz,
		       struct gridStencil *stencil);

static
void allocateGridPoints(struct grid *grid);

static
void allocateGridAssociatedElements(struct grid *grid);

static
void associateElementsWithGridPoints(struct grid *grid);

static
void associateOneTypeElementsWithGridPoints(int elementType,
					    struct grid *grid);
static
int sphereIsFullyInsideGrid(point center, 
			    double radius,
			    struct grid *grid);
static
void closestGridPoint(int *i,int *j,int *k,
		      point x,
		      struct grid *grid);
static
void gridPointPosition(int i,int j,int k,point x,
		      struct grid *grid);
static
int gridPointsAreNeighbours(int i1,int j1,int k1, /* point 1 */
			    int i2,int j2,int k2, /* point 2 */
			    struct gridStencil *stencil);
static 
void addElementsToBothTempNeighLists(
	  int iSource,
	  int iEval,
	  struct neighbourLinkedList **sourceLinkedList,
	  struct neighbourLinkedList **evalLinkedList);
static
int findGridElementNeighbours(int ip,int jp,int kp,
			      int elementType,
			      struct grid *grid,
			      struct gridStencil *stencil,
			      int createLinks);
static
void sortGridElementNeighbourLists(int elementType,
				   struct grid *grid,
				   struct gridStencil *stencil);

static 
void makeGridNeighbourIEvalList(int ip, int jp, int kp,
				int nNeighbours,
				int *iNeighbours,
				int *jNeighbours,
				int *kNeighbours,
				int **iEvalsIn,
				int **iEvalsDirectLocationIn,
				int *nEvals, int *nMaxEvals,
				struct grid *grid,
				int verbose);

static 
void polyInterp(int Type, 
		double dx, double dy, double dz, 
		int refreshGridSizeOnly, 
		struct gridStencil *stencil);

static 
double gridElementRadius(struct grid *grid);

static
void findBasicGridValues(struct grid *grid,
			 int nCols,int nRows, int numMats);
static
double gridElementListsMemoryEstimate(struct grid *grid);
static
double gridMemoryEstimate(struct grid *grid, int countLists,
			  int verbose);
static
double directMemoryEstimate(struct grid *grid,
			    int nCols,int nRows, int numMats);
static
void freeStencilForGridResize(struct gridStencil *stencil, 
			      int *nGridElements, int sourceAndEvalDiffer);

static
void freeGridForResize(struct grid *grid);

static
int countGridElementAssociations(struct grid *grid,
				 int *maxAssocElements);

static
int stencilRange(struct gridStencil *stencil);

/*  Prototypes for functions that has to do with the 
*    direct part, but does need to know how a sparse
*    matrix is defined:
*/
static
int gridPointExistsOld(int ip, int jp,int kp,struct grid *grid);

/*  Prototypes for test functions (in grid.testfunctions.c)   */
static
void testTableElementNeighbours(int elementType,
				struct grid *grid, 
				struct gridStencil *stencil);
static
void tableGridElements(int elementType,struct grid *grid);


static 
void setupNearFieldGrid2grid(struct grid *grid, int refreshGridSizeOnly);

static 
void setupNearFieldGrid2gridUnion(struct grid *grid, int refreshGridSizeOnly);

/* 
* ==================================================================== 
* This is the main routine for setting up everything which has to
* do with the grid.
* ==================================================================== */
void *setupGrid(void *sources,int nSources,
		void *evals,  int nEvals,
		int nCols, int nRows, int numMats)
{
  char fctName[] = "setupGrid";
  struct grid *grid = NULL;
  int verbose = 1;
  double elementRadius;

  /* Allocate storage for the grid */
  assert(grid == NULL);
  grid = calloc(1,sizeof(struct grid));

  /* Setup gridElements based on calls to element-routines.
   * This basically queries for the center and radius of the 
   * bounding sphere of each element. */
  if (verbose) 
    printf("%s: CALLING setupGridElements\n",fctName);
  setupGridElements(sources,nSources,evals,nEvals,grid);

  /* TEST OUTPUT */
  if (verbose > 2) { 
    int i;
    printf("Source bounding spheres\n");
    for (i=0; i<grid->nGridElements[0]; i++){
      printf("(%f,%f,%f),   %f\n",
	     grid->gridElements[0][i]->center[0],
	     grid->gridElements[0][i]->center[1],
	     grid->gridElements[0][i]->center[2],
	     grid->gridElements[0][i]->radius);
    }
    printf("Eval bounding spheres\n");
    for (i=0; i<grid->nGridElements[1]; i++){
      printf("(%f,%f,%f),   %f\n",
	     grid->gridElements[1][i]->center[0],
	     grid->gridElements[1][i]->center[1],
	     grid->gridElements[1][i]->center[2],
	     grid->gridElements[1][i]->radius);
    }
    tableGridElements(0,grid);
  }

  /* Get some statistics on the grid element sizes */
  elementRadius = gridElementRadius(grid);
  printf("%s: Average element radius is %12.4e\n",fctName,elementRadius);

  /* Find basic grid values to minimize to total memory consumption 
   * of the program. 
   * Routine also sets up element associations, finds neighbours etc.  */   
  findBasicGridValues(grid, nCols,nRows,numMats);
  return (void *) grid;
  printf("%s: CONTINUE FROM HERE!\n",fctName); 
  exit(1);

  /* The remaining part of this routine can be deleted when we are 
   * happy with the above solution */
  /* Setup basic grid values. DX, Extend, Stencil for neighbours...
   * These may require knowledge about the spatial distribution 
   * of the elements. */
  if (verbose) 
    printf("%s: CALLING setupBasicGridValues\n",fctName);
  setupBasicGridValues(grid);

  /* Associate grid elements with grid points. 
   * ("Put elements onto the grid") */
  if (verbose) 
    printf("%s: CALLING associateElementsWithGridPoints\n",fctName);
  associateElementsWithGridPoints(grid);

  /* Locate neighbouring grid elements of different types */
  if (verbose) 
    printf("%s: Setting up neighbours (direct part)\n",fctName);
  setupGridElementNeighbours(grid,grid->directStencil);
  if (verbose) 
    printf("%s: Setting up neighbours (precondition part)\n",fctName);
  setupGridElementNeighbours(grid,grid->precondStencil);

  if (verbose) 
    printf("%s: Returning pointer to grid\n",fctName);
  return (void *) grid;
} /* End of setupGrid */


/* 
* ==================================================================== 
* This routine calls a target function for each grid point element, 
* and passes along a pointer to the "real" element (not the grid version 
* of the same) as well as the (i,j,k) index of the associated grid point.
* ==================================================================== */
void elementAndGridOperate(void *gridInput,
			   char *Type,
			   void (*targetFunction)(void *gridInput,
						  void *element,
						  char *Type,
						  int i, int j, int k,
						  void *tData),
			   void *targetData)
{
  char fctName[]="elementAndGridOperate";
  int verbose=0;
  struct grid *grid;
  gridPoint ****gridPoints;
  int ip, jp, kp, isrc;
  int selectElement, nElements, *elements;

  if (strcmp(Type,"source")==0){
    selectElement = 0;
  }
  else if (strcmp(Type,"eval")==0){
    selectElement = 1;
  }
  else {
    /* This is everything else - which is unknown=bad. 
     * Report error and exit                                           */
    fprintf(stderr,"%s: ERROR. Unknown type for source and grid stencils %s\n"
	    "    Please use \"source\" or \"eval\"\n"
	    "    Exiting\n",
	    fctName,Type);
    exit(1);
  }

  if (verbose)
    fprintf(stdout,"%s: called - using %s type stencil\n",
	    fctName,Type);

  /* Cast grid to correct type */
  grid = (struct grid*) gridInput;
  gridPoints = grid->points;

  /* For every grid point                                              */
  for (ip=0; ip<grid->nx; ip++){
    for (jp=0; jp<grid->ny; jp++){
      for (kp=0; kp<grid->nz; kp++){
	/* Skip point if there are no associated elements or no 
	 * associated sources (necessary to check for "elements" first)*/
	if (gridPoints[ip][jp][kp]==NULL)
	  continue; /* Skip to next grid point */

	nElements = gridPoints[ip][jp][kp]->nAssociatedElements[selectElement];
	elements  = gridPoints[ip][jp][kp]->associatedElements[selectElement];

	/* For each associated element (source or eval!)               */
	for(isrc=0; 
	    isrc<nElements;
	    isrc++){
	  /* Call the target function                                  */
	  targetFunction(
	    gridInput, 
            (void *)grid->gridElements[selectElement][elements[isrc]]->element,
	    Type,
	    ip, jp, kp, targetData);
	}
      }
    }
  }
  return;
} /* End of routine elementAndGridOperate */

/* 
* ==================================================================== 
* Allocates grid and kernel data.
* ==================================================================== */
void grid2gridSetup(void *gridIn)
{
  struct grid *grid;

  /* Cast input to correct type */
  grid = (struct grid*) gridIn;

  /* Allocate space for the kernel. */
  grid->kernel = createFFT(grid->nx, grid->ny, grid->nz, "kernel");

  /* Fill in kernel. */
  fillKernel(grid->kernel, grid->dx, grid->dy, grid->dz);

  /* Allocate space for the data. */
  grid->griddata = createFFT(grid->nx, grid->ny, grid->nz, "data");

  return;
} /* End of routine grid2gridSetup */

/* 
* ==================================================================== 
* Deallocates grid and kernel data.
* ==================================================================== */
void grid2gridFree(void *gridIn)
{
  struct grid *grid;

  /* Cast input to correct type */
  grid = (struct grid*) gridIn;

  /* Deallocate space for the kernel. */
  freeFFT(grid->kernel);

  /* Deallocate space for the data. */
  freeFFT(grid->griddata);

  return;
} /* End of routine grid2gridFree */

/* 
* ==================================================================== 
* Precorrects the direct interaction.
* Extract the "kernel" pointer from "grid" and pass it to a routine 
* which knows what to do with it. Pass through also the matrices for 
* the direct part, the projection part and the interpolation part.
*
* dmat: Direct part matrix (to be precorrected)
* pmat: Projection matrix (used in precorrection process)
* imat: Interpolation matrix (used in precorrection process)
* ==================================================================== */
void precorrect(void *gridIn, void *dmat, void *pmat, void *imat)
{
  struct grid *grid;
  /* Cast input data to correct type */
  grid = (struct grid*) gridIn;
  /* Call routine to do the actual work */
  precorrectMat(grid->kernel,grid->griddata, dmat, pmat, imat);
  return;
}


/* 
* ==================================================================== 
* Finishes the grid to grid interaction by fourier transforming the
* kernel.
* ==================================================================== */
void grid2gridFinish(void *gridIn)
{
  struct grid *grid;

  grid = (struct grid*) gridIn;
  rfft3(grid->kernel);
  return;
} /* End of routine grid2gridFinish */

/*
* ==================================================================== 
* Make the grid-to-grid calculations, i.e. go from (e.g.) grid sources 
* to grid-potentials. 
* Just extract the needed data and pass them on to a routine that 
* knows what to do.
* ==================================================================== */ 
void grid2gridCalculate(void *gridIn)
{
  struct grid *grid;
  int dataType = -1;

  /* Cast input pointer to correct type */
  grid = (struct grid*) gridIn;
  
  /* Call a routine to do the grid-to-grid calculations                */
  convolve(grid->griddata,grid->kernel,dataType);

  return;
} /* End of routine grid2gridCalculate */

/* 
* ==================================================================== 
* For each source element find all the eval neighbours (including 
* evaluation elements reached through neighbouring grid points, i.e.
* using the "stencil"). Then pass the source and the list of neighbours
* to a routine which operates using the information. Pass along also 
* some data to be operated on.
* If (useMultipleSources) is true, then pass all the sources associated 
* with a particular grid point together with the union of their 
* neighbours in one go.
* If (maxNeighbours) is true (different from zero), then this variable
* determines how many element neighbours that may be passed along with 
* the source(s).
* addData will make the routine bundle some grid-specific data 
* with the data to the target function. (not implemented yet).
* ==================================================================== */
void sourcesAndNeighboursOperate(void *gridInput,
				 char *neighboursType,
				 int useMultipleSources,
				 int maxNeighbours,
				 int addData,
				 void (*targetFunction)(void **sources, 
							void **evals, 
							int nSources,
							int nTotNeighbours,
							void *tData),
				 void *targetData)
{
  char fctName[]="sourcesAndNeighboursOperate";
  int verbose = 0;
  
  struct gridBundle_SANO_1 bundle1[1]; /* Used to pass additional data 
				        * to a routine that requires it*/

  int ip,jp,kp,        /* Grid point indices (counters)                */
    nSources,          /* Number of sources associated with grid point 
			*  = gridPoints[ip][jp][kp]->nAssociatedElements[0] */
    isrc, iSource,     /* Element counters                             */
    ievl, iEval;

  struct grid *grid;
  struct gridStencil *stencil;

  /* The following are short-hand to the grid     */
  gridPoint ****gridPoints;
  struct gridElement **gridElements[2];
  struct gridElementNeighbours **stencilElements[2];
  int nNeighbours,*iNeighbours,*jNeighbours,*kNeighbours;

  int nNeighboursGrid, /* Number of neighbouring elements reached 
			* through the grid (stencil)                   */
    nNeighboursSrc,    /* Number of neighbouring elements reached
			* through the neighbour list of the source     */
    nNeighboursSrcKeep,/* Copy of nNeighboursSrc                       */
    nTotNeighbours;    /* Total number of neighbours                   */

  int nSrcMax,         /* Maximum nuber of sources at a grid point     */
    nNeighMax,         /* Maximum number of neighbours                 */
    nNeighMaxKeep,     /* Copy of nNeighMax                            */
    increasedBefore,   /* 1 if the storage for evals have been 
			* increased beyond the hardcoded initial limit
			* (nNeighMax). 0 otherwise.                    */
    maxEvals,          /* Size of array for generic pointers to 
			* neighbouring elements                        */
    nEvalDist;         /* Size of array for distances                  */
  void **sources,      /* For each grid point associated sources...    */
    **evals,           /*  ...a list of the neighbouring evals         */
    **sourcesPassed;   /* Pointer to the sources that are actually 
			* passed to the target routine                 */
  int *iEvals;         /* A list of the indices of the found evals.
			* This list will (if needed) be usd for the 
			* sorting process                              */
  int *iEvalsDirectLocation; /* For each eval, the neighbouring relation 
			      * numbered as in the direct stencil will
			      * be stored. This info may be bunddled to 
			      * the target routine.                    */
  int *iuseEvals;      /* Pointer to the eval-list which is actually 
			* sent to the target routine */
  double *evalDist;    /* Distance from sources to evals (normalized)  */
  int *ievalDist;      /* Array for storing temporarily the eval list 
			* for sorting. */
  struct neighbourLinkedList *neighbourLink; /* Used to access grid evals */
  int *iNeighboursSrc, /* List of indices of grid evals that are 
			* neigbours to a source (not "through the grid") */
    iNeighboursSrcMax; /* Size of array iNeighboursSrc                 */
  double veryLarge=DBL_MAX,    /* Large positive number                */
    dist;                     /* A squared distance                    */
  int allSourcesAreDone; /* Flag used in a loop over the sources.
			  * All the sources may be done after ONE time 
			  * through the loop, thus this construction.  */
  int thisSource;      /* Counter used in loop through sources         */
  int nSourcesPassed;  /* Number of sources passed to the target routine */
  void *dataPassed;    /* Points to whatever data that are actually 
			* passed to the target routine */

  /* Cast grid to correct type */
  grid = (struct grid*) gridInput;

  if (verbose)
    fprintf(stdout,"%s: called - using %s type stencil\n",
	    fctName,neighboursType);
  /* Cast grid to correct type */
    grid = (struct grid*) gridInput;
  if (strcmp(neighboursType,"direct")==0){
    /* This uses the direct-part stencil                               */
    stencil = grid->directStencil;
  }
  else if (strcmp(neighboursType,"precond")==0){
    /* This uses the preconditioning-part stencil                      */
    stencil = grid->precondStencil;
  }
  else {
    /* This is everything else - which is unknown=bad. 
     * Report error and exit                                           */
    fprintf(stderr,"%s: ERROR. Unknown type for neighbour stencils %s\n"
	    "    Please use \"direct\" or \"precond\"\n"
	    "    Exiting\n",
	    fctName,neighboursType);
    _EXIT_;
  }

  /* Allocate initial storage for list of neighbour elements.
   * If too little storage is allocated here, then 
   * it will be reallocated later. */
  nNeighMax = 1000;
  iEvals    = (int*)    calloc(nNeighMax,sizeof(int));
  iEvalsDirectLocation = (int*) calloc(nNeighMax,sizeof(int));

  nEvalDist = 1000;
  evalDist  = (double*) calloc(nNeighMax,sizeof(double));
  ievalDist = (int*)    calloc(nNeighMax,sizeof(int));

  maxEvals  = maxNeighbours;
  if(maxNeighbours==0) maxEvals = 1000;
  evals     = (void**)  calloc(maxEvals,sizeof(void*));
  
  nSrcMax   = 1000;
  sources   = (void**) calloc(nSrcMax,sizeof(void*));

  iNeighboursSrcMax = 1000;
  iNeighboursSrc = (int*) calloc(iNeighboursSrcMax,sizeof(int));
  increasedBefore = 0;

  /* Grid "shorthands" */
  gridPoints = grid->points;
  gridElements[0] = grid->gridElements[0];
  gridElements[1] = grid->gridElements[1];
  stencilElements[0] = stencil->gridElementNeighbours[0];
  stencilElements[1] = stencil->gridElementNeighbours[1];
  nNeighbours = stencil->nStencilPoints;
  iNeighbours = stencil->ipNeighbours;
  jNeighbours = stencil->jpNeighbours;
  kNeighbours = stencil->kpNeighbours;

  /* For each grid point */
  for (ip=0; ip<grid->nx; ip++){
    for (jp=0; jp<grid->ny; jp++){
      for (kp=0; kp<grid->nz; kp++){
	/* Skip point if there are no associated elements 
	 * or no associated sources (necessary to check for 
	 * "elements" first) */
	if (gridPoints[ip][jp][kp]==NULL ||
	    gridPoints[ip][jp][kp]->nAssociatedElements[0] == 0)
	  continue; /* Skip to next grid point */
	/* Number of sources on this grid point:                       */
	nSources = gridPoints[ip][jp][kp]->nAssociatedElements[0];

	/* Set the right number of sources for the target routine      */
	if (useMultipleSources) nSourcesPassed = nSources;
	else nSourcesPassed  = 1;

	/* Make sure that the total number of sources does not exceed
	 * the maximum number of neighbours allowed. At least one eval 
	 * per source needs to be submitted to the target routine. 
	 * This is done to make sure that when the sources and the evals 
	 * "are the same", then at least the sources themselves are in 
	 * the neighbour list!                                         */
	if (maxNeighbours>0 && nSources>maxNeighbours){
	  fprintf(stderr,"%s: Error.\n"
		  "     Too many sources (%d) "
		  "associated with node [%d][%d][%d]\n"
		  "     or too few neighbours (%d) allowed.\n"
		  "     Increase the value of maxNeighbours "
		  "or refine the grid.\n"
		  "     Exiting.\n",
		  fctName,nSources,ip,jp,kp,maxNeighbours);
	  _EXIT_;
	}
	/* Store (the pointers to) the sources in an array. 
	 * Make sure the array is large enough                         */
	if (nSources>nSrcMax){
	  nSrcMax = MAX(2*nSrcMax,nSources);
	  FREE(sources);
	  sources   = (void**) calloc(nNeighMax,sizeof(void*));
	}
	/* Actually store the sources                                  */
	for(isrc=0;isrc<nSources;isrc++){
	  /* Global number of this (grid) source element               */
	  iSource = gridPoints[ip][jp][kp]->associatedElements[0][isrc];
	  /* Number of neighbours to this elements adds to the total   */
	  sources[isrc]=gridElements[0][iSource]->element;
	}
	isrc=0;
	iSource = gridPoints[ip][jp][kp]->associatedElements[0][isrc];
	assert((sources[isrc])==(gridElements[0][iSource]->element));

	/* Store indices of the neighbours (evals) in an array.
	 * Keep the present array size for the neighbours (evals) around, 
	 * so that we can see if the array size have increased 
	 * after call                                                  */
	nNeighMaxKeep = nNeighMax; 
	makeGridNeighbourIEvalList(ip, jp, kp,
				   nNeighbours,
				   iNeighbours,jNeighbours,kNeighbours,
				   &(iEvals),
				   &(iEvalsDirectLocation),
				   &nNeighboursGrid, &nNeighMax,
				   grid, 
				   0); /* The last argument is verboseness.
					* it is possible to use e.g. 0 or
					* !increasedBefore             */
	if (nNeighMax!=nNeighMaxKeep){
	  /* In this case the array size has increased. 
	   * Decrease the level of verboseness in later calls          */
	  increasedBefore = 1;
	}

	/* The following block may be done either as a loop over 
	 * the sources, or in one go - all sources at once.            */
	allSourcesAreDone=0;
	thisSource = 0;
	while(!allSourcesAreDone){
	  /* For each associated source element add all associated evals 
	   * to the list. Duplicate entries may appear here, so the list 
	   * will need sorting and possibly elimination of dublicates.
	   * In order to make the sort, these evals wil first be stored 
	   * using their gridElement number - rather than pointer to the 
	   * real element. After the sort, the elements will be stored 
	   * corectly in the main "evals" array.                       */
	  nNeighboursSrc = 0;
	  /* Count the number of neighbours (including dublicates)     */
	  if (useMultipleSources){
	    for(isrc=0;isrc<nSources;isrc++){
	      /* Global number of this (grid) source element           */
	      iSource = gridPoints[ip][jp][kp]->associatedElements[0][isrc];
	      /* Number of neighbours to this elements adds to the total*/
	      nNeighboursSrc += stencil
		->gridElementNeighbours[0][iSource]
                                ->nNeighbours;
	    }
	  }
	  else{ /* Only one source here */
	    iSource = gridPoints[ip][jp][kp]->associatedElements[0][thisSource];
	    nNeighboursSrc = stencil
	                       ->gridElementNeighbours[0][iSource]
	                        ->nNeighbours;
	  }
	  /* In many cases there won't be any neighbours in this list, 
	   * so only continue with this part of the process if the 
	   * list is not empty                                         */
	  if (!nNeighboursSrc){
	    nTotNeighbours=nNeighboursGrid;
	  }
	  else {
	    /* Check if enough memory have been allocated              */
	    if (nNeighboursSrc>iNeighboursSrcMax){
	      /* Allocate more memory. Make sure that we allocate a good
	       * deal more memory, so that this won't be needed often. */
	      FREE(iNeighboursSrc);
	      iNeighboursSrcMax = MAX(2*iNeighboursSrcMax,nNeighboursSrc);
	      /* Allocate new memory for neighbour list                */
	      iNeighboursSrc = (int*)calloc(iNeighboursSrcMax,sizeof(int));
	    }
	    /* Store indices of neighbours including dublicates in array*/
	    if (useMultipleSources){
	      nNeighboursSrc = 0;
	      for(isrc=0; 
		  isrc<gridPoints[ip][jp][kp]->nAssociatedElements[0];
		  isrc++){
		/* Global number of this (grid) source element             */
		iSource = 
		  gridPoints[ip][jp][kp]->associatedElements[0][isrc];
		/* Add these neighbours to the list                        */
		neighbourLink = 
		stencilElements[0][iSource]->neighbourLinkedList;
		for (ievl=0; 
		     ievl<stencilElements[0][iSource]->nNeighbours; 
		     ievl++){
		  /* Find the (grid) eval number                           */
		  iEval = neighbourLink->iNeighbour;
		  iNeighboursSrc[nNeighboursSrc+ievl] = iEval;
		  neighbourLink = neighbourLink->nextNeighbour;
		}
		/* Number of neighbours to this elements adds to the total */
		nNeighboursSrc += 
		  stencilElements[0][iSource]->nNeighbours;
	      }
	    }
	    else{ /* Only one source */
	      nNeighboursSrc = 0;
	      iSource = 
		gridPoints[ip][jp][kp]->associatedElements[0][thisSource];
	      /* Add these neighbours to the list                      */
	      neighbourLink = 
		stencilElements[0][iSource]->neighbourLinkedList;
	      for (ievl=0; 
		   ievl<stencilElements[0][iSource]->nNeighbours; 
		   ievl++){
		/* Find the (grid) eval number                         */
		iEval = neighbourLink->iNeighbour;
		iNeighboursSrc[nNeighboursSrc+ievl] = iEval;
		neighbourLink = neighbourLink->nextNeighbour;
	      }
	      nNeighboursSrc = stencilElements[0][iSource]->nNeighbours;
	    }

	    /* Sort the list. If there is only one source associated with 
	     * this grid point, then the sorting is not needed. In this 
	     * case the entries are already sorted, and there are no 
	     * dublicate entries. Even if there ARE multiple sources,
	     * only sort if we are using multiple sources (we could be 
	     * taking them one0-by-one).                               */
	    if (useMultipleSources && nSources>1) {
	      intSort(iNeighboursSrc,nNeighboursSrc);
	      /* Eliminate dublicate entries                           */
	      nNeighboursSrcKeep = nNeighboursSrc;
	      /* Always keep the first entry, so start from src no. 1  */
	      nNeighboursSrc = 1;
	      for(isrc=1; 
		  isrc< nNeighboursSrcKeep-1;
		  isrc++){
		/* If this and the previous element is not the same, 
		 * then keep (store) this element                      */
		if(iNeighboursSrc[isrc]!=iNeighboursSrc[isrc-1]){
		  iNeighboursSrc[nNeighboursSrc] = iNeighboursSrc[isrc];
		  nNeighboursSrc++;
		}
	      }
	    } /* End of sorting part                                   */
	    /* Add the obtained neighbours to the list of neighbours 
	     * that were found "through the grid".                     */
	    /* Test if enough memory is allocated:                     */
	    nTotNeighbours = nNeighboursGrid+nNeighboursSrc;
	    if (nTotNeighbours > nNeighMax){
	      /* Make sure the storage is increased significantly 
	       * (at least a factor of two), such that this increase won't
	       * be necessary too many times (hopefully not at all).   */
	      nNeighMax = 
		MAX(2*nNeighMax,nTotNeighbours);
	      iEvals = 
		(int*) realloc(iEvals,nNeighMax*sizeof(int));
	      iEvalsDirectLocation = 
		(int*) realloc(iEvalsDirectLocation,nNeighMax*sizeof(int));
	      increasedBefore=1;
	    }
	    /* Actually add to the list                                */
	    for (ievl=0;ievl<nNeighboursSrc;ievl++){
	      iEvals[nNeighboursGrid+ievl] = iNeighboursSrc[ievl];
	      /* This neighbour was not found through the grid, so:    */
	      iEvalsDirectLocation[nNeighboursGrid+ievl] = -1; 
	    }

	  } /* End if (nNeighboursSrc)                                 */

	  /* Reset which eval-list is used: */
	  iuseEvals = iEvals;

	  /* Possibly reduce the list of neighbours to the 
	   * given maximum - by sorting.                               */
	  if (maxNeighbours>0 && nTotNeighbours>maxNeighbours){
	    /* "Too many" neighbours have been found. 
	     * Remove those neighbours from the list that are farthest 
	     * away from the sources.
	     * Note that this sort will screw up the neighbours-through-
	     * grid part, i.e. the neighbours through the grid will not 
	     * (necessarily) be the the first in the sorted array. For 
	     * this reason, the sorting is done on a temporary array.  */
	    /* Check arrays for storing distances                      */
	    if(nTotNeighbours>nEvalDist){
	      FREE(evalDist);
	      FREE(ievalDist);
	      nEvalDist = MAX(nTotNeighbours,2*nEvalDist);
	      evalDist  = (double*) calloc(nNeighMax,sizeof(double));
	      ievalDist = (int*)    calloc(nNeighMax,sizeof(int));
	    }
	    /* Initially let all distances be "ridiculously high"      */
	    for (ievl=0;ievl<nTotNeighbours;ievl++){
	      evalDist[ievl]  = veryLarge;
	    }
	    /* Make a #sources * #totNeighbours loop to calculate 
	     * minimum (over sources) normalized distances for each 
	     * neighbouring eval.
	     * It may be bit cheaper to operate on squared distances, 
	     * however, then all radii need to be squared in the 
	     * normalization...                                        */
	    for (isrc=0;isrc<nSources;isrc++){
	      /* Global number of this (grid) source element           */
	      iSource = gridPoints[ip][jp][kp]->associatedElements[0][isrc];
	      for (ievl=0;ievl<nTotNeighbours;ievl++){
		/* The normalized distance is:                         */
		dist = 
		  /* Distance between centers                          */
		  (VDIST(grid->gridElements[0][iSource]->center,
			 grid->gridElements[1][iEvals[ievl]]->center)
		   )/
		  /* Normalized by the largest radii of the two elements */
		  MAX(grid->gridElements[0][iSource]->radius,
		      grid->gridElements[1][iEvals[ievl]]->radius);
		/* Possibly update the normalized distance for this eval */
		if (dist < evalDist[ievl]) evalDist[ievl] = dist;
	      } /* Next neighbour (ievl) */
	    }   /* Next source (isrc)    */
	    /* All distances have now been found.                      */

	    /* TEST BLOCK: 
	       {
	       printf("%s: normalized distances\n",fctName);
	       for (ievl=0;ievl<nTotNeighbours;ievl++){
	       printf("    %3d  %6d  %g \n",ievl,iEvals[ievl],evalDist[ievl]);
	       }
	       }*/
	    
	    /* Store the grid-eval numbers in a temporary array, 
	     * which may be screwed up by the sorting                  */
	    for (ievl=0;ievl<nTotNeighbours;ievl++){
	      ievalDist[ievl]  = iEvals[ievl];
	    }
	    iuseEvals = ievalDist;

	    /* Sort the entries and let the indices of the evals be 
	     * sorted as well                                          */
	    doubleSort2(evalDist,ievalDist,nTotNeighbours);
	    /* Note that iEvalsDirectLocation is NOT sorted yet 
	     * (it could be if needs come!). 
	     * Thus, if iEvalsDirectLocation is needed AND there is 
	     * a reduction of the evals, then  iEvalsDirectLocation
	     * SHOULD be sorted!                                      */
	    if (addData){
	      fprintf(stderr,"%s:, ERROR: "
		      "Adding data for direct-positions is not \n"
		      " implemented with sorting. Hardcoding may be required\n",
		      fctName);
	      _EXIT_;
	    }
	    /* TEST BLOCK: 
	       {
	       printf("%s: sorted normalized distances\n",fctName);
	       for (ievl=0;ievl<nTotNeighbours;ievl++){
	       printf("    %3d  %6d  %g \n",ievl,iEvals[ievl],evalDist[ievl]);
	       }
	    }*/

	    /* The sorting is now done. 
	     * The first elements in the array are the ones with the 
	     * smallest values of the normalized distance, i.e. the 
	     * ones "closest" to the sources.                          */
	    /* Limit the number of neighbours to the maximum           */
	    nTotNeighbours = maxNeighbours;
	    /* Note that there is no need to actually remove any 
	     * elements from the list. As long as the number is limited 
	     * and the sorting is done, we are fine.                   */
	  }/* Done reducing the list. */

	/* Make a list of pointers to elements to pass to the 
	 * target routine                                              */
	/* Make sure that enough memory is allocated for the evals     */
	if (nTotNeighbours>maxEvals){
	  /* This should NOT be necessary if there is a hard limit on 
	   * the number of eval neighbours                             */
	  assert(maxNeighbours==0);
	  /* Set new value (multiply by at least two) */
	  FREE(evals);
	  maxEvals = MAX(2*maxEvals,nTotNeighbours);
	  evals = (void**)  calloc(maxEvals,sizeof(void*));
	}
	for (ievl=0;ievl<nTotNeighbours;ievl++){
	  evals[ievl] = grid->gridElements[1][iuseEvals[ievl]]->element;
	}
	
	/* Set whether all sources should be passed on, or just one:   */
	if (useMultipleSources) {
	  /* Pass all sources on this grid point                       */
	  sourcesPassed = sources; 
	}
	else {
	  /* Only pass a single source at a time.                      */
	  sourcesPassed = &(sources[thisSource]);
	}

	if(addData==1){
	  /* Add additional data to the target routine. 
	   * Conform to bundle type 1 (struct gridBundle_SANO_1)       */
	  bundle1->origData     = targetData;
	  bundle1->nearFieldG2G = grid->nearGrid2grid;
	  bundle1->g2gIndex     = iEvalsDirectLocation;
	  bundle1->numG2G       = grid->numNearGrid2grid;
	  bundle1->G2Gn         = grid->nNearGrid2grid;
	  bundle1->G2Gm         = grid->mNearGrid2grid;
	  bundle1->nearFieldUnionG2G        = grid->nearGrid2gridUnion;
	  bundle1->nRowsNearFieldUnion      = grid->nRowsNearGrid2gridUnion;
	  bundle1->nearFieldUnionG2Gindices = grid->g2gUnionRowNumbers;
	  bundle1->unionWork                = grid->nearGrid2gridUnionWork;
	  /* Cast to generic pointer type */
	  dataPassed = (void*) bundle1;
	}
	else 
	  dataPassed = targetData;

	/* Pass the sources, the neighbours, the number of each, and 
	 * the data to a routine for further handling                  */
	targetFunction(sourcesPassed,evals,
		       nSourcesPassed,nTotNeighbours,
		       dataPassed);

	/* Test if we should continue this loop: */
	/* Maybe we are doing all the sources at once: */
	if (useMultipleSources)  allSourcesAreDone=1;
	else{
	  thisSource++; /* Next source */
	  if (thisSource==nSources) allSourcesAreDone=1;
	}

	} /* Next source on this point */

      } /* Next grid point */
    }
  }

  /* Free up the memory that has been allocated                        */
  FREE(iEvals);
  FREE(iEvalsDirectLocation);
  FREE(evalDist);
  FREE(ievalDist);
  FREE(evals);
  FREE(sources);
  FREE(iNeighboursSrc);
  return;
} /* End of sourceAndNeighboursOperate */


/* 
* ================================================================= 
* Allocate storage for the grid elements and other stuff on the 
* grid which depend only on the number of elements. 
* Grid elements contain the information about the elements which the 
* grid needs. Most importantly this is center and radius of a bounding 
* sphere to associate elements with grid points and find neighbours, 
* but also a pointer associaten to the original element is needed for 
* interaction with element-specific routines.
* ================================================================= */
static
void setupGridElements(void *sources,int nSources,
		       void *evals,int nEvals,
		       struct grid *grid)
{
  /* Allocate memory.
   * This also sets the number of grid elements, such that 
   * nSources and nEvals don't have to be passed along at 
   * later calls. */
  allocateGridElements(sources,nSources,evals,nEvals,grid);

  /* Make association with source elements and evaluation 
   * elements. */
  associateGridElements(sources,evals,grid);
  
  /* Call element routines to set up (tight) bounding spheres
   * around all elements. 
   * Also keep track of the extent of the elements in physical space.  */
  findGridElementBoundingSpheres(grid);

  /* Allocate for grid associations to elements */
  allocateGridAssociatedElements(grid);

  return;
}

/*
* ==================================================================== 
* Setup of the grid size, dimensions and number of points in each 
* cardinal direction. 
* Definition of what a neighbouring grid point is (the stencil).
* Allocation of memory needed to associate elements with grid 
* points.
* ==================================================================== */
static
void setupBasicGridValues(struct grid *grid)
{
  char fctName[] = "setupBasicGridValues";
  int refreshGridSizeOnly;

  /* Setup grid size and dimensions */
  printf("%s: WARNING: Using hard-coded grid values!\n",fctName);
  setupGridDimensions(grid,1,  8,0,0,  0, 0);

  /* Setup grid stencils                                               */
  setupGridStencils(grid, refreshGridSizeOnly=0);

  /* Allocate memory for the association of elements with 
   * grid points. */
  allocateGridPoints(grid);

  return;
} /* End of routine setupBasicGridValues                               */

/* 
* ==================================================================== 
* Allocate and setup the grid stencils. Any part of the stencils  
* which depend on the actual grid size will not be set up here.
*
* stencilTypes:
*  1: All points in a cube of side length 2*stencilSize*delta 
*     centered on the "central point"(delta being the grid size).
*  2: All points in a sphere of radius 2*stencilSize*delta 
*     centered on the "central point"(delta being the grid size).
* 10: Eight points on the corners of a cube of size length closest 
*     to  2*stencilSize*delta.
* ==================================================================== */
static
void setupGridStencils(struct grid *grid, int refreshGridSizeOnly)
{
  char fctName[] = "setupGridStencils";
  int verbose = 1;
  static int calledBefore=0;
  int npts;
  /* Information for the direct-part stencil */
  int directStencilType;
  double directStencilSize, 
    directSeparationDistanceFactor; 
  /* Information for the precondition stencil */
  int preconditionStencilType;
  double preconditionStencilSize, 
    preconditionSeparationDistanceFactor;
  /* Information for the interpolation/projection stencil */
  int projectionStencilType;
  double projectionStencilSize,
    projectionSeparationDistanceFactor;

  if((!calledBefore)&&verbose){
    /* stencilType and stencilSize sould not be hardcoded, but rather 
     * supplied by the user. Issue a warning about this!               */
    fprintf(stdout,"%s: WARNING: Using hardcoded values \n"
	    "                for grid neighbours (direct+precond)!\n",
	    fctName);
  }

  /* Choose interpolation/projection scheme (stencil) */
  /* Number of points in interpolation/projection stencil. */
  npts = 27;  /* <<- Possibly change hardcoding here! */
  projectionSeparationDistanceFactor=0.0;  /* Not used. Keep zero. */
  switch (npts){
  case (1):
    /* One-point interpolation - interpolate to/project from nearest 
     * point only. Very cheap and rather inaccurate. Have uses in some 
     * cases where low-accuracy integrated measures are needed. Local 
     * quantities may be very inaccurate when this scheme is used.     */
    projectionStencilType=1; /* or 0 - it won't matter here. */
    projectionStencilSize=0.001;
    /* Possible settings of direct size and precondition size to go 
     * with this projection/interpolation scheme:                      */
    directStencilType=1;
    directStencilSize=3.0;
    directSeparationDistanceFactor= 1.0;
    preconditionStencilType=1;
    preconditionStencilSize=1.0;
    preconditionSeparationDistanceFactor=0.01;
    break;
  case (7):
    /* Seven-point projection/interpolation. 
     * Points at the origin and on the main coordinate axes at 
     * distance delta (grid size).  
     * Rather cheap scheme - linear accurately interpolation functions
     * (not extremely accurate).                                       */
    projectionStencilType=2; 
    projectionStencilSize=1.001;
    /* Possible settings of direct size and precondition size to go 
     * with this projection/interpolation scheme:                      */
    directStencilType=2;
    directStencilSize=4.0;
    directSeparationDistanceFactor= 1.5;
    preconditionStencilType=2;
    preconditionStencilSize=2.0;
    preconditionSeparationDistanceFactor=0.01;
    break;
  case (8):
    /* Eight-point projection/interpolation scheme (stencil). 
     * Points at the corners of a cube with side length 2*delta.
     * Rather cheap scheme - linear accurately interpolation functions
     * (not extremely accurate).                                       */
    projectionStencilType=10; 
    projectionStencilSize=1.001;
    /* Possible settings of direct size and precondition size to go 
     * with this projection/interpolation scheme:                      */
    directStencilType=2;
    directStencilSize=4.0;
    directSeparationDistanceFactor= 1.5;
    preconditionStencilType=2;
    preconditionStencilSize=2.0;
    preconditionSeparationDistanceFactor=0.01;
    break;
  case (27):
    /* 27-point projection/interpolation scheme (stencil). 
     * Points in a 3x3x3 cube configuration with spacing delta.
     * Quadratic accurately interpolation functions (fairly accurate). */
    projectionStencilType=1; 
    projectionStencilSize=1.001;
    /* Possible settings of direct size and precondition size to go 
     * with this projection/interpolation scheme:                      */
    directStencilType=1;
    directStencilSize=4.0; 
    directSeparationDistanceFactor= 2.; 
    preconditionStencilType=2;
    preconditionStencilSize=2.5;
    preconditionSeparationDistanceFactor=0.01;
    break;
  case (33):
    /* 33-point projection/interpolation scheme (stencil). 
     * All points in a sphere of radius 2*delta (+epsilon).
     * Cubic accurately interpolation functions (very accurate)
     * and empirically optained scaling coefficients. */
    projectionStencilType=2; 
    projectionStencilSize=2.001;
    /* Possible settings of direct size and precondition size to go 
     * with this projection/interpolation scheme:                      */
    directStencilType=2;
    directStencilSize=5.8;
    directSeparationDistanceFactor= 2.0;
    preconditionStencilType=2;
    preconditionStencilSize=2.7;
    preconditionSeparationDistanceFactor=0.1;
    break;
  default:
    fprintf(stderr,"%s: Error! What kind of projection/interpolation\n "
	    "    stencil uses %d points?",
	    fctName,npts);
    _EXIT_;
    break;
  }

  if (verbose>1)
    printf("%s: Setting up direct-part grid neighbours: \n",fctName);
  /* Setup stencil for finding neighbouring grid points for the 
   * direct calculations                                               */
  if (!refreshGridSizeOnly){
    /* Allocate memory for the stencil                                 */
    assert(grid->directStencil == NULL);
    grid->directStencil = 
      (struct gridStencil*) calloc(1,sizeof(struct gridStencil));
  }
  /* Set values                                                    */
  setupGridPointNeighbourStencil(grid,
				 grid->directStencil,
				 directStencilType,
				 directStencilSize,
				 directSeparationDistanceFactor,
				 refreshGridSizeOnly);

  /* Test output block */
  if((!calledBefore)&&verbose){
    int i;
    printf("%s: %d neighbours in direct-part stencil",
	   fctName,grid->directStencil->nStencilPoints);
    if (verbose>1){
      printf(":\n");
      for (i=0;i<grid->directStencil->nStencilPoints;i++){
	printf("   %3d:  (%2d,%2d,%2d)\n",i,
	       grid->directStencil->ipNeighbours[i],
	       grid->directStencil->jpNeighbours[i],
	       grid->directStencil->kpNeighbours[i]);
      }
    }
    else printf("\n");
  }

  if (verbose>1)
    printf("%s: Setting up preconditioning-part grid neighbours: \n",
	   fctName);
  /* Setup stencil for finding neighbouring grid points for the 
   * preconditioning step.                                             */
  if (!refreshGridSizeOnly){
    /* Allocate memory for the stencil                                 */
    assert(grid->precondStencil == NULL);
    grid->precondStencil = 
      (struct gridStencil*) calloc(1,sizeof(struct gridStencil));
  }
  /* Set values                                                    */
  setupGridPointNeighbourStencil(grid,
				 grid->precondStencil,
				 preconditionStencilType,
				 preconditionStencilSize,
				 preconditionSeparationDistanceFactor,
				 refreshGridSizeOnly);

  /* Test output block */
  if((!calledBefore)&&verbose){
    int i;
    printf("%s: %d neighbours in preconditioning-part stencil",
	   fctName,grid->precondStencil->nStencilPoints);
    if (verbose>1){
      printf(":\n");
      for (i=0;i<grid->precondStencil->nStencilPoints;i++){
	printf("   %3d:  (%2d,%2d,%2d)\n",i,
	       grid->precondStencil->ipNeighbours[i],
	       grid->precondStencil->jpNeighbours[i],
	       grid->precondStencil->kpNeighbours[i]);
      }
    }
    else printf("\n");
  }

  if (!refreshGridSizeOnly){
    /* Allocate memory for the stencil                               */
    assert(grid->elementGridStencil == NULL);
    grid->elementGridStencil = 
      (struct gridStencil*) calloc(1,sizeof(struct gridStencil));
  }
  if (verbose>1)
    printf("%s: Setting up projection/interpolation-part grid neighbours: \n",
	   fctName);
  setupGridPointNeighbourStencil(grid,
				 grid->elementGridStencil,
				 projectionStencilType+100,
				 projectionStencilSize,
				 projectionSeparationDistanceFactor,
				 refreshGridSizeOnly);
  if((!calledBefore)&&verbose){
    int i;
    printf("%s: %d neighbours in projection/interpolation stencil",
	   fctName,grid->elementGridStencil->nStencilPoints);
    if (verbose>1){
      printf(":\n");
      for (i=0;i<grid->elementGridStencil->nStencilPoints;i++){
	printf("   %3d:  (%2d,%2d,%2d)\n",i,
	       grid->elementGridStencil->ipNeighbours[i],
	       grid->elementGridStencil->jpNeighbours[i],
	       grid->elementGridStencil->kpNeighbours[i]);
      }
    }
    else printf("\n");
  }

  calledBefore = 1; /* This routine has now been used before. 
		     * (Decrease verboseness.)                         */
  return;
} /* End of routine setupGridStencils */

/* 
* ==================================================================== 
* For all grid elements find those neighbouring elements of the 
* opposite type, which can't be reached through grid point 
* neighbours.
* ==================================================================== */
static
void setupGridElementNeighbours(struct grid *grid, 
				struct gridStencil *stencil)
{
  char fctName[] = "setupGridElementNeighbours";
  int ip,jp,kp,    /* Grid point counters */

    createLinks,   /* Passed to findGridElementNeighbours */

    elementType,   /* Type of element being worked on 
		    * 0 = Source,  1 = Eval
		    * The opposite type of elem. may be found 
		    * as !elementType                             */

    iae;           /* Counter for (associated) elements           */
  int verbose = 0;

  /* Make sure that all the neighbourLinkedList's are initiated 
   * to NULL */
  for (elementType = 0; elementType<2; elementType++){
    for (iae=0; iae < grid->nGridElements[elementType]; iae++)
      assert(stencil->gridElementNeighbours[elementType][iae]->
	        neighbourLinkedList == NULL );
  }


  if (verbose) 
    fprintf(stdout,"%s: Going over all grid points\n",fctName);
  /* For all grid points: */
  for (ip=0; ip<grid->nx; ip++){
    for (jp=0; jp<grid->ny; jp++){
      for (kp=0; kp<grid->nz; kp++){

	/* Continue processing this grid number 
	 * only if there are associated elements */
	if (grid->points[ip][jp][kp] == NULL) 
	  continue; /* Next grid point */

	/* Find neighbours to sources */
	findGridElementNeighbours(ip,jp,kp,
				  elementType = 0, 
				  grid,
				  stencil,
				  createLinks=1);

	if (grid->sourceAndEvalDiffer){
	/* Find neighbours to evals only if the sources and
	 * the evals are not exactly the same thing */
	  findGridElementNeighbours(ip,jp,kp,
				    elementType = 1,
				    grid,
				    stencil,
				    createLinks=1);
	}

      } /* Next grid point (kp) */
    }   /* Next grid point (jp) */
  }     /* Next grid point (ip) */

  /* TEST OUTPUT 
printf("Before sorting neighbours:\n");
testTableElementNeighbours(0,grid,stencil);
testTableElementNeighbours(1,grid,stencil); 
  */

  /* Sort the neighbour lists and eliminate dublicate entries */
  if (verbose) 
    fprintf(stdout,"%s: Sorting neighbour lists\n",fctName);
  sortGridElementNeighbourLists(elementType = 0,grid,stencil);
  if (grid->sourceAndEvalDiffer)
    sortGridElementNeighbourLists(elementType = 1,grid,stencil);

  /* TEST OUTPUT 
printf("After sorting neighbours:\n");
testTableElementNeighbours(0,grid,stencil);
testTableElementNeighbours(1,grid,stencil);
  */
  /* Neighbour lists now exists. Set a flag to acknowledge this:       */
  stencil->gridElementNeighboursAreFound=1;
  return;
}

/* 
* ====================================================================
* For all elements of one type, find all the neighbours, 
* put them into your own neighbourlist and yourself into theirs.
*
* Returns the number of neighbours found.
* For zero createLinks only counts the number of neighbours, but 
* does not actually create any links in the neighbour lists.
* ==================================================================== */
static
int findGridElementNeighbours(int ip,int jp,int kp,
			      int elementType,
			      struct grid *grid,
			      struct gridStencil *stencil,
			      int createLinks)
{
  char fctName[] = "findGridElementNeighbours";
  int nLinks;      /* Number of links needed/created  in the list 
		    * The value of this parameter will be returned     */
  int 
    ipp,jpp,kpp,   /* Grid point counters for closeby grid points */

    ipmin,ipmax,   /* Minimum and maximum indices for  */
    jpmin,jpmax,   /* candidate grid points associated */
    kpmin,kpmax,   /* with neighbour elements.         */

    ngp,           /* Number of candidate grid points  */

    ngpmax,        /* Maximum number of candidate grid 
		      points for which storage is already 
		      allocated */

    *igplist,      /* Candidate grid */
    *jgplist,      /* points listed  */ 
    *kgplist,      /* by index       */


    iElement[2],   /* (Associated) element global numbers */

    iae,           /* Counter for (associated) elements                  */
    iao,           /* Counter for (associated) elements of opposite type */

    ne,            /* Number of candidate neighbour elements  */

    nemax,         /* Maximum number of candidate elements */
    nemaxNew,      /* for which storage is already allocated */

    *eList,        /* Candidate elements listed  by index    */
    *eListTemp;    /* Needed when/if eList is increased      */

  double r1,       /* Maximum radius of associated elements */
    maxDist,       /* Maximum distance to other grid points 
		      for them to have candidate neighbour element */
    neighbourDist, /* Maximum distance to be a neighbour element
		    * = separationDistanceFactor * radius     */
    r1TestDist;    /* Minimum value of r1*nearbyElementFactor for
		    * neighbouring elements to be found. */

  point gp, gpLoc;  /* Grid points */

  /* Start counting new links:                                         */
  nLinks=0;

  /* Find an appropriate distance to search for near 
   * grid points for the associated source elements */
  r1 = 0.0;
  for (iae=0; 
       iae < grid->points[ip][jp][kp]->nAssociatedElements[elementType];
       iae++){
    r1 = MAX(r1,
	     grid->gridElements[elementType]
	     /* The next line is just the element number: */
	     [grid->points[ip][jp][kp]->associatedElements[elementType][iae]]
	       ->radius);
  }
  /* No need to make a list if the elements have 
   * zero maximum radius. 
   * This statement may be omitted provided that the next 
   * IF statement captures the zero case... */
  if (r1 == 0.0) return nLinks; 

  r1TestDist = stencil->minNotNeighbourDist 
                 - sqrt(grid->dx*grid->dx
			+grid->dy*grid->dy
			+grid->dz*grid->dz);
  /* The test distance should not be negative                          */
  r1TestDist = MAX(r1TestDist,0.0); 
  /* r1*nearbyElementFactor has to be at least r1TestDist
   * if neighbouring elements that can't be reach through 
   * neighbouring grid points should exist.  */
  if (r1*stencil->separationDistanceFactor <= r1TestDist) return nLinks;

  /* Allocate initial storage for list of candidate grid 
   * points. If too little storage is allocated here, then 
   * it will be reallocated later. */
  ngpmax  = 1000;
  igplist = (int*) calloc(ngpmax,sizeof(int));
  jgplist = (int*) calloc(ngpmax,sizeof(int));
  kgplist = (int*) calloc(ngpmax,sizeof(int));

  /* Allocate initial storage for list of candidate elements.
   * If too little storage is allocated here, then 
   * it will be reallocated later. */
  nemax  = 1000;
  eList = (int*) calloc(nemax,sizeof(int));

  /* If we DO make a list of candidate elements, reset the list 
   * number and calculate the distance for finding grid points: */
  ne = 0;
  maxDist = 
    r1*(stencil->separationDistanceFactor+2.0)
     +sqrt(grid->dx*grid->dx
	  +grid->dy*grid->dy
	  +grid->dz*grid->dz);

  /* Grid point position: gp */
  gridPointPosition(ip,jp,kp,gp,grid);

  /* Count the number of candidate grid points: */
  ipmin = (int) floor( (gp[0]-maxDist - grid->xmin)/grid->dx );
  ipmax = (int) ceil ( (gp[0]+maxDist - grid->xmin)/grid->dx );
  jpmin = (int) floor( (gp[1]-maxDist - grid->ymin)/grid->dy );
  jpmax = (int) ceil ( (gp[1]+maxDist - grid->ymin)/grid->dy );
  kpmin = (int) floor( (gp[2]-maxDist - grid->zmin)/grid->dz );
  kpmax = (int) ceil ( (gp[2]+maxDist - grid->zmin)/grid->dz );
  /* Remove any candidate points outside of grid extent: */
  ipmin = MAX(0,ipmin);
  jpmin = MAX(0,jpmin);
  kpmin = MAX(0,kpmin);
  ipmax = MIN(grid->nx-1,ipmax);
  jpmax = MIN(grid->ny-1,jpmax);
  kpmax = MIN(grid->nz-1,kpmax);
  /* Total number of candidate grid points */
  ngp = (ipmax-ipmin+1)*(jpmax-jpmin+1)*(kpmax-kpmin+1);

  /* If there is not enough storage assigned for the list, 
   * then free memory and assign more storage. 
   * Make sure the storage is incresed significantly 
   * (at least a factor of two), such that this increase 
   * won't be necessary too many times (which might be a
   * bad situation for reclaiming memory - also it eats
   * up system time). */
  if(ngp>ngpmax){
    /* The following procedure does *not* seem to lead to 
     * memory leaks. I've tested with different settings of 
     * initial ngpmax with no differences.
    printf("%s: WARNING: Increasing array sizes!!\n"
	   "     This could be a memory leak!! %d %d\n",
	   fctName,ngp,ngpmax,MAX(2*ngpmax,ngp));*/
    FREE(igplist);
    FREE(jgplist);
    FREE(kgplist);
    ngpmax = MAX(2*ngpmax,ngp);
    igplist = (int*) calloc(ngpmax,sizeof(int));
    jgplist = (int*) calloc(ngpmax,sizeof(int));
    kgplist = (int*) calloc(ngpmax,sizeof(int));
  }
  /* Make a list of candidate grid points, 
   * i.e. grid points with candidate neighbour elements */
  ngp = 0;
  for(ipp=ipmin; ipp<=ipmax; ipp++){
    for(jpp=jpmin; jpp<=jpmax; jpp++){
      for(kpp=kpmin; kpp<=kpmax; kpp++){
	if (/* If there are no grid elements */
	    grid->points[ipp][jpp][kpp]==NULL ||
	    /* Or the grid points are neighbours */
	    gridPointsAreNeighbours(ip,jp,kp,ipp,jpp,kpp,stencil) 
	    ) continue; /* Next candidate grid point */
	
	/* Check distance from candidate grid point to 
	 * grid point of interest.
	 * First ind the cooridinates of this grid point                */
	gridPointPosition(ipp,jpp,kpp,gpLoc,grid);
	/* Then find distance to the "central" grid point               */
	if (VDIST(gp,gpLoc) <= maxDist ){
	  /* This is a point with candidate neighbour
	   * elements. */

	  /* Check for sufficient storage to put these elements 
	   * (of opposite type) into the array of candidate 
	   * neighbour elements */
	  if (ne+grid->
	       points[ipp][jpp][kpp]->
	        nAssociatedElements[!elementType]
		> nemax){

	    /* This will rarely be executed. However, I'm not too
	     * sure about the memory efficiency, so for now 
	     * I'll write a warning every time this is executed!       */
	    printf("%s: WARNING: Increasing array sizes (2)!!\n"
		   "     This could be a memory leak??\n",
		   fctName);

	    /* Make sure the storage is increased significantly 
	     * (at least a factor of two), such that this increase 
	     * won't be necessary too many times. */
	    nemaxNew = 
	      MAX(2*nemax,
		    nemax
		    +grid->
		      points[ipp][jpp][kpp]->
		       nAssociatedElements[!elementType]);
	    /* Add new storage. */
	    eListTemp = (int*) calloc(nemaxNew,sizeof(int));
	    /* Move the existing ne values to the new array. */
	    for (iae=0; iae<ne; iae++)
	      eListTemp[iae]=eList[iae];
	    /* Deallocate old array */
	    FREE(eList);
	    /* Let eList point to the new array */
	    eList = eListTemp;
	    nemax = nemaxNew;
	    /* eListTemp should no longer point here
	     * (this is not important, but anyway): */
	    eListTemp = NULL;
	  }

	  /* Add each (opposite type) element associated with 
	   * this grid point to the list of candidate elements */
	  for (iao=0;
	       iao < grid->
                      points[ipp][jpp][kpp]->
		       nAssociatedElements[!elementType];
	       iao++){
	    eList[ne++] = grid->
                           points[ipp][jpp][kpp]->
	                    associatedElements[!elementType][iao];

	  }
	}
      }
    }
  }/* Done - list of candidate elements should now exist. */

  /* For each element, examine which (if any) of the candidate neighbour 
   * elements in the list are close enough to be a neighbour.           */
  for (iae=0; /* elementType elements */
       iae<grid->points[ip][jp][kp]->
	 nAssociatedElements[elementType];
       iae++){
    /* Global element number: */
    iElement[elementType] = 
      grid->points[ip][jp][kp]->
       associatedElements[elementType][iae];

    /* Separation distance for being a neighbour of this element: 
     * The "1.0" is added for the distance from the centre to the 
     * surface of the bounding sphere of the iae-element               */
    neighbourDist = 
      (stencil->separationDistanceFactor+1.0)
      * grid->gridElements[elementType][iElement[elementType]]->radius;

    for (iao=0;iao<ne;iao++){ /* Opposite type elements 
				 (possible neighbours) */
      /* Global eval number */
      iElement[!elementType] = eList[iao];

      /* See if the two elements are close enough together */
      if (neighbourDist
	  /* Add distance from the centre to the surface of the 
	   * bounding sphere of the iao-element                        */
	  +grid->gridElements[!elementType][iElement[!elementType]]->radius > 
	  /* Compare to distance from centre to centre                 */
	  VDIST(grid->gridElements[elementType][iElement[elementType]]->center,
		grid->gridElements[!elementType][iElement[!elementType]]->center)){
	/* This is a neighbour element!  Two new links will be needed. */
	nLinks += 2;
	/* Only actually add to lists if this is not just a count      */
	if(createLinks){
	  addElementsToBothTempNeighLists(iElement[elementType],
					  iElement[!elementType],
                &stencil->
                 gridElementNeighbours[elementType][iElement[elementType]]->
                  neighbourLinkedList,
                &stencil->
                 gridElementNeighbours[!elementType][iElement[!elementType]]->
                  neighbourLinkedList);
	}
      } /* End IF NEIGHBOUR */
    }   /* Next candidate neighbour element (!elementType)   */
  }     /* Next element to find neighbours for (elementType) */

  /* Free allocated memory */
  FREE(igplist);
  FREE(jgplist);
  FREE(kgplist);
  FREE(eList); 

  return nLinks;
} /* End of routine findGridElementNeighbours                       */

/* 
* ================================================================= 
* Allocate necessary memory to store grid element information, e.g.
* bounding spheres.
* ================================================================= */
static
void allocateGridElements(void *sources,int nSources,
			  void *evals,int nEvals,
			  struct grid *grid)
{
  int ielem; 

  /* Store number of sources and evals in grid structure: */
  grid->nGridElements[0] = nSources;
  grid->nGridElements[1] = nEvals;

  /* Allocate memory for gridElements */
  fprintf(stdout,"Allocating memory for grid elements ");
  if (sources!=evals || nSources!=nEvals){
    /* Allocate two separate arrays for the gridElements; */
    fprintf(stdout,"(two separate arrays) \n");
    grid->sourceAndEvalDiffer = 1;

    /* Sources: (elementType=0) */
    assert(grid->gridElements[0] == NULL);
    grid->gridElements[0] = 
      (struct gridElement**) calloc(nSources,sizeof(struct gridElement*));
    for(ielem=0; ielem<nSources; ielem++){
      grid->gridElements[0][ielem] = 
	(struct gridElement*) calloc(1,sizeof(struct gridElement));
    }

    /* Evals: (elementType=1) */
    assert(grid->gridElements[1] == NULL);
    grid->gridElements[1] = 
      (struct gridElement**) calloc(nEvals,sizeof(struct gridElement*));
    for(ielem=0; ielem<nEvals; ielem++){
      grid->gridElements[1][ielem] = 
	(struct gridElement*) calloc(1,sizeof(struct gridElement));
    }
  }
  else{
    /* Allocate only one array for the gridElements. */
    fprintf(stdout,"(one array with two names) \n");
    grid->sourceAndEvalDiffer = 0;
    /* Just go through the sources (the evals are exactly the same): */
    assert(grid->gridElements[0] == NULL);
    grid->gridElements[0] = 
      (struct gridElement**) calloc(nSources,sizeof(struct gridElement*));
    for(ielem=0; ielem<nSources; ielem++){
      grid->gridElements[0][ielem] = 
	(struct gridElement*) calloc(1,sizeof(struct gridElement));
    }
    /* Make two pointers to the allocated array, i.e. it has 
     * both of the "names" gridElements[0] and gridElements[1]. */
    assert(grid->gridElements[1] == NULL);
    grid->gridElements[1] = grid->gridElements[0];
  }
  return;
}

/*
* ================================================================= 
* For each element passed on, set a pointer from a grid element
* of the same type (source=[0] or eval=[1]), such that element 
* routines can later be called to supply information with regard to 
* the element.
* ================================================================= */
static
void associateGridElements(void *sources,
			   void *evals,
			   struct grid *grid)
{
  int ielem; 

  if (grid->sourceAndEvalDiffer){
    /* Associte the two separate arrays for the gridElements; */
    /* Sources: (elementType=0) */
    for(ielem=0; ielem < grid->nGridElements[0]; ielem++){
      /* Is the following line OK, or must the function call be used?  */
      /* grid->gridElements[0][ielem]->element = sources[ielem];*/
      grid->gridElements[0][ielem]->element = pointerToSource(sources,ielem);
    }
    /* Evals: (elementType=1) */
    for(ielem=0; ielem < grid->nGridElements[1]; ielem++){
      /* Is the following line OK, or must the function call be used?  */
      /* grid->gridElements[1][ielem]->element = evals[ielem];*/
      grid->gridElements[1][ielem]->element = pointerToEval(evals,ielem);
    }
  }
  else{
    /* Make association for only one array for the gridElements, 
     * say, the sources (it shouldn't matter): */
    for(ielem=0; ielem < grid->nGridElements[0]; ielem++){
      /* Is the following line OK, or must the function call be used?  */
      /* grid->gridElements[0][ielem]->element = sources[ielem]; */
      grid->gridElements[0][ielem]->element = pointerToSource(sources,ielem);
    }
  }

  /* TEST BLOCK:
   * Test the pointers for duplicate entries: 
   * (N^2 test! could be sorted first to be a n log n test instead)
  {
    int i,j,n;
    printf(" Making N^2 test!\n");
    n=grid->nGridElements[0];
    for (i=0;i<n-1;i++){
      for (j=i+1;j<n;j++){
	if(grid->gridElements[0][i]->element==
	   grid->gridElements[0][j]->element)
	{
	  printf(" Grid elements %d and %d points to the same element\n",
		 i,j);
	  exit(1);
	}
      }
    }
    printf(" Test passed - all grid elements point to unique elements.\n");
  }
  */

  return;
}

/*
* ================================================================= 
* Find reasonably tight bounding spheres around all elements using 
* calls to element-specific routines.
* ================================================================= */
static
void findGridElementBoundingSpheres(struct grid *grid)
{
  int ielem;

  /* Find bounding spheres of all grid elements */
  if(grid->sourceAndEvalDiffer){
    /* Grid sources: */
    for (ielem=0; ielem < grid->nGridElements[0]; ielem++){
      findSourceBoundingSphere(grid->gridElements[0][ielem]->center,
			       &grid->gridElements[0][ielem]->radius,
			       grid->gridElements[0][ielem]->element);
      /* Set projectLevel to unity. This should be changed later, 
       * if LARGE panels are given special treatment. IF other
       * projection levels are used, then at leastthe routine 
       * precorrectDirectMat_1 (and any other routine that uses 
       * precomputed grid-to-grid information) may need some 
       * implementation changes.                                       */
      grid->gridElements[0][ielem]->projectLevel = (unsigned short int) 1;
    }
    /* Grid evals: */
    for (ielem=0; ielem < grid->nGridElements[1]; ielem++){
      findEvalBoundingSphere(grid->gridElements[1][ielem]->center,
			     &grid->gridElements[1][ielem]->radius,
			     grid->gridElements[1][ielem]->element);
      /* Set projectLevel to unity. This should be changed later, 
       * if LARGE panels are given special treatment                   */
      grid->gridElements[0][ielem]->projectLevel = (unsigned short int) 1;
    }
  }
  else{
    /* Grid sources and evals (which are the same). 
     * So, it should not matter whether findSourceBoundingSphere 
     * or findEvalBoundingSphere is called. */
    for (ielem=0; ielem < grid->nGridElements[0]; ielem++){
      findSourceBoundingSphere(grid->gridElements[0][ielem]->center,
			       &grid->gridElements[0][ielem]->radius,
			       grid->gridElements[0][ielem]->element);
      /* Set projectLevel to unity. This should be changed later, 
       * if LARGE panels are given special treatment                   */
      grid->gridElements[0][ielem]->projectLevel = (unsigned short int) 1;
    }
  }
  return;
}

/* 
* ==================================================================== 
* Setup physical dimensions and the number of grid points in each 
* direction. Also store the coordinates of the grid points in three
* vectors, one for each cardinal direction.
*
* "gridMethod" determines which routine will be called for setting 
* grid size:
*
*  gridMethod=1:  setupGridSizeNomN (using iDir and one of [p2x,p2y,p2z])
*  gridMethod=2: setupGridSizeNi    (using iDir and one of [p2x,p2y,p2z])
*  gridMethod=3: setupGridSizeNNN   (using p2x, p2y and p2z)
* ==================================================================== */
static
void setupGridDimensions(struct grid *grid,
			 int gridMethod,
			 int p2x, int p2y, int p2z,
			 int iDir,
			 int verbose)
{
  char fctName[] = "setupGridDimensions";
  int i,j,iType, m[3], stencilrange;
  double extraLayers;
  point minCoordinates, maxCoordinates; /* Stores the minimum extent of
					 * the grid in 3D space.       */

  /* Powers of 2 in each direction (not all may be specified)          */
  m[0] = p2x;
  m[1] = p2y;
  m[2] = p2z;

  /* Figure out how large a grid will be needed - so that all 
   * gridElements lie within the grid                                  */
  minCoordinates[0]=DBL_MAX;
  minCoordinates[1]=DBL_MAX;
  minCoordinates[2]=DBL_MAX;
  maxCoordinates[0]=-(DBL_MAX/4.0);
  maxCoordinates[1]=-(DBL_MAX/4.0);
  maxCoordinates[2]=-(DBL_MAX/4.0);
  for (iType=0;iType<2;iType++){
    for (i=0;i<grid->nGridElements[iType];i++){
      for (j=0;j<3;j++){
	minCoordinates[j]= 
	  MIN(minCoordinates[j],
	      grid->gridElements[iType][i]->center[j]
	      -grid->gridElements[iType][i]->radius);
	maxCoordinates[j]= 
	  MAX(maxCoordinates[j],
	      grid->gridElements[iType][i]->center[j]
	      +grid->gridElements[iType][i]->radius);
      }
    }
  }
  /* Test output: 
  printf("%s: Minimum grid extent needed: \n"
	 "  (%16.12f, %16.12f, %16.12f)\n"
	 "  (%16.12f, %16.12f, %16.12f)\n",
	 fctName,
	 minCoordinates[0],minCoordinates[1],minCoordinates[2],
	 maxCoordinates[0],maxCoordinates[1],maxCoordinates[2]);*/
  

  /* Add a (very) little slack to the grid, so that all elements are 
   * fully inside the grid, i.e. no part of any element is on the 
   * border of the grid.                                               */
  for (j=0;j<3;j++){
    minCoordinates[j] -= sqrt(DBL_EPSILON)*ABS(minCoordinates[j]);
    maxCoordinates[j] += sqrt(DBL_EPSILON)*ABS(maxCoordinates[j]);
  }

  /* Set up the grid:                                                  */
  /* When projecting onto and interpolating from the grid, grid points 
   * must exist from/to which to interpolate/project. Check the 
   * interpolation (and projection) stencils to find how large an 
   * extent these stencils have. 
   * The maximum extension (of interpolation and projection stencils) 
   * will determine how much the grid will be larger than the 
   * domain set by the elements alone.                                 */
  stencilrange = stencilRange(grid->elementGridStencil);

  /* Set the needed number of extra layers of grid points around 
   * the domain set by the elements alone.
   * If the grid size will be very strictly set, then the number 
   * of extra layers may be stencilrange-0.5 (add a little more
   * to make sure that no element will be associated with grid points 
   * in the extra layers. If a "nominal" (flexible) size is used, then 
   * extraLayers=stencilrange should probably be used.                 */
  extraLayers = 1.0 * stencilrange;

  /* For testing only: */
  if(0){
    printf("%s: WARNING! Using hardcoded extraLayers!\n",fctName);
    extraLayers = 2.0;
  }

  if (gridMethod ==1){
    setupGridSizeNomN(grid,m[iDir],iDir,
		      minCoordinates,maxCoordinates,extraLayers,
		      verbose);
  }
  else if (gridMethod ==2){
    extraLayers = MAX(extraLayers-0.50,0.0) + sqrt(DBL_EPSILON);
    setupGridSizeNi(grid,m[iDir],iDir,
		    minCoordinates,maxCoordinates,extraLayers,
		    verbose);
  }
  else if (gridMethod ==3){
    extraLayers = MAX(extraLayers-0.50,0.0) + sqrt(DBL_EPSILON);
    setupGridSizeNNN(grid, p2x, p2y, p2z,
		     minCoordinates,maxCoordinates,extraLayers,
		     verbose);
  }
  else {
    fprintf(stderr,
	    "%s: ERROR: What kind of grid generation method is this? %d \n",
	   fctName,gridMethod);
    fprintf(stderr,"%s: ERROR: This halts the program.\n",fctName);
    exit(1);
  }
  /* More routines may be implemented to make the grid generation more
   * flexible. One could consider:
   * -> a nominal number of points in a specific direction
   *      (routine: setupGridSizeNomN)
   * -> a specified number of points in a specified direction
   *      (routine: setupGridSizeNi)
   * -> a specified number of points in each direction
   *      (routine: setupGridSizeNNN)
   * - a maximum number of points
   *      (not implemented)
   * - a nominal grid size
   *      (not implemented)
   * - a maximum grid size
   *      (not implemented)
   * - a minimum grid size
   *      (not implemented)                                            */


  /* Test output: 
  printf("%s: Used grid extent: \n"
	 "  (%16.12f, %16.12f, %16.12f)\n"
	 "  (%16.12f, %16.12f, %16.12f)\n",
	 fctName,
	 grid->xmin,grid->ymin,grid->zmin,
	 grid->xmax,grid->ymax,grid->zmax);
  */

  /* Allocate memory for the coordinate vectors */
  assert(grid->x == NULL);
  assert(grid->y == NULL);
  assert(grid->z == NULL);
  grid->x = calloc(grid->nx, sizeof(double)); 
  grid->y = calloc(grid->ny, sizeof(double)); 
  grid->z = calloc(grid->nz, sizeof(double)); 

  /* Calculate and store coordinates of the grid points */
  for (i=1;i<grid->nx-1;i++)
    grid->x[i] = grid->xmin + i * grid->dx;
  grid->x[0]          = grid->xmin;  /* Avoid truncation errors */
  grid->x[grid->nx-1] = grid->xmax;  /* Avoid truncation errors */

  for (i=1;i<grid->ny-1;i++)
    grid->y[i] = grid->ymin + i * grid->dy;
  grid->y[0]          = grid->ymin;  /* Avoid truncation errors */
  grid->y[grid->ny-1] = grid->ymax;  /* Avoid truncation errors */

  for (i=1;i<grid->nz-1;i++)
    grid->z[i] = grid->zmin + i * grid->dz;
  grid->z[0]          = grid->zmin;  /* Avoid truncation errors */
  grid->z[grid->nz-1] = grid->zmax;  /* Avoid truncation errors */

  if (verbose) 
    printf("%s: Total number of grid nodes: %d \n",
	   fctName,(grid->nx)*(grid->ny)*(grid->nz));

  return;
}

/*
* ==================================================================== 
* Stencil (or "molecule") for neighbouring grid points.
* Allocate memory for the stencil and assign values to it. 
* Note that since a grid point is close to itself for ease 
* of implementation it should also be its own neighbour.
*
* If refreshGridSizeOnly!=0, then it is assumed that the 
* stencils already exist, and that only values depending on grid
* size need to be updated.
* ==================================================================== */
static
void setupGridPointNeighbourStencil(struct grid *grid,
				    struct gridStencil *stencil,
				    int stencilTypeIn,
				    double stencilSize,
				    double separationDistanceFactor,
				    int refreshGridSizeOnly)
{
  /*char fctName[] = "setupGridPointNeighbourStencil";*/
  int storeInArray, nPoints, ielem, stencilType, 
    stencilUse; /* 1: Neighbours; 2: projection/interpolation */
  /* Figure out what kind of stencil use is intended: */
  if (stencilTypeIn<100){
    /* These are stencils for finding element neighbours, i.e. for 
     * direct part and preconditioning */
    stencilType = stencilTypeIn;
    stencilUse  = 1;
  }
  else if (stencilTypeIn<200){
    /* These are stencils for proje,ction and interpolation */
    stencilType = stencilTypeIn-100;
    stencilUse  = 2;
  }

  if (refreshGridSizeOnly==0){
    /* Set the separation distance factor                              */
    stencil->separationDistanceFactor = separationDistanceFactor;
    /* Count the number of points in the stencil                       */
    nPoints = 
      findStencilIndices(stencilType,stencilSize,
			 storeInArray=0,
			 refreshGridSizeOnly,
			 grid->dx,grid->dy,grid->dz,
			 stencil);
    stencil->nStencilPoints = nPoints;

    /* Allocate needed memory 
     *  1) for the stencil indices.                                    */
    assert(stencil->ipNeighbours == NULL);
    assert(stencil->jpNeighbours == NULL);
    assert(stencil->kpNeighbours == NULL);
    stencil->ipNeighbours = (int*) calloc(nPoints,sizeof(int));
    stencil->jpNeighbours = (int*) calloc(nPoints,sizeof(int));
    stencil->kpNeighbours = (int*) calloc(nPoints,sizeof(int));
    /*  2) for storing information on neighbours that cannot be 
     *     reached using the stencil indices.                          */

    if(stencilUse == 1) { /* Neighbors (direct/precondition) stencil.  */
      /* Sources: (elementType=0) */
      assert(stencil->gridElementNeighbours[0] == NULL);
      stencil->gridElementNeighbours[0] = 
	(struct gridElementNeighbours**) 
	calloc(grid->nGridElements[0],sizeof(struct gridElementNeighbours*));
      for(ielem=0; ielem<grid->nGridElements[0]; ielem++){
	stencil->gridElementNeighbours[0][ielem] = 
	  (struct gridElementNeighbours*) 
	  calloc(1,sizeof(struct gridElementNeighbours));
      }
      if (grid->sourceAndEvalDiffer){
	/* Setup of the evals since they are not the same as the sources */
	/* Evals: (elementType=1) */
	assert(stencil->gridElementNeighbours[1] == NULL);
	stencil->gridElementNeighbours[1] = 
	  (struct gridElementNeighbours**) 
	  calloc(grid->nGridElements[1],sizeof(struct gridElementNeighbours*));
	for(ielem=0; ielem<grid->nGridElements[1]; ielem++){
	  stencil->gridElementNeighbours[1][ielem] = 
	    (struct gridElementNeighbours*) 
	    calloc(1,sizeof(struct gridElementNeighbours));
	}
     }
      else {
	/* Since the evals are the same as the sources, just let the evals 
	 * point to the sources.                                           */
	/* Make two pointers to the allocated array, i.e. it has both of the 
	 * "names" gridElementsNeighbours[0] and gridElementsNeighbours[1] */
	assert(stencil->gridElementNeighbours[1] == NULL);
	stencil->gridElementNeighbours[1] = stencil->gridElementNeighbours[0];
      }
    }
    /* Now the structure needed for the element lists are set up, 
     * but the lists themselves have not been made:                    */
    stencil->gridElementNeighboursAreFound=0;
    /* That (obviously) means that they have not been counted either)  */
    stencil->nGridElementNeighboursAreFound=0;
  }
  else {
    /* In this case the stencil should already be allocated.           */
    assert(stencil!=NULL);
  }

  /* store stencil indices in the array (check the returned int)       */
  nPoints=findStencilIndices(stencilType,stencilSize,
			     storeInArray=1,
			     refreshGridSizeOnly,
			     grid->dx,grid->dy,grid->dz,
			     stencil);
  assert(nPoints==stencil->nStencilPoints);


  if(stencilUse == 2) {  /* Interpolation/projection stencil.          */
    polyInterp(stencilTypeIn, 
	       grid->dx, grid->dy, grid->dz,
	       refreshGridSizeOnly,
	       stencil);
  }

  return;
} /* End of setupGridPointNeighbourStencil                             */

/*
* ==================================================================== 
* Find number of points in stencil and assign the stencil indices to 
* an array if necessary. 
* Note that since a grid point is close to itself for ease 
* of implementation it should also be its own neighbour.
* It is the intention that this routine will be called twice. Once for 
* counting the number of points in the stencil, and once for storing 
* the indices in an array. The calling routing will assign needed 
* memory between the calls.
*
* The stencil type (stencilType) determines the shape of the stencil:
* 1: "Cube." The number of points is (2n+1)^3 where n is integer.
*    Make a cube with side length 2*stencilSize*"gridsize" centered 
*    on the node of interest. All grid points inside the cube are
*    neighbours. A small number will be added to stencilSize
*    such that points on the boundary are counted as "inside".
*    Sizes where the number of points change for the cube are 
*    integers. 
*     Size  #points
*      0-1       1
*      1-2      27
*      2-3     125
*      3-4     343
*      4-5     729
*
*    Example (size 1 - 27 points in 3D)
*       x: normal grid point 
*       o: stencil point
*
*                  |   |   |   |   | 
*                - x - x - x - x - x -
*                  |   |   |   |   | 
*                - x - o - o - o - x -
*                  |   |   |   |   | 
*                - x - o - o - o - x -
*                  |   |   |   |   | 
*                - x - o - o - o - x -
*                  |   |   |   |   | 
*                - x - x - x - x - x -
*                  |   |   |   |   | 
*
* 2: "Sphere." Make a sphere of radius stencilSize*dx centered on 
*    the node of interest. All grid points inside the sphere are
*    neighbours. A small number will be added to stencilSize
*    such that points on the boundary are counted as "inside".
*    Note that the "sphere" never contains more stencil points 
*    than does the cube of the same size.
*    Example (size 1 - 7 points in 3D)
*       x: normal grid point 
*       o: stencil point
*
*                  |   |   |   |   | 
*                - x - x - x - x - x -
*                  |   |   |   |   | 
*                - x - x - o - x - x -
*                  |   |   |   |   | 
*                - x - o - o - o - x -
*                  |   |   |   |   | 
*                - x - x - o - x - x -
*                  |   |   |   |   | 
*                - x - x - x - x - x -
*                  |   |   |   |   | 
*
* 3: "Single Point."  This stencil is just a single point
*    NEEDS additional types to handle projection-interpolation.
* 
* The stencil size (stencilSize) determines the number of points in 
* the stencil - see above.
*
* storeInArray: 0: Do not store. Just count number of points.
*               nonzero: Count and store stencil indices (storage must
*                        be allocated before this call).
* OUTPUT:
* ipNeighbours: First  (x) index of the stencil points. 
* jpNeighbours: Second (y) index of the stencil points. 
* kpNeighbours: Third  (z) index of the stencil points. 
*   [One set of indices (i,j,k) are given for each stencil point.]
* minNotNeighbourDist: Minimum distance between grid points that 
*                      are *not* neighbours. Actually, a somewhat more 
*                      loose requirement is needed: The distance between
*                      two points that are not neighbours must never be 
*                      smaller than minNotNeighbourDist. Thus the value
*                      could be set to zero, even though this is somewhat 
*                      inefficient.
* ==================================================================== */
static 
int findStencilIndices(int stencilType,
		       double stencilSize,
		       int storeInArray,
		       int refreshGridSizeOnly,
		       double dx,
		       double dy,
		       double dz,
		       struct gridStencil *stencil)
{
  char fctName[] = "findStencilIndices";
  int nNeighbours,  /* Number of stencil points             */
    imin,imax,      /* Search limits                        */
    i,j,k;          /* Counters                             */
  double epsi,      /* Precision for finding neighbours     */
    minDist,        /* Temporary variable for minNotNeighbourDist */
    dist;           /* Distance to candidate stencil point  */
  int verbose=0;

  /* Precision for finding neighbours     */
  epsi = sqrt(DBL_EPSILON);
  nNeighbours = 0;

  /* If we are only updating the grid size, then do not store stuff 
   * in the arrays - it is already there.                              */
  if (refreshGridSizeOnly) storeInArray=0;

  /* If the grid is not uniform report error and exit                  */
  if (dx!=dy || dy!=dz){
    fprintf(stderr,"%s: Error. Non-uniform grid encountered.\n"
	    "       (dx,dy,dz) = (%g,%g,%g)\n"
	    "       This routine is not build to handle non-uniform grids\n"
	    "       This error is terminal\n",
	    fctName,dx,dy,dz);
    _EXIT_;
  }

  /* If the results are to be stored in array, then see if the array 
   * pointers do not point to NULL.                                    */
  if(storeInArray &&
     (stencil->ipNeighbours==NULL ||
      stencil->jpNeighbours==NULL ||
      stencil->kpNeighbours==NULL) ){
    /* Should store in array, but arrays are not here!
     * Write error statement and exit.                                 */
    fprintf(stderr,"%s: Error. storeInArray set, but one or more arrays\n"
	    "       not allocated.\n"
	    "       This error is terminal\n",
	    fctName);
    _EXIT_;
  }


  if(stencilType==1){ /* Cube type */
    imax = (int) floor(stencilSize);
    imin = -imax;
    for(i=imin;i<=imax;i++){
      for(j=imin;j<=imax;j++){
	for(k=imin;k<=imax;k++){
	  if(storeInArray &&
	     stencil->ipNeighbours!=NULL &&
	     stencil->jpNeighbours!=NULL &&
	     stencil->kpNeighbours!=NULL){
	    stencil->ipNeighbours[nNeighbours]=i;
	    stencil->jpNeighbours[nNeighbours]=j;
	    stencil->kpNeighbours[nNeighbours]=k;
	  }
	  nNeighbours++;
	}   /* Next k                 */
      }     /* Next j                 */
    }       /* Next i                 */
    /* For the cube the smallest distance to a point which is not a 
     * neighbour is easy:                                              */
    stencil->minNotNeighbourDist = ((double) (imax+1))*dx;

    /* For the cube we can easily test the total number of neighbours  */
    if(nNeighbours!=intpow(2*imax+1,3)){
    fprintf(stderr,"%s: %d points found in stencil, but %d expected.\n"
	    "       This error is terminal!\n",
	    fctName,nNeighbours,intpow(imax,3));
    exit(1);
    }
  }/* End if stencilType==1  */
  else if (stencilType==2){ /* Sphere type */
    imax = (int) floor(stencilSize);
    imin = -imax;
    /* For a uniform grid this is the smallest distance to a point 
     * which are not examined in the following, i.e. all points at
     * closer range will be examined.                                  */
    minDist = ((double) (imax+1))*dx;
    for(i=imin;i<=imax;i++){
      for(j=imin;j<=imax;j++){
	for(k=imin;k<=imax;k++){
	  dist = sqrt((double)(i*i+j*j+k*k));
	  if(dist <= stencilSize){
	    /* This is a neighbour */
	    if(storeInArray &&
	       stencil->ipNeighbours!=NULL &&
	       stencil->jpNeighbours!=NULL &&
	       stencil->kpNeighbours!=NULL){
	      stencil->ipNeighbours[nNeighbours]=i;
	      stencil->jpNeighbours[nNeighbours]=j;
	      stencil->kpNeighbours[nNeighbours]=k;
	    }
	    nNeighbours++;
	  }
	  else{
	    /* This is not a neighbour - so decrease minimum 
	     * distance to other points if appropriate:                */
	    minDist = MIN(minDist,dist);
	  }
	}   /* Next k                 */
      }     /* Next j                 */
    }       /* Next i                 */
    stencil->minNotNeighbourDist = minDist*dx;
  }         /* End if stencilType==2  */
  else if (stencilType==3){ 
    /* Single point stencil. This can be obtained using types 1 or 2  
     * with a distance < 1, so this *could* be eliminated. */
    if(storeInArray &&
       stencil->ipNeighbours!=NULL &&
       stencil->jpNeighbours!=NULL &&
       stencil->kpNeighbours!=NULL){
      stencil->ipNeighbours[nNeighbours]=0;
      stencil->jpNeighbours[nNeighbours]=0;
      stencil->kpNeighbours[nNeighbours]=0;
    }
    nNeighbours++;
    stencil->minNotNeighbourDist = dx;
  }
  else if (stencilType == 10){
    /* Eight-point corners of a cube stencil */
    imax = (int) dblRound(stencilSize);
    for (i=-1; i<=1; i++){
      for (j=-1; j<=1; j++){
	for (k=-1; k<=1; k++){
	  if (i*j*k != 0){
	    if (storeInArray){
	      stencil->ipNeighbours[nNeighbours]= imax*i;
	      stencil->jpNeighbours[nNeighbours]= imax*j;
	      stencil->kpNeighbours[nNeighbours]= imax*k;
	    }
	    nNeighbours++;
	  }
	}
      }
    }

  }


  else {
    /* This is bad. This stencil type is not known. 
     * Report error and exit.                                          */
    fprintf(stderr,"%s: What type of stencil is %d?\n"
	    "       This error is terminal!\n",fctName,stencilType);
    exit(1);
  }

  /* Write obtained stencil to screen */
  if(storeInArray && verbose){
    fprintf(stdout,"%s: The following %d nodes are in the stencil:\n",
	    fctName,nNeighbours);
    for (i=0;i<nNeighbours;i++){
      fprintf(stdout,"%s: (%5d,%5d,%5d)\n",
	      fctName,
	      stencil->ipNeighbours[i],
	      stencil->jpNeighbours[i],
	      stencil->kpNeighbours[i]);
    }
  }

  return nNeighbours;
}

/*
* ==================================================================== 
* Figure out an "optimum" grid setting.
* This routine is based on a nominal number of points (2**nNom) in a
* specified direction (iDir=0,1,2) of the grid.
* Make sure that roughly extraLayers set of points are added outside the 
* physical extension of the grid (which is supplied). This is to make 
* sure that stencil points never lie outside the grid. Note that 
* "extraLayers" does not need to be an integer. For example, if 
* projection/interpolation stencils of sixe 3x3x3 are used, then the 
* (grid) elements should not be associated with the points at the 
* boundary of the grid. Presumbly, this can be satisfied by 
* extraLayers=0.5. However, since the actual grid size is unknown 
* when the routine is called, extraLayers=1.0 is recommended.
* ==================================================================== */
static
void setupGridSizeNomN(struct grid *grid,int nNomPow,int iDir,
		       point minCoordinates,point maxCoordinates,
		       double extraLayers, int verbose)
{
  char fctName[] = "setupGridSizeNomN";
  int i, j, powLow, powHigh, nLow, nHigh, iBest;
  double deltaNom, delta[3], volume[3], sideLengths[3], 
    nTarget, powTarget, deltaLow, deltaHigh, vol;
  int nxyz[3][3], nPow2[3][3], nTot[3];
  long nNom;

  /* Make a rough estimate of the grid size */
  nNom = intpow(2,nNomPow);
  deltaNom = (maxCoordinates[iDir] - minCoordinates[iDir])
    /((double)(nNom-1)-2.0*extraLayers);

  /* Calculate side lengths adding extra space needed to the 
   * size of the problem.                                              */
  for (i=0;i<3;i++){
    sideLengths[i] = 
      maxCoordinates[i]-minCoordinates[i] /* Problem size */
      +2.0*extraLayers*deltaNom;          /* Extra space  */
  }
  

  nNom = intpow(2,nNomPow);
  assert(nNom>1);
  deltaNom = (sideLengths[iDir])/(double)(nNom-1);

  if (verbose)
    printf("%s: Minimum side lengths: %g %g %g\n",fctName,
	   sideLengths[0],sideLengths[1],sideLengths[2]); 
  for (i=0;i<3;i++){
    /* Fit to the side length in direction i                           */
    nTarget   = 1.0 + sideLengths[i] / deltaNom;
    powTarget = log(nTarget)/log(2.0);
    powLow    = (int) floor(powTarget); /* This should not become zero!*/
    powLow    = MAX(powLow,1);
    powHigh   = (int) ceil(powTarget);
    nLow      = intpow(2,powLow); 
    nHigh     = intpow(2,powHigh);
    deltaLow  = (sideLengths[i])/(double)(nLow-1);  /* >= deltaNom     */
    deltaHigh = (sideLengths[i])/(double)(nHigh-1); /* <= deltaNom     */

    if( deltaLow-deltaNom < deltaNom-deltaHigh){
      /* nLow is the better choise                                     */
      delta[i]= deltaLow;
    }
    else{
      /* nHigh is the better choise                                    */
      delta[i]= deltaHigh;
    }
    /* Calculate the grid volume and total number of points needed     */
    volume[i] = 1.0;
    nTot[i]   = 1;
    for (j=0;j<3;j++){
      /* Now we *must* go for the high number of n on each side:       */
      nTarget    = 1.0 + sideLengths[j] / delta[i];
      powTarget  = log(nTarget)/log(2.0);
      powHigh    = (int) ceil(powTarget);
      nHigh      = intpow(2,powHigh);
      nxyz[i][j] = nHigh;
      nPow2[i][j]= powHigh;
      volume[i]  = volume[i] * (nHigh-1)*delta[i];
      nTot[i]    = nTot[i]*nHigh;
    }
  }
  /* Choose the fit with the smallest volume (could also be related to
   * the total number of grid points).                                 */
  iBest = 0;
  vol   = volume[0];
  for(i=1;i<3;i++){
    if(volume[i]<vol){
      iBest = i;
      vol = volume[i];
    }
  }
  /*
  printf("%s: Min/max X: %g %g\n",
	 fctName,minCoordinates[0],maxCoordinates[0]);
  printf("%s: Min/max Y: %g %g\n",
	 fctName,minCoordinates[1],maxCoordinates[1]);
  printf("%s: Min/max Z: %g %g\n",
	 fctName,minCoordinates[2],maxCoordinates[2]);
  */

  /* Store obtained data in grid structure:                            */
  grid->iDirTightFit = iBest;
  grid->nx  = nxyz[iBest][0];
  grid->ny  = nxyz[iBest][1];
  grid->nz  = nxyz[iBest][2];
  grid->p2x = nPow2[iBest][0];
  grid->p2y = nPow2[iBest][1];
  grid->p2z = nPow2[iBest][2];
  grid->dx  = delta[iBest];
  grid->dy  = delta[iBest];
  grid->dz  = delta[iBest];
  /* Set grid bounds.                                                  */
  grid->xmin=minCoordinates[0]-extraLayers*grid->dx; /* Min. x-coord.  */
  grid->ymin=minCoordinates[1]-extraLayers*grid->dy; /* Min. y-coord.  */
  grid->zmin=minCoordinates[2]-extraLayers*grid->dz; /* Min. z-coord.  */
  grid->xmax=minCoordinates[0]+grid->dx*(double) (grid->nx -1);
  grid->ymax=minCoordinates[1]+grid->dy*(double) (grid->ny -1);
  grid->zmax=minCoordinates[2]+grid->dz*(double) (grid->nz -1);

  /*
  printf("%s: GRID Min/max X: %g %g\n",
	 fctName,grid->xmin,grid->xmax);
  printf("%s: GRID Min/max Y: %g %g\n",
	 fctName,grid->ymin,grid->ymax);
  printf("%s: GRID Min/max Z: %g %g\n",
	 fctName,grid->zmin,grid->zmax);
  printf("%s: Best fit direction: %d\n",
	 fctName,iBest);
  */

  printf("%s: GRID nx ny nz delta: %d %d %d %g\n",
	 fctName,grid->nx,grid->ny,grid->nz,grid->dx);

  return;
} /* End routine setupGridSizeNomN                                     */

/*
* ==================================================================== 
* Setup the grid using a specified number of grid points (given as 
* a power of two) in a  specified cardinal direction. The number of 
* points in each of the other directions will be calculated to be 
* "large enough".
* extraLayers*delta amount of padding is added around the discretization,
* the intent being that all interpolation and projection points will 
* exist. If, for instance, an element is associated with a grid point on
* the boundary of the grid, then all neighbouring grid points (as defined 
* in the interpolation/projection stencil) may not exist. The use of 
* "extraLayers" avoids this obstacle.
* ==================================================================== */
static
void setupGridSizeNi(struct grid *grid,int nNomPow,int iDir,
		     point minCoordinates, point maxCoordinates,
		     double extraLayers, int verbose)
{
  char fctName[] = "setupGridSizeNi";
  double length, sideLength,delta, nTarget,powTarget;
  int i, ni, nPow[3], n[3];

  ni     = intpow(2,nNomPow);
  length = maxCoordinates[iDir]-minCoordinates[iDir];
  delta  = length / ((double) ni - 2.0*extraLayers -1.0);

  /*
    printf("%s: Slack = %16.8e\n",fctName,
    ABS((double)ni-1.0-length/delta-2.0*extraLayers));*/

  /* Find number of gridpoints in each cardinal direction.             */
  for (i=0;i<3;i++){
    if (i==iDir) {
      nPow[i]  = nNomPow;
      n[i]     = ni;
    }
    else {
      sideLength = maxCoordinates[i]-minCoordinates[i]+
	2.0*extraLayers*delta;
      printf("%s: Side %d has length %16.4f modified to %16.4f\n",
	   fctName,i,maxCoordinates[i]-minCoordinates[i],sideLength);
      nTarget    = 1.0 + sideLength / delta;
      powTarget  = log(nTarget)/log(2.0);
      nPow[i]    = (int) ceil(powTarget);
      n[i]       = intpow(2,nPow[i]);
    }   
  }
  /* Store obtained data in grid structure:                            */
  grid->iDirTightFit = iDir;
  grid->nx  = n[0];
  grid->ny  = n[1];
  grid->nz  = n[2];
  grid->p2x = nPow[0];
  grid->p2y = nPow[1];
  grid->p2z = nPow[2];
  grid->dx  = delta;
  grid->dy  = delta;
  grid->dz  = delta;
  /* Set grid bounds.                                                  */
  grid->xmin=minCoordinates[0]-extraLayers*grid->dx; /* Min. x-coord.  */
  grid->ymin=minCoordinates[1]-extraLayers*grid->dy; /* Min. y-coord.  */
  grid->zmin=minCoordinates[2]-extraLayers*grid->dz; /* Min. z-coord.  */
  grid->xmax=minCoordinates[0]+grid->dx*(double) (grid->nx -1);
  grid->ymax=minCoordinates[1]+grid->dy*(double) (grid->ny -1);
  grid->zmax=minCoordinates[2]+grid->dz*(double) (grid->nz -1);

  if (verbose) 
    printf("%s: GRID nx ny nz delta: %d %d %d %g\n",
	   fctName,grid->nx,grid->ny,grid->nz,grid->dx);

 return;
}
/*
* ==================================================================== 
* Setup the grid using a specified number of grid points (given as 
* a power of two) in each cardinal direction. The grid size will be 
* uniform in each direction and determined to give a "tight fit" in 
* (at least) one of the cardinal directions.
* The routine returns the direction in which a tight fit is obtained.
* ==================================================================== */
static
int setupGridSizeNNN(struct grid *grid,
		     int nNomPow0,int nNomPow1,int nNomPow2,
		     point minCoordinates, point maxCoordinates,
		     double extraLayers, int verbose)
{
  char fctName[] = "setupGridSizeNNN";
  int i, nPow[3],n[3], iTight;
  double delta, thisDelta, length;

  nPow[0] = nNomPow0;
  nPow[1] = nNomPow1;
  nPow[2] = nNomPow2;
  n[0]    = intpow(2,nNomPow0);
  n[1]    = intpow(2,nNomPow1);
  n[2]    = intpow(2,nNomPow2);

  /* Find delta and direction of tight fit 
   * (largest delta must be chosen)                                    */
  for (i=0, delta=-1.0;i<3; i++){
    length = maxCoordinates[i]-minCoordinates[i];
    thisDelta  = length / ((double) n[i] - 2.0*extraLayers -1.0);
    if (thisDelta>delta){
      delta  = thisDelta;
      iTight = i;
    }
  }

  /* Store obtained data in grid structure:                            */
  grid->iDirTightFit = iTight;
  grid->nx  = n[0];
  grid->ny  = n[1];
  grid->nz  = n[2];
  grid->p2x = nPow[0];
  grid->p2y = nPow[1];
  grid->p2z = nPow[2];
  grid->dx  = delta;
  grid->dy  = delta;
  grid->dz  = delta;
  /* Set grid bounds.                                                  */
  grid->xmin=minCoordinates[0]-extraLayers*grid->dx; /* Min. x-coord.  */
  grid->ymin=minCoordinates[1]-extraLayers*grid->dy; /* Min. y-coord.  */
  grid->zmin=minCoordinates[2]-extraLayers*grid->dz; /* Min. z-coord.  */
  grid->xmax=minCoordinates[0]+grid->dx*(double) (grid->nx -1);
  grid->ymax=minCoordinates[1]+grid->dy*(double) (grid->ny -1);
  grid->zmax=minCoordinates[2]+grid->dz*(double) (grid->nz -1);

  /*printf("%s: Tight fit in direction %d.\n",fctName,iTight);*/

  if (verbose)
    printf("%s: GRID nx ny nz delta: %d %d %d %g\n",
	   fctName,grid->nx,grid->ny,grid->nz,grid->dx);

  return iTight;
}
/* 
* ==================================================================== 
* Allocate memory for the panel-on-grid-point association info.
* ==================================================================== */
static
void allocateGridPoints(struct grid *grid)
{
  int i,j,k;

  /* Allocate memory: */
  assert(grid->points == NULL);
  grid->points = (gridPoint****) calloc((grid->nx),sizeof(gridPoint***));
  for (i=0;i<(grid->nx);i++){
    grid->points[i] = (gridPoint***) calloc((grid->ny),sizeof(gridPoint**));
    for (j=0;j<(grid->ny);j++){
      grid->points[i][j]= (gridPoint**) calloc((grid->nz),sizeof(gridPoint*));
      for (k=0;k<(grid->nz);k++){
	grid->points[i][j][k] = NULL;

	/* The grid points (pointers) point[i][j][k] are default 
	 * initialized to point to NULL. However, each grid point 
	 * have memory allocated by e.g.: 
	 * grid->points[i][j][k] = (gridPoint*) calloc(1,sizeof(gridPoint));
	 */

      }
    }
  }
  return;
} /* End of routine allocateGridPoints */
/* 
* ==================================================================== 
* Allocate memory for the panel-on-grid-point association info.
* ==================================================================== */
static
void allocateGridAssociatedElements(struct grid *grid)
{
  /* Allocate storage for associating elements of both types 
   * with grid points */
  printf("Allocating storage for %d source elements on the grid\n",
	 grid->nGridElements[0]);
  assert(grid->associatedElements[0] == NULL);
  grid->associatedElements[0] = 
    (int*) calloc(grid->nGridElements[0],sizeof(int));

  if (grid->sourceAndEvalDiffer){
    printf("Allocating storage for %d evaluation elements on the grid\n",
	   grid->nGridElements[1]);
    assert(grid->associatedElements[1] == NULL);
    grid->associatedElements[1] = 
      (int*) calloc(grid->nGridElements[1],sizeof(int));
  }
  else{
    printf(" Storage for evaluation elements on the grid not needed.\n");   
    grid->associatedElements[1] = grid->associatedElements[0];
  }
  return;
} /* End of routine allocateGridAssociatedElements */

/* 
* ================================================================= 
* Associate grid elements with grid points 
* ("Put elements onto the grid")
* ================================================================= */
static
void associateElementsWithGridPoints(struct grid *grid)
{
  if (grid->sourceAndEvalDiffer){
    /* Associate grid sources with grid points */
    associateOneTypeElementsWithGridPoints(0,grid);
    /* Associate grid evals with grid points */
    associateOneTypeElementsWithGridPoints(1,grid);
  }
  else{
    /* Associate grid sources with grid points 
     * (the evals are then exactly the same) */
    associateOneTypeElementsWithGridPoints(0,grid);
  }
  return;
}

/* 
* ================================================================= 
* For each grid element (source or eval) use the center of its 
* bounding sphere to associate it with a point on the grid. 
* Count the elements at each grid point at a first pass [find the 
* length of each list]. Then distribute the needed memory among the 
* grid points [make space for each list]. Finally, use a second pass 
* to associate each elements with a grid point [fill in each list].
* ================================================================= */ 
static
void associateOneTypeElementsWithGridPoints(int elementType,
					    struct grid *grid)
{
  char fctName[] = "associateOneTypeElementsWithGridPoints";
  int iae,ip,jp,kp;
  int *istorage, nElements, verbose=0;

  nElements = grid->nGridElements[elementType];

  if (verbose)
    printf(" In %s: %d type %d elements to process \n",
	   fctName,nElements,elementType);

  /* Count the number of associated elements at each grid point: */
  for (iae=0; iae<nElements; iae++){

   /* Check if the entire bounding sphere is located within the 
    * extend of the grid. */
    if (!sphereIsFullyInsideGrid(grid->gridElements[elementType][iae]->center,
				 grid->gridElements[elementType][iae]->radius,
				 grid)){
      /* This is bad. The grid is not big enough */
      fprintf(stderr,"%s: WARNING.\n"
	      "    Bounding sphere of grid element[%d] %d"
	      " extends beyond the mesh\n",fctName,elementType,iae);
    }

    /* Find closest grid point to the center of the bounding sphere */
    closestGridPoint(&ip, &jp, &kp, 
		     grid->gridElements[elementType][iae]->center,
		     grid);

    /* Check if the found grid point has memory allocated, if not
     * then allocate the needed memory (calloc will zero the entries) */
    if (grid->points[ip][jp][kp]==NULL)
      grid->points[ip][jp][kp] = (gridPoint*) calloc(1,sizeof(gridPoint));

    /* Step the counter of sources on this grid point */
    grid->points[ip][jp][kp]->nAssociatedElements[elementType]++;

    /* Next element */
  }

  /* Assign adequate storage to hold the element numbers;
   * the full array has already been allocated, so just
   * pointer operations are needed here. */
  /*  First element of storage: */
  istorage = &grid->associatedElements[elementType][0]; 
  for (ip=0; ip<grid->nx; ip++){
    for (jp=0; jp<grid->ny; jp++){
      for (kp=0; kp<grid->nz; kp++){
	if (grid->points[ip][jp][kp] != NULL){
	  grid->points[ip][jp][kp]->associatedElements[elementType] = istorage;
	  istorage += grid->points[ip][jp][kp]->nAssociatedElements[elementType];
	  
	  /* Zeroing the number of found sources at each point 
	   * is needed for counting the second pass */
	  grid->points[ip][jp][kp]->nAssociatedElements[elementType] = 0;
	}
      }
    }
  }

  /* Second pass over sources. Generate the lists: */

  /* Most of this is the same as for the first pass. */
  for (iae=0; iae<grid->nGridElements[elementType]; iae++){

    /* Find closest grid point */
    closestGridPoint(&ip, &jp, &kp, 
		     grid->gridElements[elementType][iae]->center,
		     grid);

    /* Check if the found grid point has memory allocated, if not
     * then report an ERROR (This shouldn't occur, though) */
    assert(grid->points[ip][jp][kp] != NULL);

    /* Make sure that the source number will replace a zero 
     * entry in the storage array. This should always be true, 
     * and is mostly intended for debugging etc. Obviously, this 
     * statement doesn't check if the very first element is 
     * overwritten, since it has number zero. */
    assert(grid->points[ip][jp][kp]->associatedElements[elementType]
	   [grid->points[ip][jp][kp]->nAssociatedElements[elementType]] == 0);

    /* Store the source number and update the number of elements 
     * associated with this point */
    grid->points[ip][jp][kp]->associatedElements[elementType]
      [grid->points[ip][jp][kp]->nAssociatedElements[elementType]++] = iae;

    /* If the sources and the evals are exactly the same, 
     * then this routine will be called only to set the sources.
     * Thus let in this case the evals be equal to the sources. */
    if (!grid->sourceAndEvalDiffer){
      /* In this case the element type should be "source": */
      assert(elementType==0);

      grid->points[ip][jp][kp]->associatedElements[1] = 
	grid->points[ip][jp][kp]->associatedElements[0];
      grid->points[ip][jp][kp]->nAssociatedElements[1] = 
	grid->points[ip][jp][kp]->nAssociatedElements[0];
    }

  } /* Next element  */
  if (verbose) printf(" %s Done \n",fctName);
  return;
}

/* 
* ================================================================= 
* Given a point and a radius deduce whether or not the sphere is 
* fully withing the extend of the grid. 
* If fully within grid return 1, otherwise return 0;
* ================================================================= */ 
static
int sphereIsFullyInsideGrid(point center,
			    double radius,
			    struct grid *grid) 
{
  if (center[0]-radius < grid->xmin || /* x-direction: */
      center[0]+radius > grid->xmax ||
      center[1]-radius < grid->ymin || /* y-direction: */
      center[1]+radius > grid->ymax ||
      center[2]-radius < grid->zmin || /* z-direction: */
      center[2]+radius > grid->zmax){
      /* The sphere extends beyond the grid. */
      return 0;
    }
  return 1;
}

/* 
* ================================================================= 
* For a given point in space, find the indices of the nearest grid 
* point. Report error and halt program if the point is outside the 
* extent of the grid.
* ================================================================= */
void closestGridPoint(int *i,int *j,int *k,point x,
		      struct grid *grid)
{
  char fctName[] = "closestGridPoint";
  if (x[0] < grid->xmin || x[0] > grid->xmax || 
      x[1] < grid->ymin || x[1] > grid->ymax ||
      x[2] < grid->zmin || x[2] > grid->zmax){
    /* This may be very bad. Report error and halt program. */
    fprintf(stderr,"ERROR: In %s:\n"
	    "   Queried for grid position outside of grid range:\n"
	    "   (%g,%g,%g)\n" 
	    "   This error is terminal\n",
	    fctName,x[0],x[1],x[2]) ;
    exit(1);
  }
  *i = (int) dblRound( (x[0]-grid->xmin)/grid->dx );
  *j = (int) dblRound( (x[1]-grid->ymin)/grid->dy );
  *k = (int) dblRound( (x[2]-grid->zmin)/grid->dz );
  return;
}

/* 
* ================================================================= 
* Given (i,j,k) indices of a grid point, return the position
* of the point in space. Report error and halt program if the 
* indices does not represent a grid point. 
* ================================================================= */
static
void gridPointPosition(int i,int j,int k,point x,
		      struct grid *grid)
{
  char fctName[] = "gridPointPosition";
  if (i < 0 || i > grid->nx - 1 || 
      j < 0 || j > grid->ny - 1 ||
      k < 0 || k > grid->nz - 1){
    /* This may be very bad. Report error and halt program. */
    fprintf(stderr,"ERROR: In %s:\n"
	    "  Queried for position of non-existing grid point:\n"
	    "  [%d][%d][%d]\n"
	    "  This error is terminal.",
	    fctName,i,j,k) ;
    exit(1);
  }
  x[0] = grid->xmin + i*grid->dx;
  x[1] = grid->ymin + j*grid->dy;
  x[2] = grid->zmin + k*grid->dz;
  return;
}

/* 
* ================================================================= 
* Given the indices of two grid points determine if they are 
* neighbours based on whatever criterion we choose to use.
* Return 1 for neighbours and 0 otherwise.
* ================================================================= */ 
static
int gridPointsAreNeighbours(int i1,int j1,int k1, /* point 1 */
			    int i2,int j2,int k2, /* point 2 */
			    struct gridStencil *stencil) 
{
  int i;

  for (i=0; i<stencil->nStencilPoints; i++){
    if (i1-i2 == stencil->ipNeighbours[i] &&
	j1-j2 == stencil->jpNeighbours[i] &&
	k1-k2 == stencil->kpNeighbours[i]){
      return 1;
    }
  }
  return 0;
}

/* 
* ================================================================= 
* Given two elements of different types, add them to each others 
* temporary lists of neighbours. 
* ================================================================= */ 
static 
void addElementsToBothTempNeighLists(
	  int iSource,
	  int iEval,
	  struct neighbourLinkedList **sourceLinkedList,
	  struct neighbourLinkedList **evalLinkedList)
{ 
  struct neighbourLinkedList *newLink;

  /* Put the eval into the list of this source */
  /*   Make a new link for the list */
  newLink =
    (struct neighbourLinkedList*)
    calloc(1,sizeof(struct neighbourLinkedList));
  /*   Put data in the new link */
  newLink->iNeighbour    = iEval;  /* Neighbour */
  newLink->nextNeighbour = 
    *sourceLinkedList;             /* Old "newest link" */
  /*   Let this be the newest link in the list */
  *sourceLinkedList = newLink;
  newLink = NULL; /* This statement is not really necessary, 
		   * but more of a "play it safe" kind'a thing. */

  /* Put the source into the list of the eval: */
  /*   Make a new link for the list */
  newLink =
    (struct neighbourLinkedList*)
    calloc(1,sizeof(struct neighbourLinkedList));
  /*   Put data in the new link */
  newLink->iNeighbour    = iSource;  /* Neighbour */
  newLink->nextNeighbour = 
    *evalLinkedList;                 /* Old "newest link" */
  /*   Let this be the newest link in the list */
  *evalLinkedList = newLink;
  newLink = NULL; /* This statement is not really necessary, 
		   * but more of a "play it safe" kind'a thing. */
  return;
}

/* 
* ================================================================= 
* Sort the neighbourlists of all grid elements of type elementType
* and eliminate dublicate entries.
* ================================================================= */ 
static 
void sortGridElementNeighbourLists(int elementType,
				   struct grid *grid,
				   struct gridStencil *stencil)
{
  int ige,          /* Counter for grid elements */

    nneigh,         /* Shorthand for number of neighbours */
    in,             /* Counter over neighbours */
    oldNeigh,nneighBak,

    *sortArray,     /* Used as temporary storage when sorting the 
		       neighbour lists */
    nSortArray;     /* Size of sortArray. */

  struct neighbourLinkedList 
    *thisLink, *nextLink;     /* Used for accessing the linked lists   */

  /* Allocate initial storage for sorting array.
   * If too little storage is allocated here, then 
   * it will be reallocated later. */
  nSortArray = 1000;
  sortArray  = (int*) calloc(nSortArray,sizeof(int));

  /* For all grid elements of this type: */
  for (ige=0; ige < grid->nGridElements[elementType]; ige++){

    /* If there are no neighbours, skip to next element */
    if (stencil
	->gridElementNeighbours[elementType][ige]
	->neighbourLinkedList == NULL)
      continue;

    /* Count the numbers of neighbours including dublicates: */
    stencil->gridElementNeighbours[elementType][ige]->nNeighbours = 0;
    thisLink = 
      stencil
      ->gridElementNeighbours[elementType][ige]
      ->neighbourLinkedList;
    while(thisLink!=NULL){
      stencil->gridElementNeighbours[elementType][ige]->nNeighbours++;
      thisLink = thisLink->nextNeighbour;
    }

    /* For convenience, copy the number of neighbour into "nneigh": */
    nneigh = stencil->gridElementNeighbours[elementType][ige]->nNeighbours;

    /* Check if enough storage is allocated.
     * If there is not enough storage assigned for the list, 
     * then free memory and assign more storage. 
     * Make sure the storage is increased significantly 
     * (at least a factor of two), such that this increase 
     * won't be necessary too many times (which might be a
     * bad situation for reclaiming memory - also it eats
     * up system time). */
    if(nneigh > nSortArray){
      FREE(sortArray);
      nSortArray = MAX(2*nSortArray,nneigh);
      sortArray  =(int*) calloc(nSortArray,sizeof(int));
    }

    /* Transfer neighbours to integer array for sorting */
    thisLink = stencil
                 ->gridElementNeighbours[elementType][ige]
                   ->neighbourLinkedList;
    for (in=0; in<nneigh; in++){
      sortArray[in]=thisLink->iNeighbour;
      thisLink = thisLink->nextNeighbour;      
    }

    /* Sort array of neighbours */
    intSort(sortArray,nneigh);

    /* Copy back to neighbour linked list ignoring dublicates */
    /* First entry: */
    thisLink = stencil
                 ->gridElementNeighbours[elementType][ige]
                   ->neighbourLinkedList;
    thisLink->iNeighbour = sortArray[0];

    /* Remaining entries: */
    oldNeigh  = sortArray[0];
    nneighBak = nneigh;
    nneigh    = 1;
    for (in=1; in<nneighBak; in++){
      if (sortArray[in] != sortArray[in-1]){
	/* Include in list */
	thisLink = thisLink->nextNeighbour;
	thisLink->iNeighbour = sortArray[in];
	nneigh++;
      }
    }
    /* Store found number of distinct links */
    stencil
      ->gridElementNeighbours[elementType][ige]
        ->nNeighbours = nneigh;
    /* Store first link to be freed */
    nextLink = thisLink->nextNeighbour;
    /* Terminate the last kept link in the list */
    thisLink->nextNeighbour = NULL;
    /* Free remaining entries of the neigbour links */
    for (in=nneigh; in<nneighBak; in++){
      thisLink = nextLink;
      assert(thisLink!=NULL);
      nextLink = thisLink->nextNeighbour;
      FREE(thisLink);
    }
    assert(nextLink==NULL);
  }

  FREE(sortArray);
  return;
} /* End routine sortGridElementNeighbourLists */

/* 
* ================================================================= 
* Given a grid point (indices) contruct a list containing the 
* grid-indices of the evals which are neighbours "thorugh the
* grid", i.e. using the grid neighbour stencil.
* 
* If the supplied array is too small to store the indices, then 
* increase its size and return alos the new size of the array.
*
* The original version of this routine returned a list of 
* generic pointers (to the elements - not the grid elements).
* Thus, in the original approach, the calling grid routine could 
* not itself use the list - it could only pass it on to an element-
* routine, which in turn could use the list.
* In the present approach, the calling routine must then construct
* a list of generic pointers to the elements if needed. This may be
* more expensive (in CPU time and memory), but allows much more 
* flexible (and intelligent) handling of the lists. 
*
* Store also the position in the grid stencil with which each 
* neighbour is associated.
* ================================================================= */ 
static 
void makeGridNeighbourIEvalList(int ip, int jp, int kp,
				int nNeighbours,
				int *iNeighbours,
				int *jNeighbours,
				int *kNeighbours,
				int **iEvalsIn,
				int **iEvalsDirectLocationIn,
				int *nEvals, int *nMaxEvals,
				struct grid *grid,
				int verbose)
{
  char fctName[]="makeGridNeighbourIEvalList";
  int 
    ignp,              /* Grid neighbour point counter                 */
    ipp,jpp,kpp,       /* Grid neighbour point indices                 */
    iae,               /* Element counters                             */
    nNeighMax,         /* Maximum number of neighbours (= *nMaxEvals)  */
    nNeighboursGrid;   /* Number of neighbouring elements reached 
			* through the grid (stencil)  (= *nEvals)      */

                       /* The following is a short-hand to the grid     */
  gridPoint ****gridPoints;
  int *iEvals, *iEvalsDirectLocation;

  /* "Shorthands" */
  nNeighMax = *nMaxEvals;
  gridPoints = grid->points;
  iEvals = *iEvalsIn;
  iEvalsDirectLocation = *iEvalsDirectLocationIn;

  nNeighboursGrid=0;
  for (ignp=0; ignp<nNeighbours; ignp++){ /* For each point in the stencil */
    /* Find the number of evaluation element associated 
     * with neighbouring grid points */
    ipp = ip + iNeighbours[ignp];
    jpp = jp + jNeighbours[ignp];
    kpp = kp + kNeighbours[ignp];

    /* Make sure that:
     * 1) These indices represents a grid point (use macro)
     * 2) AND that the grid point has associated elements
     * 3) AND that some of these elements are evals.
     * If ANY of these conditions are not met, continue with the 
     * next grid negighbour point.                                     */
    if (!gridPointExists(ipp,jpp,kpp,grid) || 
	gridPoints[ipp][jpp][kpp]==NULL    ||
	gridPoints[ipp][jpp][kpp]->nAssociatedElements[1]==0)
      continue; /* Skip to next grid point neighbour */

    /* Check for sufficient storage to put these evals 
     * into the array of neighbours */
    if (nNeighboursGrid
	+gridPoints[ipp][jpp][kpp]->nAssociatedElements[1]
	> nNeighMax){
	    /* Make sure the storage is increased significantly 
	     * (at least a factor of two), such that this increase 
	     * won't be necessary too many times (hopefully not at all). */
      nNeighMax = 
	MAX(2*nNeighMax,
	    nNeighboursGrid
	    +gridPoints[ipp][jpp][kpp]->nAssociatedElements[1]);

      if (verbose)
	printf("%s: Increasing array size %d, %d, %d\n",
	       fctName,nNeighMax,nNeighboursGrid,
	       gridPoints[ipp][jpp][kpp]->nAssociatedElements[1]);

      iEvals = (int*) realloc(iEvals,nNeighMax*sizeof(int));
      iEvalsDirectLocation = 
	(int*) realloc(iEvalsDirectLocation,nNeighMax*sizeof(int));
      *iEvalsIn = iEvals;
      *iEvalsDirectLocationIn = iEvalsDirectLocation;
    }
    /* Add the indices of these evals to the list */
    /*fprintf(stderr,"%s: Adding grid point (%d,%d,%d), evals:\n",
	    fctName,ipp,jpp,kpp);*/
    for (iae=0; 
	 iae<gridPoints[ipp][jpp][kpp]->nAssociatedElements[1]; 
	 iae++){
	    /*fprintf(stderr,"%s: %d\n",
		  fctName,
		  gridPoints[ipp][jpp][kpp]->associatedElements[1][iae]);*/
      iEvals[nNeighboursGrid+iae] = 
	gridPoints[ipp][jpp][kpp]->associatedElements[1][iae];
      iEvalsDirectLocation[nNeighboursGrid+iae] = ignp;
    }
    /* Step the counter of neighbouring evals                    */
    nNeighboursGrid += 
      gridPoints[ipp][jpp][kpp]->nAssociatedElements[1];

  } /* Next grid point neighbour in stencil                      */

  /* TEST - see if the list is OK:
   * (n^2 test - but small n)                         
  {
    int i,j;
    fprintf(stdout,"%s: Testing for dublicate entries...\n",
	    fctName);
    for (i=0;i<nNeighboursGrid-1;i++){
      for (j=i+1;j<nNeighboursGrid;j++){
	if (evals[j]==evals[i]){
	  fprintf(stdout,"%s: Test failed, %d and %d are the same!\n",
		  fctName,i,j);
	  exit(0);
	}
      }
    }
    fprintf(stdout,"%s: Test passed.\n",
	    fctName);
    exit(0);
  }
  */ 
  /* Store return values: */
  *nMaxEvals = nNeighMax;
  *nEvals    = nNeighboursGrid;
  return;
} /* End routine makeGridNeighbourIEvalList                            */



/*
* ================================================================= 
* Examine whether or not (ip,jp,kp) is the indices of an existing 
* grid point (return 1), of if is beyond the extend of the grid 
* (return 0).
*
* At the present this feature is handled by a macro named 
* gridPointExists
* ================================================================= */
static
int gridPointExistsOld(int ip, int jp,int kp,struct grid *grid)
{
  if (ip>=0 && ip<grid->nx &&
      jp>=0 && jp<grid->ny &&
      kp>=0 && kp<grid->nz  ) return 1;
  return 0;
}

/*
* ==================================================================== 
* Sets up the values of each function in the interpolation(/projection) 
* basis in each of the points of the interpolation(/projection) stencil.
* Then find the (pseudo-)inverse of this matrix to relate 
* ==================================================================== */ 
static 
void polyInterp(int Type, /*int nPoints,*/
	       double dx, double dy, double dz, 
	       int refreshGridSizeOnly, 
	       struct gridStencil *stencil)
{
  char fctName[] = "polyInterp";
  int nBasis, n, i,j,k, ip,ib, *iIndex,*jIndex,*kIndex, index, rank, dist2;
  double **tempWts, x,y,z;
  int consistentPoly, polyCubed, specialPoly, scaleEquations;
  double *eqScales;
  /* For shorthands to the stencil: */
  int polyOrder, nPoints, maxMonoOrder1D;
  int *ipNeighbours, *jpNeighbours, *kpNeighbours;
  double **weights, **values, polyScale;

  assert(Type >= 100);
  assert(dx == dy);
  assert(dy == dz);

  /* If dx is zero then this is only a "test" call, and there is no 
   * point in setting up anything fancy (or anything at all, actually, 
   * since it would be garbage computations). Even so, the correct 
   * amount of memory should be allocated. Since this routine does 
   * not contribute significantly to the total memory consumption, 
   * most of the calculations will be completed in either case (however, 
   * the arising matrix will not be passed to the pinv-routine).       */

  /* Flags for "polynomials of consistent orders", "n'th-order cubed"
   * and "other" polynomials. 
   * The polynomials to be used are for convenience divided into the 
   * the following groups:
   *  1) "consistentPoly" Polynomials containing all terms up to a 
   *     specified order and no higher-order terms. For example the 
   *     following polynomial is of "consistent order 1": 1 + x + y + z
   *  2) "polyCubed" Tensor product of three one-dimensional polynomials
   *     of specified order. This type of polynomial contains "spurious"
   *     higher-order cross terms such as "x*y*z". The following 
   *     polynomial is a "polyCubed" of order 1:
   *     (1+x)*(1+y)*(1+z)= 1 + x + y + z + x*y + x*z + y*z + x*y*z
   *  3) "specialPoly" any polynomial that does not fit into any of the 
   *     above categories. Typically, hard-coding is needed for the 
   *     individual polynomials in this case.
   * Initiate flags to zero, since we don't yet know which 
   * polynomial type to use.                                           */
  consistentPoly = 0;
  polyCubed      = 0;
  specialPoly    = 0;
  /* Flag for whether or not the equations should be scaled. 
   * [affects the accuracy for overdetermined systems.]                */
  scaleEquations = 0;

  /* Setup shorthands to stencil                                       */
  nPoints = stencil->nStencilPoints;
  ipNeighbours = stencil->ipNeighbours;
  jpNeighbours = stencil->jpNeighbours;
  kpNeighbours = stencil->kpNeighbours;

  /* We consider only a few different choises of basis functions:      */
  if (nPoints==1){
    /* Nearest-point constant-basis interpolation/projection. 
     * This is cheap in storage and compute time, but errors are high
     * (unless the direct part is really large - in which case this 
     * scheme is inefficient anyway).                                  */
    nBasis = 1;
    polyCubed = 1; /* Could also be "consistentPoly = 1". Shouldn't matter */
    maxMonoOrder1D = 0;
  }
  else if (nPoints==7){
    /* Seven points positioned on the main axis directions at the 
     * origin (the central stencil point) and at positions 
     * $\Delta (\pm 1,\pm 1,\pm 1)$. 
     * Linear basis interpolation functions with some second-order 
     * components.
     * This is fairly cheap in storage and compute time, but errors 
     * are expected to be fairly high compared to more refined schemes. 
     * The basis is the monomials in 1 + x + x^2 + y + y^2 + z + z^2 
     * Thus, this scheme does not retain the cross terms (xy, xz, yz),
     * while the "polyCubed" schemes retains the cross terms, but
     * ignores other higher-order terms.                               */
    nBasis = 7;
    specialPoly    = 1;
    maxMonoOrder1D = 2;
  }
  else if (nPoints==8){
    /* Eight points positioned at the corners of a 3x3x3 cube. 
     * Linear basis interpolation functions - i.e. a "first-order cubed"
     * polynomial. 
     * This is fairly cheap in storage and compute time, but errors 
     * also fairly high compared to more refined schemes. 
     * The basis is the monomials in (1+x)*(1+y)*(1+z)                 */
    nBasis = 8;
    polyCubed = 1; 
    maxMonoOrder1D = 1;
  }
  else if (nPoints==27){
    /* 3x3x3 stencil assumed. Use a "second-order cubed" polynomial, 
     * i.e. the basis is the monomials in (1+x+x^2)*(1+y+y^2)*(1+z+z^2)*/
    nBasis = 27;
    polyCubed = 1; 
    maxMonoOrder1D = 2;
  }
  else if (nPoints==33){
    /* 33 nearest points assumed. Use a consistently 
     * third-order polynomial, i.e. the basis is all the 
     * monomials off order less than four (20 of these).               */
    nBasis = 20;
    consistentPoly = 1;
    maxMonoOrder1D = 3;
    /* For this stencil we want to scale the equations to reduce the 
     * error as much as possible. [Empirical results for 1/r kernel 
     * used for determining the scaling coefficients]                  */
    scaleEquations = 1;
    /* Zero out the scaling used for this overdetermined system
     * of equations                                                    */
     /*for (i=0; i<nEqScales; i++) eqScales[i] = 0.0;
    Set scaling coefficients. The scaling of an equation related 
     * to a (interpolation/projection) point at distance sqrt(i) is 
     * stored in eqScales[i]                                           
    eqScales[0] = 1.00;
    eqScales[1] = 0.39;
    eqScales[2] = 0.22;
    eqScales[3] = 0.22;
    eqScales[4] = 0.02;*/
  }
  else {
    /* This is unknown = bad. Complain and exit!                       */
    printf("%s: Error!\n"
	   "       What kind of interpolation/projection\n"
	   "       stencil uses %d points?\n"
	   "       This halts the program\n",fctName,nPoints);
    _EXIT_;
  }

  /* Check that one and only one polynomial type has been chosen:      */
  if (consistentPoly+polyCubed+specialPoly!=1){
    fprintf(stderr,"%s: ERROR choosing polynomail types.\n"
	    "   Coding is needed to correct this. (%d,%d,%d) \n",
	    fctName,consistentPoly,polyCubed,specialPoly);
    _EXIT_;
  }

  /* Allocate memory needed */
  n = MAX(nBasis,nPoints);
  tempWts = (double **) calloc(n, sizeof(double *));
  for (i=0;i<n;i++)
    tempWts[i] = (double *) calloc(n, sizeof(double));

  eqScales = (double *) calloc(nPoints, sizeof(double));


  if (refreshGridSizeOnly){
    /* This is "just" a refresh, using a new grid size. 
     * Memory has already been allocated.                              */
    assert(stencil->bValues!=NULL);
    assert(stencil->weights!=NULL);
    assert(stencil->imonoOrder!=NULL);
    assert(stencil->jmonoOrder!=NULL);
    assert(stencil->kmonoOrder!=NULL);
  }
  else {
    /* This is *not* just a refresh, so we need to allocate new memory */
    /* Make sure that memory has not already been allocated
     * (simple check to avoid memory leaks)                            */
    assert(stencil->bValues==NULL);
    assert(stencil->weights==NULL);
    assert(stencil->imonoOrder==NULL);
    assert(stencil->jmonoOrder==NULL);
    assert(stencil->kmonoOrder==NULL);

    /* Values are one row per point and one column per basis function  */
    stencil->bValues = (double **) calloc(nPoints, sizeof(double *));
    for (i=0;i<nPoints;i++)
      stencil->bValues[i] = (double *) calloc(nBasis, sizeof(double));
    
    /* Weights are one row per basis function and one column per point   */
    stencil->weights = (double **) calloc(nBasis, sizeof(double *));
    for (i=0;i<nBasis;i++)
      stencil->weights[i] = (double *) calloc(nPoints, sizeof(double));

    /* Powers of each interpolation/projection basis function 
     * (assuming they are monomials)                                   */
    stencil->imonoOrder = (int *) calloc(nBasis,sizeof(int));
    stencil->jmonoOrder = (int *) calloc(nBasis,sizeof(int));
    stencil->kmonoOrder = (int *) calloc(nBasis,sizeof(int));
  }

  /* Set local shorthands to stencil                                   */
  values   = stencil->bValues;
  weights  = stencil->weights;
  iIndex   = stencil->imonoOrder;
  jIndex   = stencil->jmonoOrder;
  kIndex   = stencil->kmonoOrder;

  /* Set scaling coeficients */
  if (scaleEquations){
    /* Set the scaling coefficients based on the selected scheme       */
    if (nPoints ==33 && nBasis==20){
      fprintf(stdout,"%s: Using scaling of the equations\n",fctName);
      for (ip=0; ip<nPoints; ip++){
	dist2 =   ipNeighbours[ip]*ipNeighbours[ip] 
	        + jpNeighbours[ip]*jpNeighbours[ip] 
	        + kpNeighbours[ip]*kpNeighbours[ip] ;
	switch (dist2){
	case 0: 
	  eqScales[ip] = 1.00;
	  break;
	case 1: 
	  eqScales[ip] = 0.39;
	  break;
	case 2: 
	  eqScales[ip] = 0.22;
	  break;
	case 3: 
	  eqScales[ip] = 0.22;
	  break;
	case 4: 
	  eqScales[ip] = 0.02;
	  break;
	}
      } /* Set next scaling coefficient */
    } /* End of 33-point scaling coefficients */
    /* More schemes may be added here */
    else { /* This is unknown. Complain and exit */
      fprintf(stderr,"%s: Error! Should a scheme with %d points and \n"
	      "     %d basis functions require scaling of the equations?\n",
	      fctName,nPoints,nBasis);
      _EXIT_;
    }
  }
  else{ /* Set all the coefficients to unity */
    for (i=0; i<nPoints; i++) eqScales[i] = 1.0;
  }

  /* Scale of the monomial variables:                                  */
  if (dx != 0.0) polyScale = 1.0/dx;
  else polyScale = 0.0;

  /* Setup of basis function values */
  if (polyCubed){ /* "n'th-order cubed" polynomial, i.e. tensor product
		   * of three n'th order polynomials. Total polynomial
		   * is then of order 3*n, but contains all terms only 
		   * to order n.                                       */
    polyOrder = 3*maxMonoOrder1D;
    /* Monomial indices. Note that this is a definition, which 
     * could well be made differently. It may be possible to gain 
     * some speed by choosing these indices and the panel-moment 
     * layout consistently (but it should not be necessary for the 
     * calculations to be correct!).                                   */
    index=0;
    for (i=0;i<=maxMonoOrder1D;i++){
      for (j=0;j<=maxMonoOrder1D;j++){
	for (k=0;k<=maxMonoOrder1D;k++){
	  iIndex[index]=i;
	  jIndex[index]=j;
	  kIndex[index]=k;
	  index++;
	}
      }
    }
  }
  else if (consistentPoly){ 
    /* Consistently n'th-order polynomial, Total polynomial
     * is then of order n, and contains all terms to order n.          */
    polyOrder = maxMonoOrder1D;
    /* Monomial indices. Note that this is a definition, which 
     * could well be made differently. It may be possible to gain 
     * some speed by choosing these indices and the panel-moment 
     * layout consistently (but it should not be necessary for the 
     * calculations to be correct!).                                   */
    index=0;
    for (i=0;i<=maxMonoOrder1D;i++){
      for (j=0;j<=maxMonoOrder1D-i;j++){
	for (k=0;k<=maxMonoOrder1D-i-j;k++){
	  iIndex[index]=i;
	  jIndex[index]=j;
	  kIndex[index]=k;
	  index++;
	}
      }
    }
  }
  else if (specialPoly){ 
    if (nPoints ==7 && nBasis==7 && maxMonoOrder1D==2){
      polyOrder = maxMonoOrder1D;
      /* Hard coded monomials in the following order:
       *              1, x, x^2, y, y^2, z, z^2                        */
      iIndex[0]=0; jIndex[0]=0; kIndex[0]=0;
      iIndex[1]=1; jIndex[1]=0; kIndex[1]=0;
      iIndex[2]=2; jIndex[2]=0; kIndex[2]=0;
      iIndex[3]=0; jIndex[3]=1; kIndex[3]=0;
      iIndex[4]=0; jIndex[4]=2; kIndex[4]=0;
      iIndex[5]=0; jIndex[5]=0; kIndex[5]=1;
      iIndex[6]=0; jIndex[6]=0; kIndex[6]=2;
      index=7; /* For assertion purposes */
    }
    else { /* This is unknown. Complain and exit */
      fprintf(stderr,"%s: Error! What special polynomial of order %d\n"
	      "has %d points and %d basis functions?\n",
	      fctName,maxMonoOrder1D,nPoints,nBasis);
      _EXIT_;
    }
  }
  else { /* This should never be reached. Make sure that all possible 
	  * cases are included above.                                  */
    fprintf(stderr,"%s: ERROR. What type of basis function is this?\n"
	    "   Some coding may be needed to correct this error.\n",
	    fctName);
    _EXIT_;
  }

  assert(index==nBasis);

  /* Basis function values. The monomials are chosen for the basis, 
   * so the definition of the monomial ordering is used here.          */
  /* Unfortunately, the (linpack) SVD used does not give a decent 
   * result when dx is small. Testing in Matlab shows that the matrix 
   * starts to lose rank (numerically, not analytically) for dx<0.005.
   * At this point the result of the present computations would 
   * contain errors of RELATIVE magnitude one. The same problem occurs 
   * when dx becomes large. Obviously, the problem is a poor scaling 
   * of the columns of the "values" matrix. As a consequence a scaling 
   * of the values matrix is needed to ensure that the results are 
   * decent for all values of dx.
   * At the present (c)LAPACK SVD is used rather than the older linpack 
   * version. Even so, the scaling is maintained to keep the accuracy 
   * as high as possible.                                              */
  for (ip=0; ip<nPoints; ip++){
    /* Pick a point */
    x = dx*ipNeighbours[ip]*polyScale;
    y = dy*jpNeighbours[ip]*polyScale;
    z = dz*kpNeighbours[ip]*polyScale;
    for (ib=0; ib<nBasis; ib++){
      /* Calculate monomial values at this point                       */
      tempWts[ip][ib] = 
	pow(x,(double) iIndex[ib])
	*pow(y,(double) jIndex[ib])
	*pow(z,(double) kIndex[ib]);
    }
  }

  /* Scale equations                                                   */
  for (ip=0; ip<nPoints; ip++){
    for (ib=0; ib<nBasis; ib++) tempWts[ip][ib] *= eqScales[ip];
  } /* Next equation to scale */

  /* Store the matrix in "values" for later use (if any)               */
  for (ip=0; ip<nPoints; ip++){
    for (ib=0; ib<nBasis; ib++){
      values[ip][ib] = tempWts[ip][ib];
    }
  }

  /* Transform values to weights by inverting the matrix               */
  if (dx != 0.0) {
    /*rank = pinv_old(tempWts,nPoints,nBasis);*/
    rank = pinv(tempWts,nPoints,nBasis);

    /* Test if the system has full rank. If it doesn't, then complain 
     * and exit. Maybe someone would liek to use systems with rank 
     * deficit, but that is not the case presently.                      */
    if (rank!=nBasis){
      printf("%s: ERROR. The rank is estimated to %d. It should be %d.\n"
	     "    This halts the program.\n",fctName,rank,nBasis);
      _EXIT_;
    }
  }

  /* Stuff the inverted matrix into the weights array                  */
  for (j=0;j<nPoints;j++){
    for (i=0;i<nBasis;i++){
      weights[i][j] = tempWts[i][j];
    }
  }

  /* As a test, compute the matrix product (weights * values) and 
   * compare to the indentity matrix (of size nBasis x nBasis)         */
  if (dx != 0.0 ){ 
    double eDiag,eOffdiag;
    for (j=0;j<n;j++){
      for (i=0;i<n;i++){
	tempWts[i][j] = 0.0;
      }
    }
    for (j=0;j<nBasis;j++){
      for (i=0;i<nBasis;i++){
	for (ip=0;ip<nPoints; ip++){
	  tempWts[i][j] += weights[i][ip] * values[ip][j];
	}
      }
    }
    for (j=0,eDiag=0.0,eOffdiag=0.0;j<nBasis;j++){
      for (i=0;i<nBasis;i++){
	if (i==j) eDiag = MAX(eDiag, fabs(tempWts[i][j]-1.0));
	else eOffdiag = MAX(eOffdiag,tempWts[i][j]);
      }
    }
    /*printf("%s: Estimated errors on weights %8.2e %8.2e \n",
      fctName,eDiag,eOffdiag);*/
    if (eDiag>1.0e-6 ||eOffdiag >1.0e-6){
      printf("    Too large errors. Something must be broken!\n");
      _EXIT_;
    }
  }
  
  /* Rescale equations (column scaling of weight matrix) */
  for (j=0;j<nPoints;j++){
    for (i=0;i<nBasis;i++){
      weights[i][j] *= eqScales[j];
    }
  }
  /* Get rid of weights - we do not need to scale later! */
  for (j=0;j<nPoints;j++) eqScales[j] = 1.0;
  scaleEquations = 0;
  
  /* Store values in global array for the stencil (only needed for 
   * new stuff that is not already represented by a local pointer)     */
  stencil->nBasis         = nBasis;
  stencil->polyScale      = polyScale;
  stencil->polyOrder      = polyOrder;
  stencil->maxMonoOrder1D = maxMonoOrder1D;

  /* Free temporary storage */
  for (i=0;i<n;i++){FREE(tempWts[i]);}
  FREE(tempWts);
  FREE(eqScales);

  return;
} /* End of routine polyInterp */

/*
* ==================================================================== 
* Given the grid and a set of grid point indices (i,j,k), return the 
* coordinates (x,y,z) of the grid point.
* ==================================================================== */ 
void getGridLocation(void *gridIn,
		     int i, int j, int k,
		     double *px, double *py, double *pz)
{
  struct grid *grid;
 
  grid = (struct grid*) gridIn;

  *px = grid->x[i]; /*  i * grid->dx + grid->xmin;*/
  *py = grid->y[j]; /*  j * grid->dy + grid->ymin;*/
  *pz = grid->z[k]; /*  k * grid->dz + grid->zmin;*/

  return;
}

int gridDataIndex(void *gridIn,
		int i, int j, int k)
{
  struct grid *grid;
 
  grid = (struct grid*) gridIn;

  assert(grid->griddata!=NULL);

  return (dataIndex(grid->griddata, i, j, k));
  /*return (dataIndex(gridIn, i, j, k));*/
}

/*
* ==================================================================== 
* Returns the number of points in the interpolation stencil.
* ==================================================================== */ 
int numInterpPoints(void *gridInput)
{
  /* Define grid and recast input data                                 */
  struct grid *grid  = (struct grid*) gridInput;
  /* Return number of points in the interpolation stencil              */
  return grid->elementGridStencil->nStencilPoints;
}

/*
* ==================================================================== 
* Returns information for interpolation stencil.
* Number of points, their relative location and the weights of the 
* points.
* Note that pointers to the actual location of this information in 
* "grid" are returned. Thus, knowledge of this part of the grid 
* structure is needed in the calling routine. CAVEAT: this type of 
* coding could lead to "odd errors" if the structure of e.g. 
* "stencil->weights" is changed at a later point! 
* Also, if the pointers *piv, *pjv, *pjv, or **pwts have memory 
* associated with them, then this practice could lead to memory leaks.
* ==================================================================== */ 
int getGridInterp(void *gridIn, char *Type, 
		  int **piv, int **pjv, int **pkv, double ***pwts,
		  int **imonoOrder, int **jmonoOrder, int **kmonoOrder, 
		  int *pterms, double *polyScale,
		  int *polyOrder, int *maxOrder1D)
{
  /*char fctName[]="getGridInterp";*/
  struct grid *grid;
  struct gridStencil *stencil;

  /* At the moment sources (projection) and evals (interpolation) 
   * stencils are the same. Just check that the call is made 
   * consistently. Later, "Type" may decide which stencil 
   * will be used when returning data                                  */
  assert( strcmp(Type,"source")==0 || strcmp(Type,"eval")==0);

  /* Cast input to correct type */
  grid = (struct grid*) gridIn;
  stencil = grid->elementGridStencil;

  *piv        = stencil->ipNeighbours;
  *pjv        = stencil->jpNeighbours;
  *pkv        = stencil->kpNeighbours;
  *pwts       = stencil->weights;
  *imonoOrder = stencil->imonoOrder;
  *jmonoOrder = stencil->jmonoOrder;
  *kmonoOrder = stencil->kmonoOrder;
  *pterms     = stencil->nBasis;
  *polyScale  = stencil->polyScale;
  *polyOrder  = stencil->polyOrder;
  *maxOrder1D = stencil->maxMonoOrder1D;

  return stencil->nStencilPoints;
}

/*
* ==================================================================== 
* Project a vector (of source-element values) onto the grid. 
* This routine extracts the necessary data from the grid structure and 
* passes them on (to an nfft.c routine) for further processing, since
* the grid doesn't know what the data look like.
* ==================================================================== */ 
void projectOntoGrid(void *gridIn,double *src,void *projectMat)
{
  struct grid *grid;

  /* Cast input pointer to correct type */
  grid = (struct grid*) gridIn;

  /* Call a routine to do the projection                               */
  projectOntoGridData(grid->griddata,src,projectMat);

  return;
}

/*
* ==================================================================== 
* Interpolate the griddata values to a vector (of evaluation-element 
* values). 
* This routine extracts the necessary data from the grid structure and 
* passes them on (to an nfft.c routine) for further processing, since
* the grid doesn't know what the data look like.
* ==================================================================== */ 
void interpolateFromGrid(void *gridIn,double *dest,void *interpMat)
{
  struct grid *grid;

  /* Cast input pointer to correct type */
  grid = (struct grid*) gridIn;

  /* Call a routine to do the projection                               */
  interpolateFromGridData(grid->griddata,dest,interpMat);

  return;
} /* End of routine interpolateFromGrid */

/* 
* ==================================================================== 
* Find the sum and the sum of the squares of the grid element 
* bounding sphere radii. Stores values with the grid.
* Returns the average radius.
* ==================================================================== */ 
static 
double gridElementRadius(struct grid *grid)
{
  int ielem, eType, maxType;
  double sum1[2],sum2[2];

  /* Make sure that the existing counts of radius and radius squared 
   * are zero (i.e. that this is the first call) */
  assert(grid->elementRadiusSum[0]==0.0);
  assert(grid->elementRadiusSum[1]==0.0);
  assert(grid->elementRadiusSquareSum[0]==0.0);
  assert(grid->elementRadiusSquareSum[1]==0.0);

  if (grid->sourceAndEvalDiffer) maxType = 1;
  else maxType=0;

  for (eType=0;eType<=maxType;eType++){
    for (sum1[eType]=0.0,sum2[eType]=0.0,ielem=0; 
	 ielem< grid->nGridElements[eType];
	 ielem++){
      sum1[eType] += grid->gridElements[eType][ielem]->radius;
      sum2[eType] += grid->gridElements[eType][ielem]->radius
	            *grid->gridElements[eType][ielem]->radius;
    }
    
    /* Store with grid */
    grid->elementRadiusSum[eType]       = sum1[eType];
    grid->elementRadiusSquareSum[eType] = sum2[eType];
    /* Normalize sums for output (not stored with grid) */
    sum1[eType] /= (double) grid->nGridElements[eType];
  }
  /* Return maximum average radius */
  return (MAX(sum1[0],sum1[maxType]));
} /* End of routine gridElementRadius */

/*
* ==================================================================== 
* Based on memory consumption find an "optimum" grid size and layout.
*
* This routine / method needs to be improved. At the present it is 
* using WAY too much time. 
* ==================================================================== */
static
void findBasicGridValues(struct grid *grid,
			 int nCols,int nRows, int numMats)
{
/*# include <time.h> timing routines may eventually not be needed ? */
/*# define SETTIMER(time){time = ((double) clock())/((double) CLOCKS_PER_SEC);}*/
  char fctName[] = "findBasicGridValues";
  int verbose = 1; /* Level of verboseness. Use 0 for silent. 
		    * 1 for normal and 2 for very verbose.             */
  double gridMemory, directMemory, minMemory, eleSize;
  int i, refreshGridSizeOnly, calculateNextGrid;
  int maxAssocElements[2], nAssocGridPoints;
  int p2x,p2y,p2z, p2xBest,p2yBest,p2zBest;

  double goodDeal=3.0;   /* Factor that memory consumption must 
			  * increase from a local minimum before it will 
			  * be considered to be the global minimum.
			  * The computations may also be stoped by 
			  * other criteria (see below in code)         */

  int maxAssocs = 100;   /* Maximum # elements associated with a 
			  * single grid point                          */

  int accurateEstimate=0;/* 0: Not so accurate memory estimate for the 
			  * linked-list parts (some overestimation)
			  * 1: Accurate, but requires setup of linked 
			  * list (allocation and deallocation) for every 
			  * grid examined */

  int calculateMemory;   /* Set to zero for a particular grid, if it is 
			  * estimated that the grid is not reasonable, 
			  * e.g. too coarse or too fine.               */

  double interactVolume, /* An average element will "want to" interact
			  * with other elements (of opposite type) if 
			  * they are in this volume (approximation)    */
    stencilVolume,       /* Volume of the direct stencil (dx^3 times 
			  * number of stencil points)                  */
    rad;

  double dxOld;          /* dx for last grid for which memory was 
			  * calculated (set to negative initially)     */

  double gridVolume,     /* The volume of a grid, i.e.                 */
    gridVolumeLast;      /* (nx-1)*(nx-1)*(nx-1)*(delta^3)
			  * A grid should only be considered if the 
			  * volume has decreased compared to the last 
			  * grid with large delta.                     */

  double minGridToEleRatio=0.5;
                         /* If the grid size (grid->dx,dy,dz) becomes 
			  * small compared to an (average) element, then 
			  * most panel interactions will take place 
			  * through linked lists. To avoid this scenario 
			  * a minimum grid size (relative to an avarage 
			  * panel size) can be set. The grid size will 
			  * be compared to the average diameter of the 
			  * bounding spheres of the grid elements. 
			  * Setting this parameter to zero will enable 
			  * any grid size.                             */ 

  double maxGridToEleRatio=10.0;
                         /* As minGridToEleRatio except for maximum 
			  * grids size. Set to a really large number 
			  * (not as large as DBL_MAX, though ;) to 
			  * disable this test.                         */

  double minInteractVolumeRatio=0.5;
                         /* Minimum allowed value for the ratio between 
			  * the volume of the direct stencil for the 
			  * grid and the interaction volume of an average 
			  * element based on radius and 
			  * "stencil->separationDistanceFactor". Most of 
			  * the elements should have their interaction 
			  * range to fall fully inside the direct stencil. 
			  * This factor can be used to adjust how much 
			  * smaller the interaction volume of an average 
			  * element should be relative to the direct 
			  * stencil. The factor is used to determine when 
			  * a grid is "too fine".                      */

  /* Temporarily needed: 
  FILE *dumpFile;
  double tStart,tStop;*/

  /* So far no minimum has been found. Set the minimum to be 
   * *REALLY* large (larger than what would occur normally             */
  minMemory = DBL_MAX;

  /* So far no grid has been examine. Set the volume to be large       */
  gridVolume = DBL_MAX;
  gridVolumeLast = DBL_MAX;
  /*dumpFile = fopen("memory.out","w");*/

  /*printf("%s: Writing memory use to memory.out \n",fctName);
    dumpFile = fopen("memory.out","w");*/

  /* Set minimum numbers of grid nodes - as powers of two.
   * Note that the minimum should be chosen so that the interpolation 
   * and projection stencil(s) can be contained (with still a little 
   * slack to go).                                                     */
  p2x = p2y = p2z = 2;

  /* So far no grids have been examined                                */
  dxOld = -1.0;

  /* Initiate grid stencils - they will need recalculation when the 
   * grid size has been determined. Note that the grid itself has  
   * been initialized, but not values have been found. Thus the 
   * stencil setup has to not accept e.g. dx=0.0                       */
  setupGridStencils(grid, refreshGridSizeOnly=0);

  /* Set a measure for the element size (diameter) */
  eleSize = 2.0*MAX(grid->elementRadiusSum[0]/(double)(grid->nGridElements[0]),
		    grid->elementRadiusSum[1]/(double)(grid->nGridElements[1]));
  if (eleSize <= 0.0){
    fprintf(stderr,"%s: ERROR: Non-positive maximum element size found: %g\n"
	    "   This error is terminal\n",fctName,eleSize);
    _EXIT_;
  }

  /* Now that the stencils are set up we can find the "interaction 
   * volume" for an average element. This is the volume that is 
   * covered by its own interaction range, i.e. the volume of a 
   * sphere that is somewhat larger than an average bounding sphere.   */
  rad = 0.5*eleSize*(1.0+grid->directStencil->separationDistanceFactor);
  interactVolume = (4.0/3.0)*PI*pow(rad,3.0);

  /* Initiate flag and a loop counter */
  calculateNextGrid=1; i=0;
  while(calculateNextGrid){
    calculateMemory = 1; /* By default calculate memory for this grid  */

    /* Store the volume of the last grid examined:                     */
    gridVolumeLast = gridVolume;

    /* For tests only. 
     * Hardcode grid size (as power of two in each direction)          */
    if (0){
      calculateNextGrid=0;
      p2x = 9-1;
      p2y = 5;
      p2z = 4-1;
      fprintf(stderr,
	      "%s: WARNING! Using hardcoded # grid nodes!\n"
	      "   p2x=%d, p2y=%d, p2z=%d\n",fctName,p2x,p2y,p2z);
    }

    if (verbose>1) printf("%s: Loop no. %d.\n",fctName,i++);

    /* Setup basic grid values. DX, Extend, Stencil for neighbours...
     * These may require knowledge about the spatial distribution 
     * of the elements.                                                */

    /* Setup grid size and dimensions */
    setupGridDimensions(grid, 
			3,           /* Grid generation type           */
			p2x,p2y,p2z, /* # grid points (powers of two)  */
			-7,          /* Not used for type 3            */
			MAX(0,verbose-1)); /* Verboseness of called function */         

    if (verbose) printf("%s: Testing grid (nx,ny,nz)=(%d,%d,%d), DX=%g\n",
			fctName,grid->nx,grid->ny,grid->nz,
			grid->dx);
    /* Calculate the volume of the present grid:                       */
    gridVolume = grid->dx*grid->dy*grid->dz*((grid->nx-1)*(grid->ny-1)*(grid->nz-1));


    /* Check grid vs. various conditions... */
    if(grid->dx>maxGridToEleRatio*eleSize){
      /* This grid is too coarse. Don't examine it further.            */
      calculateMemory   = 0;
      if (verbose) 
	printf("%s: This grid is too coarse: DX/eleSize = %g. Skipping to next\n",
	       fctName,grid->dx/eleSize);
    }
    else if (grid->dx<minGridToEleRatio*eleSize && dxOld>0.0){
      /* This grid is too fine. Don't examine it further. 
       * Also, do not examine finer grids!                             */
      if (verbose) 
	printf("%s: This grid is too fine: DX/eleSize = %g.\n",
	       fctName,grid->dx/eleSize);
      calculateMemory = 0;
      calculateNextGrid = 0;
    }
    else if ( (gridMemory=gridMemoryEstimate(grid, 0, MAX(0,verbose-1))) 
	      >= minMemory){
    /* The grid memory alone (not including memory used on linked lists) 
     * on this grid is larger than the total memory on the optimum grid 
     * found so far. Thus, we can stop searching. As the grids become 
     * denser, the grid memory will only increase, and so the total 
     * memory will always be larger than what we have found so far.    */
      if (verbose) 
	printf("%s: This grid is too fine. Grid memory >= %gMb\n",
	       fctName,gridMemory/pow(2.0,20.0));
      calculateMemory   = 0;
      calculateNextGrid = 0;
    }
    /* See if the grid volume is smaller than the volume of the last 
     * grid. If not, then the present grid should not be considered 
     * for computations.                                               */
    else if (gridVolume>gridVolumeLast && dxOld>0.0){
      /* This grid is larger than it has to be. This is inefficient 
       * with respect to storage of the grid data and the time 
       * required for the convolution. 
       * Skip this grid size and move on unless NO grid has been 
       * found so far. In that case we'll better examine this grid 
       * for good measure.                                             */
      if (verbose) 
	printf("%s: This grid is disregarded based on volume size\n",
	       fctName);
      calculateMemory   = 0;
    }
    else { /* Continue grid setup */
      /* Setup grid stencils                                           */
      setupGridStencils(grid, refreshGridSizeOnly=1);

      /* Allocate memory for the association of elements with 
       * grid points.                                                  */
      allocateGridPoints(grid);

      /* Associate grid elements with grid points. 
       * ("Put elements onto the grid")                                */
      if (verbose>1) 
	printf("%s: CALLING associateElementsWithGridPoints\n",fctName);
      associateElementsWithGridPoints(grid);

      /* Find max numbers of elements per grid point                   */
      if (verbose>1) 
	printf("%s: Counting associations \n",fctName);
      nAssocGridPoints = 
	countGridElementAssociations(grid,maxAssocElements);

      /* Estimate memory only if the grid is not ludicrously coarse.   */
      if(maxAssocElements[0] > maxAssocs || 
	 maxAssocElements[1] > maxAssocs) {
	calculateMemory = 0;
	if (verbose) 
	  printf("%s: This grid is too coarse (too many elements per grid point)\n",
	       fctName);
      }

      /* See if this grid is reasonable with respect to grid size, 
       * direct stencil size and panel size (interaction range). 
       * MOST of the interactions should be using the stencil. If this 
       * is violated, then the run time could be severely hampered.    */
      /* Volume of the direct stencil */
      stencilVolume  = (grid->directStencil->nStencilPoints) 
	              * pow(grid->dx,3.0);
      /* If the interaction volume for an average element is not covered 
       * by the direct interaction stencil, then the grid is too fine in
       * the sence that too many interactions will need to be obtained 
       * through an interaction list. Do not attempt to calculated needed
       * memory, and do not find next grid (it will only be even finer).*/
      if (verbose>1){
	printf("%s: interaction volume %g; stencil volume %g\n",
	       fctName,interactVolume,stencilVolume);
	printf("%s: DX= %g; avg diameter = %g\n",
	       fctName,grid->dx,eleSize);
      }
      if (stencilVolume < minInteractVolumeRatio*interactVolume){
	if (verbose) 
	  printf("%s: This grid is estimated to be too fine\n"
		 "        based on the interaction volume %g\n"
                 "        and the stencil volume %g\n",
		 fctName,interactVolume,stencilVolume);
	calculateNextGrid=0;
	/* If there has not already been found a grid, then we HAVE to 
	 * use the present grid size, even though a lot of interactions 
	 * may be calculated (not too bad, since the problem is probably 
	 * "small" in that case).                                      */
	if(dxOld>0) calculateMemory = 0;
      }
    }/* End of setting up this grid. */

    /* Continue calculations only if the grid the memory needs to be 
     * estimated (not on very fine or very coarse grids)               */
    if( calculateMemory ){
      /* If we need an accurate estimate, then actually find neighbours, 
       * store them in linked lists and eliminate duplicate entries    */
      if (accurateEstimate){
	/* For the direct part:                                        */
	if (verbose>1) 
	  printf("%s: Setting up neighbours (direct part)\n",fctName);
	setupGridElementNeighbours(grid,grid->directStencil);
	/* For the preconditioning part:                               */
	if (verbose>1) 
	  printf("%s: Setting up neighbours (precondition part)\n",fctName);
	setupGridElementNeighbours(grid,grid->precondStencil);
      }
      /* Estimate memory used on the grid (include memory used 
       * for linked lists)                                             */
      gridMemory = gridMemoryEstimate(grid,1, MAX(0,verbose-1));
      /* Estimate memory used on direct part (incl. preconditioning)   */
      directMemory = directMemoryEstimate(grid, nCols,nRows,numMats);

      if(verbose)
	printf("%s: \n"
	       "    Memory estimate for direct interactions:  %9.3e Mb\n"
	       "    Memory estimate for grid:                 %9.3e Mb\n"
	       "    Total memory estimate:                    %9.3e Mb\n",
	       fctName,
	       directMemory/pow(2.0,20.0),
	       gridMemory/pow(2.0,20.0),
	       (gridMemory+directMemory)/pow(2.0,20.0));
      
      /* Set this value of dx to the most recently calulated */
      dxOld = grid->dx;

      /* Store result if memory is lower than previously best bet      */
      if ( (gridMemory+directMemory) < minMemory){
	minMemory = gridMemory+directMemory;
	p2xBest = p2x;
	p2yBest = p2y;
	p2zBest = p2z;
      }
      /* If the memory is "a great deal more" than the minimum found 
       * so far, then the memory consumption will probably increase from 
       * here as the grid is refined. No need to go on and examine yet 
       * finer grids.
       * Also, if the memory used on the grid is larger than the memory 
       * used on the direct problem (and if this is not a minimum), then 
       * we probably won't do any better by refining the grid further.*/
      if ( (gridMemory+directMemory> goodDeal*minMemory) || 
	   (gridMemory+directMemory> 1.001*minMemory && 
	    gridMemory>directMemory)    ){
	calculateNextGrid = 0;
      }
    }
    /* Check the "tightest fit" and set grid size for next grid to try 
     * (slightly denser grid)                                          */
    if      (grid->iDirTightFit == 0) p2x++;
    else if (grid->iDirTightFit == 1) p2y++;
    else if (grid->iDirTightFit == 2) p2z++;
    else {
      printf("%s: ERROR: What kind of direction is this? %d\n",
	     fctName,grid->iDirTightFit);
      printf("%s: This halts the program!\n",fctName);
      _EXIT_;
    }

    /* Free the grid for later resize                                  */
    freeGridForResize(grid);
  }

  /* Having chosen the grid, now actually set it up: 
   * [Comments for each step are above]                                */
  setupGridDimensions(grid, 3, p2xBest,p2yBest,p2zBest, -7, 
		      MAX(0,verbose-1));
  if (verbose) {
    printf("%s:\n\n",fctName);
    printf("  Choosing grid size: (nx,ny,nz)=(%d,%d,%d), DX=%g\n",
	   grid->nx,grid->ny,grid->nz,grid->dx);
    printf("  Estimated memory for direct interactions and grid data: %.1fMb\n",
	   minMemory/pow(2.0,20.0));
    printf("\n");
  }
  setupGridStencils(grid, refreshGridSizeOnly=1);
  allocateGridPoints(grid);
  associateElementsWithGridPoints(grid);
  setupGridElementNeighbours(grid,grid->directStencil);
  setupGridElementNeighbours(grid,grid->precondStencil);
  /* Setup near-field grid2grid values to speed up precorrection 
   * calculations. (Only one of these two should eventually be used.)  */
  printf("   Precorrection setup I:\n");
  setupNearFieldGrid2grid(grid, refreshGridSizeOnly=0);
  printf("   Precorrection setup II:\n");
  setupNearFieldGrid2gridUnion(grid, refreshGridSizeOnly=0);

  return;
}
/*
* ==================================================================== 
* Estimate the amount of memory that will be required to realize the 
* present grid.
* Double precision is used to represent the memory requirement 
* in bytes (this is to avoid problems if the total estimated memory 
* exceeds the maximum value of an integer 4G for a 4b unsigned).
*
* countLists: 0: Do not add memory for the linked lists of grid 
*                elements neighbours (may be expensive to count).
*             1: Estimate all memory use.
* ==================================================================== */
static
double gridMemoryEstimate(struct grid *grid, int countLists, int verbose)
{
  char fctName[]="gridMemoryEstimate";
  int nodes;
  double memCount, thisMem;

  nodes = grid->nx * grid->ny * grid->nz;
  memCount = 0.0;

  /* Memory needed for grid data: */
  thisMem   = griddataMemoryEstimate(grid->nx,grid->ny,grid->nz, MAX(0,verbose-1));
  memCount += thisMem;
  if (verbose) printf("%s: Estimated memory for grid data  %10.3e Mb\n",
		      fctName, thisMem/pow(2.0,20.0));

  /* Memory used by linked lists of grid elements */  
  if (countLists){
    thisMem   = gridElementListsMemoryEstimate(grid);
    memCount += thisMem;
    if (verbose) printf("%s: Estimated memory for link-lists %10.3e Mb\n",
			fctName, thisMem/pow(2.0,20.0));
  }
  return memCount;
}
/*
* ==================================================================== 
* Estimate the amount of memory that will be required to realize the 
* direct part (incl. preconditioner, projection and interpolation) 
* using the present grid.
* ==================================================================== */
static
double directMemoryEstimate(struct grid *grid,
			    int nCols,int nRows, int numMats)
{
  char fctName[] = "directMemoryEstimate";
  int nodes, iele, nUnknowns, nEntries, nCols2, numMats2;
  int verbose = 0;
  double memCount, thisMem;

  nodes = grid->nx * grid->ny * grid->nz;
  memCount = 0.0;
  /* These needs to know #entries in sparse matrices and #bytes/entry  */
  /* Memory needed for the direct part                                 */

  thisMem = directPartMemoryEstimate((void *) grid, 
				     nCols, nRows, numMats,
				     "direct");
  memCount += thisMem;
  if (verbose) printf("%s: Estimated memory for direct part  %10.3e Mb\n",
		      fctName, thisMem/pow(2.0,20.0));
  /* Memory needed for the preconditioner part                         */
  thisMem = directPartMemoryEstimate((void *) grid, 
				     nCols, nRows, numMats,
				     "precond");
  memCount += thisMem;
  if (verbose) printf("%s: Estimated memory for precondition %10.3e Mb\n",
		      fctName, thisMem/pow(2.0,20.0));

  /* Memory needed for the projection part                             */
  /* Count total number of unknowns (sources)                          */
  for (nUnknowns=0,iele=0;iele<grid->nGridElements[0];iele++){
    nUnknowns += numSourceUnknowns(grid->gridElements[0][iele]->element);
  }
  /* Multiply by the number of projection basis functions              */
  nEntries = nUnknowns*grid->elementGridStencil->nBasis;
  nCols2   = nUnknowns;
  thisMem  = spMatMemoryEstimate(nCols2, nEntries, numMats2 = 1);
  memCount += thisMem;
  if (verbose) printf("%s: Projection estimated to    %10.3e Mb\n",
		      fctName,thisMem/pow(2.0,20.0));
  /* Memory needed for the interpolation part                          */
  /* Count total number of unknowns (evals)                            */
  for (nUnknowns=0,iele=0;iele<grid->nGridElements[1];iele++){
    nUnknowns += numEvalUnknowns(grid->gridElements[1][iele]->element);
  }
  /* Multiply by the number of interpolation basis functions           */
  nEntries = nUnknowns*grid->elementGridStencil->nBasis;
  nCols2   = nUnknowns;
  thisMem  = spMatMemoryEstimate(nCols2, nEntries, numMats2 = 1);
  memCount += thisMem;
  if (verbose) printf("%s: Interpolation estimated to %10.3e Mb\n",
		      fctName,thisMem/pow(2.0,20.0));

  return memCount;
}

/* 
*===================================================================== 
* Count the number of linked lists used in the grid, and figure out 
* how much memory this is using.
*===================================================================== */ 
static
double gridElementListsMemoryEstimate(struct grid *grid)
{
  char fctName[] = "gridElementListsMemoryEstimate";
  struct gridStencil *stencil;
  int i, iele,nele, nlink, ntotLinks, returnZeroIfAssociated;
  char *Type, directType[]="direct", precondType[]="precond";
  int verbose = 1;

  ntotLinks=0;
  /* Count number of lists                                             */
  nele=grid->nGridElements[0];
  for (i=0;i<2;i++){
    if (i==0){
      stencil = grid->directStencil;
      Type    = directType;
    }
    else if (i==1){
      stencil = grid->precondStencil;
      Type    = precondType;
    }

    /* Count links:                                                    */
    if (stencil->gridElementNeighboursAreFound!=0) {
      if(stencil->nGridElementNeighboursAreFound ==0 ||
	 stencil->nGridElementNeighboursAreFound ==1){
	/* Count the links - they have not been counted properly yet   */
	for (iele=0,nlink=0; iele<nele; iele++){
	  if(stencil->gridElementNeighbours[0][iele]!=NULL){
	    nlink += stencil->gridElementNeighbours[0][iele]->nNeighbours;
	    /*printf(" %d has %d linked neighbours, now found %d\n",iele,
	      stencil->gridElementNeighbours[0][iele]->nNeighbours,nlink);*/
	  }
	}
	/* Store number of links in the stencil and set flag to show 
	 * that they have been properly counted                        */
	stencil->nGridElementNeighbours = nlink;
	stencil->nGridElementNeighboursAreFound = 2;
	if (verbose)
	  printf("%s: Found %d links in the %s part\n",fctName,nlink,Type);
      }
      else if (stencil->nGridElementNeighboursAreFound ==2){ 
	/* The links have already been counted. 
	 * Just go get the number from the previous count              */
	nlink = stencil->nGridElementNeighbours;
      }
      else { /* This is unknown = bad. Report error and exit! */
	printf("%s: ERROR! What kind of flag is this? \n"
	       "  nGridElementNeighboursAreFoundFound = %d in %s stencil\n"
	       "  (In file: %s, line # %d)\n"
	       "  This halts the program!\n",
	       fctName,stencil->nGridElementNeighboursAreFound, Type,
	       __FILE__,__LINE__);
	exit(1);
      }
    }
    else { /* stencil->gridElementNeighboursAreFound==0:
	    * Grid elementy linked lists have not been set up.         */
      if (stencil->nGridElementNeighboursAreFound==0){
	/* Count (estimate) # links */
	nlink = countPossibleNeighbourLinks((void*) grid,
					    returnZeroIfAssociated=0,
					    Type);
	/* [The number of links and the flag will be set in the stencil 
	 * by the routine that actually does the counting / estimating,
	 * in this case "countPossibleNeighbourLinks"]                 */
	if (verbose)
	  printf("%s: Estimated %d links in the %s part\n",
		 fctName,nlink,Type);
      }
      else if (stencil->nGridElementNeighboursAreFound==1){
	/* Get # links from previous estimate */
	nlink = stencil->nGridElementNeighbours;
      }
      else if (stencil->nGridElementNeighboursAreFound==2){
	/* This should never occur! Print error and exit */
	printf("%s: ERROR! How can the links in the %s stencil be \n"
	       "counted when they are not even set up?\n"
	       "  (In file: %s, line # %d)\n"
	       "  This halts the program!\n",
	       fctName,Type,__FILE__,__LINE__);
	exit(1);
      }
      else { /* This is unknown = bad. Report error and exit           */
	printf("%s: ERROR! What kind of flag is this? \n"
	       "  nGridElementNeighboursAreFoundFound = %d in %s stencil\n"
	       "  (In file: %s, line # %d)\n"
	       "  This halts the program!\n",
	       fctName,stencil->nGridElementNeighboursAreFound,Type,
	       __FILE__,__LINE__);
	exit(1);
      }
    }
    ntotLinks += nlink;
  }

  return (((double)ntotLinks) 
	  * ((double)sizeof(struct neighbourLinkedList)) );
}

/* 
* ==================================================================== 
* Free all parts of the grid that depends on the choise of grid size.
* It is assumed that some parts of the grid structure may not be 
* allocated. If so, then the respective pointers must point to NULL.
* ==================================================================== */
static
void freeGridForResize(struct grid *grid)
{
  int ip,jp,kp, iele, itype,typemax;

  FREE(grid->x);
  FREE(grid->y);
  FREE(grid->z);

  if (grid->points!=NULL){
      /* For every grid point                                          */
    for (ip=0; ip<grid->nx; ip++){
      for (jp=0; jp<grid->ny; jp++){
	for (kp=0; kp<grid->nz; kp++){
	  /* Free point if it is allocated */
	  IFFREE(grid->points[ip][jp][kp]);
	}
	FREE(grid->points[ip][jp]);
      }
      FREE(grid->points[ip]);
    }
    FREE(grid->points);
  }
  /* Clear grid size dependent data from stencils:                     */
  if (grid->directStencil!=NULL){
    freeStencilForGridResize(grid->directStencil,grid->nGridElements,
			     grid->sourceAndEvalDiffer);
  }
  if (grid->precondStencil!=NULL){
    freeStencilForGridResize(grid->precondStencil,grid->nGridElements,
			     grid->sourceAndEvalDiffer);
  }

  /* Zero out grid point references to grid elements (these references 
   * will change for a different grid size, but the size of this array 
   * itself won't change.)                                             */
  if (grid->sourceAndEvalDiffer) typemax=1;
  else typemax=0;
  for (itype=0; itype<=typemax; itype++){
    for (iele=0; iele<grid->nGridElements[itype]; iele++){
      grid->associatedElements[itype][iele]=0;
    }
  }
  /* Clear memory for grid data and grid kernel                        */
  if (grid->griddata != NULL){
    freeFFT(grid->griddata);
    assert(grid->griddata==NULL);
  }
  if (grid->kernel != NULL){
    freeFFT(grid->kernel);
    assert(grid->kernel==NULL);
  }
  return;
}
/* 
* ==================================================================== 
* Free all parts of a grid stencil that depends on the choise of grid 
* size. It is assumed that some parts of the stencil structure may not 
* be allocated. If so, then the respective pointers must point to NULL.
* ==================================================================== */
static
void freeStencilForGridResize(struct gridStencil *stencil, 
			      int *nGridElements, int sourceAndEvalDiffer)
{
  char fctName[] = "freeStencilForGridResize";
  int iele, ilink, itype,typemax;
  struct neighbourLinkedList *nextLink, *thisLink;
  double freedMem,totlinksFreed;
  int verbose = 0;

  if (sourceAndEvalDiffer) typemax=1;
  else typemax=0;
  /* Free linked lists and clear gridElementNeighbours array           */
  for (itype=0,totlinksFreed=0; itype<=typemax; itype++){
    for (iele=0; iele<nGridElements[itype]; iele++){
      nextLink = stencil->gridElementNeighbours[itype][iele]->neighbourLinkedList;
      totlinksFreed += stencil->gridElementNeighbours[itype][iele]->nNeighbours;
      for (ilink=0; 
	   ilink<stencil->gridElementNeighbours[itype][iele]->nNeighbours;
	   ilink++){
	
	thisLink = nextLink;
	nextLink = thisLink->nextNeighbour;
	FREE(thisLink);
      }
      /* Set the first pointer to point to NULL (it doesn't point to 
       * anywhere logical now anyway, and it is used as flag to tell 
       * other parts of the program that the pointer is not associated 
       * with memory). Also set the number of neighbouring elements to zero. */
      stencil->gridElementNeighbours[itype][iele]->neighbourLinkedList=NULL;
      stencil->gridElementNeighbours[itype][iele]->nNeighbours = 0;
    }
  }
  /* Neighbours are no longer found, nor counted or numbers estimated. 
   * Toggle flags to reflect this:                                     */
  stencil->gridElementNeighboursAreFound  = 0;
  stencil->nGridElementNeighboursAreFound = 0;
  /* For good measure, set the link count to zero                      */
  stencil->nGridElementNeighbours = 0;

  freedMem = ((double) totlinksFreed)
	* ((double) sizeof(struct neighbourLinkedList));
  if (verbose) 
    printf("%s: %12.3eMb memory freed from linked lists\n",
	   fctName,freedMem/pow(2.0,20.0));
  return;
}

/* 
* ==================================================================== 
* Count the number links needed for the direct / precorrect 
* interactions without removing doubly counted entries. 
* The count makes it possible to allocate all the needed links in 
* one contiguous part of memory, rather than allocating on the fly 
* as they are needed. The reason for abandoning this latter procedure,
* which was used earlier, was the possible fragmentation of memory 
* if the memory for the links are later deallocated and needs to be 
* reclaimed.
* Further, the present routine will make it possible to estimate the 
* memory use without actually allocating the linked lists.
*
* Input parameters: 
*   gridIn: The grid cast to a generic pointer.
*   returnZeroIfAssociated: Flag. If non-zero causes the 
* returnZeroIfAssociated is an input parameter 
* ==================================================================== */
int countPossibleNeighbourLinks(void *gridIn,
				int returnZeroIfAssociated, 
				char *Type)
{
  char fctName[] = "countPossibleNeighbourLinks";
  int ip,jp,kp;   /* Grid point counters  */
  int nLinks;     /* Number of links needed to store neighbours */
  int createLinks,elementType;
  struct gridStencil *stencil;
  struct grid *grid;
  int verbose = 0;

  /* Check that we have a valid call                                   */
  assert(strcmp(Type,"direct")==0 || strcmp(Type,"precond")==0);

  /* Cast grid to correct type:                                        */
  grid = (struct grid*) gridIn;

  /* Choose correct stencil to examine                                 */
  if(strcmp(Type,"direct")==0)  stencil = grid->directStencil;
  if(strcmp(Type,"precond")==0) stencil = grid->precondStencil;

  /* Return zero without counting if that is called for                */
  if (returnZeroIfAssociated!=0 &&
      stencil->gridElementNeighboursAreFound!=0) {
    nLinks = 0;
  }
  /* If the links have already been counted, then just return the 
   * previous result (estimated or (perhaps) actually counted)         */
  else if (stencil->nGridElementNeighboursAreFound){
    nLinks = stencil->nGridElementNeighbours;
  }
  else { /* Actually count */
    /* For every grid point find the number of needed links in the 
     * neighbour lists using "findGridElementNeighbours" */
    nLinks = 0;
    for (ip=0; ip<grid->nx; ip++){
      for (jp=0; jp<grid->ny; jp++){
	for (kp=0; kp<grid->nz; kp++){
	  if (grid->points[ip][jp][kp] != NULL){
	    nLinks += findGridElementNeighbours(ip,jp,kp,
						elementType=0,
						grid,
						stencil,
						createLinks=0);
	  }
	}
      }
    }
    /* Store result in the stencil and set flag to show that the 
     * the number of links has been estimated.                         */
    stencil->nGridElementNeighbours = nLinks;
    stencil->nGridElementNeighboursAreFound = 1;
  }
  if(verbose) printf("%s: %d links needed\n",fctName,nLinks);
  /* Return the total number of links */
  return nLinks;
}
/* 
* ==================================================================== 
* Find the maximum number of elements associated with a grid point 
* and the number of grid points which have elements associated with 
* them (which may yield the average number of elements per non-NULL 
* grid point).
* This will make it possible to estimate the number of counting 
* in the direct part. If, for instance all elements are associated 
* with the same grid node, then an O(N^2) effort will be spent just 
* in counting how many direct-part entries there will be. Admittedly, 
* the factor in front of the N^2 is small, but it is still O(N^2). 
* The present routine can be used to discard a grid, before an estimate 
* of the direct part is made, solely on the basis of, say, the maximum
* number of elements associated with a particular grid point.
*
* The maximum numbers of asociated elements (sources and evals) are 
* returned thorugh the argument list, while the number of grid points 
* visited are returned as the function value.
* ==================================================================== */
static
int countGridElementAssociations(struct grid *grid,
				 int *maxAssocElements)
{
  char fctName[] = "countGridElementAssociations";
  int verbose = 0;
  int ip,jp,kp;     /* Grid point counters  */
  int nVisited;     /* Number of grid points with associated elements  */
  int maxSrcs, maxEvals; /* Counters of maximum values                 */
  gridPoint ****points;  /* Shorthand                                  */

  points = grid->points;

  /* For every grid point find the number of associated elements, 
   * keep only the maximum.                                            */
  maxSrcs  = 0;
  maxEvals = 0;
  nVisited = 0;
  for (ip=0; ip<grid->nx; ip++){
    for (jp=0; jp<grid->ny; jp++){
      for (kp=0; kp<grid->nz; kp++){
	if (grid->points[ip][jp][kp] != NULL){
	  nVisited++;
	  maxSrcs  = MAX(points[ip][jp][kp]->nAssociatedElements[0],maxSrcs);
	  maxEvals = MAX(points[ip][jp][kp]->nAssociatedElements[0],maxEvals);
	}
      }
    }
  }

  if (verbose)
    fprintf(stdout,"%s: MAX %d sources and %d evals associated for %d points\n",
	    fctName,maxSrcs,maxEvals,nVisited);

  maxAssocElements[0] = maxSrcs;
  maxAssocElements[1] = maxEvals;

  return nVisited;
}
/* 
* ==================================================================== 
* Return the maximum extend of a stencil in any of the three directions 
* of the grid.
* If the stencil is not allocated (it points to NULL), then return zero.
* ==================================================================== */
static
int stencilRange(struct gridStencil *stencil)
{
  int range, ip;

  range = 0;
  if (stencil != NULL){
    for (ip=0; ip<stencil->nStencilPoints; ip++){
      range = MAX(range, ABS(stencil->ipNeighbours[ip]));
      range = MAX(range, ABS(stencil->jpNeighbours[ip]));
      range = MAX(range, ABS(stencil->kpNeighbours[ip]));
    }
  }

  return range;
} /* End of routine stencilRange */

/* 
* ==================================================================== 
* Given lists of projection and interpolation points 
* return kernel values for all combinations of point pairs.
* This assumes that the kernel has not been convolved.
* ==================================================================== */
void supplyKernelValues(void *gridInput, 
			int nProjectPoints, int *projectPoints,
			int nInterpPoints, int *interpPoints,
			double **kernelValuesArray){
  struct grid *grid;
  /* Cast input data to correct type                                 */
  grid  = (struct grid*) gridInput;

  /* Extract data and call grid-data routine */
  supplyKernelDataValues(grid->kernel,grid->griddata, 
			 nProjectPoints, projectPoints,
			 nInterpPoints, interpPoints,
			 kernelValuesArray);
  return;
} /* End of routine supplyKernelValues */

/* 
* ==================================================================== 
* Grid-based precorrection scheme. Precorrect together all columns 
* which are associated with the same grid point. 
* [Really: precorrect together all columns for which the source element 
* have identical points in the projection stencil.]
* ==================================================================== */
void precorrectThroughGridOld(void *directMatIn, 
			      void *projectMatIn, 
			      void *interpMatIn,
			      void *gridIn){
  char fctName[] = "precorrectThroughGrid";
  struct grid *grid;
  gridPoint ****gridPoints;
  struct gridElement **gridElements[2];
  int ip,jp,kp, ic, nSources, iSource,iele, nCols, localColNumber, 
    *colIndices;
  int maxSources=0, *nColumns=NULL, maxCols=0, *columns=NULL;
  int verbose = 0;

  /* Cast input data to correct type                                 */
  grid  = (struct grid*) gridIn;

  /* Grid "shorthands" */
  gridPoints = grid->points;
  gridElements[0] = grid->gridElements[0];

  /* For each grid point */
  for (ip=0; ip<grid->nx; ip++){
    for (jp=0; jp<grid->ny; jp++){
      for (kp=0; kp<grid->nz; kp++){
	/* Skip point if there are no associated elements 
	 * or no associated sources (necessary to check for 
	 * "elements" first) */
	if (gridPoints[ip][jp][kp]==NULL ||
	    (nSources = gridPoints[ip][jp][kp]->nAssociatedElements[0]) == 0)
	  continue; /* Skip to next grid point */

	/* Allocate memory if necessary                                */
	if (nSources>maxSources){
	  if (verbose) printf("%s: Allocating memory, part I\n",fctName);
	  /* Free old memory if allocated */
	  IFFREE(nColumns);
	  /* Set new size */
	  maxSources = MAX(2*maxSources,nSources);
	  /* Allocate new memory */
	  nColumns = (int *) calloc(maxSources, sizeof(int));
	}

	/* CAVEAT: If a "multi-level" projection scheme is employed, 
	 * then only sources that are projectet at the same level 
	 * should be precorrected together. 
	 * It might be a good strategy to precorret all columns 
	 * corresponding to "weird" levels (>base) here, and then 
	 * precorrect the remaining columns below.                     */
	for (iele=0; iele<nSources; iele++){
	  iSource = gridPoints[ip][jp][kp]->associatedElements[0][iele];
	  assert( (int)gridElements[0][iele]->projectLevel == 1 );
	}

	/* Find out how many distinct columns in the direct matrix 
	 * that are correspond to these sources                        */
	for (iele=0, nCols=0; iele<nSources; iele++){
	  /* Number of this source (as grid-element)                   */
	  iSource = gridPoints[ip][jp][kp]->associatedElements[0][iele];
	  /* Number of distinct columns from this source:              */
	  nColumns[iele] = 
	    numSourceUnknowns(gridElements[0][iSource]->element);
	  /* Add to total                                              */
	  nCols   += nColumns[iele];
	}

	/* Allocate memory if necessary                                */
	if (nCols>maxCols){
	  if (verbose) printf("%s: Allocating memory, part II\n",fctName);
	  /* Free old memory if allocated */
	  IFFREE(columns);
	  /* Set new size */
	  maxCols = MAX(2*maxCols,nCols);
	  /* Allocate new memory */
	  columns = (int *) calloc(maxCols, sizeof(int));
	}
	/* Store column entries in in array */
	for (iele=0, localColNumber=0; iele<nSources; iele++){
	  iSource = gridPoints[ip][jp][kp]->associatedElements[0][iele];
	  /* Get column entries for this source */
	  colIndices = sourceIndices(gridElements[0][iSource]->element);
	  for (ic=0; ic<nColumns[iele]; ic++,localColNumber++){
	    columns[localColNumber]=colIndices[ic];
	  }
	}
	/* Test that we have counted correctly */
	assert(localColNumber==nCols);
	
	/* Precorrect these columns together */
	precorrectColumnsOld(nCols, columns,
			     directMatIn, projectMatIn, interpMatIn, 
			     gridIn);

      }
    }
  }
  /* Free memory used here: */
  IFFREE(nColumns);
  IFFREE(columns);

  return;
} /* End of routine precorrectThroughGridOld */

/* 
* ==================================================================== 
* Grid-based precorrection scheme. Precorrect together all columns 
* which are associated with closely spaced grid points, i.e. grid points 
* lying in a nxnxn cube, n being a natural number.
* This approach allows compression of projection points, BUT it will 
* make some unneccessary computations. 
* This approach does at the present not seem to do improve the speed, 
* but it does not decrease efficiancy either (as long as n is kept low, 
* e.g. 1 or 2). Maybe "n" should depend on the average panel size to 
* grid size ratio? At present it is made an input parameter.
* The use of n=1 is recommended, but feel free to play around with 
* other values if you think it may benefit your particular problem.
* If the projection stencil is "sparse", as the eight corners of a 
* 3x3x3 cube, then the approach taken in this routine (with n>1) may 
* actually increase the CPU time needed to the precorrection.
* ==================================================================== */
void precorrectThroughGrid(void *directMatIn, 
			   void *projectMatIn, 
			   void *interpMatIn,
			   void *gridIn,
			   int n){
  char fctName[] = "precorrectThroughGrid";
  struct grid *grid;
  gridPoint ****gridPoints;
  struct gridElement **gridElements[2];
  int ip,jp,kp, ipg,jpg,kpg, ipl,jpl,kpl, 
    ic, nTotSources, iSource,iele,ieleTot, ipoint, 
    nCols, localColNumber, *colIndices,
    nTotPoints, iProgress;
  char mark = '#';
  int maxSources, *nColumns, *iSources, maxCols, *columns;
  int nPoints, *nSources;
  int nCalls;
  int verbose = 1;

  /* Cast input data to correct type                                   */
  grid  = (struct grid*) gridIn;

  /* Grid "shorthands" */
  gridPoints = grid->points;
  gridElements[0] = grid->gridElements[0];

  /* Initialize memory and counters */
  maxSources=0;
  nColumns=NULL;
  iSources=NULL;
  maxCols=0;
  columns=NULL;

  /* Set size of cube to precorrect together: */
  nPoints = n*n*n;
  /* Number of sources for each point considered */
  nSources = (int *) calloc(nPoints,sizeof(int));

  nCalls = 0;

  /* Set up stuff to mark progress */
  nTotPoints = grid->nx*grid->ny*grid->nz;
  iProgress  = 0;
  if (verbose) plotProgressMark(0, 1, 0, nTotPoints,1,mark);
  /* For each grid point */
  for (ipg=0; ipg<grid->nx; ipg+=n){
    for (jpg=0; jpg<grid->ny; jpg+=n){
      for (kpg=0; kpg<grid->nz; kpg+=n){
	/* Plot progress */
	iProgress += nPoints;
	if (verbose) plotProgressMark(iProgress, 0, 0, nTotPoints,0,mark);
	
	/* Loop over points in this (little) cube */
	for (ipoint = 0, nTotSources=0, ipl=0; ipl<n; ipl++){
	  for (jpl=0; jpl<n; jpl++){
	    for (kpl=0; kpl<n; kpl++, ipoint++){
	      ip = ipg + ipl; /* "Global" plus "local" */
	      jp = jpg + jpl;
	      kp = kpg + kpl;
	      if ( !gridPointExists(ip,jp,kp,grid) || 
		   gridPoints[ip][jp][kp]==NULL){
		/* If the point does not exists or if it has no 
		 * elements associated with it, then there is 
		 * certainly no associated sources.                    */
		nSources[ipoint] = 0;
	      }
	      else {
		/* The point exists and has associated elements. Still 
		 * there may be no sources, but at least we can now 
		 * count them.                                         */
		nSources[ipoint] = 
		  gridPoints[ip][jp][kp]->nAssociatedElements[0];
		nTotSources += nSources[ipoint];
	      }
	    } /* Next local k (kpl) */
	  }   /* Next local j (jpl) */
	}     /* Next local i (ipl) */

	/* Only continue if there are sources on some of these points  */
	if (nTotSources==0) continue;

	/* Allocate memory if necessary                                */
	if (nTotSources>maxSources){
	  if (verbose>1) printf("%s: Allocating memory, part I\n",fctName);
	  /* Free old memory if allocated */
	  IFFREE(nColumns);
	  IFFREE(iSources);
	  /* Set new size */
	  maxSources = MAX(2*maxSources,nTotSources);
	  /* Allocate new memory */
	  nColumns = (int *) calloc(maxSources, sizeof(int));
	  iSources = (int *) calloc(maxSources, sizeof(int));
	}

	/* Find out how many (distinct) columns in the direct matrix 
	 * that correspond to these sources                            */
	for (ipoint = 0, ieleTot=0, ipl=0, nCols=0; ipl<n; ipl++){
	  for (jpl=0; jpl<n; jpl++){
	    for (kpl=0; kpl<n; kpl++, ipoint++){
	      ip = ipg + ipl; /* "Global" plus "local" */
	      jp = jpg + jpl;
	      kp = kpg + kpl;
	      if (nSources[ipoint]>0){
		for (iele=0; iele<nSources[ipoint]; iele++){
		  /* Number of this source (as grid-element)           */
		  iSource = 
		    gridPoints[ip][jp][kp]->associatedElements[0][iele];
		  /* Store number of this source (to avoid triple loop 
		   * later)                                            */
		  iSources[iele+ieleTot] = iSource;
		  /* Number of distinct columns from this source:      */
		  nColumns[iele+ieleTot] = 
		    numSourceUnknowns(gridElements[0][iSource]->element);
		  /* Add to total                                      */
		  nCols   += nColumns[iele];
		}
		/* Update how many sources have been examined:         */
		ieleTot += nSources[ipoint];
	      } /* End if (nSources[]>0) */
	    } /* Next local k (kpl) */
	  }   /* Next local j (jpl) */
	}     /* Next local i (ipl) */
	/* Test that this count was made correctly: */
	assert(ieleTot==nTotSources);

	/* Allocate memory if necessary                                */
	if (nCols>maxCols){
	  if (verbose>1) printf("%s: Allocating memory, part II\n",fctName);
	  /* Free old memory if allocated */
	  IFFREE(columns);
	  /* Set new size */
	  maxCols = MAX(2*maxCols,nCols);
	  /* Allocate new memory */
	  columns = (int *) calloc(maxCols, sizeof(int));
	}

	/* Store column entries in in array */
	for (iele=0, localColNumber=0; iele<nTotSources; iele++){
	  /* Number of this source stored from earlier                 */
	  iSource = iSources[iele];
	  /* Get column entries for this source */
	  colIndices = sourceIndices(gridElements[0][iSource]->element);
	  for (ic=0; ic<nColumns[iele]; ic++,localColNumber++){
	    columns[localColNumber]=colIndices[ic];
	  }
	}
	/* Test that we have counted correctly */
	assert(localColNumber==nCols);
	
	/* Precorrect these columns together                           */
	precorrectColumns(nCols, columns,
			   directMatIn, projectMatIn, interpMatIn, 
			   gridIn);
	/* Update the number of calls to make, 
	 * for statistical purposes only.                              */
	nCalls++;
      }
    }
  }
  /* End progress marking. */
  if (verbose) plotProgressMark(iProgress, -1, 0, nTotPoints,0,mark);

  if (verbose) printf("%s: Precorrection through %d calls\n",
		      fctName,nCalls);
 

  /* Free memory used here: */
  IFFREE(nSources);
  IFFREE(nColumns);
  IFFREE(columns);

  return;
} /* End of routine precorrectThroughGrid */

/* 
* ==================================================================== 
* This routine is NOT finished. Calling it will graciously 
* halt the program! The routine is intended for deciding the grid 
* level based on the size of the elements. 
* The following comments is intended for the implementation 
* of the routine:
*
* This routine counts the number (of a specified type) of elements, 
* which are "large enough" to require special interpolation/projection 
* scheme (i.e. "go up grid levels". If many (by some measure) elements 
* require special treatment, then the grid is probably to dense and a 
* coarser discretization should be taken.
*
* Input parameters:
*
*   grid:     The grid structure
*
*   elementType: Type of grid-elements to operate on
*              0: Sources; 
*              1: Evals;
*
*   nearDistanceSphere: If the maximum distance fromt the stencil 
*              origin to a point on a gridElement exceeds this value, 
*              then a go up a grid level.
*
*   nearDistanceCube: Similar to "nearDistanceSphere", but used as 
*              maximum distance in one spatial dimension, i.e. giving 
*              (half the side length of) a cube, rather than (the radius 
*              of) a sphere, in which the grid-element must lie to be 
*              projected/interpolated to/from the finest grid level.
*
*              nearDistanceSphere and nearDistanceCube are both given 
*              normalized by the grid size (dx). A negative value may be 
*              given to ommit the specific test, but at least one of 
*              the two must be positive.
*              Example: To ensure that elements do not "hang out" of 
*              cubes (at all!) nearDistanceCube=1.0 may be used in 
*              conjunction with a negative value of nearDistanceSpere.
*              In general it is better to allow elements to "hang out" 
*              a little rather than to increase the level at which they 
*              are projected/interpolated. However, nearDistanceCube=1.0 
*              may still be used to get an idea of the general size of 
*              the elements compared to the grid...
*
*   countOnly: If countOnly=1, then the number of elements that need 
*              projection (interpolation) on a level different from the 
*              finest is found, but the level is not stored with each 
*              grid element.
*
* The return value (int) of this routine is the number of elements 
* that need projection (interpolation) on a level different from 
* the finest one.
* ==================================================================== */
static
int findElementStencilLevels(struct grid *grid, 
			     int elementType,
			     double nearDistanceCube,
			     double nearDistanceSphere,
			     int countOnly)
{
  char fctName[] = "findElementStencilLevels";
  struct gridElement **gridElements; 

  /* Test input data */
  if (elementType!=0 && elementType!=1){
    /* This is unknown = bad. Report error and exit */
    fprintf(stderr,
	    "%s: ERROR. Unknown element type: %d. Should be 0 or 1.\n",
	    fctName,elementType);
    _EXIT_;
  }

  /* Shorthand: */
  gridElements = grid->gridElements[elementType];

  /* CONTINUE FROM HERE!! */
  _EXIT_;
  return 0;
} /* End of routine findElementStencilLevels */


/* Test functions and stuff that are not really needed (but
 * which need access to the variables in this file) can go
 * in the following file: */
#include "grid.testfunctions.c"
