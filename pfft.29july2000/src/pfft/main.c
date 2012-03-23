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
/* Example on use of the pfft routines. */

#include "pfft.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <time.h>

/* 
* ==================================================================== 
* MACROS
* ==================================================================== */ 
/* The _EXIT_ macro is convenient when execution is stopped due to error
 * encounters and so forth. Also nice for debugging / code developing: */
#define _EXIT_ {\
 printf("Exited from line # %d in '%s'\n",__LINE__,__FILE__); \
 exit(1);}

#define TRUE 1
#define FALSE 0

#define TYPESTR_MAX 100 

/* Macros for easy timing. For long calculations on some systems clock 
 * (which counts the CPU time) may overflow. In particular, on linux 
 * boxes, time overflows after 2.148secs due to a high value (10^6) of 
 * CLOCKS_PER_SEC. If you are compiling on a UNIX (or Linux) system, 
 * then you may want to use SETTIMER_UNIX, which use the NONSTANDARD 
 * times() function (in the <times.h> headers), see "man times". Comment
 * out the next line to used standard timings (leave in for UNIX ;).   */
#define USE_UNIX_TIMERS
/* ANSI-C implementation (doesn't work very well on Linux systems)     */
#define SETTIMER_ANSI(settimer,setWall){settimer = clock();setWall=time(NULL);}
#define SHOWTIME_ANSI(showtimer,showwall,str){\
        clock_t newTimer;       \
        time_t newWall;       \
                  SETTIMER(newTimer,newWall);\
                  printf("%s: Time for %s was %.2f secs \n"\
			 "          (%.0f secs wall clock time)\n",   \
			 fctName,str,\
                         ((double)(newTimer-showtimer)/(double)CLOCKS_PER_SEC),\
                         difftime(newWall,showwall));                \
                  SETTIMER(showtimer,showwall);\
                        }
/* UNIX / Linux implementation (not as portable as the above!) 
 * Comment out the inclusion of <sys/times.h> and the SETTIMER_UNIX 
 * if you get error messages. In that case you probably have to 
 * stick to the standard timing routines.                              */
#ifdef USE_UNIX_TIMERS
# include <sys/times.h>
/* The following few lines is a hack. (I got errors of 
 * missing CLK_TCK when compiling a debug version of the code.         */
#ifndef CLK_TCK
# define CLK_TCK 100
#endif
# define SETTIMER_UNIX(settimer,setWall){\
  struct tms temp_clock;\
  times(&temp_clock);\
  settimer = temp_clock.tms_utime;\
  setWall=time(NULL);\
  }
# define SHOWTIME_UNIX(showtimer,showwall,str){\
        clock_t newTimer;       \
        time_t newWall;       \
                  SETTIMER(newTimer,newWall);\
                  printf("%s: Time for %s was %.2f secs \n"\
			 "          (%.0f secs wall clock time)\n",   \
			 fctName,str,\
                         ((double)(newTimer-showtimer)/(double)CLK_TCK),\
                         difftime(newWall,showwall));                \
                  SETTIMER(showtimer,showwall);\
                        }
# define SETTIMER SETTIMER_UNIX
# define SHOWTIME SHOWTIME_UNIX
#else
# define SETTIMER SETTIMER_ANSI
# define SHOWTIME SHOWTIME_ANSI
#endif
static
void parseCommandLine(int argc, char *argv[],
		      char *sourceInputFile,
		      char *evalInputFile,
		      int *inputFileArguments);

int main(int argc, char *argv[])
{
  char fctName[]="main";  /* Name of this routine                      */
  int  nSources,          /* Total number of source elements           */
    nEvals,               /* Total number of evaluation elements       */
                          /* Note that nSources and nEvals are the total 
			   * number of ELEMENTs, not the total number of 
			   * unknowns. The total number of unknowns are 
			   * stored in numColsDirect and numRowsDirect, 
			   * respectively, for the sources and the 
			   * evals.                                    */
    numDirectMats;        /* Number of (sparse) direct matrices        */
  void *sources,          /* Pointer to source elements                */
    *evals,               /* Pointer to evaluation elements            */
    *grid,                /* Pointer to the grid                       */
    *directMat,           /* Pointer to the (sparse) direct matrices   */
    *precondMat,          /* Pointer to the (sparse) preconditioning 
			   * matrices                                  */
    *projectMat,	  /* Pointer to sparse grid projection matrix. */
    *interpMat, 	  /* Pointer to sparse grid interpolation matrix.  */
    *bconds,              /* Pointer to boundary condition values      */
    *rightHandSide,       /* Pointer to the right-hand-side vector     */
    *solutionVector;      /* Pointer to vector for storing solution    */
  char sourceInputFile[FILENAME_MAX],/* File names of input files for  */
    evalInputFile[FILENAME_MAX];     /* elements (sources and evals)   */
  int numColsDirect,                 /* Number of columns in the 
			              * (sparse) direct matrix. This 
			              * corresponds to the degrees of 
				      * freedom for the sources        */
    numRowsDirect,                   /* Number of "rows" in the (sparse) 
			              * direct matrix. Note, that the 
			              * number of rows is not needed 
			              * for allocating the sparse matrix. 
			              * This corresponds to the degrees 
				      * of freedom for the evals       */
    maxPrecNeighbours;               /* Maximum numbers of neighbouring 
			              * elements to use in the 
			              * preconditioning step           */
  int isConductors;                  /* !=0 if we solve Ax=b 
				      * (as opposed to Ax=By)          */
  int maxBoundaryPart;               /* Largest value for boundary parts 
				      * (e.g. conductor numbers) in the 
				      * element discretization. This value 
				      * is used for setting up boundary 
				      * conditions and integrating the 
				      * total solution.                */
  /* Arguments to be passed to other parts of the program              */
  int inputFileTypes[2];
  char bcType[TYPESTR_MAX];
  int  bcIdata;

  int swapNeumann=0; /* Set to 1 to swap sign on Neuman boundary 
		      * conditions (as to negate the normal vector)    */
  /* Variables for timing :*/
  clock_t timer, totTimer;
  time_t  wallClock,totWallClock;
  
  /* Define start timer */
  SETTIMER(totTimer,totWallClock);

  /* This sets the boundary condition type and values that are in use: 
   * Choose from (see bconds.c)
   *  bcType :                 bcIdata
   *    "conductor i"           *ignored*
   *    "unity"                 0,1
   *    "uniform"               0,1,*
   *    "constant"              0,1,*
   *    "translating sphere"    *ignored*
   *    "sphere charge"         *ignored*                              */
  isConductors = 1;
  strcpy(bcType,"unity"); /* Start finding TOTAL SELF CAPACITANCE      */
  /*printf(" >> SOLVING sphere charge NONUNIFORM PROBLEM <<\n");
    strcpy(bcType,"sphere charge");*/
  bcIdata = 1;            /* Dirichlet conditions everywhere           */

  SETTIMER(timer,wallClock);
  /* Parse command line to obtain input files                          */
  parseCommandLine(argc,argv,
		   sourceInputFile,
		   evalInputFile,
		   inputFileTypes);
  SHOWTIME(timer,wallClock,"command line parsing");

  /* Setup of source elements (reading file):                          */
  sources = setupSources(&nSources,
			 &maxBoundaryPart,
			 sourceInputFile,
			 inputFileTypes,
			 bcType,
			 bcIdata);
  SHOWTIME(timer,wallClock,"source read and setup");

  /* Possibly write elements to a tecplot file 
  dumpElementsToTecplot(sources,nSources, 0,NULL);
  SHOWTIME(timer,wallClock,"writing elements to file in TecPlot format");
  _EXIT_;*/

  /* Check source elements: 
  tableElements(sources,nSources); 
  */

  /* Check if more than one input file is specified                    */
  if (!strcmp(sourceInputFile,evalInputFile)){
    /* Sources and evals are the same. */
    evals = sources;
    nEvals = nSources;
  }
  else{
    /* Setup of evaluation elements (reading file).*/
    evals = setupEvals(&nEvals,evalInputFile,inputFileTypes);
    SHOWTIME(timer,wallClock,"eval read and setup");
  }

  /* Assign elements with indices in the sparse direct matrix. 
   * Sources are associated with columns and evals with rows.          
   * Return the number of columns needed for the direct part, and a 
   * vector for storing the number of nonzero rows in each column      */
  printf("%s: Assigning indices for direct part to elements\n",fctName);
  assignElementDirectIndices(sources, nSources,
			     evals, nEvals,
			     &numColsDirect,&numRowsDirect);
  SHOWTIME(timer,wallClock,"element indices assignments");

  /* Setup the grid elements and other grid values. 
   * Find the neighbour lists...                                       */
  printf("%s: Setting up the grid\n",fctName);
  grid = setupGrid(sources,nSources,evals,nEvals,
		   numColsDirect,numRowsDirect,numDirectMats=1); 
  SHOWTIME(timer,wallClock,"grid setup");

  if (!isConductors){ /* "Ax = By" form: "Direct formulation". */
    /* === RIGHT-HAND-SIDE SETUP === */
    /* Set up the grid2grid data. 
     * Note: Some of this is needed when building the projection 
     * and interpolation matrices                                        */
    printf("%s: Setting up grid-to-grid part\n"
	   "  --==> Some work is needed for the right-hand-side? <==-- \n",
	   fctName);
    grid2gridSetup(grid);
    SHOWTIME(timer,wallClock,"grid-to-grid setup");

    /* Build the structure needed for the sparse projection matrix.      */
    printf("%s: Setting up the projection matrix\n",fctName);
    projectMat = callocElementGridSparse(grid,"source", nSources, (int) 1);
    /* calculate the projection (sparse) matrix.                         */
    printf("%s: Computing values of the projection matrix\n",fctName);
    calcElementGridSparse(grid, "source", "rhs", projectMat);
    SHOWTIME(timer,wallClock,"projection setup");

    /* Build the structure needed for the sparse interpolation matrix.   */
    printf("%s: Setting up the interpolation matrix\n",fctName);
    interpMat = callocElementGridSparse(grid,"eval", nEvals, (int) 1);
    /* calculate the interpolation (sparse) matrix.                      */
    printf("%s: Computing values of the interpolation matrix\n",fctName);
    calcElementGridSparse(grid, "eval", "rhs", interpMat);
    SHOWTIME(timer,wallClock,"interpolation setup");

    /* Build the structure needed for the sparse direct matrix.          */
    printf("%s: Setting up the direct part structure\n",fctName);
    directMat =callocDirectSparse(grid,numColsDirect,numRowsDirect,
				numDirectMats=1);
    /* Build the right-hand-side matrix */
    printf("%s: Setting up the direct part (for right-hand-side)\n",fctName);
    calcDirectSparse(grid,directMat,nSources,"rhs");
    SHOWTIME(timer,wallClock,"direct part for rhs");
    /* Precorrect the right-hand-side matrix 
     * Precorrect all sources at a cluster of grid points (nxnxn) 
     * together (choose n=1)                                             
    printf("%s: Precorrecting direct part\n",fctName);
    precorrectThroughGrid(directMat, projectMat, interpMat, grid,  1);
    SHOWTIME(timer,wallClock,"precorrection calculations");*/
    /* Build one large (skinny) near-field grid-to-grid matrix */
    printf("Precorrecting direct part (one large)\n");
    calcPrecorrectSparse(directMat,projectMat,interpMat,grid,nSources,"one large");
    SHOWTIME(timer,wallClock,"precorrection calculations");

    /* Complete the grid setup by transforming the kernel.               */
    printf("%s: Transforming the kernel\n",fctName);
    grid2gridFinish(grid);
    SHOWTIME(timer,wallClock,"kernel transformation");

    /* Allocate memory for right-hand-side vector */
    bconds = allocateSystemVector(numColsDirect);
    /* Fill right-hand-side vector with boundary condition values        */
    setupElementBoundaryConditionValues(sources, nSources,bconds,
					bcType,bcIdata);
    /* Possibly swap sign on some entries in the boundary condition vector */
    if(swapNeumann) swapSignVectorElement(sources, nSources,
					  (double*) bconds, "neumann");
    
    /* Compute right-hand-side, i.e. multiply matrix by vector of 
     * boundary conditions                                               */
    rightHandSide = allocateSystemVector(numColsDirect);
    /* [vector should be zeroed out in case it has been used before] */

    /* TEST SYSTEM. Let the first entry be unity */
    if (0) { 
      double *vect;
      vect = (double*) bconds;
      vect[0] = 1.0;
      fprintf(stderr,
	      ">>>>WARNING: MAIN IS ALTERING THE BOUNDARY CONDITIONS!!  <<<<");
    }
    /* The following is "bad form" - change it later (main knows too much) */
    systemMatrixMultiply(directMat,interpMat, projectMat, 
			 NULL,   /* precondMat  */ 
			 "none", /* precondType */
			 grid,
			 (double *) bconds, (double  *)rightHandSide,
			 0, /* matNum */
			 numColsDirect, numRowsDirect);

    /* Clear the direct projection and interpolation matrices            */
    spMatReinit(directMat, 0);
    spMatReinit(projectMat, 0);/* Mat should possibly be totally deleted?*/
    spMatReinit(interpMat, 0); /* Mat should possibly be totally deleted?*/
    /* Clear the grid data and kernel setup.                             */
    grid2gridFree(grid);
  } 
  /* End "B*y"-type right hand side */ 
  else {
    /* This is for e.g. capacitance calculations 
     *          "Ax = y" form  - easy right hand side! */
    /* Allocate system matrices */
    projectMat = callocElementGridSparse(grid,"source", nSources, (int) 1);
    interpMat = callocElementGridSparse(grid,"eval", nEvals, (int) 1);
    directMat =callocDirectSparse(grid,numColsDirect,numRowsDirect,
				  numDirectMats=1);
    /* Allocate memory for right-hand-side vector */
    bconds = allocateSystemVector(numColsDirect);
    /* Fill right-hand-side vector with boundary condition values        */
    setupElementBoundaryConditionValues(sources, nSources,bconds,
					bcType,bcIdata);
    /* Easy right hand side - equal to the boundary conditions: */
    rightHandSide = bconds;
  }
  /* === LEFT-HAND-SIDE SETUP === */
  /* Setup grid-to-grid interactions for the left-hand-side            */
  /* Set up the grid2grid data. 
   * Note: Some of this is needed when building the projection 
   * and interpolation matrices                                        */
  printf("%s: Setting up grid-to-grid part\n",fctName);
  grid2gridSetup(grid);
  SHOWTIME(timer,wallClock,"grid-to-grid setup");

  /* calculate the projection (sparse) matrix.                         */
  printf("%s: Computing values of the projection matrix\n",fctName);
  calcElementGridSparse(grid, "source", "lhs", projectMat);
  SHOWTIME(timer,wallClock,"projection setup");
  /* calculate the interpolation (sparse) matrix.                      */
  printf("%s: Computing values of the interpolation matrix\n",fctName);
  calcElementGridSparse(grid, "eval", "lhs", interpMat);
  SHOWTIME(timer,wallClock,"interpolation setup");

  /* Setup the direct (sparse) matrix for solving the system 
   * (left-hand-side matrix).                                          */
  printf("%s: Setting up the direct part\n",fctName);
  calcDirectSparse(grid,directMat,nSources,"lhs");  /* MUST be "lhs" */
  SHOWTIME(timer,wallClock,"direct part setup");

  /* Build preconditioner if it is needed. 
   * It is strongly recommended to use the preconditioner.             */
  if (1) {
    /* If direct interactions are to be recycled for the preconditioner, 
     * then the preconditioner must be build before the direct matrix is 
     * precorrected.                                                   */
    /* Build the structure needed for the preconditioner.              */
    printf("%s: Setting up the preconditioner structure\n",fctName);
    precondMat = callocPrecondSparse(grid,directMat,
				     maxPrecNeighbours=100,
				     numColsDirect,numRowsDirect,
				     numDirectMats);
    /* Setup the preconditioning (sparse) matrix.                      */
    printf("%s: Setting up the values of the preconditioner\n",fctName);
    calcPrecondSparse(grid,directMat,precondMat,maxPrecNeighbours);
    SHOWTIME(timer,wallClock,"precondition setup");
  }
  

  /* Precorrection must be made before the kernel is transformed       */
  /* Precorrect the direct interaction. Several different 
   * implementations of this presently exist. All should give 
   * the same result, but get there by slightly different means 
   * (the order of the computations may differ, etc.) At the time 
   * of writing "calcPrecorrectSparse" is fastest for the problems 
   * tested. The naive implementetion at the bottom is about five 
   * times slower (rates may vary).                                    */
  /* Build one large (skinny) near-field grid-to-grid matrix */
  printf("Precorrecting direct part (one large)\n");
  calcPrecorrectSparse(directMat,projectMat,interpMat,grid,nSources,"one large");
  SHOWTIME(timer,wallClock,"precorrection calculations");

  /* Build many small near-field grid-to-grid matrices 
  printf("Precorrecting direct part\n");
  calcPrecorrectSparse(directMat,projectMat,interpMat,grid,nSources,"many small");
  SHOWTIME(timer,wallClock,"precorrection calculations");*/

  /* Precorrect all sources at a cluster of grid points (nxnxn) together
   * n=1 recommended 
  {
    int nloc=2; char line[60];
    printf("%s: Precorrecting direct part (%dx%dx%d points at a time) \n",
	   fctName,nloc,nloc,nloc);
    precorrectThroughGrid(directMat, projectMat, interpMat, grid,  1);
    sprintf(line,"%dx%dx%d points precorrection",2,2,2);
    SHOWTIME(timer,wallClock,line);
  }*/

  /* Precorrect all sources at a grid point together 
  printf("%s: Precorrecting direct part (point by point)\n",fctName);
  precorrectThroughGridOld(directMat, projectMat, interpMat, grid);
  SHOWTIME(timer,wallClock,"point-by-point precorrection");*/

  /* Precorrect one column at a time 
  printf("%s: Precorrecting direct part (column by column)\n",fctName);
  colPrec(directMat, projectMat, interpMat, grid);
  SHOWTIME(timer,wallClock,"column-by-column precorrection");*/

  /* Precorret one element at a time 
  printf("%s: Precorrecting direct part (naive implementation)\n",fctName);
  precorrect(grid, directMat, projectMat, interpMat);
  SHOWTIME(timer,wallClock,"naive precorrection calculations");*/

  /*dumpSparseMatrixTranspose(directMat, 0);
    _EXIT_;*/

  /* Complete the setup by transforming the kernel. */
  printf("%s: Transforming the kernel\n",fctName);
  grid2gridFinish(grid);
  SHOWTIME(timer,wallClock,"kernel transformation");

  /* TEST GRID OPERATIONS: 
  testProject2grid(grid, projectMat, interpMat, numColsDirect);
  testGrid2grid(grid);
  exit(1);*/

  /* Allocate solution vector */
  solutionVector = allocateSystemVector(numColsDirect);

  /* Solve the system with a dummy right-hand-side  */
  /* Use "none", "right" or "left" preconditioner.   */
if(0){
  printf("%s: Solving the system.\n",fctName);
  solveFull(directMat, interpMat, projectMat,
	    precondMat, "right", grid,
	    (double *) rightHandSide, (double  *) solutionVector,
	    0, /* matNum */
 	    numColsDirect, numRowsDirect);
  SHOWTIME(timer,wallClock,"right-preconditioned solve");
  /* Possibly swap sign on some entries in the solution vector */
  if(swapNeumann) swapSignVectorElement(evals, nEvals,
					(double*) solutionVector, "neumann");

  if (isConductors) 
    fprintf(stdout,"%s: Integrated solution: %16.6e\n",
	    fctName,
	    integrateOverElements(sources, nSources,evals, nEvals, 
				  (double*) solutionVector,
				  "whole",
				  NULL,-1) /* Last two arguments not used */
	    );
  
  /* Possibly test the solution if both Diriclet and Neumann 
   * conditions are known: */
  if (strcmp(bcType,"uniform")==0            || 
      strcmp(bcType,"constant")==0           || 
      strcmp(bcType,"translating sphere")==0 || 
      strcmp(bcType,"sphere charge")==0 ) {
    testSolution(evals, nEvals, (double*) solutionVector,
		 bcType,bcIdata);
  }
}
/*_EXIT_;*/
  /* Possibly find the capacitance matrix (conductor problem only)     */
  if (isConductors && maxBoundaryPart>0){
    double *capvector, **capmat, *solv, *rhs;
    int i,j;
    SETTIMER(timer,wallClock);
    solv = (double*) solutionVector;
    rhs  = (double*) rightHandSide;
    capvector = (double*) calloc(maxBoundaryPart+1,sizeof(double));
    capmat = (double**) calloc(maxBoundaryPart,sizeof(double*));
    for (j=0; j<maxBoundaryPart; j++) 
      capmat[j] = (double*) calloc(maxBoundaryPart,sizeof(double*));
    for (i=1; i<=maxBoundaryPart; i++){
      /* Zero out solution vector */
      for (j=0; j<numColsDirect; j++)    solv[j] = 0.0;
      for (j=0; j<numColsDirect; j++)    rhs[j]  = 0.0;
      for (j=0; j<=maxBoundaryPart; j++) capvector[j]  = 0.0;
      /* Set boundary conditions */
      setupElementBoundaryConditionValues(sources, nSources, rightHandSide,
					  "conductor i",i);
      /* Solve system (left-hand-side is unchanged) */
      solveFull(directMat, interpMat, projectMat,
		precondMat, "right", grid,
		(double *) rightHandSide, (double  *) solutionVector,
		0, /* matNum */
		numColsDirect, numRowsDirect);
      /* Integrate solution over each part of the boundary */ 
      integrateOverElements(sources, nSources,evals, nEvals, 
			    (double *)solutionVector,
			    "all parts",
			    capvector,
			    -1); /* Last argument not used */
      /* Display information */
      fprintf(stdout,"%s: CAPACITANCE #%d:",fctName,i);
      for (j=1; j<=maxBoundaryPart; j++) 
	fprintf(stdout," %14.6e",capvector[j]);
      fprintf(stdout,"\n");
    }
    /* Free memory */
    free(capvector);
    SHOWTIME(timer,wallClock,"capacitance matrix");
  }

  /* Write total run time to screen */
  SHOWTIME(totTimer,totWallClock,">>> TOTAL RUN <<<");
  _EXIT_;
  /* Temporarily write solution to stdout (swap back bcond values!)*/
  {
    int i, ip;
    double *sv, *bc;
    bc = (double *) bconds;
    sv = (double *) solutionVector;
    if(swapNeumann) swapSignVectorElement(sources, nSources,
					  bc, "neumann");
    ip = 4*nSources/6;
    /*printf(" val[%d] = %16.8e\n",ip,sv[ip]);
      _EXIT_;*/
    printf("BCs:\n");
    for (i=0;i<numColsDirect;i++) printf("%d: %16.8e %16.8e\n",
					 i,bc[i], sv[i]);
  }

_EXIT_;
  printf("%s: Solving the system. This needs some work\n",fctName);
  solveFull(directMat, interpMat, projectMat,
	    precondMat, "none", grid,
	    (double *) rightHandSide, (double  *) solutionVector,
	    0, /* matNum */
    numColsDirect, numRowsDirect);
  SHOWTIME(timer,wallClock,"no-preconditioned solve");
  solveFull(directMat, interpMat, projectMat,
	    precondMat, "left", grid,
	    (double *) rightHandSide, (double  *) solutionVector,
	    0, /* matNum */
	    numColsDirect, numRowsDirect);
  SHOWTIME(timer,wallClock,"left-preconditioned solve");

  /* Setup boundary conditions and right-hand-side.
   * Probably, this should be done before setting up the direct matrix.
   * I.e. set up a direct matrix for this part, then calculate the RHS, 
   * Evaluate a new direct matrix (with swapped Diriclet and Neumann 
   * boundary conditions), setup the preconditioner, and finally solve 
   * the direct part!                                                  */

  /*dumpSparseMatrixTranspose(directMat, 0);*/
  /* Test output */
  plotSizeOf();
  tableElements(sources,nSources);
  tableElements(evals,nEvals);
  exit(1);
  return 1;
}

/*
* ====================================================================
* Command line parsing routine.
*
* This routine was originially written by Tom Korsmeyer, 
* though now fairly heavily modified. 
* ==================================================================== */
static
void parseCommandLine(int argc, char *argv[],
		      char *sourceInputFile,
		      char *evalInputFile,
		      int *inputFileArguments)
{
  int cmderr, i, filesSoFar, inext, ielem, iformat;
  char **chkp, *chk, *elemType;
  char sourceType[]="source elements", 
    evalType[]="evaluation elements",
    allType[]="all elements";
  long strtol();
  FILE *testInputFile;

  cmderr = FALSE;
  chkp = &chk;			/* pointers for error checking */

  /* Set defaults (may be changed below by command line arguments      */
  inputFileArguments[0]=1;
  inputFileArguments[1]=1;

  filesSoFar=0; /* Number of input file names found */

  /* Look at input args and figure out what they are. */
  for(i = 1; i < argc && cmderr == FALSE; i++) {
    if(argv[i][0] == '-') {
      if (argv[i][1] == 'h') {
	/* Someone needs help. Put verbatim stuff here. */
	fprintf(stdout, "%s: The following help is presently available:\n",
		argv[0]);
	fprintf(stdout," Options:\n");
	/* Help options: */
	fprintf(stdout,
		"   -h :  Help (this message)\n\n");
	/* Input file format options: */
	fprintf(stdout,
		"   -i :  Input file format. \n"
                "         Presently, the only valid formats are \"1\" and \"2\",\n"
		"         but more formats may be added. An \"s\" (source elements)\n"
		"         or an \"e\" (evaluation elements) to denote that the \n"
		"         format is to be used only when reading a particular \n"
		"         element type. If only one input file is specified\n"
		"         (see below) then the format of the evaluation element\n"
		"         input file (-ie#) is not used.\n" 
		"         The default value is -is%d -ie%d\n"
		"         \n",
		inputFileArguments[0], /* For -i[s] */
		inputFileArguments[1]  /* For -i[e] */
		);
	fprintf(stdout,
		"         Formats:\n"
		"          1: For constant-type elements (panels).\n"
		"             Initial lines may contain an identifying header.\n"
		"             The first character of a header line must be '0'.\n"
		"             The first character of each subsequent line must be \n"
		"             in identifyer for the line. If information of an\n"
		"             element is given, then the identifyer must be '3'\n"
		"             or '4' denoting the number of corners of the element.\n"
		"             'T' or 't' (triangular) may be used instead of '3'\n"
		"             and 'Q' or 'q' (quadrilateral) may be substituded \n"
		"             for '4'.\n"
		"             After the element identifier an integer should indicate \n"
		"             what piece of the boundary the present element belongs to.\n"
		"             This identifyer must be a positive integer. Preferably\n"
		"             if n conductors are considered, then the identifiers should\n"
                "             be given numbers 1,...,n. The ordering does not matter.\n"
		"             Subsequently the corners of the elements are given as \n"
		"             x1 y1 z1  x2 y2 z2  x3 y3 z3  [x4 y4 z4]\n"
		"             Comment lines are accepted any time after the initial\n"
		"             header lines. The first character of a comment line\n"
		"             must be '#', '*' or '%%'.\n"
		"             For compability with fastcap '*' is preferred for comments.\n"
		"             This file type is typically named \"*.qui\".\n"
                "          2: As \"1\", except that the integer indicating boundary\n"
		"             part identification is omitted.\n"
		"             This file type is typically named \"*.dat\".\n"
		" Input files:\n"
		"     One or two file names must be given for reading element \n"
		"     data. The file containing the source elements (e.g. panels)\n"
		"     must be given before the file containing evaluation \n"
		"     elements (e.g. collocation points).\n"
		"     If only one file name is given, then it is assumed that\n"
		"     the evals are defined from the evals, e.g. as collocation\n"
		"     points at the panel centroids.\n\n"
		);
      }
      else if (argv[i][1] == 'i') {
	/* This is for setting the expected format of the input 
	 * files (for sources and evals                         */
	if (argv[i][2] == 's') {
	  /* Input format for the sources only */
	  ielem = 0; /* Sources */
	  inext = 3; /* Already used [i][2] */
	  elemType = sourceType;
	}
	else if (argv[i][2] == 'e') {
	  /* Input format for the evals only */
	  ielem = 1; /* Evals */
	  inext = 3; /* Already used [i][2] */
	  elemType = evalType;
	}
	else {
	  ielem = -1; /* All element types (sources and evals) */
	  inext =  2; /* Still need to use [i][2] */
	  elemType = allType;
	}
	/* Figure out what file format is used */
	if (argv[i][inext] == '1') iformat = 1;      /* File format #1 */
	else if (argv[i][inext] == '2') iformat = 2; /* File format #2 */
	/* Other formats? */
	else { /* This is unknown = bad. Report error and set flag */
	  fprintf(stderr, "%s: What file format is this \"%s\"\n"
		  " Acceptable values (with the -i argument) is: 1 and 2\n", 
		  argv[0], &(argv[i][inext]) );
	  cmderr = TRUE;
	}
	if (!cmderr) {
	  if (ielem == -1) { /* All element types (sources and evals) */
	    inputFileArguments[0]=iformat;
	    inputFileArguments[1]=iformat;
	  }
	  else if (ielem == 0 || ielem == 1) inputFileArguments[ielem]=iformat;
	  /* Write statement to show that the argument has been accepted */
	  fprintf(stdout,"%s: Using file format %d when reading %s\n",
		  argv[0],iformat,elemType);
	}
      }
      else if (FALSE) { /* (argv[i][1] == 'd') { */
	/* Might want to set the level
	*numLev = (int) strtol(&(argv[i][2]), chkp, 10);
	if(*chkp == &(argv[i][2]) || *numLev < 0) {
	  fprintf(stderr, "%s: bad partitioning depth `%s'\n", 
		  argv[0], &argv[i][2]);
	  cmderr = TRUE;
	  break;
	}
	else *autlev = OFF;*/
      }
      /* Might want grid order.
      else if (argv[i][1] == 'G') {
        if(sscanf(&(argv[i][2]), "%d", gsize) != 1) cmderr = TRUE;
        else if (*gsize < 2 || *gsize > 8) cmderr = TRUE;
      }*/
      /* Might want wavenumber or frequency.
      else if(argv[i][1] == 'W') {
	*wavenum = (double) strtod(&(argv[i][2]), chkp);
      }*/
      /* Who knows what else... 
      else if(argv[i][1] == 'n') {
	*pan_per_body = (int) strtol(&(argv[i][2]), chkp, 10);
      }*/
      else {
	fprintf(stderr, "%s: illegal option -- %s\n", argv[0], &(argv[i][1]));
	cmderr = TRUE;
	break;
      }
    }
    /* This isn't an option, so it must be an input file. */
    else {
      if (filesSoFar==0){
	filesSoFar++;
	strcpy(sourceInputFile,delcr(argv[i]));
	fprintf(stdout, "%s: Will read source elements from file %s\n", 
		argv[0],sourceInputFile);
      }
      else if(filesSoFar==1){
	filesSoFar++;
	strcpy(evalInputFile,delcr(argv[i]));
	fprintf(stdout, "%s: Will read evaluation elements from file %s\n", 
		argv[0],evalInputFile);
      }
      else {
	fprintf(stderr, "%s: illegal input file: %s\n"
		"   Maximum number of input files exceeded.\n"
		"   Already included files %s and %s.\n", 
		argv[0], argv[i],
		sourceInputFile,
		evalInputFile);
	cmderr = TRUE;
	break;
      }
    }
  }

  /* If no errors were reported check the input files                  */
  if (!cmderr){
    /* At least one input file MUST be specified. 
     * Only print this error message if there are no other errors.     */
    if (filesSoFar < 1){
      fprintf(stderr, "%s: No input file(s) specified!\n",argv[0]);
      cmderr = TRUE;
    } 
    /* If only one input file is specified assume that sources and 
     * evals are specified by the same file                            */
    else if (filesSoFar == 1){
      fprintf(stdout, "%s: Evaluation elements and source elements "
	      "will be assumed equal,\n",argv[0]);
      strcpy(evalInputFile,sourceInputFile);
      /* Check that the file exists */
      testInputFile = fopen(sourceInputFile,"r");
      if (testInputFile == NULL){
	fprintf(stderr, "%s: Could not open input file %s\n",
		argv[0],sourceInputFile);
	cmderr = TRUE;
      }
      else fclose(testInputFile);
    }
    else if (filesSoFar == 2){
      /* Check that the input files exist                              */
      testInputFile = fopen(sourceInputFile,"r");
      if (testInputFile == NULL){
	fprintf(stderr, "%s: Could not open source input file %s\n",
		argv[0],sourceInputFile);
	cmderr = TRUE;
      }
      else {fclose(testInputFile);
      }
      testInputFile = fopen(evalInputFile,"r");
      if (testInputFile == NULL){
	fprintf(stderr, "%s: Could not open eval input file %s\n",
		argv[0],evalInputFile);
	cmderr = TRUE;
      }
      else {fclose(testInputFile);
      }
    }
  }
  /* Put a dumpConfig here to tell the user what legal 
   * input looks like.                                                 */
  if(cmderr == TRUE) {
    fprintf(stderr, "\n"
	    "A valid command line has the format:\n"
	    "  %s [-h] [-i[s/e]#] <src.qui> [<eval.qui>]\n\n"
	    "   (Presently only ONE file accepted - "
	    "do not specify <eval.qui>.)\n"
	    "Use:    \"%s -h\"   for further help.\n\n",
	    argv[0],argv[0]);
    exit(0);
  }
} /* End of routine parseCommandLine */

