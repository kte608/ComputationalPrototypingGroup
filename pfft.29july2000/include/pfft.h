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
/* FILE: pfft.h
 *
 * Header files for routines that are intended to be callable 
 * from external programs.
 * If you plane to use the existing code with very few modifications, 
 * then call these routines directly.
 * Most likely you should modify the routines in bconds.c as well.
 */

/* Only include this file if it hasn't been included before: */
#if ! defined pfftDotHIsIncluded
#define pfftDotHIsIncluded

/* BASIC routines (in stdmath.c) */
char *delcr(char *str);
void plotSizeOf(void);

/* ELEMENT routines (in elements.c) */
void *setupSources(int *nSourcesReturn,
		   int *maxBoundaryPart, 
		   char *sourceInputFileName,
		   int *inputFileTypes,
		   char *conditionType,
		   int idata);

void *setupEvals(int *nEvalsReturn,
		 char *evalInputFileName,
		 int *inputFileTypes); /* NOT IMPLEMENTED! */

void assignElementDirectIndices(void *sources, int nSources,
				void *evals,   int nEvals,
				int *nColsDirect, int *nRowsDirect);

void setupElementBoundaryConditionValues(void *sourcesIn,
					 int nSources,
					 void *rightHandSideIn,
					 char *conditionType,
					 int idata);

void swapSignVectorElement(void *elementsIn, int nElements,
			   double *vector, char *type);

double integrateOverElements(void *sourcesIn, int nSources, 
			     void *evalsIn, int nEvals, 
			     double *solution,
			     char *integrationType,
			     double *integrationOnPart,
			     int ipart);

/* GRID routines (in grid.c) */
void *setupGrid(void *sources,int nSources,
		void *evals,  int nEvals,
		int nCols, int nRows, int numMats);

void grid2gridSetup(void *gridIn);

void precorrectThroughGrid(void *directMatIn, 
			   void *projectMatIn, 
			   void *interpMatIn,
			   void *gridIn,
			   int n);

void grid2gridFinish(void *gridIn);

void grid2gridFree(void *gridIn);

/* Grid Data Routines (in gridData.c). */
void *callocElementGridSparse(void *grid, char *Type, 
			      int n, int numMats);

void calcElementGridSparse(void *grid, char *elementType, 
			   char *systemType, void *mat);

/* SPARSE matrix routines (in sparse.c) */
void spMatReinit(void *matIn, int matNum);

/* DIRECT part routines (in direct.c) */
void *callocDirectSparse(void *grid, 
			int nCols, int nRows, 
			int numDirectMats);

void calcDirectSparse(void *grid, 
		      void *directMat,
		      int nSources,
		      char *type);

void calcPrecorrectSparse(void *directMat,
			  void *projectMat,
			  void *interpMat,
			  void *grid,
			  int nSources,
			  char *type);
void *callocPrecondSparse(void *grid, 
			  void *directMat, 
			  int maxNeighbours,
			  int nCols, int nRows, 
			  int numPrecondMats);

void calcPrecondSparse(void *grid, 
		       void *directMat, 
		       void *precondMat,
		       int maxNeighbours);


/* SOLVE part routines (in solve.c) */
void *allocateSystemVector(int nEntries);

void systemMatrixMultiply(void *directMat, void *interpMat, 
			  void *projectMat, void *precondMat, 
			  char *precondType, void *grid,
			  double *sourceVector, double  *destVector,
			  int matNum, 
			  int numColsDirect, int numRowsDirect);

void solveFull(void *directMat, void *interpMat, void *projectMat,
	       void *precondMat, char *precondType, void *grid,
	       double *sourceVector, double  *destVector,
	       int matNum, 
	       int numColsDirect, int numRowsDirect);

/* ==== TEST FUNCTIONS =====
 * (not really needed - or may need reimplementation) */ 
void testSolution(void *evalsIn, int nEvals, double *solution,
		  char *bcType, int bcIdata);

void tableElements(void **sources,int nSources);

void dumpElementsToTecplot(void *elementsIn, 
			   int nElements,
			   int nParts,
			   int *parts);


#endif
/* ====== DONT PUT ANYTHING BELOW THIS LINE ====== */

