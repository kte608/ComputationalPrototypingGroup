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
/* FILE: Global.h
*
* Globally available macros, type definitions and structures.
* These definitions and headers should be available everywhere 
* internally in the pfft code.
*
*/

/* Only include this file if it hasn't been included before: */
#if ! defined globalDotHIsIncluded
#define globalDotHIsIncluded

#include <assert.h>

/* Include headers in pfft.h: */
#include "pfft.h"

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

# define MAX(a,b) ( ((a)>(b)) ? (a) : (b) )
# define MIN(a,b) ( ((a)>(b)) ? (b) : (a) )

#define ABS(A) ( ( (A) > 0 ) ? (A) : (-(A)) )

#define FREE(p) { assert((p)!=NULL); free(p); (p)=NULL;}

#define IFFREE(p) { if( (p)!=NULL ) FREE(p) }

/* For doubles of length 3 - vectors in real space. 
 * These doubles of length 3 are defined as type "point"               */
# define DOTPROD(V1,V2) (V1[0]*V2[0]+V1[1]*V2[1]+V1[2]*V2[2])

# define VNORM(V1) sqrt(V1[0]*V1[0]+V1[1]*V1[1]+V1[2]*V1[2])

# define VDIST2(V1,V2) \
      (V1[0]-V2[0])*(V1[0]-V2[0]) \
     +(V1[1]-V2[1])*(V1[1]-V2[1]) \
     +(V1[2]-V2[2])*(V1[2]-V2[2])

# define VDIST(V1,V2) \
sqrt( (V1[0]-V2[0])*(V1[0]-V2[0]) \
     +(V1[1]-V2[1])*(V1[1]-V2[1]) \
     +(V1[2]-V2[2])*(V1[2]-V2[2]))
/* VDIST(V1,V2) = sqrt(VDIST2(V1,V2)) */

#define VOP(V,V1,op,V2){ \
V[0] = V1[0] op V2[0];\
V[1] = V1[1] op V2[1];\
V[2] = V1[2] op V2[2];}

#define XPROD(V,V1,V2){ \
V[0] = V1[1]*V2[2]-V1[2]*V2[1]; \
V[1] = V1[2]*V2[0]-V1[0]*V2[2]; \
V[2] = V1[0]*V2[1]-V1[1]*V2[0];}

#define VASSIGN(V, V1){ \
V[0] = V1[0]; \
V[1] = V1[1]; \
V[2] = V1[2];}

#define PI 3.141592653589793


/* 
* ====================================================================
* TYPE DEFINITIONS (structures)
* ==================================================================== */ 
typedef double point[3];
/* The accuracy of different parts of the code could be set here.      */


/* 
====================================================================== 
  PROTOTYPES
   for functions available globally within this program
====================================================================== */ 

/* BASIC math (and other) routines (in stdmath.c) */
int intpow (int base, int n);
double dblRound (double);

void intSort(int* Array, int n);
void intQuickSort(int* Array, int nStart, int nEnd);

void intSort2(int* Array, int* Array2, int n);
void intQuickSort2(int* Array, int* Array2, int nStart, int nEnd);

void doubleSort2(double* Array, int* Array2, int n);
void doubleQuickSort2(double* Array, int* Array2, int nStart, int nEnd);

int setupCubatureRuleSquare(int degree, int irule, /* NOT TESTED */
			    double *x, double *y, double *w);

int setupCubatureRuleTriangle(int degree, int irule, 
			      double *a, double *b, double *c, 
			      double *w);

void setupGaussLegendreRule(int order,  /* NOT TESTED */
			    double *points, double *weights);

void plotProgressMark(int iProgress, int iReset,
		      int iFirst,int iLast,
		      int iMarkLast,char mark);


/* ELEMENT routines (in elements.c) */
void calcNumNonzeroDirectEntries(void **sourceIn,
				 void **evalsIn, 
				 int nSources, 
				 int nTotNeighbours,
				 void *colsData);
void calcNumNonzeroPrecondEntries(void **sourceIn,
				  void **evalsIn, 
				  int nSources,
				  int nTotNeighbours,
				  void *colsData);
void *pointerToSource(void *sourcesIn,
		      int iSource);
void *pointerToEval(void *evalsIn,
		    int iEval);
void findSourceBoundingSphere(point center,
			      double *radius,
			      void *sourceIn);
void findEvalBoundingSphere(point center,
			    double *radius,
			    void *evalIn);
void calcElementMoments(void *elementIn, 
			char *elementType,
			char *systemType,
			double xc, double yc, double zc,
			double polyScale, 
			double **moments,
			int terms, int polyOrder, int maxMonoOrder1D,
			int *imonoOrder,int *jmonoOrder,int *kmonoOrder);

int numEvalUnknowns(void *elementIn);

int *evalIndices(void *elementIn);

int numSourceUnknowns(void *elementIn);

int *sourceIndices(void *elementIn);

void addSelfTermConstant(void **sourcesIn, int nSources,
			 void **evalsIn,   int nEvals,
			 void *mat,        int matNum,
			 char *type,       double value);

/* GRID routines (in grid.c) */
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
				 void *targetData);

void elementAndGridOperate(void *gridInput,
				char *Type,
				void (*targetFunction)(void *grid,
						       void *element,
						       char *Type,
						       int i, int j, int k,
						       void *tData),
				void *targetData);

int getGridInterp(void *gridIn, 
		  char *Type, 
		  int **piv, int **pjv, int **pkv,
		  double ***pwts,
		  int **imonoOrder, int **jmonoOrder, int **kmonoOrder, 
		  int *pterms,
		  double *polyScale,
		  int *polyOrder, int *maxOrder1D);

int numInterpPoints(void *gridInput);

void getGridLocation(void *gridIn,
		     int i, int j, int k,
		     double *px, double *py, double *pz);

int gridDataIndex(void *gridIn,
		  int i, int j, int k);

void grid2gridCalculate(void *gridIn);

void precorrect(void *gridIn, void *dmat, void *pmat, void *imat); /* OBSOLETE */

void projectOntoGrid(void *gridIn,double *src,void *projectMat);

void interpolateFromGrid(void *gridIn,double *dest,void *interpMat);

int countPossibleNeighbourLinks(void *gridIn,
				int returnZeroIfAssociated, 
				char *Type);
void supplyKernelValues(void *gridInput, 
			int nProjectPoints, int *projectPoints,
			int nInterpPoints, int *interpPoints,
			double **kernelValuesArray);

void precorrectThroughGridOld(void *directMatIn, /* OBSOLETE */
			      void *projectMatIn, 
			      void *interpMatIn,
			      void *gridIn);

/* Grid Data Routines (in gridData.c). */

/* FFT Routines. (in nfft.c) */

double precorrectVal(void *kernelIn, void *griddataIn, /* OBSOLETE */
		   double *projectVals, int *projectGpoints, int projectPoints,
		     double *interpVals, int *interpGpoints, int interpPoints); 

void supplyKernelDataValues(void *kernelIn, void *griddataIn,
			    int nProjectPoints, int *projectPoints,
			    int nInterpPoints, int *interpPoints,
			    double **kernelValuesArray);

void rfft3(void *fftdata);

double griddataMemoryEstimate(int nx, int ny, int nz, int verbose);

/* SPARSE matrix routines (in sparse.c) */
void *spMatAlloc(int n, 
		 int *colSizes, 
		 int numMats);
void spMatColumnStore(void *matIn, int matNum, int iCol, 
		      int nrows, int *rows, double *vals, 
		      char *rowStructure,
		      char *oldValues);
void spMatColumnGet(void *matIn, int matNum, int iCol, 
		    int *nrows, int **rows, double **vals, 
		    char *rowStructure);
void spMatEntryAddto(void *matIn, int matNum, int iCol, int iRow, 
		     double val);
void spMatMultiply(void *mat, 
		   double *src, 
		   double *dest, 
		   int matNum);
void spTransposeMatMultiply(void *mat, 
			    double *src, 
			    double *dest, 
			    int matNum);

double spMatMemoryEstimate(int n, int nEntries, int numMats);

void precorrectColumns(int nCols, int *dCols,
		       void *directMatIn, 
		       void *projectMatIn, 
		       void *interpMatIn,
		       void *grid);

void precorrectColumnsOld(int nCols, int *dCols, 
			  void *directMatIn, 
			  void *projectMatIn, 
			  void *interpMatIn,
			  void *grid);/* OBSOLETE */
void precorrectPartialColumns(int nCols, int *dCols,
			      int nRows, int *dRows,
			      double **precVals,
			      void *directMatIn, 
			      void *projectMatIn, 
			      void *interpMatIn,
			      void *grid);
void colPrec(void *directMatIn,  
	     void *projectMatIn, 
	     void *interpMatIn,
	     void *grid);/* OBSOLETE */

void precorrectColumn(int dCol,  
		      void *directMatIn, 
		      void *projectMatIn, 
		      void *interpMatIn,
		      void *grid);/* OBSOLETE */

void dumpSparseMatrixTranspose(void *matIn,int matNum); /* (For testing) */


/* DENSE matrix routines (in dense.c) */
void LUdecompNoPivoting(double *mat,int n);
void LUbacksolveNoPivoting(double *mat,int n,double* b);
int pinv(double **A, int rows, int cols);
void dumpDenseMatrix(double *mat,int n); /* (For testing) */

/* DSVDC (singular value decomposition) routines */
int dsvdc_(double *x, int *ldx, int *n, int *p, 
	   double *s, double *e, double *u, int *ldu, 
	   double *v, int *ldv, double *work, int *job, int *info);

/* DIRECT part routines (in direct.c) */
void solveDirectSparse(void *directMat, 
		       void *precondMat, 
		       char *precondType,
		       int nColsDirect);

void directMatrixMultiplyNoPrec(double *src,
				double *dest,
				int size,
				void *inputData);
void directMatrixMultiplyRightPrec(double *src,
				double *dest,
				int size,
				void *inputData);
void directMatrixMultiplyLeftPrec(double *src,
				  double *dest,
				  int size,
				  void *inputData);
void calcDirectMat_1(void **sources, 
		     void **evals, 
		     int nSources,
		     int nTotNeighbours,
		     void *bundleData);
void precorrectDirectMat_1(void **sourcesIn, 
			   void **evalsIn, 
			   int nSources,
			   int nTotNeighbours,
			   void *bundleIn);
void precorrectDirectMat_2(void **sourcesIn, 
			   void **evalsIn, 
			   int nSources,
			   int nTotNeighbours,
			   void *bundleIn);
void calcPrecondMat_1(void **sourcesIn,
		      void **evalsIn, 
		      int nSources,
		      int nTotNeighbours,
		      void *twoMatsIn);
void precorrectMat(void *kernel, void *griddata,
		   void *directmat, void *projectmat, 
		   void *interpmat); /* OBSOLETE */

double directPartMemoryEstimate(void *grid, 
				int nCols, int nRows, int numMats,
				char *Type);

/* GMRES routines (in gmres.c) */
int gmres(double *solutionVector, double *rightHandSide,
	  int size, double tol, int maxiter,
	  void (*matmulFunction)(double *src,
				 double *dest,
				 int size,
				 void *mData),
	  void *matmulData);

/* SOLVE part routines (in solve.c) */
void fullMatrixMultiply(double *src,double *dest,
			int size,
			void *inputData);
void freeSystemVector(void *vector);

/* BOUNDARY CONDITION part routines (in bconds.c) */
int boundaryConditionType(point coordinate, 
			  int iboundary,
			  char *conditionType,
			  int idata);
double boundaryConditionValue(int bcType, 
			      point coordinate, double *normal, 
			      int iboundary,
			      char *conditionType,
			      int idata);

/* ==== TEST FUNCTIONS (not really needed) ==== */ 
void testGrid2grid(void *gridIn);

void testProject2grid(void *gridIn, void *projectMat, 
		      void *interpMat, int ncols);

void testLocalMoments(void);

void testMoments(void);

void testPinv(void);

/* This is the last statement to make sure that 
   all of the above is only included once */
#endif
/* ====== DONT PUT ANYTHING BELOW THIS LINE ====== */
