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
/* FILE: sparse.c
* 
* Routines for definitions of and operations on/with sparse 
* matrices. 
*
*/

#include <stdlib.h>
#include <stdio.h>

#include "Global.h"
#include "Sparse.h"

/* 
* ==================================================================== 
* Allocates a sparse matrix of size n, with colSizes nonzeros per column.
* Note that a sparse matrix can have several sets of values for a given
* structure by specifying numMats to be larger than 1.
* The sparse matrix is set up by columns to increase the accesibility
* of the structure on a source-by-source basis. Recall that each 
* source element will correspond to one or more columns of the 
* direct part matrix, and that each of the neighbouring evals will
* correspond to one or more (nonzero) rows in these columns.
* Obviously, the sparse matrix structure may be used not only for 
* the sparse direct part.
* Note that the routine returns a generic type pointer (void*),
* rather than a pointer of type sparseMat.
* Further, note that the indices of the nonzero rows in each column
* does not (necessary) have to be sorted. 
* ==================================================================== */
void *spMatAlloc(int n, int *colSizes, int numMats)
{
  char fctName[]="spMatAlloc";  /* Name of this routine                */
  sparseMat *mat;
  int j, im, nonzeros, *colIndex;

  /* Allocate a sparse matrix with n columns. An extra column entry is 
   * used to denote the end (size+1) of the array of values (matVals). */
  mat = (sparseMat *) calloc(1, sizeof(sparseMat));
  mat->n = n;
  mat->colIndex = (int *) calloc(n+1, sizeof(int));
  colIndex = mat->colIndex ; 

  /* For each column set the first entry number in the array of values 
   * (matVals). Also count the total number of nonzero entries: */
  for(colIndex[0] = 0, nonzeros = 0, j = 0; j < n; j++) {
    colIndex[j+1] = colIndex[j] + colSizes[j];
    nonzeros += colSizes[j]; 
  } 
  /* Note that (hopefully) colIndex[n+1]=nonzeros */
  

  /* Allocate array to store the row number for each entry 
   * in each column */
  mat->rowIndex = (int *) calloc(nonzeros, sizeof(int));
  /* Allocate array to store the values of each matrix for each entry 
   * in each column */
  mat->nMats = numMats;
  mat->matVals = (double **) calloc(numMats, sizeof(double *));
  /* Allocate arrays to store the values of each entry 
   * in each matrix */
  for(im = 0; im < numMats; im++) {
    mat->matVals[im] = (double *) calloc(nonzeros, sizeof(double));
  }
  fprintf(stdout,"%s: Allocated %d integers and %d reals\n",
	  fctName,n+nonzeros,numMats*nonzeros);
  return (void*) mat;
}

/* 
* ==================================================================== 
* Re-initialize (zero out) the entries of a sparse matrix, but keep 
* the structure.
* ==================================================================== */
void spMatReinit(void *matIn, int matNum)
{
  char fctName[]="spMatReinit";  /* Name of this routine               */
  sparseMat *mat;
  double *vals;
  int i, istart, istop;

  /* Cast input data to correct type */
  mat = (sparseMat*) matIn;

  /* Check validity of number of matrix values */
  if (matNum<0 || matNum>= mat->nMats){
    printf("%s: This set of matrix values does not exist: %d\n",
	   fctName,matNum);
    _EXIT_;
  }

  /* Shorthand notation */
  vals = mat->matVals[matNum];

  /* Zero out this set of matrix values */
  istart = mat->colIndex[0];
  istop  = mat->colIndex[mat->n];
  for (i=istart; i<istop; i++) vals[i]=0.0;

  return;
} /* End routine spMatReinit */

/* 
* ==================================================================== 
* Store a column in a sparse matrix. 
* The matrix must already be initialized. 
* 
* matIn:   Matrix to operate on.
* matNum:  Sets of matrix values to operate on, i.e. mat->matVals[matNum].
* iCol:    Column number to store data in (0, ... , mat->n-1)
* rows:    Row numbers for rows with (possibly) non-zero values.
* vals:    Values of the entries of the rows, i.e. vals[i] contains the 
*          matrix entry (icol(+1), rows[i](+1)) [The "+1" notation to 
*          signify that entries are numbered from 0 to n-1, rather than 
*          from 1 to n.
* rowStructure: 
*          Determines what to do with existing matrix structure for 
*          this column.
*          "overwrite": Overwrite existing row numbers uncritically.
*          "keep":      Do not overwrite row numbers. Keep the old.
*          "check":     Check that the new and the old row numbers
*                       match (and are sorted in the same way).
* oldValues:
*          Determines what to do with existing matrix values for 
*          this column.
*          "discard": Overwrite existing values uncritically.
*          "addTo":   Add the new values to the old.
*          "subtractFrom":  Subtract the new values from the old.
* ==================================================================== */
void spMatColumnStore(void *matIn, int matNum, int iCol, 
		      int nrows, int *rows, double *vals, 
		      char *rowStructure,
		      char *oldValues)
{
  char fctName[] = "spMatColumnStore";
  int i, iColStart, iColStop, *rowIndex;
  double *matVals;
  sparseMat *mat;

  /* Cast input matrix to correct type                                 */
  mat = (sparseMat *) matIn;

  /* Short hands:                                                      */
  iColStart = mat->colIndex[iCol];
  iColStop  = mat->colIndex[iCol+1];
  rowIndex  = mat->rowIndex;
  matVals   = mat->matVals[matNum];

  /* Check that the matrix value number is in the correct range        */
  assert(matNum>=0 && matNum<mat->nMats);
  /* Check that the matrix column number is in the correct range       */
  assert(iCol>=0 && iCol<mat->n);
  /* Check that the total number of row entries is not too large       */
  assert(nrows>=0 && nrows<=iColStop-iColStart);
  /* Check that the "structure" parameter is of a recognized format    */
  if (  strcmp(rowStructure,"overwrite") &&
        strcmp(rowStructure,"keep") &&
        strcmp(rowStructure,"check") ){
    printf("%s: Unrecognized value of parameter rowStructure: \"%s\"\n",
	   fctName,rowStructure);
    printf("     Please use \"overwrite\", \"keep\", or \"check\".\n");
    _EXIT_;
  }
  if (  strcmp(oldValues,"discard")&&
        strcmp(oldValues,"addTo")&&
        strcmp(oldValues,"subtractFrom") ){
    printf("%s: Unrecognized value of parameter rowStructure: \"%s\"\n",
	   fctName,rowStructure);
    printf("     Please use \"discard\", \"addTo\", or \"subtractFrom\".\n");
    _EXIT_;
  }

  /* Verification part */
  if (strcmp(rowStructure,"check")==0){
    assert(nrows==iColStop-iColStart);
    for (i=0; i<nrows; i++) assert(rowIndex[iColStart+i]==rows[i]);
  }

  /* Store new row index structure if necessary                        */
  if (strcmp(rowStructure,"overwrite")==0){
    for (i=0;i<nrows;i++) 
      rowIndex[iColStart+i] = rows[i];
  }
  /* Store new values                                                  */
  if (strcmp(rowStructure,"overwrite")==0){
    for (i=0;i<nrows;i++) 
      matVals[iColStart+i] = vals[i];
  }
  else if (strcmp(rowStructure,"addTo")==0){ 
    for (i=0;i<nrows;i++) 
      matVals[iColStart+i] += vals[i];
  }
  else { /* strcmp(rowStructure,"subtractFrom")==0 */
    for (i=0;i<nrows;i++) 
      matVals[iColStart+i] -= vals[i];
  }

  return;
} /* End of routine spMatColumnStore */

/* 
* ==================================================================== 
* Return a column from a sparse matrix. 
* The matrix must already be initialized. 
* 
*Input:
* matIn:   Matrix to operate on.
* matNum:  Sets of matrix values to operate on, i.e. mat->matVals[matNum].
* iCol:    Column number to get data from (0, ... , mat->n-1)
*
*Output:
* rows:    Pointer to row numbers for rows with (possibly) non-zero values.
* vals:    Pointer to values of the entries of the rows, i.e. vals[i] 
*          contains the matrix entry (icol(+1), rows[i](+1)) 
*          [The "+1" notation to signify that entries are numbered from 
*          0 to n-1, rather than from 1 to n.]
*          
* rowStructure: 
*          Presently not used. See e.g. routine spMatColumnStore for 
*          possible future use.
* ==================================================================== */
void spMatColumnGet(void *matIn, int matNum, int iCol, 
		    int *nrows, int **rows, double **vals, 
		    char *rowStructure)
{
  /* char fctName[] = "spMatColumnStore"; */
  int iColStart, iColStop, *rowIndex;
  double *matVals;
  sparseMat *mat;

  /* Cast input matrix to correct type                                 */
  mat = (sparseMat *) matIn;

  /* Short hands:                                                      */
  iColStart = mat->colIndex[iCol];
  iColStop  = mat->colIndex[iCol+1];
  rowIndex  = mat->rowIndex;
  matVals   = mat->matVals[matNum];

  /* Check that the matrix value number is in the correct range        */
  assert(matNum>=0 && matNum<mat->nMats);
  /* Check that the matrix column number is in the correct range       */
  assert(iCol>=0 && iCol<mat->n);

  /* Set the number of row entries for this column */
  *nrows = iColStop-iColStart;
  /* Pointer to the row indices */
  *rows  = &(rowIndex[iColStart]);
  /* Pointer to row values */
  *vals  = &(matVals[iColStart]);

  return;
} /* End of routine spMatColumnGet */

/* 
* ==================================================================== 
* Add a specified value to a specified matrix entry given as 
* (row,column). This requires a search for the row number in the matrix 
* column. If the matrix columns are known to be sorted, then a binary 
* search could be imagined. However, since no knowledge a priori is 
* given, the row entries are scanned from one end until the correct 
* row entry is met. 
* If the row does not exist, then nothing is done. A call to assert()
* is used to signal problems in this case.
*
* Note: the present routine should be used sparingly since it is a 
* costly way to assign values to matrix entries.
* ==================================================================== */
void spMatEntryAddto(void *matIn, int matNum, int iCol, int iRow, 
		     double val)
{
  /*char fctName[]="spMatEntryAddto";*/  /* Name of this routine           */
  sparseMat *mat;
  int i, istart, istop, *rowIndex, entryNotFound;

  /* Cast input data to correct type */
  mat = (sparseMat*) matIn;


  /* Check validity of number of matrix values */
  assert(matNum>=0 && matNum< mat->nMats);

  /* Check that the column exists */
  assert(iCol>=0 && iCol < mat->n);

  /* Set flag to signal that the correct entry is not yet found.       */
  entryNotFound = 1;

  /* Short hands */
  rowIndex  = mat->rowIndex;

  istart = mat->colIndex[iCol];
  istop  = mat->colIndex[iCol+1];
  for(i = istart;  i < istop; i++) {
    if(rowIndex[i]==iRow){
      /* This is the correct entry.                                    */
      /* Add value to entry                                            */
      mat->matVals[matNum][i] += val;
      /* Set flag to notify that the entry has been found and the 
       * value added.                                                  */
      entryNotFound = 0;
      break; /* Skip remaining entries of the column.                  */
    }
  }
  /* Make sure that the entry was found */
  assert(entryNotFound==0);

  return;
} /* End of routine spMatEntryAddto */

/* 
* ================================================================= 
* Performs a sparse matrix vector multiply.
* NOTE: forms dest = dest + mat * src.
*
* Please note that the sparse matrix is passed as a generic type 
* pointer, rather than a pointer of type sparseMat. This enables 
* routines to require a matrix-vector multiply, without them 
* "knowing" what a sparse matrix "is".
* ================================================================= */
void
spMatMultiply(void *passedMat, double *src, double *dest, int matNum)
{
  double acc, *vals;
  int i, j, n, startRow, endRow, *rowIndex, *colIndex;
  sparseMat *mat;

  mat = (sparseMat*) passedMat;

  n = mat->n;
  rowIndex = mat->rowIndex;
  colIndex = mat->colIndex;
  vals = mat->matVals[matNum];

  for(j = 0; j < n; j++) {
    startRow = colIndex[j];
    endRow   = colIndex[j+1]; 
    acc = src[j];
    for(i = startRow;  i < endRow; i++) {
      dest[rowIndex[i]] += vals[i] * acc;
    }
  }
  return;
} /* End routine spMatMultiply */

/* 
* ================================================================= 
* Performs a sparse matrix transpose vector multiply.
* NOTE: forms dest = dest + mat' * src.
*
* Please note that the sparse matrix is passed as a generic type 
* pointer, rather than a pointer of type sparseMat. This enables 
* routines to require a matrix-vector multiply, without them 
* "knowing" what a sparse matrix "is".
* ================================================================= */
void
spTransposeMatMultiply(void *passedMat, double *src, double *dest, int matNum)
{
  double acc, *vals;
  int i, j, n, startRow, endRow, *rowIndex, *colIndex;
  sparseMat *mat;

  mat = (sparseMat*) passedMat;

  n = mat->n;
  rowIndex = mat->rowIndex;
  colIndex = mat->colIndex;
  vals = mat->matVals[matNum];

  for(j = 0; j < n; j++) {
    startRow = colIndex[j];
    endRow   = colIndex[j+1];

    acc = 0.0;
    for(i = startRow;  i < endRow; i++) {
      acc += vals[i] * src[rowIndex[i]];
    }
    dest[j] += acc;
  }
  return;
} /* End routine spTransposeMatMultiply */

/* 
* ==================================================================== 
* Precorrect the direct matrix based on the projection and 
* interpolation parts.
* ==================================================================== */
void precorrectMat(void *kernel,  void *griddata, 
		   void *directmatIn, void *projectmatIn, 
		   void *interpmatIn)
{
  char fctName[] = "precorrectMat";
  sparseMat *dmat, *pmat, *imat;
  double coeff;
  int dCol, dStart, dStop, pStart, pStop, i, iStart, iStop;
  int ihash, nhash=60;
  int verbose = 1;

  /* Cast input matrices to correct type */
  dmat = (sparseMat*) directmatIn;
  pmat = (sparseMat*) projectmatIn;
  imat = (sparseMat*) interpmatIn;

  if(dmat->nMats > 1) {
    printf("%s: Multiple mats not yet supported.",fctName);
    _EXIT_;  /* Just one mat for now. */
  }

  /* Go through the direct matrix, column by column 
   * (unknown by unknown).                                             */
  ihash = 0;
  if (verbose) {
    /* Show how many hash marks to expect */
    for (i=0;i<nhash;i++) fprintf(stdout," ");
    fprintf(stdout,"v\n");
  }
  for(dCol = 0; dCol < dmat->n; dCol++) {
    /* Print a line of hash marks to show progress                     */
    if (   (double)dCol /(double)(dmat->n -1)
	>= (double)ihash/(double)nhash) {
      fprintf(stdout,"#");
      fflush(stdout);
      ihash++;
    }
    /* Get the column in the direct matrix associated with source dCol.*/
    dStart = dmat->colIndex[dCol];
    dStop = dmat->colIndex[dCol+1]; 

    /* Get the projection column associated with this source.          */
    pStart = pmat->colIndex[dCol];
    pStop = pmat->colIndex[dCol+1];

    /* For each (nearby) eval in the direct column, 
     * get precorrect coeff.                                           */
    for(i = dStart; i < dStop; i++) {
      iStart = imat->colIndex[dmat->rowIndex[i]];
      iStop = imat->colIndex[dmat->rowIndex[i] + 1];
      /* Here is the tricky part, we need the grid to grid matrix that 
       * relates the projection grid points to the interpolation 
       * gridpoints for this source-eval pair.  But, we don't want to 
       * store it, we just need the single result.  Since the info 
       * needed is stored in the untransformed kernel, we pass that 
       * along with the projection and interpolation vectors associated 
       * with the source-eval pair.
       */
      coeff = precorrectVal(kernel, griddata, 
			    &(pmat->matVals[0][pStart]), 
			    &(pmat->rowIndex[pStart]), 
			    pStop-pStart,
			    &(imat->matVals[0][iStart]), 
			    &(imat->rowIndex[iStart]), 
			    iStop-iStart);
      /* Finally, we update the direct matrix entry. */
      dmat->matVals[0][i] -= coeff;
    }
			    
  }
  printf("\n");
} /* End routine precorrectMat */ 

/* 
* ==================================================================== 
* Dump a matrix to a file.
* ==================================================================== */
void dumpSparseMatrixTranspose(void *matIn,int matNum)
{
  int i,iCol,iRow,iRowSp,icolstart, iStop;
  int *rowIndex, *colIndex;
  double value, *vals;
  sparseMat *mat;
  static int icall=0;
  FILE *dumpFile;
  char *dumpFileName;

  dumpFileName = (char*) calloc(10,sizeof(char));
  sprintf(dumpFileName, "dmp%d.out", icall);

  fprintf(stdout," Dumping a transpose of a sparse matrix to file %s\n",
	  dumpFileName);
  dumpFile = fopen(dumpFileName,"w");
  mat = (sparseMat*) matIn;

  /* The following dumps all entries - including zero ones! */
  if (1==1){
    for (iCol=0;iCol<mat->n;iCol++){
      icolstart=mat->colIndex[iCol];
      for (iRow=0,iRowSp=0;iRow<mat->n;iRow++){
	if (iRow==mat->rowIndex[icolstart+iRowSp]){
	  value = mat->matVals[matNum][icolstart+iRowSp];
	  iRowSp++;
	}
	else value = 0.0;
	fprintf(dumpFile," %.20e",value);
      }
      fprintf(dumpFile,"\n");
    }
  }
  else{
  /* This dumps only the non-zero entries (no information on 
   * position given - only values:                             */
    /* Changed routine to one that gives rows and cols.*/

    rowIndex = mat->rowIndex;
    colIndex = mat->colIndex;
    vals = mat->matVals[matNum];

    for(iCol = 0; iCol < mat->n; iCol++) {
      fprintf(dumpFile, "printing column %d\n", iCol);
      icolstart=mat->colIndex[iCol];
      iStop   = colIndex[iCol+1]; 
      for(i = icolstart;  i < iStop; i++) {
	fprintf(dumpFile, "%d %16.8e   ", rowIndex[i], vals[i]);
      }
      fprintf(dumpFile, "\n");
    }
  }
  fclose(dumpFile);
  icall=icall+1;
  return;
} /* End of routine dumpSparseMatrixTranspose */

/* 
* ==================================================================== 
* Estimate the memory consumption for setting up a sparse matrix 
* with n columns and nEntries nonzero entries in each of numMats
* matrices.
* ==================================================================== */
double spMatMemoryEstimate(int n, int nEntries, int numMats)
{
  double memCount;
  memCount = 0.0;
  /* Memory for storing column start and end */
  memCount += (double) ((n+1) * sizeof(int));
  /* Memory for storing the row number of each entry */
  memCount += ((double) nEntries) * ((double) sizeof(int));
  /* Memory for storing values of each entry */
  memCount += ((double) numMats * sizeof(double)) * ((double) nEntries) ;
  return memCount;
}

/* 
* ==================================================================== 
* Precorrects the nCols columns in the direct matrix.
* The efficiency of the present routine relies on the corresponding 
* columns in the projection matrix to be structurally similar, if not 
* structurally identical. The VALUES of the projection entries may 
* differ significantly from column to column.
* ==================================================================== */
void precorrectColumns(int nCols, int *dCols,
		       void *directMatIn, 
		       void *projectMatIn, 
		       void *interpMatIn,
		       void *grid){
  char fctName[]="precorrectColumns";  /* Name of this routine        */
  sparseMat *dmat, *pmat, *imat;
  int dStart, dStop, pStart, pStop, iStart,iStop,
    i,j,ir,ic, np, 
    nInterpPoints, nInterpPointsPerCol, nDistinctInterpPoints, 
    nProjPoints,   nProjPointsPerCol,   nDistinctProjPoints, 
    iRow, iipoint, ippoint, iiLoc, ipLoc, nRows, 
    dCol;
  double pval, ival, cval;
  static int 
    maxInterpPoints=0, *interpPoints=NULL, *interpPoints2=NULL, 
    *localInterpPoints=NULL,
    maxDistinctInterpPoints=0, 
    maxProjPoints=0, *projPoints=NULL, *projPoints2=NULL, 
    *localProjPoints=NULL,
    maxDistinctProjPoints=0;
  static double **localKernel, *pTimesG2G;
  int verbose = 0;

  /* If negative number of columns is supplied then 
   * simply free memory and return                                     */
  if (nCols<0){
    IFFREE(projPoints);
    IFFREE(projPoints2);
    IFFREE(localProjPoints);
    IFFREE(interpPoints);
    IFFREE(interpPoints2);
    IFFREE(localInterpPoints);
   if (localKernel!=NULL){
      for (i=0; i<maxProjPoints; i++) FREE(localKernel[i]);
      FREE(localKernel);
    }
    IFFREE(pTimesG2G);
    return;
  }

  if (verbose>1) {
    printf("%s: Precorrecting columns:",fctName);
    for (ic=0; ic<nCols; ic++)  printf(" %d",dCols[ic]);
    printf("\n");
  }
  /* Cast input matrices to correct type */
  dmat = (sparseMat*) directMatIn;
  pmat = (sparseMat*) projectMatIn;
  imat = (sparseMat*) interpMatIn;

  if(dmat->nMats > 1) {
    printf("%s: Multiple mats not yet supported.",fctName);
    _EXIT_;  /* Just one mat for now. */
  }

  /* COMPRESSION OF PROJECTION POINTS */
  /* Count the total number of interpolation points 
   * (including duplicates) and set up memory.                         */
  nProjPoints = 0;     /* Total number of projection points here       */
  nProjPointsPerCol=0; /* Max # proj. points per column (sources)      */
  for (ic=0; ic<nCols; ic++){
    /* Source number (column number in the direct matrix)              */
    dCol = dCols[ic];
    /* Get number of projection points for this column number in the 
     * projection matrix                                               */
    np    = pmat->colIndex[dCol+1]-imat->colIndex[dCol];
    nProjPoints += np;
    nProjPointsPerCol = MAX(nProjPointsPerCol,np);
  }
  if (verbose>1){
    printf("%s: %d projection points to be treated.\n",
	   fctName,nProjPoints);
  }

  /* Setup memory if needed */
  if (nProjPoints > maxProjPoints){
    IFFREE(projPoints);
    IFFREE(projPoints2);
    IFFREE(localProjPoints);
    if (verbose) printf("%s: Allocating memory, part Ia\n",fctName);
    /* Memory for sorting projection points                            */
    maxProjPoints   = MAX(2*maxProjPoints,nProjPoints);
    projPoints      = (int *) calloc(maxProjPoints,sizeof(int));
    projPoints2     = (int *) calloc(maxProjPoints,sizeof(int));
    localProjPoints = (int *) calloc(maxProjPoints,sizeof(int));
  }

  /* Find (few?) related grid points from projection matrix            */
  /* Store numbers of interpolation points in array for sorting        */
  ippoint = 0;
  for (ic=0; ic<nCols; ic++){
    /* Source number (column number in the direct matrix)              */
    dCol = dCols[ic];
    /* Get projection points for this column number in the 
     * projection matrix                                               */
    pStart = pmat->colIndex[dCol];
    pStop  = pmat->colIndex[dCol+1];
    for (j=pStart; j<pStop; j++, ippoint++){
      projPoints[ippoint]  = pmat->rowIndex[j];
      projPoints2[ippoint] = ippoint;
      /*printf("(%d,%d): %d  %d\n",i,j,ipoint,interpPoints[ipoint]);*/
    }
  }
  assert(ippoint==nProjPoints); 

  /* Sort the projection points (keep track of original relations)     */
  intSort2(projPoints, projPoints2, nProjPoints);
  /* Eliminate duplicate points while keeping track of relations to 
   * original entries. I know this is a bit messy, but here is a little 
   * help:
   *   "projPoints" is the sorted array. Duplicate entries will 
   *                  be eliminated here.
   *   "projPoints2" contains the indices of "interpPoints" relating
   *                  to before the sorting. Thus, interpPoints2 is used 
   *                  to keep track of the sorting process.
   *   "localProjPoints" is used for storing the local indices (into 
   *                  interpPoints), but ordered in the same way as 
   *                  interpPoints was BEFORE the sorting.
   * Thus, if the original index (before sorting) is "i", then 
   * localInterpPoints[i] contains the index of the point in the sorted 
   * and eliminated array. (I hope).                                   */
  /* Initialize array to capture errors (this should eventually 
   * not be necessary) */
  for (i=1; i<nProjPoints; i++) localProjPoints[i] = -1;
  /* Always keep the first point:                                      */
  ippoint = 0;
  localProjPoints[projPoints2[0]] = 0;
  for (i=1; i<nProjPoints; i++){
    if (projPoints[i]!=projPoints[i-1]){
      assert(projPoints[i]>projPoints[i-1]);
      ippoint++;
      projPoints[ippoint] = projPoints[i];
    }
    localProjPoints[projPoints2[i]] = ippoint;
  }
  nDistinctProjPoints = ippoint+1;
  /* Projection points are now compressed                              */

 /* Make sure that all points have been visited (this should 
   * eventually not be neccesary) */
  for (i=1; i<nProjPoints; i++) assert(localProjPoints[i]>=0);
  if (verbose>1) 
    printf("Of these %d are distinct points:\n",nDistinctProjPoints);

  /* COMPRESSION OF INTERPOLATION POINTS */
  /* Count the total number of interpolation points 
   * (including duplicates) and set up memory.                         */
  nInterpPoints = 0;     /* Total number of interpolation points here  */
  nInterpPointsPerCol=0; /* # interp. points per column (eval)         */
  for (ic=0; ic<nCols; ic++){
    dCol = dCols[ic];
    /* Get the column in the direct matrix associated with source dCol.*/
    dStart = dmat->colIndex[dCol];
    dStop  = dmat->colIndex[dCol+1]; 

    for (i=dStart; i<dStop; i++){
      iRow  = dmat->rowIndex[i];
      np    = imat->colIndex[iRow+1]-imat->colIndex[iRow];
      nInterpPoints += np;
      nInterpPointsPerCol = MAX(nInterpPointsPerCol,np);
    }
  }
  if (verbose>1){
    printf("%s: %d interpolation points to be treated\n",
	   fctName,nInterpPoints);
  }
  /* Setup memory if needed */
  if (nInterpPoints > maxInterpPoints){
    IFFREE(interpPoints);
    IFFREE(interpPoints2);
    IFFREE(localInterpPoints);
    if (verbose) printf("%s: Allocating memory, part I\n",fctName);
    /* Memory for sorting interpolation points                         */
    maxInterpPoints = MAX(2*maxInterpPoints,nInterpPoints);
    interpPoints    = (int *) calloc(maxInterpPoints,sizeof(int));
    interpPoints2   = (int *) calloc(maxInterpPoints,sizeof(int));
    localInterpPoints = (int *) calloc(maxInterpPoints,sizeof(int));
  }

  /* Find (many) related grid points from interpolation matrix         */
  /* Store numbers of interpolation points in array for sorting        */
  iipoint = 0;
  for (ic=0; ic<nCols; ic++){
    dCol = dCols[ic];
    /* Get the column in the direct matrix associated with source dCol.*/
    dStart = dmat->colIndex[dCol];
    dStop  = dmat->colIndex[dCol+1]; 

    for (ir=dStart; ir<dStop; ir++){
      iRow = dmat->rowIndex[ir];
      iStart = imat->colIndex[iRow];
      iStop  = imat->colIndex[iRow+1];
      for (j=iStart; j<iStop; j++, iipoint++){
	interpPoints[iipoint]  = imat->rowIndex[j];
	interpPoints2[iipoint] = iipoint;
	/*printf("(%d,%d): %d  %d\n",i,j,iipoint,interpPoints[iipoint]);*/
      }
    }
  }
  assert(iipoint==nInterpPoints); 
  /* Sort the grid points (keep track of original relations)           */
  intSort2(interpPoints, interpPoints2, nInterpPoints);
  /* Eliminate duplicate points while keeping track of relations to 
   * original entries. I know this is a bit messy, but here is a little 
   * help:
   *   "interpPoints" is the sorted array. Duplicate entries will 
   *                  be eliminated here.
   *   "interpPoints2" contains the indices of "interpPoints" relating
   *                  to before the sorting. Thus, interpPoints2 is used 
   *                  to keep track of the sorting process.
   *   "localInterpPoints" is used for storing the local indices (into 
   *                  interpPoints), but ordered in the same way as 
   *                  interpPoints was BEFORE the sorting.
   *   "interpVals"   Contains the values of the 
   * Thus, if the original index (before sorting) is "i", then 
   * localInterpPoints[i] contains the index of the point in the sorted 
   * and eliminated array. (I hope).                                   */
  /* Initialize array to capture errors (this should eventually 
   * not be necessary) 
  for (i=1; i<nInterpPoints; i++) localInterpPoints[i] = -1;*/
  /* Always keep the first point:                                      */
  iipoint = 0;
  localInterpPoints[interpPoints2[0]] = 0;
  for (i=1; i<nInterpPoints; i++){
    if (interpPoints[i]!=interpPoints[i-1]){
      assert(interpPoints[i]>interpPoints[i-1]);
      iipoint++;
      interpPoints[iipoint] = interpPoints[i];
    }
    localInterpPoints[interpPoints2[i]] = iipoint;
  }
  nDistinctInterpPoints = iipoint+1;
  /* Interpolation points are now compressed */

 /* Make sure that all points have been visited (this should 
   * eventually not be neccesary) 
  for (i=1; i<nInterpPoints; i++) assert(localInterpPoints[i]>=0);*/

  if (verbose>1) 
    printf("Of these %d are distinct points.\n",nDistinctInterpPoints);
  /* Test output 
  for (i=0; i<nDistinctInterpPoints; i++) printf("%d\n",interpPoints[i]);
  */

 /* Test to see if the sorting and indexing are OK:                   
  i=0;
  for (ic=0; ic<nCols; ic++){
    dCol = dCols[ic];
    dStart = dmat->colIndex[dCol];
    dStop  = dmat->colIndex[dCol+1]; 
    for (ir=dStart; ir<dStop; ir++){
      iRow = dmat->rowIndex[ir];
      iStart = imat->colIndex[iRow];
      iStop  = imat->colIndex[iRow+1];
      for (j=iStart; j<iStop; j++, i++){
	assert((imat->rowIndex[j])==(interpPoints[localInterpPoints[i]]));
      }
    }
  }*/
  /* Simple test to catch stupid errors: 
  for (i=0; i<nInterpPoints; i++)
      assert(localInterpPoints[i]<nDistinctInterpPoints && iiLoc>=0);*/

  /* Setup memory for storing grid-to-grid (kernel) operations.
   * The kernel values (grid-to-grid interactions) are stored in 
   * "localKernel".                                                    */
  if (nDistinctInterpPoints>maxDistinctInterpPoints || 
      nDistinctProjPoints>maxDistinctProjPoints){
    /* Possibly free memory */
    if (localKernel!=NULL){
      for (i=0; i<maxDistinctProjPoints; i++) FREE(localKernel[i]);
      FREE(localKernel);
      FREE(pTimesG2G);
    }
    if (verbose) printf("%s: Allocating memory, part III\n",fctName);
    if (nDistinctProjPoints>maxDistinctProjPoints)
      maxDistinctProjPoints = MAX(nDistinctProjPoints,2*maxDistinctProjPoints);
    if (nDistinctInterpPoints>maxDistinctInterpPoints)
      maxDistinctInterpPoints = MAX(nDistinctInterpPoints,2*maxDistinctInterpPoints);

    localKernel = (double **) calloc(maxDistinctProjPoints, sizeof(double *));
    for (i=0; i<maxDistinctProjPoints; i++) 
      localKernel[i] = 
	(double *) calloc(maxDistinctInterpPoints, sizeof(double));

    /* For storing the sum of PROJ*G2G (projection times grid-to-grid) 
     * for each distinct interpolation point                           */
    pTimesG2G = (double *) calloc(maxDistinctInterpPoints, sizeof(double));
  }

  /* Setup kernel values (grid-to-grid interactions). These are found 
   * for all distinct pairs of projection and interpolation points.    */
  supplyKernelValues(grid, 
		     nDistinctProjPoints, projPoints,
		     nDistinctInterpPoints, interpPoints,
		     localKernel);

  /* The precorrection consists of calculating PROJ*G2G*INTERP and 
   * subtracting the sum (over projection and interpolation points) 
   * from each matrix entry. Basically this may be written as
   *
   *        \sum_i \sum_j P_i G_{ij} I_j  = P*G*I
   *
   * or, equivalently as:
   *
   *        \sum_j (I_j ( \sum_i P_i G_{ij} ) ) 
   *
   * The latter aproach is taken here.                                 */

  /* From here-on each column must be treated separately.              */
  ippoint = 0; /* Counters starts out here! */
  iipoint = 0; 
  for (ic=0; ic<nCols; ic++){
    dCol = dCols[ic];
    dStart = dmat->colIndex[dCol];
    dStop  = dmat->colIndex[dCol+1]; 
    nRows  = dStop-dStart;

    pStart = pmat->colIndex[dCol];
    pStop = pmat->colIndex[dCol+1];
    nProjPoints = pStop-pStart;

    /* Zero out the working array                                      */
    for (j=0; j<nDistinctInterpPoints; j++) pTimesG2G[j] = 0.0;
    /* Calculate  SUM(PROJ*G2G)  for all distinct interpolation points */
    /* Loop over projection points                                     */
    for (i=0; i<nProjPoints; i++, ippoint++){
      /* Local value of this projection entry */
      ipLoc = localProjPoints[ippoint]; /* Projection point number 
					 * in the local sorted arrays  */
      /* PROJ value for this projection point                          */
      pval = pmat->matVals[0][i+pStart];
      /* Loop over interpolation points                                */
      for (j=0; j<nDistinctInterpPoints; j++){
	pTimesG2G[j] += localKernel[ipLoc][j] * pval;
      }
    }

    /* For each matrix entry in the direct matrix, loop over the 
     * equivalent column in the  interpolation matrix and subtract 
     * (PROJ*G2G)*INTERP from the direct matrix entry                  */
    for (ir=dStart; ir<dStop; ir++){
      iRow = dmat->rowIndex[ir];
      iStart = imat->colIndex[iRow];
      iStop  = imat->colIndex[iRow+1];
      for (cval=0.0, j=iStart; j<iStop; j++, iipoint++){
	iiLoc = localInterpPoints[iipoint]; /* Interpolation point number 
					   * in the local sorted arrays*/
	assert(iiLoc<nDistinctInterpPoints && iiLoc>=0);
	/* INTERP value */
	ival = imat->matVals[0][j];
	/* Add this value of PROJ*G2G*INTERP the the correction value  */
	cval += ival*pTimesG2G[iiLoc];
      }

      {/* Test vs. existing routine 
#include "Grid.h"
      struct grid *g;
      double cval2, ee;

      g=(struct grid*) grid;
      cval2 = precorrectVal(g->kernel,g->griddata,
			    &(pmat->matVals[0][pStart]), 
			    &(pmat->rowIndex[pStart]), 
			    pStop-pStart,
			    &(imat->matVals[0][iStart]), 
			    &(imat->rowIndex[iStart]), 
			    iStop-iStart);

      ee = ABS(cval-cval2);

      if (ee>1e-6){
	printf("DISCREPANCY ");
	printf("(%d,%d-%d): %16.8e %16.8e\n",
	       ir,iStart,iStop-1,cval,cval2);
      }
      assert(ee<1e-6);*/
      }

      /* Subtract the calculated corretion value from the matrix entry */
      dmat->matVals[0][ir] -= cval;
    }
  }
  /* Check that we have been through all the interpolation points.     */
  assert(iipoint == nInterpPoints);

  /*for (i=0; i<iipoint; i++) printf(" %d %d\n",i,interpPoints[i]);*/

  return;
} /* End of routine precorrectColumns */

/* 
* ==================================================================== 
* Precorrects the nCols columns in the direct matrix.
* Each of the corresponding columns IN THE PROJECTION MATRIX *must* 
* be structurally identical, i.e. must have non-zero entries in exactly 
* the same places. The VALUES of the entries may differ from column to 
* column.
* ==================================================================== */
void precorrectColumnsOld(int nCols, int *dCols,
			  void *directMatIn, 
			  void *projectMatIn, 
			  void *interpMatIn,
			  void *grid){
  char fctName[]="precorrectColumnsOld";  /* Name of this routine         */
  sparseMat *dmat, *pmat, *imat;
  int dStart, dStop, pStart, pStop, iStart,iStop,
    i,j,ir,ic, np, nInterpPoints, nProjectPoints, 
    iRow, ipoint, nInterpPointsPerCol, iLoc,
    nRows, nDistinctInterpPoints,
    dCol, iRowTest;
  double pval, ival, cval;
  static int maxInterpPoints=0, *interpPoints=NULL, *interpPoints2=NULL, 
    *localInterpPoints=NULL,
    maxDistinctInterpPoints=0, maxProjPoints=0;
  static double **localKernel, *pTimesG2G;
  int verbose = 0;

  /* If negative number of columns is supplied then 
   * simply free memory and return                                     */
  if (nCols<0){
    IFFREE(interpPoints);
    IFFREE(interpPoints2);
    IFFREE(localInterpPoints);
    if (localKernel!=NULL){
      for (i=0; i<maxProjPoints; i++) FREE(localKernel[i]);
      FREE(localKernel);
    }
    IFFREE(pTimesG2G);
    return;
  }

  /* Cast input matrices to correct type */
  dmat = (sparseMat*) directMatIn;
  pmat = (sparseMat*) projectMatIn;
  imat = (sparseMat*) interpMatIn;

  if(dmat->nMats > 1) {
    printf("%s: Multiple mats not yet supported.",fctName);
    _EXIT_;  /* Just one mat for now. */
  }

  /* Start out using the initial column (all columns MUST have 
   * identical structure, i.e. have nonzero entries in exactly the 
   * same places. The VALUES of the entries may vary, though.          */
  dCol = dCols[0];
  /* Get the column in the direct matrix associated with source dCol.  */
  dStart = dmat->colIndex[dCol];
  dStop  = dmat->colIndex[dCol+1]; 
  nRows  = dStop-dStart;
  /* Find (few) related grid points from projection matrix             */
  /* Get the projection column associated with this source.            */
  pStart = pmat->colIndex[dCol];
  pStop = pmat->colIndex[dCol+1];
  nProjectPoints = pStop-pStart;

  /* Test that columns are structurally identical                      */
  for (i=0; i<nProjectPoints; i++){
    /* Row number of the test column */
    iRowTest = pmat->rowIndex[pStart+i];
    /* Test this vs. the same row in all other supplied columns        */
    for (ic=1; ic<nCols; ic++){
      assert(pmat->rowIndex[pmat->colIndex[dCols[ic]]+i] == iRowTest);
    }
  }

  /* Count the total number of grid points (including duplicates) 
   * and set up memory.                                                */
  nInterpPoints = 0;     /* Total number of interpolation points here  */
  nInterpPointsPerCol=0; /* # interp. points per column (eval)         */
  for (ic=0; ic<nCols; ic++){
    dCol = dCols[ic];
    /* Get the column in the direct matrix associated with source dCol.*/
    dStart = dmat->colIndex[dCol];
    dStop  = dmat->colIndex[dCol+1]; 

    for (i=dStart; i<dStop; i++){
      iRow  = dmat->rowIndex[i];
      np    = imat->colIndex[iRow+1]-imat->colIndex[iRow];
      nInterpPoints += np;
      nInterpPointsPerCol = MAX(nInterpPointsPerCol,np);
    }
  }
  if (verbose>1){
    printf("%s: %d interpolation points to be treated for columns:",
	   fctName,nInterpPoints);
    for (ic=0; ic<nCols; ic++)  printf(" %d",dCols[ic]);
    printf("\n");
  }
  /* Setup memory if needed */
  if (nInterpPoints > maxInterpPoints){
    IFFREE(interpPoints);
    IFFREE(interpPoints2);
    IFFREE(localInterpPoints);
    if (verbose) printf("%s: Allocating memory, part I\n",fctName);
    /* Memory for sorting interpolation points                         */
    maxInterpPoints = MAX(2*maxInterpPoints,nInterpPoints);
    interpPoints    = (int *) calloc(maxInterpPoints,sizeof(int));
    interpPoints2   = (int *) calloc(maxInterpPoints,sizeof(int));
    localInterpPoints = (int *) calloc(maxInterpPoints,sizeof(int));
  }

  /* Find (many) related grid points from interpolation matrix         */
  /* Store numbers of interpolation points in array for sorting        */
  ipoint = 0;
  for (ic=0; ic<nCols; ic++){
    dCol = dCols[ic];
    /* Get the column in the direct matrix associated with source dCol.*/
    dStart = dmat->colIndex[dCol];
    dStop  = dmat->colIndex[dCol+1]; 

    for (ir=dStart; ir<dStop; ir++){
      iRow = dmat->rowIndex[ir];
      iStart = imat->colIndex[iRow];
      iStop  = imat->colIndex[iRow+1];
      for (j=iStart; j<iStop; j++, ipoint++){
	interpPoints[ipoint]  = imat->rowIndex[j];
	interpPoints2[ipoint] = ipoint;
	/*printf("(%d,%d): %d  %d\n",i,j,ipoint,interpPoints[ipoint]);*/
      }
    }
  }
  assert(ipoint==nInterpPoints); 

  /* Compression of interpolation points:                              */
  /* Sort the grid points (keep track of original relations)           */
  /*intSort(interpPoints, maxInterpPoints);*/
  intSort2(interpPoints, interpPoints2, nInterpPoints);
  /* Eliminate duplicate points while keeping track of relations to 
   * original entries. I know this is a bit messy, but here is a little 
   * help:
   *   "interpPoints" is the sorted array. Duplicate entries will 
   *                  be eliminated here.
   *   "interpPoints2" contains the indices of "interpPoints" relating
   *                  to before the sorting. Thus, interpPoints2 is used 
   *                  to keep track of the sorting process.
   *   "localInterpPoints" is used for storing the local indices (into 
   *                  interpPoints), but ordered in the same way as 
   *                  interpPoints was BEFORE the sorting.
   *   "interpVals"   Contains the values of the 
   * Thus, if the original index (before sorting) is "i", then 
   * localInterpPoints[i] contains the index of the point in the sorted 
   * and eliminated array. (I hope).                                   */
  /* Initialize array to capture errors (this should eventually 
   * not be necessary) 
  for (i=1; i<nInterpPoints; i++) localInterpPoints[i] = -1;*/
  /* Always keep the first point:                                      */
  ipoint = 0;
  localInterpPoints[interpPoints2[0]] = 0;
  for (i=1; i<nInterpPoints; i++){
    if (interpPoints[i]!=interpPoints[i-1]){
      assert(interpPoints[i]>interpPoints[i-1]);
      ipoint++;
      interpPoints[ipoint] = interpPoints[i];
    }
    localInterpPoints[interpPoints2[i]] = ipoint;
  }
  nDistinctInterpPoints = ipoint+1;
  /* Interpolation points are now compressed                           */

 /* Make sure that all points have been visited (this should 
   * eventually not be neccesary) 
  for (i=1; i<nInterpPoints; i++) assert(localInterpPoints[i]>=0);*/

  if (verbose>1) 
    printf("Of these %d are distinct points:\n",nDistinctInterpPoints);
  /* Test output 
  for (i=0; i<nDistinctInterpPoints; i++) printf("%d\n",interpPoints[i]);
  */

  /* Test to see if the sorting and indexing are OK:                   
  i=0;
  for (ic=0; ic<nCols; ic++){
    dCol = dCols[ic];
    dStart = dmat->colIndex[dCol];
    dStop  = dmat->colIndex[dCol+1]; 
    for (ir=dStart; ir<dStop; ir++){
      iRow = dmat->rowIndex[ir];
      iStart = imat->colIndex[iRow];
      iStop  = imat->colIndex[iRow+1];
      for (j=iStart; j<iStop; j++, i++){
	assert((imat->rowIndex[j])==(interpPoints[localInterpPoints[i]]));
      }
    }
  }*/
  /* Simple test to catch stupid errors: 
  for (i=0; i<nInterpPoints; i++)
      assert(localInterpPoints[i]<nDistinctInterpPoints && iLoc>=0);*/

  /* Setup memory for storing grid-to-grid (kernel) operations.
   * The kernel values (grid-to-grid interactions) are stored in 
   * "localKernel".                                                    */
  if (nDistinctInterpPoints>maxDistinctInterpPoints || 
      nProjectPoints>maxProjPoints){
    /* Possibly free memory */
    if (localKernel!=NULL){
      for (i=0; i<maxProjPoints; i++) FREE(localKernel[i]);
      FREE(localKernel);
      FREE(pTimesG2G);
    }
    if (verbose) printf("%s: Allocating memory, part III\n",fctName);
    if (nProjectPoints>maxProjPoints)
      maxProjPoints = MAX(nProjectPoints,2*maxProjPoints);
    if (nDistinctInterpPoints>maxDistinctInterpPoints)
      maxDistinctInterpPoints = MAX(nDistinctInterpPoints,2*maxDistinctInterpPoints);

    localKernel = (double **) calloc(maxProjPoints, sizeof(double *));
    for (i=0; i<maxProjPoints; i++) 
      localKernel[i] = 
	(double *) calloc(maxDistinctInterpPoints, sizeof(double));

    /* For storing the sum of PROJ*G2G (projection times grid-to-grid) 
     * for each distinct interpolation point                           */
    pTimesG2G = (double *) calloc(maxDistinctInterpPoints, sizeof(double));
  }

  /* Setup kernel values (grid-to-grid interactions). These are found 
   * for all distinct pairs of projection and interpolation points. 
   * Note that all the supplied values will be used, since all pairs 
   * of porjction and interpolation points will be used in at least one
   * entry in the direct matrix.                                       */
  supplyKernelValues(grid, 
		     nProjectPoints, &(pmat->rowIndex[pStart]),
		     nDistinctInterpPoints, interpPoints,
		     localKernel);

  /* The precorrection consists of calculating INTERP*G2G*PROJ and 
   * subtracting the sum (over projection and interpolation points) 
   * from each matrix entry. Basically this may be written as
   *
   *        \sum_i \sum_j P_i G_{ij} I_j  = P*G*I
   *
   * or, equivalently as:
   *
   *        \sum_j (I_j ( \sum_i P_i G_{ij} ) ) 
   *
   * The latter aproach is taken here.                                 */

  /* From here-on each column must be treated separately.              */
  ipoint = 0; /* Counter starts out here! */
  for (ic=0; ic<nCols; ic++){
    dCol = dCols[ic];
    dStart = dmat->colIndex[dCol];
    dStop  = dmat->colIndex[dCol+1]; 
    nRows  = dStop-dStart;

    pStart = pmat->colIndex[dCol];
    pStop = pmat->colIndex[dCol+1];
    nProjectPoints = pStop-pStart;

    /* Zero out the working array                                      */
    for (j=0; j<nDistinctInterpPoints; j++) pTimesG2G[j] = 0.0;
    /* Calculate  SUM(PROJ*G2G)  for all distinct interpolation points */
    /* Loop over projection points                                     */
    for (i=0; i<nProjectPoints; i++){
      /* PROJ value for this projection point                          */
      pval = pmat->matVals[0][i+pStart];
      /* Loop over interpolation points                                */
      for (j=0; j<nDistinctInterpPoints; j++){
	pTimesG2G[j] += localKernel[i][j] * pval;
      }
    }

    /* For each matrix entry in the direct matrix, loop over the 
     * equivalent column in the  interpolation matrix and subtract 
     * (PROJ*G2G)*INTERP from the direct matrix entry                  */
    for (ir=dStart; ir<dStop; ir++){
      iRow = dmat->rowIndex[ir];
      iStart = imat->colIndex[iRow];
      iStop  = imat->colIndex[iRow+1];
      for (cval=0.0, j=iStart; j<iStop; j++, ipoint++){
	iLoc = localInterpPoints[ipoint]; /* Interpolation point number 
					   * in the local sorted arrays*/
	assert(iLoc<nDistinctInterpPoints && iLoc>=0);
	/* INTERP value */
	ival = imat->matVals[0][j];
	/* Add this value of PROJ*G2G*INTERP the the correction value  */
	cval += ival*pTimesG2G[iLoc];
      }

      {/* Test vs. existing routine 
#include "Grid.h"
      struct grid *g;
      double cval2, ee;

      g=(struct grid*) grid;
      cval2 = precorrectVal(g->kernel,g->griddata,
			    &(pmat->matVals[0][pStart]), 
			    &(pmat->rowIndex[pStart]), 
			    pStop-pStart,
			    &(imat->matVals[0][iStart]), 
			    &(imat->rowIndex[iStart]), 
			    iStop-iStart);

      ee = ABS(cval-cval2);

      if (ee>1e-6){
	printf("DISCREPANCY ");
	printf("(%d,%d-%d): %16.8e %16.8e\n",
	       ir,iStart,iStop-1,cval,cval2);
      }
      assert(ee<1e-6);*/
      }

      /* Subtract the calculated corretion value from the matrix entry */
      dmat->matVals[0][ir] -= cval;
    }
  }
  /* Check that we have been through all the interpolation points.     */
  assert(ipoint == nInterpPoints);

  /*for (i=0; i<ipoint; i++) printf(" %d %d\n",i,interpPoints[i]);*/

  return;
} /* End of routine precorrectColumnsOld */
/* 
* ==================================================================== 
* Calculates precorrection for nCols columns in the direct matrix with 
* respect to some given rows.
* Each of the corresponding columns IN THE PROJECTION MATRIX *must* 
* be structurally identical, i.e. must have non-zero entries in exactly 
* the same places (this means that the corresponding elements - or 
* basis functions - must use the sam projection scheme centered on the 
* same grid point). The VALUES of the entries may differ from column to 
* column.
* ==================================================================== */
void precorrectPartialColumns(int nCols, int *dCols,
			      int nRows, int *dRows,
			      double **precVals,
			      void *directMatIn, 
			      void *projectMatIn, 
			      void *interpMatIn,
			      void *grid){
  char fctName[]="precorrectPartialColumns";  /* Name of this routine  */
  sparseMat *dmat, *pmat, *imat;
  int dStart, dStop, pStart, pStop, iStart,iStop,
    i,j,ir,ic, np, nInterpPoints, nProjectPoints, 
    iRow, ipoint, nInterpPointsPerCol, iLoc,
    nDistinctInterpPoints,
    dCol, iRowTest;
  double pval, ival, cval;
  static int maxInterpPoints=0, *interpPoints=NULL, *interpPoints2=NULL, 
    *localInterpPoints=NULL,
    maxDistinctInterpPoints=0, maxProjPoints=0;
  static double **localKernel, *pTimesG2G;
  int verbose = 0;

  /* If negative number of columns is supplied then 
   * simply free memory and return                                     */
  if (nCols<0){
    IFFREE(interpPoints);
    IFFREE(interpPoints2);
    IFFREE(localInterpPoints);
    if (localKernel!=NULL){
      for (i=0; i<maxProjPoints; i++) FREE(localKernel[i]);
      FREE(localKernel);
    }
    IFFREE(pTimesG2G);
    return;
  }

  /* Cast input matrices to correct type */
  dmat = (sparseMat*) directMatIn;
  pmat = (sparseMat*) projectMatIn;
  imat = (sparseMat*) interpMatIn;

  if(dmat->nMats > 1) {
    printf("%s: Multiple mats not yet supported.",fctName);
    _EXIT_;  /* Just one mat for now. */
  }

  /* Start out using the initial column (all columns MUST have 
   * identical structure, i.e. have nonzero entries in exactly the 
   * same places. The VALUES of the entries may vary, though.          */
  dCol = dCols[0];
  /* Find (few) related grid points from projection matrix             */
  /* Get the projection column associated with this source.            */
  pStart = pmat->colIndex[dCol];
  pStop = pmat->colIndex[dCol+1];
  nProjectPoints = pStop-pStart;

  /* Test that columns are structurally identical                      */
  for (i=0; i<nProjectPoints; i++){
    /* Row number of the test column */
    iRowTest = pmat->rowIndex[pStart+i];
    /* Test this vs. the same row in all other supplied columns        */
    for (ic=1; ic<nCols; ic++){
      assert(pmat->rowIndex[pmat->colIndex[dCols[ic]]+i] == iRowTest);
    }
  }

  /* Count the total number of interpolation grid points (including 
   * duplicates) and set up memory.                                    */
  nInterpPoints = 0;     /* Total number of interpolation points here  */
  nInterpPointsPerCol=0; /* # interp. points per column (eval)         */
  for (ir=0; ir<nRows; ir++){
    iRow  = dRows[ir];
    np    = imat->colIndex[iRow+1]-imat->colIndex[iRow];
    nInterpPoints += np;
    nInterpPointsPerCol = MAX(nInterpPointsPerCol,np);
  }
  if (verbose>1){
    printf("%s: %d interpolation points to be treated for columns:",
	   fctName,nInterpPoints);
    for (ic=0; ic<nCols; ic++)  printf(" %d",dCols[ic]);
    printf("\n");
  }
  /* Setup memory if needed */
  if (nInterpPoints > maxInterpPoints){
    IFFREE(interpPoints);
    IFFREE(interpPoints2);
    IFFREE(localInterpPoints);
    if (verbose) printf("%s: Allocating memory, part I\n",fctName);
    /* Memory for sorting interpolation points                         */
    maxInterpPoints = MAX(2*maxInterpPoints,nInterpPoints);
    interpPoints    = (int *) calloc(maxInterpPoints,sizeof(int));
    interpPoints2   = (int *) calloc(maxInterpPoints,sizeof(int));
    localInterpPoints = (int *) calloc(maxInterpPoints,sizeof(int));
  }

  /* Find (many) related grid points from interpolation matrix         */
  /* Store numbers of interpolation points in array for sorting        */
  ipoint = 0;
  for (ir=0; ir<nRows; ir++){
    iRow  = dRows[ir];
    iStart = imat->colIndex[iRow];
    iStop  = imat->colIndex[iRow+1];
    for (j=iStart; j<iStop; j++, ipoint++){
      interpPoints[ipoint]  = imat->rowIndex[j];
      interpPoints2[ipoint] = ipoint;
      /*printf("(%d,%d): %d  %d\n",iRow,j,ipoint,interpPoints[ipoint]);*/
    }
  }

  assert(ipoint==nInterpPoints); 

  /* Compression of interpolation points:                              */
  /* Sort the grid points (keep track of original relations)           */
  /*intSort(interpPoints, maxInterpPoints);*/
  intSort2(interpPoints, interpPoints2, nInterpPoints);
  /* Eliminate duplicate points while keeping track of relations to 
   * original entries. I know this is a bit messy, but here is a little 
   * help:
   *   "interpPoints" is the sorted array. Duplicate entries will 
   *                  be eliminated here.
   *   "interpPoints2" contains the indices of "interpPoints" relating
   *                  to before the sorting. Thus, interpPoints2 is used 
   *                  to keep track of the sorting process.
   *   "localInterpPoints" is used for storing the local indices (into 
   *                  interpPoints), but ordered in the same way as 
   *                  interpPoints was BEFORE the sorting.
   *   "interpVals"   Contains the values of the 
   * Thus, if the original index (before sorting) is "i", then 
   * localInterpPoints[i] contains the index of the point in the sorted 
   * and eliminated array. (I hope).                                   */
  /* Initialize array to capture errors (this should eventually 
   * not be necessary) 
  for (i=1; i<nInterpPoints; i++) localInterpPoints[i] = -1;*/
  /* Always keep the first point:                                      */
  ipoint = 0;
  localInterpPoints[interpPoints2[0]] = 0;
  for (i=1; i<nInterpPoints; i++){
    if (interpPoints[i]!=interpPoints[i-1]){
      assert(interpPoints[i]>interpPoints[i-1]);
      ipoint++;
      interpPoints[ipoint] = interpPoints[i];
    }
    localInterpPoints[interpPoints2[i]] = ipoint;
  }
  nDistinctInterpPoints = ipoint+1;
  /* Interpolation points are now compressed                           */

 /* Make sure that all points have been visited (this should 
   * eventually not be neccesary) 
  for (i=1; i<nInterpPoints; i++) assert(localInterpPoints[i]>=0);*/

  if (verbose>1) 
    printf("Of these %d are distinct points.\n",nDistinctInterpPoints);
  /* Test output 
  for (i=0; i<nDistinctInterpPoints; i++) printf("%d\n",interpPoints[i]);
  */

  /* Test to see if the sorting and indexing are OK:                   
  i=0;
  for (ir=0; ir<nRows; ir++){
    iRow  = dRows[ir];
    iStart = imat->colIndex[iRow];
    iStop  = imat->colIndex[iRow+1];
    for (j=iStart; j<iStop; j++, i++){
      assert((imat->rowIndex[j])==(interpPoints[localInterpPoints[i]]));
    }
  }*/


  /* Simple test to catch stupid errors: 
  for (i=0; i<nInterpPoints; i++)
      assert(localInterpPoints[i]<nDistinctInterpPoints && iLoc>=0);*/

  /* Setup memory for storing grid-to-grid (kernel) operations.
   * The kernel values (grid-to-grid interactions) are stored in 
   * "localKernel".                                                    */
  if (nDistinctInterpPoints>maxDistinctInterpPoints || 
      nProjectPoints>maxProjPoints){
    /* Possibly free memory */
    if (localKernel!=NULL){
      for (i=0; i<maxProjPoints; i++) FREE(localKernel[i]);
      FREE(localKernel);
      FREE(pTimesG2G);
    }
    if (verbose) printf("%s: Allocating memory, part III\n",fctName);
    if (nProjectPoints>maxProjPoints)
      maxProjPoints = MAX(nProjectPoints,2*maxProjPoints);
    if (nDistinctInterpPoints>maxDistinctInterpPoints)
      maxDistinctInterpPoints = MAX(nDistinctInterpPoints,2*maxDistinctInterpPoints);

    localKernel = (double **) calloc(maxProjPoints, sizeof(double *));
    for (i=0; i<maxProjPoints; i++) 
      localKernel[i] = 
	(double *) calloc(maxDistinctInterpPoints, sizeof(double));

    /* For storing the sum of PROJ*G2G (projection times grid-to-grid) 
     * for each distinct interpolation point                           */
    pTimesG2G = (double *) calloc(maxDistinctInterpPoints, sizeof(double));
  }

  /* Setup kernel values (grid-to-grid interactions). These are found 
   * for all distinct pairs of projection and interpolation points. 
   * Note that all the supplied values will be used, since all pairs 
   * of porjction and interpolation points will be used in at least one
   * entry in the direct matrix.                                       */
  supplyKernelValues(grid, 
		     nProjectPoints, &(pmat->rowIndex[pStart]),
		     nDistinctInterpPoints, interpPoints,
		     localKernel);

 /* The precorrection consists of calculating PROJ*G2G*INTERP and 
   * subtracting the sum (over projection and interpolation points) 
   * from each matrix entry. Basically this may be written as
   *
   *        \sum_i \sum_j P_i G_{ij} I_j  = P*G*I
   *
   * or, equivalently as:
   *
   *        \sum_j (I_j ( \sum_i P_i G_{ij} ) ) 
   *
   * The latter aproach is taken here.                                 */

  /* From here-on each column must be treated separately.              */
  ipoint = 0; /* Counter starts out here! */
  /* For each source scan through all the given evals */
  for (ic=0; ic<nCols; ic++){
    dCol = dCols[ic];

    pStart = pmat->colIndex[dCol];
    pStop = pmat->colIndex[dCol+1];
    nProjectPoints = pStop-pStart;

    /* Zero out the working array                                      */
    for (j=0; j<nDistinctInterpPoints; j++) pTimesG2G[j] = 0.0;
    /* Calculate  SUM(PROJ*G2G)  for all distinct interpolation points */
    /* Loop over projection points                                     */
    for (i=0; i<nProjectPoints; i++){
      /* PROJ value for this projection point                          */
      pval = pmat->matVals[0][i+pStart];
      /* Loop over interpolation points                                */
      for (j=0; j<nDistinctInterpPoints; j++){
	pTimesG2G[j] += localKernel[i][j] * pval;
      }
    }

    /* For each matrix entry in the direct matrix, loop over the 
     * equivalent column in the  interpolation matrix and subtract 
     * (PROJ*G2G)*INTERP from the direct matrix entry                  */

    for (ir=0; ir<nRows; ir++){
      iRow  = dRows[ir];
      iStart = imat->colIndex[iRow];
      iStop  = imat->colIndex[iRow+1];
      for (cval=0.0, j=iStart; j<iStop; j++, ipoint++){
	iLoc = localInterpPoints[ipoint]; /* Interpolation point number 
					   * in the local sorted arrays*/
	assert(iLoc<nDistinctInterpPoints && iLoc>=0);
	/* INTERP value 
	ival = imat->matVals[0][j];*/
	/* Add this value of PROJ*G2G*INTERP the the correction value  */
	/*cval += ival*pTimesG2G[iLoc];*/
	cval += (imat->matVals[0][j])*pTimesG2G[iLoc];
      }

      if (0) {/* Test vs. existing routine */
#include "Grid.h"
	struct grid *g;
	double cval2, ee;
	
	g=(struct grid*) grid;
	cval2 = precorrectVal(g->kernel,g->griddata,
			      &(pmat->matVals[0][pStart]), 
			      &(pmat->rowIndex[pStart]), 
			      pStop-pStart,
			      &(imat->matVals[0][iStart]), 
			      &(imat->rowIndex[iStart]), 
			      iStop-iStart);
	
	ee = ABS(cval-cval2);

	if (ee>1e-6){
	  printf("DISCREPANCY ");
	  printf("(%d,%d-%d): %16.8e %16.8e\n",
		 ir,iStart,iStop-1,cval,cval2);
	}
	assert(ee<1e-6);
      }

      /* Store this value in the return array */
      precVals[ic][ir] = cval;
    }
  }
  /* Check that we have been through all the interpolation points.     */
  assert(ipoint == nInterpPoints);

  /*for (i=0; i<ipoint; i++) printf(" %d %d\n",i,interpPoints[i]);*/

  return;
} /* End of routine precorrectPartialColumns */

/* 
* ==================================================================== 
* Precorrects the dCol'th column in the direct matrix.
* ==================================================================== */
void precorrectColumn(int dCol,
		      void *directMatIn, 
		      void *projectMatIn, 
		      void *interpMatIn,
		      void *grid){
  char fctName[]="precorrectColumn";  /* Name of this routine          */
  sparseMat *dmat, *pmat, *imat;
  int dStart, dStop, pStart, pStop, iStart,iStop,
    i,j,ir, np, nInterpPoints, nProjectPoints, 
    iRow, ipoint, nInterpPointsPerCol, iLoc,
    nRows, nDistinctInterpPoints;
  double pval, ival, cval;
  static int maxInterpPoints=0, *interpPoints=NULL, *interpPoints2=NULL, 
    *localInterpPoints=NULL,
    maxDistinctInterpPoints=0, maxProjPoints=0;
  static double **localKernel, *pTimesG2G;
  int verbose = 0;

  /* If negative column number is supplied then simply free memory 
   * and return */
  if (dCol<0){
    IFFREE(interpPoints);
    IFFREE(interpPoints2);
    IFFREE(localInterpPoints);
    if (localKernel!=NULL){
      for (i=0; i<maxProjPoints; i++) FREE(localKernel[i]);
      FREE(localKernel);
    }
    IFFREE(pTimesG2G);
    return;
  }

  /* Cast input matrices to correct type */
  dmat = (sparseMat*) directMatIn;
  pmat = (sparseMat*) projectMatIn;
  imat = (sparseMat*) interpMatIn;

  if(dmat->nMats > 1) {
    printf("%s: Multiple mats not yet supported.",fctName);
    _EXIT_;  /* Just one mat for now. */
  }

  /* Get the column in the direct matrix associated with source dCol.  */
  dStart = dmat->colIndex[dCol];
  dStop  = dmat->colIndex[dCol+1]; 
  nRows  = dStop-dStart;
  /* Find (few) related grid points from projection matrix             */
  /* Get the projection column associated with this source.            */
  pStart = pmat->colIndex[dCol];
  pStop = pmat->colIndex[dCol+1];
  nProjectPoints = pStop-pStart;

  /* Count the total number of grid points (including duplicates) 
   * and set up memory.                                                */
  nInterpPoints = 0;     /* Total number of interpolation points here  */
  nInterpPointsPerCol=0; /* # interp. points per column (eval)         */
  for (i=dStart; i<dStop; i++){
    iRow  = dmat->rowIndex[i];
    np    = imat->colIndex[iRow+1]-imat->colIndex[iRow];
    nInterpPoints += np;
    nInterpPointsPerCol = MAX(nInterpPointsPerCol,np);
  }
  if (verbose>1)
    printf("%s: %d interpolation points to be treated for column %d\n",
	   fctName,nInterpPoints,dCol);
  /* Setup memory if needed */
  if (nInterpPoints > maxInterpPoints){
    IFFREE(interpPoints);
    IFFREE(interpPoints2);
    IFFREE(localInterpPoints);
    if (verbose) printf("%s: Allocating memory, part I\n",fctName);
    /* Memory for sorting interpolation points                         */
    maxInterpPoints = MAX(2*maxInterpPoints,nInterpPoints);
    interpPoints    = (int *) calloc(maxInterpPoints,sizeof(int));
    interpPoints2   = (int *) calloc(maxInterpPoints,sizeof(int));
    localInterpPoints = (int *) calloc(maxInterpPoints,sizeof(int));
  }

  /* Find (many) related grid points from interpolation matrix         */
  /* Store numbers of interpolation points in array for sorting        */
  ipoint = 0;
  for (ir=dStart; ir<dStop; ir++){
    iRow = dmat->rowIndex[ir];
    iStart = imat->colIndex[iRow];
    iStop  = imat->colIndex[iRow+1];
    for (j=iStart; j<iStop; j++, ipoint++){
      interpPoints[ipoint]  = imat->rowIndex[j];
      interpPoints2[ipoint] = ipoint;
      /*printf("(%d,%d): %d  %d\n",i,j,ipoint,interpPoints[ipoint]);*/
    }
  }
  assert(ipoint==nInterpPoints); /*maxInterpPoints*/
  /* Sort the grid points (keep track of original relations)           */
  /*intSort(interpPoints, maxInterpPoints);*/
  intSort2(interpPoints, interpPoints2, nInterpPoints);

  /* Eliminate duplicate points while keeping track of relations to 
   * original entries. I know this is a bit messy, but here is a little 
   * help:
   *   "interpPoints" is the sorted array. Duplicate entries will 
   *                  be eliminated here.
   *   "interpPoints2" contains the indices of "interpPoints" relating
   *                  to before the sorting. Thus, interpPoints2 is used 
   *                  to keep track of the sorting process.
   *   "localInterpPoints" is used for storing the local indices (into 
   *                  interpPoints), but ordered in the same way as 
   *                  interpPoints was BEFORE the sorting.
   *   "interpVals"   Contains the values of the 
   * Thus, if the original index (before sorting) is "i", then 
   * localInterpPoints[i] contains the index of the point in the sorted 
   * and eliminated array. (I hope).                                   */
  /* Initialize array to capture errors (this should eventually 
   * not be necessary) 
  for (i=1; i<nInterpPoints; i++) localInterpPoints[i] = -1;*/
  /* Always keep the first point:                                      */
  ipoint = 0;
  localInterpPoints[interpPoints2[0]] = 0;
  for (i=1; i<nInterpPoints; i++){
    if (interpPoints[i]!=interpPoints[i-1]){
      assert(interpPoints[i]>interpPoints[i-1]);
      ipoint++;
      interpPoints[ipoint] = interpPoints[i];
    }
    localInterpPoints[interpPoints2[i]] = ipoint;
  }
  nDistinctInterpPoints = ipoint+1;

  /* Make sure that all points have been visited (this should 
   * eventually not be neccesary) 
  for (i=1; i<nInterpPoints; i++) assert(localInterpPoints[i]>=0);*/

  if (verbose>1) 
    printf("Of these %d are distinct points:\n",nDistinctInterpPoints);
  /* Test output 
  for (i=0; i<nDistinctInterpPoints; i++) printf("%d\n",interpPoints[i]);
  */

  /* Test to see if the sorting and indexing are OK:                   */
  i=0;
  for (ir=dStart; ir<dStop; ir++){
    iRow = dmat->rowIndex[ir];
    iStart = imat->colIndex[iRow];
    iStop  = imat->colIndex[iRow+1];
    for (j=iStart; j<iStop; j++, i++){
      assert((imat->rowIndex[j])==(interpPoints[localInterpPoints[i]]));
    }
  }
  /* Simple test to catch stupid errors: */
  for (i=0; i<nInterpPoints; i++)
      assert(localInterpPoints[i]<nDistinctInterpPoints && 
	     localInterpPoints[i]>=0);

  /* Setup memory for storing grid-to-grid (kernel) operations.
   * The kernel values (grid-to-grid interactions) are stored in 
   * "localKernel".                                                    */
  if (nDistinctInterpPoints>maxDistinctInterpPoints || 
      nProjectPoints>maxProjPoints){
    /* Possibly free memory */
    if (localKernel!=NULL){
      for (i=0; i<maxProjPoints; i++) FREE(localKernel[i]);
      FREE(localKernel);
      FREE(pTimesG2G);
    }
    if (verbose) printf("%s: Allocating memory, part III\n",fctName);
    if (nProjectPoints>maxProjPoints)
      maxProjPoints = MAX(nProjectPoints,2*maxProjPoints);
    if (nDistinctInterpPoints>maxDistinctInterpPoints)
      maxDistinctInterpPoints = MAX(nDistinctInterpPoints,2*maxDistinctInterpPoints);

    localKernel = (double **) calloc(maxProjPoints, sizeof(double *));
    for (i=0; i<maxProjPoints; i++) 
      localKernel[i] = 
	(double *) calloc(maxDistinctInterpPoints, sizeof(double));

    /* For storing the sum of PROJ*G2G (projection times grid-to-grid) 
     * for each distinct interpolation point                           */
    pTimesG2G = (double *) calloc(maxDistinctInterpPoints, sizeof(double));
  }
  /* Setup kernel values (grid-to-grid interactions). These are found 
   * for all distinct pairs of projection and interpolation points.    */

  supplyKernelValues(grid, 
		     nProjectPoints, &(pmat->rowIndex[pStart]),
		     nDistinctInterpPoints, interpPoints,
		     localKernel);

  /* The precorrection consists of calculating PROJ*G2G*INTERP and 
   * subtracting the sum (over projection and interpolation points) 
   * from each matrix entry. Basically this may be written as
   *
   *        \sum_i \sum_j P_i G_{ij} I_j  = P*G*I
   *
   * or, equivalently as:
   *
   *        \sum_j (I_j ( \sum_i P_i G_{ij} ) ) 
   *
   * The latter aproach is taken here.                                 */
  /* Zero out the working array                                        */
  for (j=0; j<nDistinctInterpPoints; j++) pTimesG2G[j] = 0.0;
  /* Calculate  SUM(PROJ*G2G)  for all distinct interpolation points   */
  /* Loop over projection points                                       */
  for (i=0; i<nProjectPoints; i++){
    /* PROJ value for this projection point                            */
    pval = pmat->matVals[0][i+pStart];
    /* Loop over interpolation points                                  */
    for (j=0; j<nDistinctInterpPoints; j++){
      pTimesG2G[j] += localKernel[i][j] * pval;
    }
  }

  /* For each matrix entry in the direct matrix, loop over the 
   * equivalent column in the  interpolation matrix and subtract 
   * (PROJ*G2G)*INTERP from the direct matrix entry                    */
  ipoint = 0;
  for (ir=dStart; ir<dStop; ir++){
    iRow = dmat->rowIndex[ir];
    iStart = imat->colIndex[iRow];
    iStop  = imat->colIndex[iRow+1];
    for (cval=0.0, j=iStart; j<iStop; j++, ipoint++){
      iLoc = localInterpPoints[ipoint]; /* Interpolation point number 
					 * in the local sorted arrays  */
      assert(iLoc<nDistinctInterpPoints && iLoc>=0);
      /* INTERP value */
      ival = imat->matVals[0][j];
      /* Add this value of PROJ*G2G*INTERP the the correction value    */
      cval += ival*pTimesG2G[iLoc];
    }

    {/* Test vs. existing routine 
#include "Grid.h"
      struct grid *g;
      double cval2, ee;

      g=(struct grid*) grid;
      cval2 = precorrectVal(g->kernel,g->griddata,
			    &(pmat->matVals[0][pStart]), 
			    &(pmat->rowIndex[pStart]), 
			    pStop-pStart,
			    &(imat->matVals[0][iStart]), 
			    &(imat->rowIndex[iStart]), 
			    iStop-iStart);

      ee = ABS(cval-cval2);

      if (ee>1e-6){
	printf("DISCREPANCY ");
	printf("(%d,%d-%d): %16.8e %16.8e\n",
	       ir,iStart,iStop-1,cval,cval2);
      }
      assert(ee<1e-6);*/
    }

    /* Subtract the calculated corretion value from the matrix entry   */
    dmat->matVals[0][ir] -= cval;
  }
  assert(ipoint == nInterpPoints);

  /*for (i=0; i<ipoint; i++) printf(" %d %d\n",i,interpPoints[i]);*/

  return;
} /* End of routine precorrectColumn */


/*
* ==================================================================== 
* Precorrection direct matrix column-by-column.
* ==================================================================== */
void colPrec(void *directMatIn, 
	     void *projectMatIn, 
	     void *interpMatIn,
	     void *grid){
  sparseMat *dmat;
  int i;
  /* Cast input matrix to correct type */
  dmat = (sparseMat*) directMatIn;
  /* Precorrect one column at a time: */
  for (i=0; i<dmat->n; i++)
    precorrectColumn(i,directMatIn, projectMatIn, interpMatIn, grid);

  /* Free memory */
  precorrectColumn((int) -1,directMatIn, projectMatIn, interpMatIn, grid);
  return;
}


#include "sparse.testfunctions.c"
