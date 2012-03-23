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
/* FILE: elements.c
*
*  This file contains routines that operate on elements, 
*  i.e "sources" and "evals".
*
*/

/* Max length of a line in an input file: */
# define maxLineLength 1000

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "Global.h"

#include "Elements.h"

#include "Calcp.h"

/* ======= Prototypes for static routines ======= */
static 
int identifyElementLine1(char *inputLine, char **restOfLine);
static 
int identifyBoundaryLine1(char *inputLine, char **restOfLine);
static 
void readCorners1(char *inputLine,int nCorners,int iLine,point *corners);

static 
void readCorners(char *inputLine,int nCorners,int iLine,point *corners);

static
int quadNeedsCracking(point *corners, double crackTol);

static
void crackElement(int iElement,
		  int iElement2,
		  int iDiag,
		  struct constantElement **elements);

static
void setupElementCentroids(struct constantElement **elements,
			   int nElements);
static
void setupElementCentroid(struct constantElement *element);

static
void setupElementBoundaryConditionTypes(struct constantElement **elements,
					int nElements,
					char *conditionType,
					int idata);

static
double elementNormal(struct constantElement *element, double *normal);

static
int assignSourceDirectIndices(struct constantElement **sources, 
			      int nSources);
static
int assignEvalDirectIndices(struct constantElement **evals, 
			      int nEvals);

static void localElementMoments(int m, int n, double **moments);

static
void addTriangularElementMoments(point *corners, 
				 double ***moments,
				 double polyScale,
				 int polyOrder, 
				 int maxOrder1D,
				 double *direction);
static
void addMonomialValues(point point, 
		       double ***moments, 
		       double polyScale,
		       int polyOrder, int maxOrder1D,
		       double *direction);

/* The following should eventually *not* be static */
static double ***ComputeMomentsGlobal(double *xp, double *yp, double *zp, 
				      int order, double ***moments);

/* 
* ==================================================================== 
* Setup of source elements. Read from data file.
* ==================================================================== */
/* The following macros are used to ease the */
void *setupSources(int *nSourcesReturn, 
		   int *maxBoundaryPart, 
		   char *sourceInputFileName,
		   int *inputFileTypes,
		   char *bcType,
		   int bcIdata)
{
  char fctName[]="setupSources";  /* Name of this routine              */

  /* The following two could be specified from the outside(?)          */
  int performCracking=1; /* Crack skew quadrilateral elements? 
			  * 1=yes, 0=no                                */
  double crackTol=1.E-6; /* Crack criterion. Crack if the normalized 
			  * normal vectors after cracking deviate more 
			  * than this from each other.                 */
  /* Structure for the sources. List of pointers each pointing 
   * to a source as e.g. sources[4]->substructure                      */
  struct constantElement **sources;

  int i, nCorners, iElement, iDiag, readNextLine, inputType;
  int nLines, nHeaderLines, nElementLines, nCommentLines;
  int nSources, /* Number of sources as read from file                 */
      nCracked; /* Number of elements needed to be cracked.            */
  FILE *sourceInputFile;
  point corners[4];      /* Corner points of element                   */
  char *flagString, *inputLine, *restOfLine;
  fpos_t endOfTitle;
  int nParts;

  /* Write what we do:                                                 */
  fprintf(stdout,"%s: Reading the input file named %s\n",
	  fctName,sourceInputFileName);
  inputType = inputFileTypes[0];
  fprintf(stdout,"%s: Using input file format %d\n",
	  fctName,inputType);

  /* Allocate room for input line */
  inputLine = (char*) calloc(maxLineLength,sizeof(char));

  sourceInputFile = fopen(sourceInputFileName,"r");
  /* Capture any errors */
  if (sourceInputFile==NULL){
    fprintf(stderr,"%s: ERROR: Couldn't open file %s\nExiting!\n",
	    fctName,sourceInputFileName);
    _EXIT_;
  }

  /* First scan the file to figure out the total number of elements.
   * (This could easily be avoided by inluding the number of elements 
   * in the begining of the input file.)                               */

  /* Set initial position of the file (in case there are no header)    */
  if( fgetpos(sourceInputFile,&endOfTitle) ){
    fprintf(stdout,
	    "%s: ERROR! Could not get position in file '%s'\n",
	    fctName,sourceInputFileName);
    _EXIT_;
  }
  /* Read first line: The title line.                                  */
  nHeaderLines=0;
  nLines=0;
  readNextLine=1;
  while (readNextLine){
    nLines++;
    flagString = fgets(inputLine,maxLineLength,sourceInputFile);
    /* Capture any errors */
    if (flagString==NULL){
      fprintf(stderr,"%s: ERROR: Couldn't read line %d of %s\nExiting!\n",
	      fctName,nLines,sourceInputFileName);
      _EXIT_;
    }
    /* Verify that the header line starts with '0'. If not, then assume 
     * that it is *not* a header line after all.                         */
    if (inputLine[0]=='0'){
      nHeaderLines++;
      /* If this is the first header line, then write a comment on 
       * what this output is all about:                                  */
      if(nHeaderLines==1)
	fprintf(stdout,
		"%s: Title of source element input file is:\n\n",fctName);
      /* Write the title line to standard out. 
       * It might be advantageous later to store the title, such that it 
       * can be written with the final solution data as identification. 
       * Do not write the initial '0' of the title line.                 */
      fprintf(stdout,"%s",&(inputLine[1]));
      /* Set this position to be the end of the title (it may be changed 
       * again later if more title lines are found). 
       * Make sure no error is encountered!                              */
      if( fgetpos(sourceInputFile,&endOfTitle) ){
	fprintf(stdout,
		"%s: ERROR! Could not get position in file '%s'\n",
		fctName,sourceInputFileName);
	_EXIT_;
      }
    }
    else { /* This is probably not a header line.                        */
      readNextLine=0;
    }
  } /* No more header lines */
  /* If no header lines were found, then write a statement 
   * to that effect                                                      */
  if (nHeaderLines==0)
    fprintf(stdout,
	    "%s: No title lines found in file '%s'.\n",
	    fctName,sourceInputFileName);
  else fprintf(stdout,"\n");

  /* Go to the end of the header lines. Capture errors.                */
  if ( fsetpos(sourceInputFile,&endOfTitle) ){
    fprintf(stdout,
	    "%s: ERROR! Could not go to position in file '%s'\n",
	    fctName,sourceInputFileName);
    _EXIT_;
  }
  /* Scan the remaining part of the input file to count the number of 
   * source elements to store.                                        */
  fprintf(stdout,
	  "%s: Scanning file %s (pass 1)\n"
	  "              [counting the number of source elements]\n",
	  fctName,sourceInputFileName);

  readNextLine=1;
  nSources=0;
  nCracked=0;
  nCorners=0;
  nElementLines=0;
  nCommentLines=0;
  while(readNextLine){
    /* Read a line                                                     */
    flagString = fgets(inputLine,maxLineLength,sourceInputFile);
    if (flagString==NULL){
      readNextLine = 0;
      continue;
    }
    nLines++; /* A line has been read */
    /* Indentify the data on the line */
    nCorners = identifyElementLine1(inputLine,&restOfLine);

    switch(nCorners){
    case 0:
      /* Line is invalid */
      readNextLine=0;
      break;

    case 3: case 4:
      /* This is a source element (and a source element line)          */
      nSources++;
      nElementLines++;
      break;

    case -1:
      /* Comment line                                                  */
      nCommentLines++;
      break;

    default:
      /* This is the bad stuff (error capture of return data 
       * from identifyElementLine1)                                    */
      fprintf(stdout,
	      "%s: ERROR. %d corners on line %d???\n",
	      fctName,nCorners,nLines);
      _EXIT_;
    } /* End switch */

    /* Get boundary number for this element. Not used in first pass, but 
     * needs to be stripped out of the linea before the read continues */
    if( (inputType==1) && (nCorners==3 || nCorners==4)){
      identifyBoundaryLine1(restOfLine,&restOfLine);
    }
    /* Test if a quadrilateral panel needs cracking. If a quadrilateral 
     * panel is non-planer (skew) then it will be "cracked" (as a tile) 
     * into two plane triangular panels.                               */
    if (performCracking && nCorners==4){
      /* Read the input line - store in variable "corners".
       * Strip the first character off the line to avoid dealing with 
       * the format descriptor a second time                           */
      readCorners1(restOfLine,nCorners,nLines,corners);
      /*      if( VNORM(vec[0]) > crackTol){*/
      /* Test whether cracking is needed.                              */
      if ( quadNeedsCracking(corners,crackTol) ){
	/* This element will need cracking                             */
	nCracked++;
	nSources++;
      }
    }   /* End cracking test                                           */
  } /* Read next line */
  /* Write the findings to standard out.                               */
  fprintf(stdout,"%s: Read %d source elements lines and %d comment lines\n",
	  fctName,nElementLines,nCommentLines);
  if (nCracked==0){
    fprintf(stdout,
	    "%s: All quadrilateral elements are planar \n"
	    "           (specified tolerance: %g)\n",
	    fctName,crackTol);
  }
  else {
    fprintf(stdout,
	    "%s: %d quadrilateral elements will need cracking.\n"
	    "           (specified tolerance: %g)\n",
	    fctName,nCracked,crackTol);
  }
  fprintf(stdout,"%s: %d total source elements\n",fctName,nSources);

  /* Allocate memory for the source elements                           */
  sources = 
    (struct constantElement**) 
    calloc(nSources,sizeof(struct constantElement*));

  for (i=0;i<nSources;i++){
    sources[i] = 
      (struct constantElement*) 
      calloc(1,sizeof(struct constantElement));
  }

  /* Make a second pass through the input file, this time reading the 
   * corner data of each element.                                      */
  fprintf(stdout,
	  "%s: Scanning file %s (pass 2)\n"
	  "              [reading source element data - %d elements]\n",
	  fctName,sourceInputFileName,nSources);

  /* Go to the end of the header lines. Capture errors.                */
  if ( fsetpos(sourceInputFile,&endOfTitle) ){
    fprintf(stdout,
	    "%s: ERROR! Could not go to position in file '%s'\n",
	    fctName,sourceInputFileName);
    _EXIT_;
  }
  nLines = nHeaderLines;

  /* Second pass. Store elements in the allocated storage              */
  if (nCracked == 0){
    /* If no elements need cracking, then we dont need to check them 
     * again. Thus - make a fast pass if no elements need              */
     /* Read the corner data for each element                          */
    for (iElement=0;iElement<nSources;){
      /* Read a line                                                   */
      fgets(inputLine,maxLineLength,sourceInputFile);
      nLines++;
      /* Get the number of corners                                     */
      nCorners = identifyElementLine1(inputLine,&restOfLine);
      if (nCorners==3 || nCorners==4){ /* If this is an element        */
	if (inputType==1){
	  /* Get boundary number for this element.                     */
	  sources[iElement]->boundaryNumber = 
	    (short int) identifyBoundaryLine1(restOfLine,&restOfLine);
	}
	/* Store the numbers of corners                                */
	sources[iElement]->shape = nCorners;
	/* Read the remainder of the line and store in element         */
	readCorners1(restOfLine,nCorners,nLines,
		     sources[iElement]->corners);
	sources[iElement]->crackedNeighbourIndex = iElement;
	/* Update the number of elements                               */
	iElement++;
      } /* End if element */
    }  /* Read next line                                               */
  }    /* End if (no crakced elements)                                 */
  else {
    /* In this case some of the element need cracking. They all need to
     * be examined again to figure out which ones it is (since we didn't)
     * store any of this information.                                  */
    for (iElement=0;iElement<nSources;){
      /* Read a line                                                   */
      fgets(inputLine,maxLineLength,sourceInputFile);
      nLines++;
      /* Get the number of corners                                     */
      nCorners = identifyElementLine1(inputLine,&restOfLine);
      if (nCorners==3 || nCorners==4){ /* If this is an element        */
	/* Store the numbers of corners                                */
	sources[iElement]->shape = nCorners;
	if (inputType==1) {
	  /* Get and store the boundary number                         */
	  sources[iElement]->boundaryNumber = 
	    (short int) identifyBoundaryLine1(restOfLine,&restOfLine);
	}
	/* Read the element as if no cracking is necessary.            */
	readCorners1(restOfLine,nCorners,nLines,
		     sources[iElement]->corners);
	sources[iElement]->crackedNeighbourIndex = iElement;

	/* For quadrilateral panels, check if they need cracking:      */
	if (nCorners==4 &&
	    (iDiag=quadNeedsCracking(sources[iElement]->corners,crackTol)) 
	    ){
	  /* Crack the element                                         */
	  crackElement(iElement,iElement+1,iDiag,sources);
	  iElement++;
	}
	/* Update the number of elements                               */
	iElement++;
      }
    }  /* Read next line                                               */
  }    /* End else (if cracked)                                        */

  /* Close input file - capture errors */
  if ( fclose(sourceInputFile) )
    fprintf(stderr,
	    "%s: WARNING: Could not close file %s\n"
	    "             Continuing.\n",
	    fctName,sourceInputFileName);

  /* For each element set up the stuff that is still needed            */
  /* Calculate element centroids                                       */
  setupElementCentroids(sources,nSources);

  /* Find element boundary types and values                            */
  setupElementBoundaryConditionTypes(sources,nSources,bcType,bcIdata);

  /* Find maximum value for the number of different boundary parts:    */
  nParts=-1;
  for (iElement=0; iElement<nSources; iElement++) 
    nParts = MAX(nParts,(int) sources[iElement]->boundaryNumber);
  assert(nParts>=0);

  /* Set appropriate return values pointers in the list 
   * of dummy arguments                                                */
  *nSourcesReturn  = nSources;
  *maxBoundaryPart = nParts;

  /* Return the pointer to the sources, recast as void */
  return (void *) sources;
} /* End of setupSources */

/* 
* ==================================================================== 
* Setup of evaluation elements. Read from data file.
*
* This routine is basically a copy of setupSources - only with "eval"
* substituted for "source" everywhere (roughly so anyway).
* ==================================================================== */
void *setupEvals(int *nEvalsReturn,
		 char *evalInputFileName, 
		 int *inputFileTypes)
{
  char fctName[]="setupEvals";  /* Name of this routine                */
  struct constantElement **evals;


  evals = (struct constantElement **) NULL;
  fprintf(stderr,"%s: ERROR. This routine is not set up yet.\n"
	  "       Use only one input file with this implementation\n"
	  "       (constant type elements).\n"
	  "       The implementation can be fixed without too great\n"
	  "       problems. (Copy from routine \"setupSources\".)\n"
	  "       This error is terminal.\n",fctName);
  _EXIT_;
  /* The following is simply to avoid warnings: */
  *nEvalsReturn = 0;
  fprintf(stderr,"THIS SHOULD NOT BE PRINTED! %s %d\n",
	  evalInputFileName,inputFileTypes[0]);
  return (void *) evals;
} /* End of setupEvals */

/*
* ==================================================================== 
* This routine reads a line and interprets it based on its first 
* character (a line format descriptor) as:
*  1) a triangular element
*  2) a quadrilateral element
*  3) a comment
*  4) garbage (unidentifiable)
* For elements the number of corners is returned.
* For a comment -1 is returned. 
* For garbage 0 is returned.
* 
* The format descriptors accepted are (no quotes should be in the 
* actual input files:
* '3' 't' 'T' : Triangular element (nCorners=3)
* '4' 'q' 'Q' : Quadrilateral element (nCorners=4)
* '*' '#' '%' : Comment line (nCorners=-1)
*
* Every other character (including a space!) will make the line 
* being interpreted as "garbage", and nCorners=0 will be returned.
*
* ==================================================================== */
static 
int identifyElementLine1(char *inputLine, char **restOfLine)
{
  char oneChar; /* One character  */
  int lineType;
  lineType=0;

  /* Read the initial character on the line */
  if (sscanf(inputLine, "%c",&oneChar)!=1){
    /* If one character cannot be read, then this line certainly 
     * qualifies as "garbage"                                         */
    lineType=0;
  }
  /*else if( !strcmp(oneChar,"T") || !strcmp(oneChar,"t") 
    || !strcmp(oneChar,"3")  ){*/
  else if( oneChar=='T' || oneChar=='t' || oneChar=='3'){
    /* This is a triangular element */
    lineType=3;
  }
  else if( oneChar=='Q' || oneChar=='q' || oneChar=='4'){
    /* This is a quadrilateral element */
    lineType=4;
  }
  else if( oneChar=='*' || oneChar=='#' || oneChar=='%'){
    lineType=-1;
  }
  else { /* Everything else is garbage */
    lineType=0;    
  }
  /* Strip initial character from line remainder */
  *restOfLine = &(inputLine[1]);

  return lineType;
} /* End of routine identifyElementLine1 */

/*
* ==================================================================== 
* This routine reads an integer (couold be changed to "name" later?) 
* from a given line (string) and returns the integer as the number 
* of the surface (boundary part) that the corresponding element is 
* on (fastcap generic file input format, second column). 
* ==================================================================== */
static 
int identifyBoundaryLine1(char *inputLine, char **restOfLine)
{
  char fctName[]="identifyBoundaryLine1";  /* Name of this routine     */
  char *strippedLine, *returnLine;
  int iBoundary, ic;
  /* Find first non-zero character [' ' (space) or '\t' (tab)]         */
  for (ic=0; inputLine[ic]==' '||inputLine[ic]=='\t';ic++);
  strippedLine = &(inputLine[ic]);
  /* Read the integer (if we can) */
  if (sscanf(strippedLine,"%d",&(iBoundary))!=1){
    fprintf(stderr,"%s: ERROR. Could not read integer of line:\n",fctName);
    fprintf(stderr,"%s\n",inputLine);
    _EXIT_;
  }
  /* Go to next zero character [' ' (space) or '\t' (tab)]             */
  for (ic=0; strippedLine[ic]!=' '&&strippedLine[ic]!='\t';ic++);
  returnLine = &(strippedLine[ic]);
  *restOfLine = returnLine;
  return iBoundary;
}

/*
* ==================================================================== 
* This routine reads the corners given in an input line for a 
* "constant element" (format type 1).
* The line format descriptor (the first character on the line) has 
* already been read and is not passed to this routine (only the 
* remainder of the string should be passed along). The number of corners 
* corresponding to the format descriptor is supplied as an argument. 
* iLine is the input line number, and is supplied for error handling only.
* The corner coordinates are stored in corners, which is of the form 
* "point corners[4]".
* The string (input line) passed to this routine must contain only 
* the coordinates of the corners of the panel to be read on the form:
*    x1 y1 z1  x2 y2 z2  x3 y3 z3  [x4 y4 z4]
* The last set of coordinates is ommitted if nCorners==3.
* ==================================================================== */
static 
void readCorners1(char *inputLine,int nCorners,int iLine,point *corners)
{
  char fctName[]="readCorners1";  /* Name of this routine              */
  int nRead, nExcessNodes;
  double dummy;

  nExcessNodes = 0; /* Change if a wrong number of nodes is obtained    */

  if (nCorners==3){
    /* Make sure that the line has the correct format. It should 
     * contain only nine reals. Complain if there are too 
     * many reals (this could indicate an input error). However, 
     * accept any string following the nine reals (not starting by 
     * a valid real).                                                  */
    nRead = sscanf(inputLine,
		   "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		   &(corners[0][0]),&(corners[0][1]),&(corners[0][2]),
		   &(corners[1][0]),&(corners[1][1]),&(corners[1][2]),
		   &(corners[2][0]),&(corners[2][1]),&(corners[2][2]),
		   &dummy );
    nExcessNodes = nRead-9;
  }
  else if (nCorners==4){
    nRead = sscanf(inputLine,
		   "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		   &(corners[0][0]),&(corners[0][1]),&(corners[0][2]),
		   &(corners[1][0]),&(corners[1][1]),&(corners[1][2]),
		   &(corners[2][0]),&(corners[2][1]),&(corners[2][2]),
		   &(corners[3][0]),&(corners[3][1]),&(corners[3][2]),
		   &dummy );
    nExcessNodes = nRead-12;
  }
  else { /* This is a weird panel (unknown type). Report error and 
	  * exit (this error should never occur).                      */
    fprintf(stderr,
	    "%s: ERROR: What kind of panel has %d corners?\n"
	    "  This halts the program.\n",
	      fctName,iLine);
    _EXIT_;
  }
  /* Deal with any unexpected read results                             */
  if (nExcessNodes<0){
    /* Too few nodes cannot be accepted.                               */
    fprintf(stderr,
	    "%s: ERROR: Too few corner nodes found on line %d:\n"
	    "%s"
	    "         This error is terminal!!\n",
	    fctName,iLine,inputLine);
    _EXIT_;
  }
  else if (nExcessNodes>0){
    /* Too many nodes can be accepted - by discarding the excess ones.
     * However, this is considered extremely bad practice, so issue a 
     * warning (for every input line with this property)               */
    fprintf(stderr,
	    "%s: WARNING: Too many corner nodes found on line %d:\n\n"
	    "%s\n"
	    "    This *could* indicate a problem!! (continuing)\n",
	    fctName,iLine,inputLine);
  }
  return;
} /* End of routine readCorners1                                       */

/*
* ==================================================================== 
* This routine reads the corners given in an input line for a 
* "constant element". nCorners (the first integer on the line) has 
* already been read, and is supplied as an argument. 
* iLine is the input line number, and is supplied for error handling only.
* The corner coordinates are stored in corners, which is of the form 
* "point corners[4]".
*
* Presently, this routine is not used and it can probably be deleted.
* ==================================================================== */
static 
void readCorners(char *inputLine,int nCorners,int iLine,point *corners)
{
  char fctName[]="readCorners";  /* Name of this routine               */
  int idum;
  double dummy;
  /* A number of different input formats may be accepted here. 
   * At the present time, it is not quite clear how to chose between 
   * each input format, but that can't really be important at the 
   * present stage of implementation. 
   * For now the format of "FastWamit" will be adapted - in particular 
   * the input file "unit25.dat" for the MOB is used as example.
   * The first line is a title line for later identification of the 
   * calculations. 
   * Every subsequent line starts by an integer (3 or 4) denoting 
   * wether the element is a triangle (3) or a quadrilateral (4). 
   * The integer "k" is followed by 3*k reals giving the corners of the
   * element as: x1 y1 z1  x2 y2 z2  x3 y3 z3  [x4 y4 z4]
   * 
   * Another format is the following (p.t. not supported):
   * Each line will start by a "Q", "q", "T" or "t" to denote
   * a quadrilateral or triangular element. Then follows an integer 
   * (which is usually 1 - I guess this denotes what part of the 
   * boundary the panel is part of) and nine or twelve reals 
   * giving the coordinates of the corners of the element.
   * These files typically are denoted *.in                            */

  if(nCorners==3){
    /* Make sure that the line has the correct format. It should 
     * contain one integer and nine reals. Complain if there are too 
     * many reals (this could indicate an input error). However, 
     * accept any string following the nine reals (not starting by 
     * a valid real).                                                  */
    if (sscanf(inputLine,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	       &idum,
	       &dummy,&dummy,&dummy,
	       &dummy,&dummy,&dummy,
	       &dummy,&dummy,&dummy,
	       &dummy)==11){
      fprintf(stderr,
	      "%s: WARNING: Too many corner nodes found on line %d:\n"
	      "%s"
	      "         This *could* indicate a problem!! (continuing)\n",
	      fctName,iLine,inputLine);
    }
    if (sscanf(inputLine,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	       &idum,
	       &(corners[0][0]),&(corners[0][1]),&(corners[0][2]),
	       &(corners[1][0]),&(corners[1][1]),&(corners[1][2]),
	       &(corners[2][0]),&(corners[2][1]),&(corners[2][2]))!=10){
      fprintf(stderr,
	      "%s: ERROR: Too few corner nodes found on line %d:\n"
	      "%s"
	      "         This error is terminal!!\n",
	      fctName,iLine,inputLine);
      _EXIT_;
    }
  }
  else if (nCorners==4){
    /* Make sure that the line has the correct format. It should 
     * contain one integer and nine reals. Complain if there are too 
     * many reals (this could indicate an input error). However, 
     * accept any string following the nine reals (not starting by 
     * a valid real).                                                */
    if (sscanf(inputLine,"%d "
	       "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	       &idum,
	       &dummy,&dummy,&dummy,
	       &dummy,&dummy,&dummy,
	       &dummy,&dummy,&dummy,
	       &dummy,&dummy,&dummy,
	       &dummy)==14){
      fprintf(stderr,
	      "%s: WARNING: Too many corner nodes found on line %d:\n"
	      "%s"
	      "         This *could* indicate a problem!! (continuing)\n",
	      fctName,iLine,inputLine);
    }
    if (sscanf(inputLine,"%d "
	       "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	       &idum,
	       &(corners[0][0]),&(corners[0][1]),&(corners[0][2]),
	       &(corners[1][0]),&(corners[1][1]),&(corners[1][2]),
	       &(corners[2][0]),&(corners[2][1]),&(corners[2][2]),
	       &(corners[3][0]),&(corners[3][1]),&(corners[3][2]))!=13){
      fprintf(stderr,
	      "%s: ERROR: Too few corner nodes found on line %d:\n"
	      "%s"
	      "         This error is terminal!!\n",
	      fctName,iLine,inputLine);
      _EXIT_;
    }
  }
  else{ 
    /* This position should never be reached - even if the input file 
     * contains errors. The errors should be captured earlier.     */
    fprintf(stderr,
	    "%s: ERROR: %d corners in element??\n"
	    "       This error is terminal!\n",
	    fctName,nCorners);
    _EXIT_;
  }
  return;
} /* End of routine readCorners                                        */

/*
* ==================================================================== 
* This routine determines wheter or not a quadrilateral element 
* (defined by its corner points) needs cracking, i.e. being split
* into two triangular elements. 
* If the the panel does not need cracking, then the return value is 
* zero.
* If the element needs cracking, then the return status is the 
* number of the diagonal along which cracking is to be performed:
*  diagonal 1 connects corners 0 and 2
*  diagonal 2 connects corners 1 and 3
* 
* Presently cracking may be caused by one of two different effects:
*
* 1) If the element is not planar within tolerance, it will be cracked.
*    The two new triangles wil each be planar. The new element edge will 
*    be along the shortest diagonal of the original quad, since this 
*    gives the best (lowest) aspect ratios (and smallest bounding 
*    spheres) of the triangular panels.
* 
* 2) If a quadrilateral element has one corner angle larger than
*    pi (180 deg), then the centroid may not be inside the panel. 
*    Also there may be other problems with this kind of panel.
*    Presently, this kind of panel is not accepted: The panel will 
*    be cracked into two triangles. Cracking should be performed along 
*    that diagonal which makes the two new triangles have normals that 
*    point in the same direction. THIS PART NOT YET IMPLEMENTED!
*
* CAVEAT: If an element is both skew AND has a angle larger than pi, 
* then it may not be cracked along the "most reasonable" diagonal. 
* This may result in the two new elements having normal that are *almost* 
* in opposite directions. Possibly, the normals of the new elements to 
* see along what direction cracking is reasonable. (make sure n1 dot n2 
* is positive!)
* ==================================================================== */
static
int quadNeedsCracking(point *corners, double crackTol)
{
  /* char fctName[]="quadNeedsCracking";  Name of this routine         */
  int shortDiag;
  point vec[2],normal[2];  /* Work vectors                             */
  double ldiag[2], invNorm;
  int trigPts[2][3];       /* Indices of corners of cracked element    */

  /* Finding the shortest diagonal (to crack element by if skew)       */
  ldiag[0] = VDIST2(corners[0],corners[2]);
  ldiag[1] = VDIST2(corners[1],corners[3]);
  if (ldiag[0]<=ldiag[1]){
    shortDiag=1;   /* Return value if cracking is needed               */
    /* First triangle: */
    trigPts[0][0]=0;
    trigPts[0][1]=1;
    trigPts[0][2]=2;
    /* Second triangle: */
    trigPts[1][0]=0;
    trigPts[1][1]=2;
    trigPts[1][2]=3;
  }
  else {
    shortDiag=2;   /* Return value if cracking is needed               */
    /* First triangle: */
    trigPts[0][0]=0;
    trigPts[0][1]=1;
    trigPts[0][2]=3;
    /* Second triangle: */
    trigPts[1][0]=1;
    trigPts[1][1]=2;
    trigPts[1][2]=3;
  }
  /* Find the normal vector on each of the two element parts           */
  /*  Vectors along triangle sides:                                    */
  VOP(vec[0],corners[trigPts[0][1]],-,corners[trigPts[0][0]]);
  VOP(vec[1],corners[trigPts[0][2]],-,corners[trigPts[0][0]]);
  /*  First normal vector                                              */
  XPROD(normal[0],vec[0],vec[1]);
  /*  Vectors along triangle sides:                                    */
  VOP(vec[0],corners[trigPts[1][1]],-,corners[trigPts[1][0]]);
  VOP(vec[1],corners[trigPts[1][2]],-,corners[trigPts[1][0]]);
  /*  Second normal vector                                             */
  XPROD(normal[1],vec[0],vec[1]);
  /* Normalize the normal vectors                                      */
  invNorm = 1.0/VNORM(normal[0]);
  normal[0][0] = invNorm*normal[0][0];
  normal[0][1] = invNorm*normal[0][1];
  normal[0][2] = invNorm*normal[0][2];
  invNorm = 1.0/VNORM(normal[1]);
  normal[1][0] = invNorm*normal[1][0];
  normal[1][1] = invNorm*normal[1][1];
  normal[1][2] = invNorm*normal[1][2];
  /* Compare the normal vectors to see if the points 
   * on the quadrilateral are all in the same plane                    */
  VOP(vec[0],normal[0],-,normal[1]);

  if( VNORM(vec[0]) <= crackTol){
    /* In this case we don't need to crack this element.               */
    return 0;
  }
  return shortDiag;
} /* End of routine quadNeedsCracking                                  */

/*
* ==================================================================== 
* Crack one quadrilateral element into two triangular elements.
* Input:
*  iElement:  Number of elment to crack.
*  iElement2: Number of element to store the new cracked element.
*  iDiag:     Number of the diagonal to crack by.
*  elements:  Storage of the elements.
* ==================================================================== */
static
void crackElement(int iElement,
		  int iElement2,
		  int iDiag,
		  struct constantElement **elements){
  char fctName[]="crackElement";  /* Name of this routine             */
  /* Set up the *new* element (first) to avoid overwriting data       */
  if (iDiag==1){
    /* The diagonal between corners 0 and 2 are the shortest          */
    /* The new elements has corners (0,2,3)                           */
    elements[iElement2]->shape = 3;
    elements[iElement2]->corners[0][0]=elements[iElement]->corners[0][0];
    elements[iElement2]->corners[0][1]=elements[iElement]->corners[0][1];
    elements[iElement2]->corners[0][2]=elements[iElement]->corners[0][2];
    elements[iElement2]->corners[1][0]=elements[iElement]->corners[2][0];
    elements[iElement2]->corners[1][1]=elements[iElement]->corners[2][1];
    elements[iElement2]->corners[1][2]=elements[iElement]->corners[2][2];
    elements[iElement2]->corners[2][0]=elements[iElement]->corners[3][0];
    elements[iElement2]->corners[2][1]=elements[iElement]->corners[3][1];
    elements[iElement2]->corners[2][2]=elements[iElement]->corners[3][2];
    /* The "old" elements has corners (0,1,2), i.e. it does not need 
     * change                                                          */
    elements[iElement]->shape = 3;
  }
  else if(iDiag==2){
    /* The diagonal between corners 1 and 3 are the shortest           */
    /* The new elements has corners (1,2,3)                            */
    elements[iElement2]->shape = 3;
    elements[iElement2]->corners[0][0]=elements[iElement]->corners[1][0];
    elements[iElement2]->corners[0][1]=elements[iElement]->corners[1][1];
    elements[iElement2]->corners[0][2]=elements[iElement]->corners[1][2];
    elements[iElement2]->corners[1][0]=elements[iElement]->corners[2][0];
    elements[iElement2]->corners[1][1]=elements[iElement]->corners[2][1];
    elements[iElement2]->corners[1][2]=elements[iElement]->corners[2][2];
    elements[iElement2]->corners[2][0]=elements[iElement]->corners[3][0];
    elements[iElement2]->corners[2][1]=elements[iElement]->corners[3][1];
    elements[iElement2]->corners[2][2]=elements[iElement]->corners[3][2];
    /* The "old" elements has corners (0,1,3) [only the last changes]  */
    elements[iElement]->shape = 3;
    elements[iElement]->corners[2][0]=elements[iElement]->corners[3][0];
    elements[iElement]->corners[2][1]=elements[iElement]->corners[3][1];
    elements[iElement]->corners[2][2]=elements[iElement]->corners[3][2];
  }
  else {
    /* This is bad. The shortest diagonal should be numbered 1 or 2.  */
    fprintf(stderr,"%s: What kind of diagonal is this? %d ... Fix me!\n"
	    "        This error is terminal.\n",fctName,iDiag);
    _EXIT_;
  }
  /* For good measure, reset the now redundant fourth corner point 
   * in the original element.                                          */
  elements[iElement]->corners[3][2]=
    elements[iElement]->corners[3][1]=
    elements[iElement]->corners[3][0]= 0.0;
  /* Copy various information from the old to the new element.         */
  elements[iElement2]->boundaryNumber = 
    elements[iElement]->boundaryNumber;

  /* Recall to cross-reference the indices of the cracked pair         */
  elements[iElement]->crackedNeighbourIndex = iElement2;
  elements[iElement2]->crackedNeighbourIndex = iElement;
  return;
} /* End of routine crackElement                                       */

/* 
* ==================================================================== 
* Calculate the centroid of each element supplied.
* ==================================================================== */
static
void setupElementCentroids(struct constantElement **elements,
			   int nElements){
  int i;
  /* These calculations assume that the coordinates of the 
   * fourth corner point of a triangular panel is set to zero.         */
  for (i=0;i<nElements;i++){
    /* Calculate centroid */
    setupElementCentroid(elements[i]);
  }
  return;
} /* End of routine setupElementCentroids                              */

/*
* ==================================================================== 
* Calculate the centroid of one constant-element.
* ==================================================================== */ 
static
void setupElementCentroid(struct constantElement *element)
{
  char fctName[]="elementCentroid";  /* Name of this routine           */
  double cent0[3], cent1[3], *centroids[2], side0[3], side1[3], 
    normal0[3], normal1[3], *normals[2], areas[2];
  double area2, areaNorm;
  int i, it, triag0[3],triag1[3], *triangles[2];
  double ONE_THIRD = 0.3333333333333333;
  point *corners;
  double *centroid;
  /* Shorthand notation: */
  corners = element->corners;
  centroid= element->centroid;

  if (element->shape==3){
    /* This is a triangle, so the centroid is the average
     * of the corner points:                                           */
    for (i=0; i<3; i++){
      centroid[i] = ONE_THIRD*
	               (corners[0][i]+corners[1][i]+corners[2][i]);
    }
  }
  else if (element->shape==4){
    /* This is a quadrilateral. The centroid is a little bit more 
     * complicated in this setting. 
     * Set up the quadrilateral as the sum (or possibly difference!?!) 
     * of two triangles. Calculate the centroid and area of each 
     * triangle. The centroid of the quadrilateral is then obtained 
     * as the average of the two triangle centroids weighted by the 
     * triangle areas, i.e.
     *   C_q = (a0*C0+a1*C1)/(a0+a1)
     * To capture oddities areas are computed including sign.          */
    /* Set up arrays to use (this eliminates need for static storage, 
     * and/or allocation/reallocation in the present routine)          */
    centroids[0] = cent0;
    centroids[1] = cent1;
    normals[0]   = normal0;
    normals[1]   = normal1;
    triangles[0] = triag0;
    triangles[1] = triag1;
    /* Set two corner points of the two triangles to consider:         */
    triangles[0][0] = 0;
    triangles[0][1] = 1;
    triangles[0][2] = 2;
    triangles[1][0] = 0;
    triangles[1][1] = 2;
    triangles[1][2] = 3;
    /* For each triangle calculate centroid                            */
    for (it=0; it<2; it++){
      for (i=0; i<3; i++){
	centroids[it][i] = ONE_THIRD*
	  ( corners[triangles[it][0]][i]
           +corners[triangles[it][1]][i]
           +corners[triangles[it][2]][i]);
	/* Vector along two sides: */
	side0[i] = corners[triangles[it][1]][i]
	          -corners[triangles[it][0]][i];
	side1[i] = corners[triangles[it][2]][i]
	          -corners[triangles[it][0]][i];
      }
      /* Normal (times area times two):                               */
      XPROD(normals[it],side1,side0);
      /* And area (times two)                                         */
      areas[it] = VNORM(normals[it]);
    }
    /* Test if the normals point in the same or in opposite directions:*/
    area2 = DOTPROD(normals[0],normals[1]);

    /* Check if the dot product is OK (i.e. if the normals are collinear)*/
    assert(ABS(ABS(area2)-areas[0]*areas[1])<1.0e-6);

    if( area2 < 0.0){
      /* Change sign on ONE of the two areas (which one doesn't matter */
      /* NOTE: It is probably BAD if the normals point in different 
       * directions. It means that the panels "look weird". Also, it 
       * means that the centroid of the panel may not be located at 
       * the panel, which is definitely bad if the centroid is used as 
       * a collocation point! */
      areas[1] *= -1.0;
    }
   
    /* Calculate quad centroid as weighted average of the centroids of 
     * the two triangles:                                              */
    areaNorm = 1.0/(areas[0]+areas[1]);
    for (i=0; i<3; i++){
      centroid[i] = areaNorm*
	            (areas[0]*centroids[0][i]+areas[1]*centroids[1][i]);
    }
  }
  else{ /* Weird panel shape? */
    fprintf(stderr,"%s: Weird element shape: %d\n",
	    fctName,element->shape);
    _EXIT_;
  }
  return;
} /* End of routine setupElementCentroid */

/* 
* ==================================================================== 
* Setup values for the boundary condition types.
* ==================================================================== */
static
void setupElementBoundaryConditionTypes(struct constantElement **elements,
					int nElements,
					char *conditionType,
					int idata){
  int i;
  for (i=0;i<nElements;i++){
    /* Setup type for one element based on the coordinate at its 
     * centroid.                                                       */
    elements[i]->bcType = 
      (short int) boundaryConditionType(elements[i]->centroid,
					(int) elements[i]->boundaryNumber,
					conditionType,
					idata);
  }
  return;
} /* End of routine setupElementBoundaryConditionTypes                 */

/* 
* ==================================================================== 
* Setup boundary condition values for each element.
* ==================================================================== */
void setupElementBoundaryConditionValues(void *sourcesIn,
					 int nSources,
					 void *rightHandSideIn,
					 char *conditionType,
					 int idata){
  struct constantElement **elements;
  double *rightHandSide, normal[3], *direction;
  int ie;

  /* Cast input data to correct type */
  elements      = (struct constantElement**) sourcesIn;
  rightHandSide = (double*) rightHandSideIn;

  for (ie=0;ie<nSources;ie++){
    /* Setup type for one element based on the coordinate 
     * at its centroid.                                                */
    if ((int)elements[ie]->bcType == 0){ /* Neumann boundary condition */
      elementNormal(elements[ie], normal);
      direction=normal;
    }
    else direction=NULL; /* The "direction" pointer is not strictly 
			  * needed, but it makes sure that a wrong 
			  * normal direction is not used when setting 
			  * boundary conditions.                      */
    
    rightHandSide[elements[ie]->directIndex] = 
      boundaryConditionValue((int)elements[ie]->bcType, 
			     elements[ie]->centroid,
			     direction, 
			     (int) elements[ie]->boundaryNumber,
			     conditionType,
			     idata);
  }
  return;
} /* End of routine setupElementBoundaryConditionValues                */

/*
* ==================================================================== 
* Calculate the normal of one constant-element. Since, during the 
* process of finding the normal, the panel area can be obtained almost 
* for free, the panel area is returned as the function value.
* Ideas have been taken from "calcp". 
*
* By definition(!) the normal direction is set to follow left-hand rule 
* (clockwise ordered points has normal pointing up), just as in calcp.
* ==================================================================== */ 
static
double elementNormal(struct constantElement *element, double *normal)
{
  double X[3],Y[3];
  double nnorm;
  int i;
  point *corners;
  /* Shorthand notation: */
  corners = element->corners;

  /* Use the "diagonals" to find the normal direction.                 */
  /* Find the vectors along the "diagonals", calculating  
   * partially the panel coordinate system.                            */
  VOP(X, corners[2], -, corners[0]);
  if(element->shape == 3) {
    VOP(Y, corners[1], -, corners[0]); 
  } else {
    VOP(Y, corners[1], -, corners[3]); 
  }
  /* Form the vector product (cross-product) of these vectors 
   * (this is the normal direction) */
  XPROD(normal, X, Y);
  /* Normalize the vector */
  nnorm = VNORM(normal);
  for (i=0;i<3;i++)  normal[i] /= nnorm;
  /* Return the element area (half of the length of the normal 
   * vector before normalization)                                      */
  return 0.5*nnorm;
}/* End of routine elementNormal */


/* 
* ==================================================================== 
* Assigns each source element with one or more columns and each 
* evaluation element with one or more rows of the of the direct matrix.
* Stores the associations with the elements.
*
* If the spurces and the evals are the same, then the association with
* rows and columns will be the same (this is a choise). I.e. in this 
* case, if element i as a source is associated with column j, then 
* element i as an eval will be associated with row j.
*
* Allocates the memory for numNonzeroInCols in which the number of 
* nonzero entries in each column later will be stored.
* 
* Returns the total number of columns needed in the direct part.
* Note that the number of rows needed is found, but not used nor 
* returned to the calling routine.
* ==================================================================== */ 
void assignElementDirectIndices(void *sourcesIn, int nSources,
				void *evalsIn,   int nEvals,
				int *nColsDirect, int *nRowsDirect)
{
  /* char fctName[]="assignElementDirectIndices"; */
  int nCols, nRows;
  struct constantElement **sources, **evals;

  /* Recast input pointers                                             */
  sources = (struct constantElement**) sourcesIn;
  evals   = (struct constantElement**) evalsIn;

  /* Count and index the source elements                               */
  nCols = assignSourceDirectIndices(sources, nSources);
  
  /* Count and index the evaluation elements.
   * This is not needed if the sources and the evals are equal.        */
  if (sources != evals || nSources != nEvals)
    nRows = assignEvalDirectIndices(evals, nEvals);
  else
    nRows = nCols;

  /* Return the total number of unknowns connected with 
   * sources and evals.                                                */
  *nColsDirect = nCols;
  *nRowsDirect = nRows;

  return;
}

/* 
* ==================================================================== 
* Find the number of nonzero entries in each column of the 
* (sparse) direct matrix:
*
* Given a particular source element and  list of its neighbours (evals)
* fill in the number of needed rows in the entries corresponding 
* to the columns used by the source in the direct (sparse) matrix. 
* This information is needed to be able to initialize the direct matrix.
* ==================================================================== */ 
void calcNumNonzeroDirectEntries(void **sourcesIn,
				 void **evalsIn, 
				 int nSources,
				 int nTotNeighbours,
				 void *colsData)
{
  char fctName[]="calcNumNonzeroDirectEntries"; 
  struct constantElement *source,**evals;
  int *numNonzeroInColsDirect;

  assert(nSources==1); /* One source only, please! */

  source = (struct constantElement*)  sourcesIn[0];
  evals  = (struct constantElement**) evalsIn;
  numNonzeroInColsDirect = (int*) colsData;

  /* TEST OUTPUT: 
  {
    int i;
    fprintf(stdout,"%s: For src (column) %d I have been "
	    "handed the neighbourlist:\n",
	    fctName,source->directIndex);
    for (i=0;i<nTotNeighbours;i++)
      fprintf(stdout," %d",evals[i]->directIndex);
    fprintf(stdout,"\n");
  }*/

  /* The following test assumes that calloc has been used for 
   * allocation (or that the elements of numNonzeroInColsDirect has 
   * been zeroed initially by other means. 
   * An entry, which is to be changed, should be zero. If not so, then
   * it has been changed before, i.e. this routine may be overwriting 
   * important information.                                           */
  assert(numNonzeroInColsDirect[source->directIndex]==0);

  /* For the present case each source needs only one column, 
   * and the number of nonzero entries equals the number of 
   * neighvbouring evals (i.e. the actual INDICES of the evals 
   * are presently not used, but it may well be used later!)          */
  numNonzeroInColsDirect[source->directIndex] = nTotNeighbours;

  return;
} /* End of calcNumNonzeroDirectEntries */

/* 
* ==================================================================== 
* Find the number of nonzero entries in each column of the 
* (sparse) matrix for preconditioning the direct matrix:
*
* Given a list source elements and list of their neighbours (evals)
* count the number of evaluation elements in the list and store
* this number at the column for each source element passed.
*
* This information is needed to be able to initialize the 
* preconditioning matrix.
*
* NOTE: This procedure is used for the preconditioning part and should 
* only be used when sources and evals "are the same". If the sources and 
* the evals are defined independently, then it may be difficult to set
* up a square system for the preconditioning.
*
* ==================================================================== */ 
void calcNumNonzeroPrecondEntries(void **sourcesIn,
				 void **evalsIn, 
				 int nSources,
				 int nTotNeighbours,
				 void *colsData)
{
#if defined(MYDEBUGFLAG)
  char fctName[]="calcNumNonzeroPrecondEntries";
  int verbose=0;
#endif
  struct constantElement **sources,**evals;
  int *numNonzeroInColsDirect;
  int i;

  sources = (struct constantElement**) sourcesIn;
  evals   = (struct constantElement**) evalsIn;
  numNonzeroInColsDirect = (int*) colsData;

  /* Output block for debugging runs. May further be toggled using 
   * the "verbose" flag.                                               */
#if defined(MYDEBUGFLAG)
  if (verbose){
    printf("%s: Passed along sources (%d):\n",fctName,nSources);  
    for (i=0;i<nSources;i++){
      printf("  %d\n",sources[i]->directIndex);
    }
    printf("%s: Passed along evals (%d):\n",fctName,nTotNeighbours);  
    for (i=0;i<nTotNeighbours;i++){
      printf("  %d\n",evals[i]->directIndex);
    }
  }
#endif

  /* For each of the sources listed, the preconditioner will have 
   * nTotNeighbours non-zero entries.                                  */
  for (i=0;i<nSources;i++){
    numNonzeroInColsDirect[sources[i]->directIndex]=nTotNeighbours;
  }
  return;
} /* End of calcNumNonzeroPrecondEntries                               */

/* 
* ==================================================================== 
* Associate the source elements with columns in the direct matrix.
* Store the association(s) with each element. Note that each element 
* may well be associated with more than one column, as long as the 
* definition of the elements support this.
* Make sure that each column is associated with one and only one 
* source element.
* Return the total number of columns needed.
* ==================================================================== */ 
static
int assignSourceDirectIndices(struct constantElement **sources, 
			      int nSources)
{
  int isrc, ncols;
  /* For now, just let the column order be the same as the source order,
   * and assign one column per source                                  */
  for (isrc=0; isrc<nSources; isrc++){
    sources[isrc]->directIndex = isrc;
  }
  ncols = nSources;
  return ncols;
} /* End of assignSourceDirectIndices                                  */

/* 
* ==================================================================== 
* Associate the evaluation elements with rows in the direct matrix.
* Store the association(s) with each element. Note that each element 
* may well be associated with more than one row, as long as the 
* definition of the elements support this.
* Make sure that each row is associated with one and only one 
* evaluation element.
* Return the total number of rows needed. At the moment the total number
* of rows is not used for anything, but the chosen layout emphasizes the
* similarity between the source and the eval associations to the 
* direct matrix. 
* ==================================================================== */ 
static
int assignEvalDirectIndices(struct constantElement **evals, 
			      int nEvals)
{
  int ieval, nrows;
  /* For now, just let the row order be the same as the eval order,
   * and assign one row per eval                                       */
  for (ieval=0; ieval<nEvals; ieval++){
    evals[ieval]->directIndex = ieval;
  }
  nrows = nEvals;
  return nrows;
} /* End of assignEvalDirectIndices                                    */

/*
* ==================================================================== 
* Given the source elements and a element number return a pointer
* to this particular element.
* ==================================================================== */ 
void *pointerToSource(void *sourcesIn,int iSource)
{
  struct constantElement **sources;
  sources = (struct constantElement**) sourcesIn;
  return (void*) sources[iSource];
}

/*
* ==================================================================== 
* Given the evaluation elements and a element number return a pointer
* to this particular element.
* ==================================================================== */ 
void *pointerToEval(void *evalsIn,int iEval)
{
  struct constantElement **evals;
  evals = (struct constantElement**) evalsIn;
  return (void*) evals[iEval];
}

/* 
* ==================================================================== 
* Given a source element calculate the center and radius of its 
* bounding sphere.
* ==================================================================== */ 
void findSourceBoundingSphere(point center,
			      double *radius,
			      void *sourceIn)
{
  int i;
  /* At the moment this routine is the same as 
   * findEvalBoundingSphere, but it need not to be that way.
   * To exploit knowledge of the differences between sources and 
   * evals two separate routines may be used. */
  struct constantElement *source;

  /* Cast input pointer to the element type */
  source = (struct constantElement *) sourceIn;

  /* For the constant element use the centroid as the center of the 
   * bounding sphere. This will not always lead to the smallest radius 
   * of the bounding sphere, but is gives a good estimate of the 
   * location of the element, and the bounding sphere can be made 
   * relatively tight.                                                 */
  center[0] = source->centroid[0];
  center[1] = source->centroid[1];
  center[2] = source->centroid[2];

  /* Find largest distance from centroid to corner point:              */
  *radius = 0.0;
  for (i=0;i<source->shape;i++){
    *radius = MAX(*radius, VDIST2(source->centroid,source->corners[i]));
  }
  *radius = sqrt(*radius);
  return;
} /* End of findSourceBoundingSphere */

/* 
* ==================================================================== 
* Given an evaluation element calculate the center and radius of 
* its bounding sphere.
* ==================================================================== */ 
void findEvalBoundingSphere(point center,
			    double *radius,
			    void *evalIn)
{
  /* At the moment this routine is the same as 
   * findSourceBoundingSphere, but it need not to be that way.
   * To exploit knowledge of the differences between sources and 
   * evals two separate routines may be used. */
  struct dummyElement *eval;
  
  /* Cast input pointer to the element type */
  eval = (struct dummyElement *) evalIn;

  /* For the dummy element type the values are directly 
   * accesible: */
  center[0] = eval->center[0];
  center[1] = eval->center[1];
  center[2] = eval->center[2];
  *radius   = eval->radius;

  return;
} /* End of findEvalBoundingSphere */

/*
* ==================================================================== 
* Given an source element pointer, a Type, and a center, returns 
* integrals of polynomials over the element basis functions.
* elementIn   - pointer to the element
* elementType - "eval" or "source"
* sideType    - "rhs" or "lhs" for right-hand-side or left-hand-side
*               (system matrix) build.
* xc,yc,zc - coordinate system origin.
* moments
*
* NOTE: This routine does not free the memory that it allocates.
* The memory is stored and reused, and freed only if more memory 
* is needed. This should eventually be changed, so that the memory 
* can be freed after the last call.
*
* TODO: Use TYPE and a second string (to be made) such as "rhs" or "lhs"
* to determine whether monopoles or derivatives should be used.
* Most probably, the boundary condition type specifier can be used, but 
* is it enough?
* Evals are used for evaluation, i.e. interpolation is made to evals. 
* The interpolation does not depend on rhs/lhs concepts, it is 
* always the same(?)
* Source elements are projected onto the grid, either on the form 
* INT(G phi_n) or INT(G_n phi). Thus moments of the element are needed 
* either as moments in the normal sence, or calculated by integrating the 
* normal derivative of the the basis functions (monomials) over the element
* (depending on the G vs. G_n part).
* For "source": DIRICLET boundary condition states that phi is known, 
* project monopoles. NEUMANN boundary condition specifies that phi_n 
* is known, i.e. project derivatives!
* NO! This is no good. MUST think if we are setting up the right-hand-side 
* or the system matrix itself. We DO need this information, since it 
* will reverse the whole concept!
* ==================================================================== */ 
void calcElementMoments(void *elementIn, 
			char *elementType,
			char *systemType,
			double xc, double yc, double zc,
			double polyScale, 
			double **moments, 
			int nterms, int polyOrder, int maxOrder1D,
			int *imonoOrder,int *jmonoOrder,int *kmonoOrder)
{
  char fctName[]="calcElementMoments";
  struct constantElement *element;
  int i,j;
  point collocationPoint;
  double normal[3], *direction;
  static point corners[3];    /* Used when nterms==27                  */
  static double ***monoMoments=NULL, *monoMomentsStore=NULL;
  static int iSize=-1;

  element = (struct constantElement *) elementIn;

  /* Setup local storage if not of the correct size:                   */
  if (iSize != maxOrder1D+1) { /* Could be "if(maxOrder1D+1>iSize)" */
    /* Free memory if necessary */
    if (monoMoments!=NULL){
      for (i=0;i<iSize;i++) FREE(monoMoments[i]);
      FREE(monoMoments);
      FREE(monoMomentsStore);
    }
    /* Setup local storage                                             */
    /* Size of moment arrays etc:                                      */
    iSize = maxOrder1D+1;
    /* Allocate memory for integration of monomials over element       */
    monoMomentsStore = (double *) calloc(iSize*iSize*iSize,
					   sizeof(double));
    monoMoments  = (double ***) calloc(iSize,sizeof(double **));
    for (i=0;i<iSize;i++){
      monoMoments[i] = (double **) calloc(iSize,sizeof(double *));
      for (j=0;j<iSize;j++){
	monoMoments[i][j] = &(monoMomentsStore[i*iSize*iSize +j*iSize]);
      }
    }
  }/* Local storage is now in place                                    */
  /* Zero out the (local) moment storage                               */
  for (i=0;i<iSize*iSize*iSize;i++){
    monoMomentsStore[i]=0.0;
  }

  if(strcmp(elementType,"eval") == 0){
    /* Assumes collocation by using delta function in moment           */
    /* Use the centroid as collocation point                           */
    VASSIGN(collocationPoint,element->centroid);
    /* Translate to specified origin                                   */
    collocationPoint[0] = collocationPoint[0] - xc;
    collocationPoint[1] = collocationPoint[1] - yc;
    collocationPoint[2] = collocationPoint[2] - zc;
    /* Calculate monomial values (use direction=NULL)*/
    addMonomialValues(collocationPoint,monoMoments,polyScale,
		      polyOrder,maxOrder1D,
		      NULL);
  }
  else if(strcmp(elementType,"source") == 0) {
    /* Determine if we are to calculate moments using monopoles or 
     * derivatives of monopoles                                        */
    if ( (strcmp(systemType,"rhs")==0 && (int)element->bcType==0)||
	 (strcmp(systemType,"lhs")==0 && (int)element->bcType==1)){
      /* This is a INT(G phi_) term (right hand side for a Neumann  
       * condition or left hand side for a Diriclet condition).
       * Thus we need to project the "monopole", i.e. find moments 
       * based on the monomials.                                      */
      direction = NULL;
    }
    else if ( (strcmp(systemType,"rhs")==0 && (int)element->bcType==1)||
	      (strcmp(systemType,"lhs")==0 && (int)element->bcType==0)){
      /* This is a INT(G_n phi) term (right hand side for a Diriclet 
       * condition or left hand side for a Neumann condition).
       * Thus we need to project the "dipole" - INT(G_n), i.e. find 
       * moments based on normal derivatives of monomials              */
      /* Find normal direction of this element                         */
      elementNormal(element, normal);
      /* TEST:
      for (i=0;i<3;i++) normal[i] *= -1.0;*/
      direction = normal;
    }
    else { /* This is unknown = bad. Complain and exit!                */
      fprintf(stderr,
	      "%s: What kind of system type and boundary condition is this?\n"
	      "   systemType=\"%s\", bcType=%d\n"
	      "   This halts the program.\n",
	      fctName,systemType,(int)element->bcType);
    }
   /* Setup corners for the panel.                                    */
    for (i=0;i<3;i++){
      /* Get data from element structure                               */
      VASSIGN(corners[i],element->corners[i]);
      /* Translate to specified origin */
      corners[i][0] -= xc;
      corners[i][1] -= yc;
      corners[i][2] -= zc;
    }
    /* Calculate moments (integrals of monomials; use direction=NULL)  */
    addTriangularElementMoments(corners,monoMoments,polyScale,
				polyOrder,maxOrder1D,
				direction);

    /* For quadrilateral elements repeat the procedure to add the 
     * moments of the second half of the quad.                       */
    if (element->shape==4) {
      /* Repeat the setup for the second part of the quadrilateral   */
      for (i=0;i<3;i++){
	/* Really, we want the corners no. 0,2,3 so:                 */
	if (i==0) j=0;
	else j=i+1;
	/* Get data from element structure                           */
	VASSIGN(corners[i],element->corners[j]);
	/* Translate to specified origin */
	corners[i][0] -= xc;
	corners[i][1] -= yc;
	corners[i][2] -= zc;
      }
      /* Calculate moments (integrals of monomials; use direction=NULL)*/
      addTriangularElementMoments(corners,monoMoments,polyScale,
				  polyOrder,maxOrder1D,
				  direction);
    }
  }
  else {
    /* This is everything else - which is unknown=bad. 
     * Report error and exit                                           */
    fprintf(stderr,"%s: ERROR. Unknown type for element moments %s\n"
	    "    Please use \"source\" or \"eval\"\n"
	    "    Exiting\n",
	    fctName,elementType);
    _EXIT_;
  } /* End if elementType */

  /* Transfer monomial moments to moments array.
   * Presently, the monomials are basis functions, so no 
   * calculations are needed here, just copy the stuff to the 
   * "moments" array in the right order (as given by "?monoOrder")   */
  for (i=0;i<nterms;i++)
    moments[0][i]= 
      monoMoments[imonoOrder[i]][jmonoOrder[i]][kmonoOrder[i]];

  for (i=0;i<nterms;i++)
  /* Test output 
  if (direction!=NULL && strcmp(elementType,"source") == 0){
    printf("%s: Analysis for element #%d; polyScale=%16.8e\n",
	   fctName,element->directIndex,polyScale);
    printf("%s: Direction (%16.6e %16.6e %16.6e)\n",
	   fctName,direction[0],direction[1],direction[2]);
     printf("%s: Coll pnt (%16.6e %16.6e %16.6e)\n",
	   fctName,xc,yc,zc);
   for (i=0;i<nterms;i++)
      printf("%s: Moment[%d][%d][%d] =  %16.6e\n",
	     fctName,imonoOrder[i],jmonoOrder[i],kmonoOrder[i],
	     monoMoments[imonoOrder[i]][jmonoOrder[i]][kmonoOrder[i]]);
  }*/

  return;

} /* End of routine calcElementMoments */


/*
* ==================================================================== 
* Given an eval element pointer, return number of unknowns (i.e. number 
* of equations/collocation points associated with this eval).
* ==================================================================== */ 
int numEvalUnknowns(void *elementIn)
{
  struct constantElement *element;
  element = (struct constantElement *) elementIn;
  return (int) 1;
}

/*
* ==================================================================== 
* Given an eval element pointer, returns an array of indices.
* ==================================================================== */ 
int *evalIndices(void *elementIn)
{
  struct constantElement *element;
  element = (struct constantElement *) elementIn;
  return &(element->directIndex);
}


/*
* ==================================================================== 
* Given an source element pointer, return number of unknowns
* ==================================================================== */ 
int numSourceUnknowns(void *elementIn)
{
  struct constantElement *element;
  element = (struct constantElement *) elementIn;
  return (int) 1;
}

/*
* ==================================================================== 
* Given an source element pointer, returns an array of indices.
* ==================================================================== */ 
int *sourceIndices(void *elementIn)
{
  struct constantElement *element;
  element = (struct constantElement *) elementIn;
  return &(element->directIndex);
}


/*
* ==================================================================== 
* Calculate the values of the monomials of order [i][j][k] for 
* given terms, number of terms, maximum orders etc. 
* For the most commonly used cases the FOR loops are written out 
* in full to gain maximum speed.
* If the pointer "direction" is not pointing to NULL, then it must 
* be a direction vector, and the derivatives of the monomials in the 
* specified direction, rather than the monomials themselves, are 
* calculated and returned. Note: The direction vector is not being 
* normalized in the present routine, so the length of "direction" 
* will affect the output. A vector of unit length gives the "right" 
* result, i.e. the spatial derivatives in the direction of the vector.
* ==================================================================== */
static
void addMonomialValues(point point, 
		       double ***moments, 
		       double polyScale,
		       int polyOrder, int maxOrder1D,
		       double *direction)
{
  char fctName[]="addMonomialValues";
  double x1,y1,z1, x2,y2,z2, x3,y3,z3, 
    xpow,ypow,zpow, xpowx,ypowy,zpowz;
  int i,j,k, jSubtract, kSubtract, ifSubtract;

  /* For now this routine is set up to handle generic orders only. 
   * If this routine at some point takes up a lot of CPU time, then
   * the more sophisticated last part of the routine could be used 
   * instead.                                                          */
  if (polyOrder == maxOrder1D) ifSubtract=1; /* Stop calculations at 
					      * "consistent" order.    */
  else ifSubtract=0; /* Calculate all terms to given order in each dir.*/

  x1 = polyScale*point[0];
  y1 = polyScale*point[1];
  z1 = polyScale*point[2];

  if (direction==NULL) {
    /* Calculate monomial values */
    for (i=0, xpow=1.0; i<=maxOrder1D; i++){
      jSubtract = ifSubtract*i;
      for (j=0, ypow=1.0; j<=maxOrder1D-jSubtract; j++){
	kSubtract = ifSubtract*(i+j);
	for (k=0, zpow=1.0; k<=maxOrder1D-kSubtract; k++){
	  moments[i][j][k] = xpow * ypow * zpow;
	  zpow *= z1;
	}
	ypow *= y1;
      }
      xpow *= x1;
    }
  }
  else { /* direction!=NULL */
    /* Calculate spatial derivatives in direction given by "direction" */
    for (i=0, xpow=1.0, xpowx=0.0; i<=maxOrder1D; i++){
      jSubtract = ifSubtract*i;
      for (j=0, ypow=1.0, ypowy=0.0; j<=maxOrder1D-jSubtract; j++){
	kSubtract = ifSubtract*(i+j);
	for (k=0, zpow=1.0, zpowz=0.0; k<=maxOrder1D-kSubtract; k++){
	  moments[i][j][k] = 
	     direction[0]*xpowx * ypow * zpow
	    +direction[1]*ypowy * xpow * zpow
	    +direction[2]*zpowz * xpow * ypow;
	  /* New values and derivatives. Note how the scaling of the 
	   * variables affects the derivatives:                        */
	  zpowz = polyScale*(k+1)*zpow;
	  zpow *= z1;
	}
	ypowy = polyScale*(j+1)*ypow;
	ypow *= y1;
      }
      xpowx = polyScale*(i+1)*xpow;
      xpow *= x1;
    }
  } 
  return; /* <<<< NOTE RETURN HERE! (always used!) */

  /* The remaining part of this routine is a (tested) version that 
   * may be slightly faster than the above [for normal==NULL only. 
   * A version for the derivatives is *not* implemented at the moment]. 
   * However, the code is also 
   * less readable (to humans), and so is not used at the moment (the 
   * CPU time spent in this routine is, at present, negligible.
   * Timings suggest that the lower version is not a great deal 
   * faster... (no difference in the interpolation setup time could be 
   * seen)             */

  switch (100*polyOrder + maxOrder1D) {
  case 0:
    /* Constant needed only */
    moments[0][0][0] = 1.0     ;
    break;

  case 301:
    /* Max quadratic in each direction needed */
    x1 = polyScale*point[0];
    y1 = polyScale*point[1];
    z1 = polyScale*point[2];
    moments[0][0][0] = 1.0     ;
    moments[0][0][1] =       z1;
    moments[0][1][0] =    y1   ;
    moments[0][1][1] =    y1*z1;
    
    moments[1][0][0] = x1      ;
    moments[1][0][1] = x1*   z1;
    moments[1][1][0] = x1*y1   ;
    moments[1][1][1] = x1*y1*z1;
    break;

  case 602:
    /* Max quadratic in each direction needed */
    x1 = polyScale*point[0];
    y1 = polyScale*point[1];
    z1 = polyScale*point[2];
    x2 = x1*x1;
    y2 = y1*y1;
    z2 = z1*z1;
    moments[0][0][0] = 1.0     ;
    moments[0][0][1] =       z1;
    moments[0][0][2] =       z2;
    moments[0][1][0] =    y1   ;
    moments[0][1][1] =    y1*z1;
    moments[0][1][2] =    y1*z2;
    moments[0][2][0] =    y2   ;
    moments[0][2][1] =    y2*z1;
    moments[0][2][2] =    y2*z2;
    
    moments[1][0][0] = x1      ;
    moments[1][0][1] = x1*   z1;
    moments[1][0][2] = x1*   z2;
    moments[1][1][0] = x1*y1   ;
    moments[1][1][1] = x1*y1*z1;
    moments[1][1][2] = x1*y1*z2;
    moments[1][2][0] = x1*y2   ;
    moments[1][2][1] = x1*y2*z1;
    moments[1][2][2] = x1*y2*z2;

    moments[2][0][0] = x2      ;
    moments[2][0][1] = x2*   z1;
    moments[2][0][2] = x2*   z2;
    moments[2][1][0] = x2*y1   ;
    moments[2][1][1] = x2*y1*z1;
    moments[2][1][2] = x2*y1*z2;
    moments[2][2][0] = x2*y2   ;
    moments[2][2][1] = x2*y2*z1;
    moments[2][2][2] = x2*y2*z2;
    break;
   
  case 202:
    /* All monomials up to second order (10 of these) */
    x1 = polyScale*point[0];
    y1 = polyScale*point[1];
    z1 = polyScale*point[2];
    moments[0][0][0] = 1.0     ;
    moments[0][0][1] =       z1;
    moments[0][0][2] =       z1*z1;
    moments[0][1][0] =    y1   ;
    moments[0][1][1] =    y1*z1;
    moments[0][2][0] =    y1*y1   ;
    
    moments[1][0][0] = x1      ;
    moments[1][0][1] = x1*   z1;
    moments[1][1][0] = x1*y1   ;

    moments[2][0][0] = x1*x1      ;
    break;

  case 303:
    /* All monomials up to third order (20 of these)*/
    x1 = polyScale*point[0];
    y1 = polyScale*point[1];
    z1 = polyScale*point[2];
    x2 = x1*x1;
    y2 = y1*y1;
    z2 = z1*z1;
    moments[0][0][0] = 1.0     ;
    moments[0][0][1] =       z1;
    moments[0][0][2] =       z2;
    moments[0][0][3] =       z2*z1;
    moments[0][1][0] =    y1   ;
    moments[0][1][1] =    y1*z1;
    moments[0][1][2] =    y1*z2;
    moments[0][2][0] =    y2   ;
    moments[0][2][1] =    y2*z1;
    moments[0][3][0] =    y2*y1   ;

    moments[1][0][0] = x1      ;
    moments[1][0][1] = x1*   z1;
    moments[1][0][2] = x1*   z2;
    moments[1][1][0] = x1*y1   ;
    moments[1][1][1] = x1*y1*z1;
    moments[1][2][0] = x1*y2   ;

    moments[2][0][0] = x2      ;
    moments[2][0][1] = x2*   z1;
    moments[2][1][0] = x2*y1   ;

    moments[3][0][0] = x2*x1      ;
    break;

  case 404:
    /* All monomials up to third order (35 of these)*/
    x1 = polyScale*point[0];
    y1 = polyScale*point[1];
    z1 = polyScale*point[2];
    x2 = x1*x1;
    y2 = y1*y1;
    z2 = z1*z1;
    x3 = x1*x1*x1;
    y3 = y1*y1*y1;
    z3 = z1*z1*z1;
    moments[0][0][0] = 1.0     ;
    moments[0][0][1] =       z1;
    moments[0][0][2] =       z2;
    moments[0][0][3] =       z3;
    moments[0][0][4] =       z3*z1;
    moments[0][1][0] =    y1   ;
    moments[0][1][1] =    y1*z1;
    moments[0][1][2] =    y1*z2;
    moments[0][1][3] =    y1*z3;
    moments[0][2][0] =    y2   ;
    moments[0][2][1] =    y2*z1;
    moments[0][2][2] =    y2*z2;
    moments[0][3][0] =    y3   ;
    moments[0][3][1] =    y3*z1;
    moments[0][4][0] =    y3*y1   ;

    moments[1][0][0] = x1      ;
    moments[1][0][1] = x1*   z1;
    moments[1][0][2] = x1*   z2;
    moments[1][0][3] = x1*   z3;
    moments[1][1][0] = x1*y1   ;
    moments[1][1][1] = x1*y1*z1;
    moments[1][1][2] = x1*y1*z2;
    moments[1][2][0] = x1*y2   ;
    moments[1][2][1] = x1*y2*z1;
    moments[1][3][0] = x1*y3   ;

    moments[2][0][0] = x2      ;
    moments[2][0][1] = x2*   z1;
    moments[2][0][2] = x2*   z2;
    moments[2][1][0] = x2*y1   ;
    moments[2][1][1] = x2*y1*z1;
    moments[2][2][0] = x2*y2   ;

    moments[3][0][0] = x3      ;
    moments[3][0][1] = x3*   z1;
    moments[3][1][0] = x3*y1   ;
    
    moments[4][0][0] = x3*x1      ;
    break;

    /* Other hardcoding needed ? */
  default:
    /* Any orders could go here */
    /* General purpose - "Consistent-order" polynomial */
    if (polyOrder == maxOrder1D) {
      for (i=0; i<=polyOrder; i++){
	xpow = pow(polyScale*point[0],(double) i);
	for (j=0; j<=polyOrder-i; j++){
	  ypow = pow(polyScale*point[1],(double) j);
	  for (k=0; k<=polyOrder-i-j; k++){
	    moments[i][j][k] = 
	      xpow * ypow * pow(polyScale*point[2],(double) k);
	  }
	}
      }
    }
    /* General purpose go up to max order in each direction */
    else { 
      for (i=0; i<=maxOrder1D; i++){
	xpow = pow(polyScale*point[0],(double) i);
	for (j=0; j<=maxOrder1D; j++){
	  ypow = pow(polyScale*point[1],(double) j);
	  for (k=0; k<=maxOrder1D; k++){
	    moments[i][j][k] = 
	      xpow * ypow * pow(polyScale*point[2],(double) k);
	  }
	}
      }
    }

  } /* End switch */

  return;

  /* Test section. When testing is done, then put "return" statement 
   * above this line.
   * The test was completed successfully with maxdiff < 10^-16
   * for the following cases (see above): 0, 301, 602, 202, 303, 404
   * + default "1" ("case 505") and default "2" ("case 1505") 
   * each for several different elements as input. 
   * Note: Only tested like this for monomials (not derivatives)       */
  {
    int i,j,k;
    double diff, maxdiff;
    maxdiff=0.0;
    for (i=0;i<maxOrder1D;i++){
      for (j=0;j<maxOrder1D;j++){
	for (k=0;k<maxOrder1D;k++){
	  if (i+j+k<=polyOrder){
	    diff = fabs(moments[i][j][k]
			-pow(polyScale*point[0],(double) i)
			*pow(polyScale*point[1],(double) j)
			*pow(polyScale*point[2],(double) k));
	    assert(diff<1.e-16);
	    maxdiff=MAX(diff,maxdiff);
	  }
	}
      }
    }
    printf("%s: Test OK! order=%d, each=%d\n"
	   "Max diff = %16.6e\n",
	   fctName,polyOrder,maxOrder1D,maxdiff);
  }
  return;
} /* End of routine addMonomialValues */


/*
* ==================================================================== 
* Calculate the moments of a triangular flat element with regard to 
* the origin. If moments with regard to another point is needed, then
* simply subtract the coordinates of this point from each corner 
* coordinate before calling the present routine.
* This routine EITHER calculates the moments of order [i][j][k] 
* specified ranges of i, j and k, OR calculates the moments using 
* the spatial derivative of the above monomials in a specified 
* direction.
*
* corners: Element corners, corners[i][j] contains the value of the 
*          j'th coordinate (x,y,z) of corner #i (i=0,1,2).
*
* moments: On output moments[i][j][k] contains the moment M_ijk, 
*          i.e. (x-xzero)^i * (y-yzero)^j * (z-zzero)^k integrated 
*          over the element. moments must have size (at least) [3][3][3]
*          to contain the 27 moments on return.
*
* direction: 
* If the pointer "direction" is not pointing to NULL, then it must 
* be a direction vector, and the derivatives of the monomials in the 
* specified direction, rather than the monomials themselves, are 
* calculated and returned. Note: The direction vector is not being 
* normalized in the present routine, so the length of "direction" 
* will affect the output. A vector of unit length gives the "right" 
* result, i.e. the spatial derivatives in the direction of the vector.
*
* Note: This routine does not clean up memory after last call. It doesn't 
* use a lot of memory, but still not cleaning up might be considered 
* bad form.
* ==================================================================== */
static
void addTriangularElementMoments(point *corners, 
				 double ***moments,
				 double polyScale,
				 int polyOrder, 
				 int maxOrder1D,
				 double *direction)
{
  char fctName[]="addTriangularElementMoments";
  int i,j,k, jSubtract, kSubtract, ifSubtract, ic;
  double area_factor, X[3], Y[3], Z[3], x1,y1,z1,x2,y2,z2, dummy;
  double xpow, ypow, zpow, xpowx, ypowy, zpowz;
  /* For cubatures: */
  static double 
    **cubatureABC=NULL, /* Barycentric coordinates                     */
    *cubatureW=NULL;    /* Weights of cubature scheme                  */
  static int npointsCubature, degreeCubature, iruleCubature,
    lastPolyOrder=-1;

  /* Figure out what degree cubature rule is needed                    */
  /* Generally, the order of the cubature rule must be (at least) 
   * as large as the order of the polynomial being integrated  
   * (over the flat surface). If this is not what is set, then 
   * setup the cubature rule. 
   * If this routine is called altenately with different requirements, 
   * then the used approach here should be cahnged to store more than 
   * one cubature rule or only change the rule if the order is not high 
   * enough (i.e. don't change if the order is too high).              */
  if (polyOrder!=lastPolyOrder){ /* could be changed from "!="to ">"   */
    lastPolyOrder = polyOrder;
    /* Deallocate memory from any existing cubature scheme             */
    if (cubatureABC!=NULL){
      for (i=0;i<3;i++) FREE(cubatureABC[i]);
      FREE(cubatureABC);
      FREE(cubatureW);
    }
    /* The cubature rule should be (at least) of the same order as 
     * the polynomial to integrate                                     */
    degreeCubature = polyOrder;
    /* In general pick the first rule of a given order                 */
    iruleCubature  = 1;
    /* In some cases we may want to change the order due to the 
     * specific layout of the cubature rules                           */
    switch (polyOrder){
    case 0: /* Use first-order rule (1-point) */
      degreeCubature = 1;
      iruleCubature  = 1; 
      break;
      
    case 3:
      /* Rule # 1 has only four points, but it contains negative weights 
       * (could lead to large relative errors due to subtraction).
       * Rule # 2 has six points and no negative weights.
       * However, there is also a fourth-order rule with six points,
       * so why not use that?
       * Note, however, that some speed-up may be gained by choosing 
       * rule #1 of third order instead.                               */
      degreeCubature = 4;
      iruleCubature  = 1; 
      break;
      
    case 6: 
      /* The sixth-order rules has 12 points. So has the implemented 
       * seventh-order rule, so why not use that?                      */
      degreeCubature = 7;
      iruleCubature  = 1; 
      break;

      /* The cubature routine will complain if too high order is 
       * asked for. We do not need to do that here.                    */
    }

    /* Setup cubature points and weights.                              */
    printf("%s: Setting up cubature scheme of degree %d\n",
	   fctName,degreeCubature);
    npointsCubature=setupCubatureRuleTriangle(degreeCubature,
					      -iruleCubature,
					      NULL,NULL,NULL,NULL);
    /* Allocate memory */
    cubatureW  = (double *) calloc(npointsCubature,sizeof(double));
    cubatureABC = (double **) calloc(3,sizeof(double));
    for (i=0;i<3;i++){
      cubatureABC[i] = (double *) calloc(npointsCubature,sizeof(double));
    }
    /* Setup the points (in barycentric coordinates) and weights 
     * of the cubature scheme */
    i=setupCubatureRuleTriangle(degreeCubature,iruleCubature,
				cubatureABC[0],cubatureABC[1],
				cubatureABC[2],cubatureW);
    assert(i==npointsCubature);
  } /* The cubature points and weights now exist for a rule of 
     * sufficient order.                                               */

  /* The cubature rules are set up for a triangle with area 1/2. 
   * Thus, the Jacobian of the transfomation from this "unit" 
   * triangle the the present element is equal to TWICE the 
   * area of the element.                                              */
  VOP(X, corners[2], -, corners[0]);
  VOP(Y, corners[1], -, corners[0]); 
  XPROD(Z, X, Y);
  area_factor= VNORM(Z);
  /* This factor should be multiplied to EITHER all the cubature 
   * weights, OR to the resulting integrals. Which is better depends 
   * on the number of cubature points relative to the number of 
   * integrals. However, the first solution avoids a (triple) loop 
   * through the moments, so this is preferred, even though an extra 
   * array is needed to store the weights. */

  /* Calculate the moments using the curbature formula.
   * Note that the formula is based on the highest degree of the 
   * polynomial being integrated.  */

  if (polyOrder == maxOrder1D) ifSubtract=1; /* Stop calculations at 
					      * "consistent" order.    */
  else ifSubtract=0; /* Calculate all terms to given order in each dir.*/

  if (direction==NULL) {
    /* Integrate monomial values */
    for (ic=0;ic<npointsCubature;ic++){
      /* Setup Cartesian coordinates for this cubature node            */
      x1 = cubatureABC[0][ic]*corners[0][0]
          +cubatureABC[1][ic]*corners[1][0]
          +cubatureABC[2][ic]*corners[2][0];
      y1 = cubatureABC[0][ic]*corners[0][1]
          +cubatureABC[1][ic]*corners[1][1]
          +cubatureABC[2][ic]*corners[2][1];
      z1 = cubatureABC[0][ic]*corners[0][2]
          +cubatureABC[1][ic]*corners[1][2]
          +cubatureABC[2][ic]*corners[2][2];
      /* Scale the value coordinates to get scaled monomials           */
      x1 *= polyScale;
      y1 *= polyScale;
      z1 *= polyScale;
      
      /* Include area size and cubature weights with the x-power 
       * (to do this multiplication as few times as possible)          */
      for (i=0, xpow= area_factor*cubatureW[ic] ; i<=maxOrder1D; i++){
	jSubtract = ifSubtract*i;
	for (j=0, ypow=1.0; j<=maxOrder1D-jSubtract; j++){
	  kSubtract = ifSubtract*(i+j);
	  for (k=0, zpow=1.0; k<=maxOrder1D-kSubtract; k++){
	    moments[i][j][k] += xpow*ypow*zpow;
	    zpow *= z1;
	  }
	  ypow *= y1;
	}
	xpow *= x1;
      }
    } /* Next cubature point */
  } /* End if integrate monomials */
  else { /* direction!=0 */
    /* Integrate spatial derivatives in direction given by "direction" */
    for (ic=0;ic<npointsCubature;ic++){
      /* Setup Cartesian coordinates for this cubature node            */
      x1 = cubatureABC[0][ic]*corners[0][0]
          +cubatureABC[1][ic]*corners[1][0]
          +cubatureABC[2][ic]*corners[2][0];
      y1 = cubatureABC[0][ic]*corners[0][1]
          +cubatureABC[1][ic]*corners[1][1]
          +cubatureABC[2][ic]*corners[2][1];
      z1 = cubatureABC[0][ic]*corners[0][2]
          +cubatureABC[1][ic]*corners[1][2]
          +cubatureABC[2][ic]*corners[2][2];
      /* Scale the value coordinates to get scaled monomials           */
      x1 *= polyScale;
      y1 *= polyScale;
      z1 *= polyScale;

      /* It is possible to include area size and cubature weights with 
       * e.g. a temporary normal vector to decrease the number of 
       * multiplications made. 
       * [Do this only if the present routine getes expensive:
       *  the present formulation is simpler and less error prone]     */
      for (i=0, xpow=1.0, xpowx=0.0 ; i<=maxOrder1D; i++){
	jSubtract = ifSubtract*i;
	for (j=0, ypow=1.0, ypowy=0.0 ; j<=maxOrder1D-jSubtract; j++){
	  kSubtract = ifSubtract*(i+j);
	  for (k=0, zpow=1.0, zpowz=0.0 ; k<=maxOrder1D-kSubtract; k++){
	    moments[i][j][k] +=  
	      area_factor*cubatureW[ic]
	      *( direction[0]*xpowx * ypow * zpow
		+direction[1]*ypowy * xpow * zpow
		+direction[2]*zpowz * xpow * ypow);
	  /* New values and derivatives: */
	    zpowz = polyScale*(k+1)*zpow;
	    zpow *= z1;
	  }
	  ypowy = polyScale*(j+1)*ypow;
	  ypow *= y1;
	}
	xpowx = polyScale*(i+1)*xpow;
	xpow *= x1;
      }
    } /* Next cubature point */
  } /* End if integrate monomial derivatives */

  return; /* <<<< NOTE RETURN HERE! (always used!) */

  /* If a speed-up of this part of the computation is needed, then it 
   * may be considered to explicitly write out the computations made 
   * in each case (explicit unrolling the tripple loop) and using a 
   * SWITCH to choose between them. Alternatively, to avoid the SWITCH,
   * separate routines could be made for different orders.
   * See code below for a (tested) version of code that deals with 
   * the case maxOrder1D=2, polyOrder=6, directions=NULL. In this 
   * case using a switch here (based on input parameters) to select 
   * the hardcoded version rather  than the above results in an 
   * increase in speed of around 15% for  the projection part. 
   * Calling a  single routine specifically written for this case 
   * gives the same level of speed-up.                                 */
  /*switch(100*polyOrder + maxOrder1D)
   *...
   *case 602:*/
  for (ic=0;ic<npointsCubature;ic++){
    /* Setup Cartesian coordinates for this cubature node              */
    x1 = cubatureABC[0][ic]*corners[0][0]
        +cubatureABC[1][ic]*corners[1][0]
        +cubatureABC[2][ic]*corners[2][0];
    y1 = cubatureABC[0][ic]*corners[0][1]
        +cubatureABC[1][ic]*corners[1][1]
        +cubatureABC[2][ic]*corners[2][1];
    z1 = cubatureABC[0][ic]*corners[0][2]
        +cubatureABC[1][ic]*corners[1][2]
        +cubatureABC[2][ic]*corners[2][2];
    /* Scale the value coordinates to get scaled monomials             */
    x1 *= polyScale;
    y1 *= polyScale;
    z1 *= polyScale;
    /* Squared values of the coordinates */
    x2=x1*x1;
    y2=y1*y1;
    z2=z1*z1;

    /* Add to the values of the moments */
    dummy = area_factor*cubatureW[ic];
    moments[0][0][0] += dummy         ;
    moments[0][0][1] += dummy      *z1;
    moments[0][0][2] += dummy      *z2;
    moments[0][1][0] += dummy   *y1   ;
    moments[0][1][1] += dummy   *y1*z1;
    moments[0][1][2] += dummy   *y1*z2;
    moments[0][2][0] += dummy   *y2   ;
    moments[0][2][1] += dummy   *y2*z1;
    moments[0][2][2] += dummy   *y2*z2;

    dummy = area_factor*cubatureW[ic]*x1;
    moments[1][0][0] += dummy      ;
    moments[1][0][1] += dummy   *z1;
    moments[1][0][2] += dummy   *z2;
    moments[1][1][0] += dummy*y1   ;
    moments[1][1][1] += dummy*y1*z1;
    moments[1][1][2] += dummy*y1*z2;
    moments[1][2][0] += dummy*y2   ;
    moments[1][2][1] += dummy*y2*z1;
    moments[1][2][2] += dummy*y2*z2;

    dummy = area_factor*cubatureW[ic]*x2;
    moments[2][0][0] += dummy      ;
    moments[2][0][1] += dummy   *z1;
    moments[2][0][2] += dummy   *z2;
    moments[2][1][0] += dummy*y1   ;
    moments[2][1][1] += dummy*y1*z1;
    moments[2][1][2] += dummy*y1*z2;
    moments[2][2][0] += dummy*y2   ;
    moments[2][2][1] += dummy*y2*z1;
    moments[2][2][2] += dummy*y2*z2;
  }
  return;
} /* End of routine addTriangularElementMoments */


/*
* ==================================================================== 
* Add a constant to each matrix entry corresponding to a self-term.
* For constant type elements and "consistently" numbering this 
* is the diagonal entries of the matrix.
*
* The string 'type' is used to determine which entries (columns) to 
* modify. Allowed values are "all", "neumann" and "diriclet". The values 
* relate to the boundary condition type of the element corresponding to 
* each column.
*
* The idea behind the present routine is to be able to solve problems 
* using a direct formulation like
*  2pi phi = INT( G phi_n - G_n phi)
* In this case, the present routine may be used to add the "2pi" 
* term to the appropriate entries of a matrix if a integration routine 
* is used, which returns zero for the self integration of Gn.
* If a version of, e.g. calcp, is used which returns 2pi for the singular 
* integration of Gn, then the present routine should *not* be needed.
* HOWEVER, if you - after having obtained the solution - want to 
* evaluate the potential at some interior points, then you may want 
* to use the present routine to add 4*pi to the self term.
* ==================================================================== */
void addSelfTermConstant(void **sourcesIn, int nSources,
			 void **evalsIn,   int nEvals,
			 void *mat,        int matNum,
			 char *type,       double value)
{
  char fctName[]="addSelfTermConstant";  /* Name of this routine       */
  struct constantElement **sources, **evals;
  int isrc;

  /* Cast input parameters to correct type                             */
  sources = (struct constantElement **) sourcesIn;
  evals   = (struct constantElement **) evalsIn;

  /* At present, it is required that the sources and the evals 
   * "are the same". If they are not, then some method must be made 
   * which relates each source with it's corresponding eval.           */
  if (sources!=evals || nSources!=nEvals){
    fprintf(stderr,"%s: Sources and evals must be the same.\n"
	    "  If they are deliberately not so, then some "
	    "implementation is needed.\n",fctName);
    _EXIT_;
  }

  /* Check that the parameter 'type' is of a recognized format         */ 
  if (strcmp(type,"all")!=0 && 
      strcmp(type,"diriclet")!=0 && 
      strcmp(type,"neumann")!=0 ){
    fprintf(stderr,"%s: Parameter 'type' not recognized: \"%s\".\n"
	    "  Please use one of \"all\", \"diriclet\" and \"neumann\"\n",
	    fctName,type);
    _EXIT_;
  }

  /* Only add self term under given conditions: */
  if (strcmp(type,"all")==0){
    /* Add to ALL entries */
    for (isrc=0; isrc<nSources; isrc++){
      spMatEntryAddto(mat, matNum, 
		      sources[isrc]->directIndex, 
		      evals[isrc]->directIndex, 
		      value);
    }
  }
  else if (strcmp(type,"neumann")==0){
    /* Add self term to entries with Neumann boundary conditions */
    for (isrc=0; isrc<nSources; isrc++){
      if ((int)sources[isrc]->bcType==0)
	spMatEntryAddto(mat, matNum, 
			sources[isrc]->directIndex, 
			evals[isrc]->directIndex, 
			value);
    }
  }
  else {
    /* Add self term to entries with Diriclet boundary conditions */
    for (isrc=0; isrc<nSources; isrc++){
      if ((int)sources[isrc]->bcType==1)
	spMatEntryAddto(mat, matNum, 
			sources[isrc]->directIndex, 
			evals[isrc]->directIndex, 
			value);
    }
  } /* Next source */

  return;
} /* End of routine addSelfTermConstant */

/*
* ==================================================================== 
* To avoid keeping track of the sign of changes in the system matrices 
* as G or Gn are used on right-hand-side or left-hand-side in problems 
* with mixed boundary conditions, allow a change of sign in the values 
* based on the boundary condition type of the corresponding element.
*
* For example, a change of sign of all Neumann boundary condition values 
* combined with a change of sign of the final solution for evaluation 
* points on the Neumann boundaries should work. 
* This approach is only needed for problems with mixed boundary 
* conditions (some Diriclet and some Neumann), but I guess it won't do 
* any harm to swap signs on the Neumann conditions and corresponding 
* solutions always (it will take a little CPU time, but the result 
* should be OK.
* ==================================================================== */
void swapSignVectorElement(void *elementsIn, int nElements,
		     double *vector, char *type)
{
  char fctName[]="swapSignVectorElement";  /* Name of this routine     */
  struct constantElement **elements;
  int iele;

  /* Cast input parameters to correct type                             */
  elements = (struct constantElement **) elementsIn;
  
  /* Check that the parameter 'type' is of a recognized format         */ 
  if (strcmp(type,"all")!=0 && 
      strcmp(type,"diriclet")!=0 && 
      strcmp(type,"neumann")!=0 ){
    fprintf(stderr,"%s: Parameter 'type' not recognized: \"%s\".\n"
	    "  Please use one of \"all\", \"diriclet\" and \"neumann\"\n",
	    fctName,type);
    _EXIT_;
  }

  /* Only add self term under given conditions: */
  if (strcmp(type,"all")==0){
    /* Change sign on ALL entries */
    for (iele=0; iele<nElements; iele++){
      vector[elements[iele]->directIndex] *= -1.0;
    }
  }
  else if (strcmp(type,"neumann")==0){
    /* Change sign on entries with Neumann boundary conditions */
    for (iele=0; iele<nElements; iele++){
      if ((int) elements[iele]->bcType==0)
	vector[elements[iele]->directIndex] *= -1.0;
    }
  }
  else {
    /* Change sign on entries with Diriclet boundary conditions */
    for (iele=0; iele<nElements; iele++){
      if ((int) elements[iele]->bcType==1)
	vector[elements[iele]->directIndex] *= -1.0;
    }
  } /* Next element */

  return;
} /* End of routine swapSignVectorElement */

/* 
* ====================================================================
* Integrate a solution over the boundary.
* Integration may be over:
* 1 "whole" The entire boundary without discremination.
* 2 "all parts" The entire boundary keeping different parts of the 
*       boundary (e.g. conducters) separated in "integrationOnPart".
*       Total result is also calculated. 
*       Note that each part will be ADDED TO whatever is in the 
*       integrationOnPart vector.
* 3 "part" Integrate only over boundary part ipart. Result is returned 
*       as the value of the routine.
* ==================================================================== */ 
double integrateOverElements(void *sourcesIn, int nSources, 
			     void *evalsIn, int nEvals, 
			     double *solution,
			     char *integrationType,
			     double *integrationOnPart,
			     int ipart){
  char fctName[] = "integrateOverElements";
  struct constantElement **sources;
  int iele;
  double value, area, normal[3];

  /* If the sources and the evals are not numbered consistently, then 
   * we have no way of knowing (?) which area piece (source) should go 
   * with each solution value (eval). */
  assert(sourcesIn == evalsIn);
  assert(nSources  == nEvals);

  /* Cast input data to correct type */
  sources = (struct constantElement **) sourcesIn;

  if (strcmp(integrationType,"whole")==0){
    /* Integrate over entire boundary only */
    for (value=0,iele=0; iele<nSources; iele++){
      /* Get element area */
      area = elementNormal(sources[iele],normal);
      /* Add area*solution to the total value */
      value += area*solution[sources[iele]->directIndex];
    }
  }
  else if(strcmp(integrationType,"all parts")==0){
    /* Integrate over entire boundary keeping parts separate */
    for (value=0,iele=0; iele<nSources; iele++){
      /* Get element area */
      area = elementNormal(sources[iele],normal);
      /* Add area*solution to the local value */
      integrationOnPart[(int)sources[iele]->boundaryNumber] 
	+= area*solution[sources[iele]->directIndex];
      /* Add area*solution to the total value */
      value += area*solution[sources[iele]->directIndex];
    }
  }
  else if(strcmp(integrationType,"part")==0){
    /* Integrate over a single part of the boundary only */
    for (value=0,iele=0; iele<nSources; iele++){
      if ((int)sources[iele]->boundaryNumber == ipart){
	/* Get element area */
	area = elementNormal(sources[iele],normal);
	/* Add area*solution to the solution */
	value += area*solution[sources[iele]->directIndex];
      }
    }
  }
  else { /* This is unknown = bad. Complain and exit! */
    fprintf(stderr,"%s: Unknown integration type: \"%s\"\n",
	    fctName,integrationType);
    _EXIT_;
  }

  return value;
} /* End of routine integrateOverElements */



#include "elements.testfunctions.c"


