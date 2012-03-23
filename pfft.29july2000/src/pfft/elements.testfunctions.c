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
/* ==== TEST FUNCTIONS (not really needed) ==== */ 
#include <time.h>


/* 
* ====================================================================
* Temporary test routine. Test a solution by examining the 
* calculated values with values obtained directly from the 
* analytic solution (boundaryConditionValue).
* ==================================================================== */ 
void testSolution(void *evalsIn, int nEvals, double *solution,
		  char *bcType, int bcIdata){
  struct constantElement **evals;
  int iele, id, ibcType, ie2[2], itype;
  double err[2], err2[2], locErr, val, normal[3], intErr[2],
    area, fact, fact2, testarea[2], testcharge[2];
  char *types[2], typep[]="phi", typepn[]="phi_n";
  FILE *outFile;
  /* Cast input to correct type */
  evals = (struct constantElement **) evalsIn;

  /* Initialize arrays for storing error terms: */
  err[0]    = 0.0;
  err[1]    = 0.0;
  err2[0]   = 0.0;
  err2[1]   = 0.0;
  intErr[0] = 0.0;
  intErr[1] = 0.0;
  testarea[0] = 0.0;
  testarea[1] = 0.0;
  testcharge[0] = 0.0;
  testcharge[1] = 0.0;
  ie2[0]    = 0;
  ie2[1]    = 0; 
  /* Examine error on each element. Store 2-norm and inf-norm */
  for (iele=0; iele<nEvals; iele++){
    id = evals[iele]->directIndex;
    /* Get the oposite type as what the boundary condition type 
     * indicates, i.e. if phi_n is known as a boundary condition, 
     * then get phi. */
    if ((int) (evals[iele]->bcType)==0) ibcType = 1;
    else ibcType = 0;
    /* Normal direction for this element */
    area = elementNormal(evals[iele], normal);
    /* Get exact value for this element */
    val = boundaryConditionValue(ibcType, 
				 evals[iele]->centroid,
				 normal,
				 (int) evals[iele]->boundaryNumber,
				 bcType,
				 bcIdata);
    /* Compare to calculated value, and store error terms */
    locErr = ABS(solution[id]-val);

    /* Inf-norm error: */
    err[ibcType]  =  MAX(err[ibcType],locErr);
    /* RMS-error */
    err2[ibcType] +=  locErr*locErr;
    /* Integrated error (1-norm error weighted by element area) */
    intErr[ibcType] +=  locErr*area;
    testarea[ibcType] += area;
    if (evals[iele]->centroid[2]>0.0) 
      testcharge[ibcType] += area*solution[id];

    ie2[ibcType]++;
  }
  /* Calculate RMS error: */

  for (itype=0; itype<2; itype++){
    if (ie2[itype]>0){
      err2[itype] /= (double) (ie2[itype]);
      err2[itype] =  sqrt(err2[itype]);
    }
  }
  /* Hardcoding for normalization. */
  fact = (8.0*PI)/3.0;
  printf(" Normalizing inf-norm and rms errors with %g\n",
	 1./fact);
  fact2 = 8.0/6.0;
    printf(" Normalizing int(|err|) by Q=int(|q|)=2*3/8=%g\n",
	   1./fact2);
  for (itype=0; itype<2; itype++){
    err[itype]    *= fact;
    err2[itype]   *= fact;
    intErr[itype] *= fact2;
  }

  types[0] = typepn;
  types[1] = typep;
  /* Write calculated errors to screen: */
  for (itype=0; itype<2; itype++){
    if(ie2[itype] == 0) continue; /* Next type (this type not encountered) */
    printf(" Errors on caluclated %s, (%d elements):\n",
	   types[itype],nEvals);
    printf("  Int|error|=%11.6e, Max error=%11.6e, RMS error=%11.6e\n",
	   intErr[itype],err[itype], err2[itype]);
    printf("  Area of %s-part integrated to %g\n",
	   types[itype],testarea[itype]);
    printf("  Charge on half sphere (%s-part) integrated to %g\n",
	   types[itype],testcharge[itype]);
  }

  /* See if output file exists */
  outFile = fopen("errors.dat","r");
  if (outFile==NULL){
    printf("Writing header to error file \"errors.dat\"\n");
    outFile = fopen("errors.dat","w");
    fprintf(outFile,"#elements   int|errors|     inf|error|"
	    "      RMS-error\n");
  }
  fclose(outFile);

  /* Write errors to file: */
  if (ie2[1]>0) itype=1; /* Examine calculated PHI data only    */
  else itype = 0;        /* Examine calculated PHI_n data only  */
  printf(" Writing errors on %s to \"errors.dat\"\n",types[itype]);
  outFile = fopen("errors.dat","a"); /* Append to existing data */
 
  fprintf(outFile,"%7d  %12.6e  %12.6e  %12.6e\n",
	  nEvals,intErr[ibcType],err[itype],err2[itype]);
  fclose(outFile);
  

  return;
}


/* Dump elements to tecplot for print Since we do not know 
* how the panel neighbours are arranged, four distinct points 
* are used for each panel (if the panel is a triangle, then 
* one of the points will be written twice).  */
void dumpElementsToTecplot(void *elementsIn, int nElements,
			   int nParts, int *parts){
  int iele,i, iplot, ip, nWrote, writeMe;
  struct constantElement **elements, *element;
  FILE *tecplotOutputFile;
  /* Cast input date to correct type */
  elements = (struct constantElement **) elementsIn;

  /* Open output file */
  tecplotOutputFile = fopen("panels_tp7.out","w");

  /* Write tecplot header */
  fprintf(tecplotOutputFile,"VARIABLES = \"X\", \"Y\", \"Z\"\n");
  fprintf(tecplotOutputFile,
	  "ZONE  N = %d  E = %d, F=FEPOINT, ET=QUADRILATERAL\n",
	  4*nElements, nElements);

  /* Write point data */
  nWrote = 0;
  for (iele=0; iele<nElements; iele++){
    writeMe=0;

    element = elements[iele];
    /* Only write if we want this element out: */
    if (nParts==0){ /* Write all elements */
      writeMe=1;
    }
    else {
      for (i=0; i<nParts; i++) 
	if (parts[nParts]==(int)element->boundaryNumber)
	  writeMe=1;
    }
    if (writeMe){
      for (i=0; i<4; i++){
	if (element->shape==3 && i==3) iplot = 2;
	else iplot = i;
	fprintf(tecplotOutputFile,
		"%12.6e %12.6e %12.6e\n",
		element->corners[iplot][0],
		element->corners[iplot][1],
		element->corners[iplot][2]);
	nWrote++;
      }
    }
  }

  /* Write elements as combinations of points */
  for (iele=0, ip=0; iele<nWrote; iele++, ip+=4){
    fprintf(tecplotOutputFile,
	      "%d %d %d %d\n",
	      ip+1,ip+2,ip+3,ip+4);
  }


  /* Close output file */
  fclose(tecplotOutputFile);
}



/* Test the calculation of moments over flat panels */
void testMoments(void)
{
  double xp[3], yp[3], zp[3];
  point corners[3], xzero;
  double ***moments1,***moments2,***moments3,*moments3B, maxerr, t1,t2;
  int i,j,k, many=1;

  /* NOTE: On the final timing on 1000000=10^6 triangles the 
   * new implementation used 7.36 secs while the old implementation 
   * used 355.39 secs. This corresponds to a speedup of a factor 
   * of 48 (and that ain't bad!). */

printf("TESTMOMENTS:\n");

  moments1 = (double ***) calloc(3,sizeof(double **));
  moments2 = (double ***) calloc(3,sizeof(double **));
  for (i=0;i<3;i++){
    moments1[i] = (double **) calloc(3,sizeof(double *));
    moments2[i] = (double **) calloc(3,sizeof(double *));
    for (j=0;j<3;j++){
      moments1[i][j] = (double *) calloc(3,sizeof(double));
      moments2[i][j] = (double *) calloc(3,sizeof(double));
    }
  }

  moments3B = (double *)   calloc(3*3*3,sizeof(double));
  moments3  = (double ***) calloc(3,sizeof(double **));
  for (i=0;i<3;i++){
    moments3[i] = (double **) calloc(3,sizeof(double *));
    for (j=0;j<3;j++){
      moments3[i][j] = &(moments3B[i*3*3+j*3]);
    }
  }

  /* Set up panel coordinates */
  corners[0][0]=0.0;
  corners[0][1]=0.0;
  corners[0][2]=-1.0;
  corners[1][0]=1.0;
  corners[1][1]=0.0;
  corners[1][2]=0.0;
  corners[2][0]=2.0;
  corners[2][1]=2.0;
  corners[2][2]=3.0;
  for (i=0;i<3;i++){
    xp[i] = corners[i][0];
    yp[i] = corners[i][1];
    zp[i] = corners[i][2];
  }

  xzero[0]=0.0;
  xzero[1]=0.0;
  xzero[2]=0.0;
  printf("Calculating a lot of *new* moments \n");
  t1= (double) clock()/CLOCKS_PER_SEC;
  for (i=0;i<many;i++){
    /* Zero out moments */
    for (j=0;j<3*3*3;j++) moments3B[j]=0.0;
    /*addTriangularElementMomentsOrder2(corners,moments3,(double) 1.0);*/
  }
  t2= (double) clock()/CLOCKS_PER_SEC;
  printf(" New method used %10.2f seconds\n",t2-t1);
  /*add one more set of moments */
  /*addTriangularElementMomentsOrder2(corners,moments3,(double) 1.0);*/

  printf("Calculating a lot of *old* moments \n");
  t1= (double) clock()/CLOCKS_PER_SEC;
  for (i=0;i<many;i++){
    ComputeMomentsGlobal(xp, yp, zp, 2, moments2);
  }
  t2= (double) clock()/CLOCKS_PER_SEC;
  printf(" Old method used %10.2f seconds\n",t2-t1);

  for (i=0,maxerr=0.0;i<3;i++){
    for (j=0;j<3;j++){
      for (k=0;k<3;k++){
	/* printf(" M_{%d,%d,%d} = %12.8f %12.8f   %12.6e\n",i,j,k,
	       moments1[i][j][k],moments2[i][j][k],
	       fabs(moments1[i][j][k]-moments2[i][j][k]) );*/
	maxerr=MAX(fabs(moments3[i][j][k]-2.0*moments2[i][j][k])/fabs(moments3[i][j][k]),maxerr);
      }
    }
  }
  printf(" Maximum difference on moments is %12.6e\n",maxerr);
  
  return;
}

void testLocalMoments(void)
{
  int i,j,ii,m=2,n=2, order;
  double xp[3], yp[3], zp[3];
  double **momentsLoc, ***momentsGlob;

  if (0) {
    momentsLoc = (double **) calloc(m+1,sizeof(double *));
    for (i=0;i<=m;i++){
      momentsLoc[i]= (double *) calloc(n+1,sizeof(double));
    }
    for (i=0;i<1000000;i++){
      localElementMoments(m,n,momentsLoc);
    }
    printf("Local moments:\n");
    for (i=0;i<=MIN(10,m);i++){
      for (j=0;j<=MIN(6,n);j++){
	printf(" %10.4e",momentsLoc[i][j]);
      }
      printf(" \n");
    }
  }
  else if (1) {
    order=2;
    xp[0]=0.0; xp[1]=0.0; xp[2]=1.0;
    yp[0]=0.0; yp[1]=1.0; yp[2]=0.0;
    zp[0]=0.0; zp[1]=0.0; zp[2]=0.0;
    for (i=0;i<100000;i++){
      /*      momentsGlob = ComputeMomentsGlobal(xp, yp, zp, order);*/
      /* FREE memory */
      for (ii=0; ii<=order; ii++) {
	for (j=0; j<=order; j++) {
	  FREE(momentsGlob[ii][j]);
	}
	FREE(momentsGlob[ii]);
      }
      FREE(momentsGlob);
    }
    /*    momentsGlob = ComputeMomentsGlobal(xp, yp, zp, order);
    for (i=0;i<=order;i++){
      for (j=0;j<=order;j++){
	for (k=0;k<=order;k++){
	  printf("(i,j,k)=(%d,%d,%d), moment=%10.4e\n",i,j,k,momentsGlob[i][j][k]);
	}
      }
    }
    */
    
  }
    
}

#define CALLOC(PNTR, NUM, TYPE, FLAG, MTYP) \
         (PNTR)=(TYPE*)calloc(NUM, sizeof(TYPE));

/* WENJING's version. This needs modifications! 
 * Eventually give this routine a new name. 
 * The original implementation leaked a lot, about 97Mb 
 * wasted when 100k triangular panels were used (order=2).
 * Double that number for quad panels, i.e. roughly 1-2kb
 * per panel are wasted here. 
 * Roughly half of the leaking is due to the repeated allocations 
 * in what used to be "ComputeMomentsLocal" the other half 
 * is due to  the allocations in "ComputeMomentsGlobal". 
 * None of the temporary arrays are freed, and they are always 
 * reallocated on call. The first part may be acceptable, but 
 * the last part certainly isn't!
 * This should probably be taken care of in the Stokes code as well. 
 */
static double ***ComputeMomentsGlobal(double *xp, double *yp, double *zp, 
				      int order,double ***Iin)
{
  int i, j, k, m, n, ii, jj, kk, mm, nn, l, iie, jje, kke;
  double g1, g2, g3, J;
  double TC[3][3];
  double ***I, **cml;
  static double **momentsLoc=NULL, **C=NULL;
  static int mloc=-1,nloc=-1,orderC=-1;

  /* Allocate C - only if not already large enough */
  if (order>orderC){
    if (C!=NULL){
      /* C should be freed here. */
    }
    orderC=order;
    CALLOC(C, 3, double*, ON, AQ2P);
    for (i=0; i <= 2; i++)
      CALLOC(C[i], (order+1)*(order+2)/2, double, ON, AQ2P);
  }
  
  /* Allocate moments matrix, Iijk. 
  CALLOC(I, order+1, double**, ON, AQ2P);
  for (i=0; i<=order; i++) {
    CALLOC(I[i], order+1, double*, ON, AQ2P);
    for (j=0; j<=order; j++) CALLOC(I[i][j], order+1, double, ON, AQ2P);
  }*/
  I=Iin;

  /* Initialize the moments matrix. */
  for (i = 0; i <= order; i++) {
    for (j = 0; j <= order; j++) {
      for (k = 0; k <= order; k++) {
        I[i][j][k] = 0.0;
      }
    }
  }

  /* Allocate temporary storage and initialize arrays. */
  for (i = 0; i <= 2; i++) {
    for (j = 0; j <= 2; j++) {
      TC[i][j] = 0.0;
    }
  }

  /* Compute the Jacobian of the panel. */
  g1 = (yp[0] - yp[2])*(zp[1] - zp[2]) - (zp[0] - zp[2])*(yp[1] - yp[2]);
  g2 = (zp[0] - zp[2])*(xp[1] - xp[2]) - (xp[0] - xp[2])*(zp[1] - zp[2]);
  g3 = (xp[0] - xp[2])*(yp[1] - yp[2]) - (yp[0] - yp[2])*(xp[1] - xp[2]);
  J = sqrt(g1 * g1 + g2 * g2 + g3 * g3);

  /* Compute coefficiences of transformation of (x,y,z) to the local system */
  TC[0][0] = xp[2];
  TC[0][1] = xp[0] - xp[2];
  TC[0][2] = xp[1] - xp[2];
  TC[1][0] = yp[2];
  TC[1][1] = yp[0] - yp[2];
  TC[1][2] = yp[1] - yp[2];
  TC[2][0] = zp[2];
  TC[2][1] = zp[0] - zp[2];
  TC[2][2] = zp[1] - zp[2];

  /* Get local moments. */
  if (mloc<3*order || nloc<3*order){
    if (momentsLoc!=NULL) FREE(momentsLoc);
    mloc=nloc=3*order;
    momentsLoc = (double **) calloc(mloc+1,sizeof(double *));
    for (i=0;i<=mloc;i++){
      momentsLoc[i]= (double *) calloc(nloc+1,sizeof(double));
    }
    localElementMoments(mloc, nloc,momentsLoc);
  }
  cml = momentsLoc;

  /* Compute expanding coefficiences of each components (x^i, y^j, z^k) */
  for (i = 0; i <= order; i++) {
    for (j = 0; j <= order; j++) {
      for (k = 0; k <= order; k++) {
        for (mm = 0; mm <= 2; mm++) {
          for (nn = 0; nn <= 2; nn++) C[mm][nn] = TC[mm][nn];
        }
        if (i == 0) {
          C[0][0] = 1.0;
          for (l=1; l<=2; l++) C[0][l] = 0.0;
        }
        if (j == 0) {
          C[1][0] = 1.0;
          for (l=1; l<=2; l++) C[1][l] = 0.0;
        }
        if (k == 0) {
          C[2][0] = 1.0;
          for (l=1; l<=2; l++) C[2][l] = 0.0;
        }
        if (i == 2) {
          C[0][0] = TC[0][0] * TC[0][0];
          C[0][1] = 2.0 * TC[0][0] * TC[0][1];
          C[0][2] = 2.0 * TC[0][0] * TC[0][2];
          C[0][3] = TC[0][1] * TC[0][1];
          C[0][4] = 2.0 * TC[0][1] * TC[0][2];
          C[0][5] = TC[0][2] * TC[0][2];
        }
        if (j == 2) {
          C[1][0] = TC[1][0] * TC[1][0];
          C[1][1] = 2.0 * TC[1][0] * TC[1][1];
          C[1][2] = 2.0 * TC[1][0] * TC[1][2];
          C[1][3] = TC[1][1] * TC[1][1];
          C[1][4] = 2.0 * TC[1][1] * TC[1][2];
          C[1][5] = TC[1][2] * TC[1][2];
        }
        if (k == 2) {
          C[2][0] = TC[2][0] * TC[2][0];
          C[2][1] = 2.0 * TC[2][0] * TC[2][1];
          C[2][2] = 2.0 * TC[2][0] * TC[2][2];
          C[2][3] = TC[2][1] * TC[2][1];
          C[2][4] = 2.0 * TC[2][1] * TC[2][2];
          C[2][5] = TC[2][2] * TC[2][2];
        }
        if (i == 3) {
          C[0][0] = TC[0][0] * TC[0][0] * TC[0][0];
          C[0][1] = 3.0 * TC[0][0] * TC[0][0] * TC[0][1];
          C[0][2] = 3.0 * TC[0][0] * TC[0][0] * TC[0][2];
          C[0][3] = 3.0 * TC[0][0] * TC[0][1] * TC[0][1];
          C[0][4] = 6.0 * TC[0][0] * TC[0][1] * TC[0][2];
          C[0][5] = 3.0 * TC[0][0] * TC[0][2] * TC[0][2];
          C[0][6] = TC[0][1] * TC[0][1] * TC[0][1];
          C[0][7] = 3.0 * TC[0][1] * TC[0][1] * TC[0][2];
          C[0][8] = 3.0 * TC[0][1] * TC[0][2] * TC[0][2];
          C[0][9] = TC[0][2] * TC[0][2] * TC[0][2];
        }
        if (j == 3) {
          C[1][0] = TC[1][0] * TC[1][0] * TC[1][0];
          C[1][1] = 3.0 * TC[1][0] * TC[1][0] * TC[1][1];
          C[1][2] = 3.0 * TC[1][0] * TC[1][0] * TC[1][2];
          C[1][3] = 3.0 * TC[1][0] * TC[1][1] * TC[1][1];
          C[1][4] = 6.0 * TC[1][0] * TC[1][1] * TC[1][2];
          C[1][5] = 3.0 * TC[1][0] * TC[1][2] * TC[1][2];
          C[1][6] = TC[0][1] * TC[0][1] * TC[0][1];
          C[1][7] = 3.0 * TC[1][1] * TC[1][1] * TC[1][2];
          C[1][8] = 3.0 * TC[1][1] * TC[1][2] * TC[1][2];
          C[1][9] = TC[1][2] * TC[1][2] * TC[1][2];
        }
        if (k == 3) {
          C[2][0] = TC[2][0] * TC[2][0] * TC[2][0];
          C[2][1] = 3.0 * TC[2][0] * TC[2][0] * TC[2][1];
          C[2][2] = 3.0 * TC[2][0] * TC[2][0] * TC[2][2];
          C[2][3] = 3.0 * TC[2][0] * TC[2][1] * TC[2][1];
          C[2][4] = 6.0 * TC[2][0] * TC[2][1] * TC[2][2];
          C[2][5] = 3.0 * TC[2][0] * TC[2][2] * TC[2][2];
          C[2][6] = TC[2][1] * TC[2][1] * TC[2][1];
          C[2][7] = 3.0 * TC[2][1] * TC[2][1] * TC[2][2];
          C[2][8] = 3.0 * TC[2][1] * TC[2][2] * TC[2][2];
          C[2][9] = TC[2][2] * TC[2][2] * TC[2][2];
        }
        if (i == 4) {
          C[0][0] = TC[0][0] * TC[0][0] * TC[0][0] * TC[0][0];
          C[0][1] = 4.0 * TC[0][0] * TC[0][0] * TC[0][0] * TC[0][1];
          C[0][2] = 4.0 * TC[0][0] * TC[0][0] * TC[0][0] * TC[0][2];
          C[0][3] = 6.0 * TC[0][0] * TC[0][0] * TC[0][1] * TC[0][1];
          C[0][4] = 12.0 * TC[0][0] * TC[0][0] * TC[0][1] * TC[0][2];
          C[0][5] = 6.0 * TC[0][0] * TC[0][0] * TC[0][2] * TC[0][2];
          C[0][6] = 4.0 * TC[0][0] * TC[0][1] * TC[0][1] * TC[0][1];
          C[0][7] = 12.0 * TC[0][0] * TC[0][1] * TC[0][1] * TC[0][2];
          C[0][8] = 12.0 * TC[0][0] * TC[0][1] * TC[0][2] * TC[0][2];
          C[0][9] = 4.0 * TC[0][0] * TC[0][2] * TC[0][2] * TC[0][2];
          C[0][10] = TC[0][1] * TC[0][1] * TC[0][1] * TC[0][1];
          C[0][11] = 4.0 * TC[0][1] * TC[0][1] * TC[0][1] * TC[0][2];
          C[0][12] = 6.0 * TC[0][1] * TC[0][1] * TC[0][2] * TC[0][2];
          C[0][13] = 4.0 * TC[0][1] * TC[0][2] * TC[0][2] * TC[0][2];
          C[0][14] = TC[0][2] * TC[0][2] * TC[0][2] * TC[0][2];
        }
        if (j == 4) {
          C[1][0] = TC[1][0] * TC[1][0] * TC[1][0] * TC[1][0];
          C[1][1] = 4.0 * TC[1][0] * TC[1][0] * TC[1][0] * TC[1][1];
          C[1][2] = 4.0 * TC[1][0] * TC[1][0] * TC[1][0] * TC[1][2];
          C[1][3] = 6.0 * TC[1][0] * TC[1][0] * TC[1][1] * TC[1][1];
          C[1][4] = 12.0 * TC[1][0] * TC[1][0] * TC[1][1] * TC[1][2];
          C[1][5] = 6.0 * TC[1][0] * TC[1][0] * TC[1][2] * TC[1][2];
          C[1][6] = 4.0 * TC[1][0] * TC[1][1] * TC[1][1] * TC[1][1];
          C[1][7] = 12.0 * TC[1][0] * TC[1][1] * TC[1][1] * TC[1][2];
          C[1][8] = 12.0 * TC[1][0] * TC[1][1] * TC[1][2] * TC[1][2];
          C[1][9] = 4.0 * TC[1][0] * TC[1][2] * TC[1][2] * TC[1][2];
          C[1][10] = TC[1][1] * TC[1][1] * TC[1][1] * TC[1][1];
          C[1][11] = 4.0 * TC[1][1] * TC[1][1] * TC[1][1] * TC[1][2];
          C[1][12] = 6.0 * TC[1][1] * TC[1][1] * TC[1][2] * TC[1][2];
          C[1][13] = 4.0 * TC[1][1] * TC[1][2] * TC[1][2] * TC[1][2];
          C[1][14] = TC[1][2] * TC[1][2] * TC[1][2] * TC[1][2];
        }
        if (k == 4) {
          C[2][0] = TC[2][0] * TC[2][0] * TC[2][0] * TC[2][0];
          C[2][1] = 4.0 * TC[2][0] * TC[2][0] * TC[2][0] * TC[2][1];
          C[0][2] = 4.0 * TC[2][0] * TC[2][0] * TC[2][0] * TC[2][2];
          C[2][3] = 6.0 * TC[2][0] * TC[2][0] * TC[2][1] * TC[2][1];
          C[2][4] = 12.0 * TC[2][0] * TC[2][0] * TC[2][1] * TC[2][2];
          C[2][5] = 6.0 * TC[2][0] * TC[2][0] * TC[2][2] * TC[2][2];
          C[2][6] = 4.0 * TC[2][0] * TC[2][1] * TC[2][1] * TC[2][1];
          C[2][7] = 12.0 * TC[2][0] * TC[2][1] * TC[2][1] * TC[2][2];
          C[2][8] = 12.0 * TC[2][0] * TC[2][1] * TC[2][2] * TC[2][2];
          C[2][9] = 4.0 * TC[2][0] * TC[2][2] * TC[2][2] * TC[2][2];
          C[2][10] = TC[2][1] * TC[2][1] * TC[2][1] * TC[2][1];
          C[2][11] = 4.0 * TC[2][1] * TC[2][1] * TC[2][1] * TC[2][2];
          C[2][12] = 6.0 * TC[2][1] * TC[2][1] * TC[2][2] * TC[2][2];
          C[2][13] = 4.0 * TC[2][1] * TC[2][2] * TC[2][2] * TC[2][2];
          C[2][14] = TC[2][2] * TC[2][2] * TC[2][2] * TC[2][2];
        }
        iie = (i+1)*(i+2)/2;
        jje = (j+1)*(j+2)/2;
        kke = (k+1)*(k+2)/2;
        for (ii = 0; ii < iie; ii++) {
          for (jj = 0; jj < jje; jj++) {
            for (kk = 0; kk < kke; kk++) {
              m = 0;
              n = 0;
              if (ii == 1) m = m + 1;
              if (jj == 1) m = m + 1;
              if (kk == 1) m = m + 1;
              if (ii == 2) n = n + 1;
              if (jj == 2) n = n + 1;
              if (kk == 2) n = n + 1;
              if (ii == 3) m = m + 2;
              if (jj == 3) m = m + 2;
              if (kk == 3) m = m + 2;
              if (ii == 4) {
                m = m + 1;
                n = n + 1;
	      }
              if (jj == 4) {
                m = m + 1;
                n = n + 1;
	      }
              if (kk == 4) {
                m = m + 1;
                n = n + 1;
	      }
              if (ii == 5) n = n + 2;
              if (jj == 5) n = n + 2;
              if (kk == 5) n = n + 2;
              if (ii == 6) m = m + 3;
              if (jj == 6) m = m + 3;
              if (kk == 6) m = m + 3;
              if (ii == 7) {
                m = m + 2;
                n = n + 1;
	      }
              if (jj == 7) {
                m = m + 2;
                n = n + 1;
	      } 
              if (kk == 7) {           
                m = m + 2;
                n = n + 1;
	      }
              if (ii == 8) {
                m = m + 1;
                n = n + 2;
	      }
              if (jj == 8) {
                m = m + 1;
                n = n + 2;
	      } 
              if (kk == 8) {           
                m = m + 1;
                n = n + 2;
	      }
              if (ii == 9) n = n + 3;
              if (jj == 9) n = n + 3;
              if (kk == 9) n = n + 3;
              if (ii == 10) m = m + 4;
              if (jj == 10) m = m + 4;
              if (kk == 10) m = m + 4;
              if (ii == 11) {
                m = m + 3;
                n = n + 1;
	      }
              if (jj == 11) {
                m = m + 3;
                n = n + 1;
	      }
              if (kk == 11) {
                m = m + 3;
                n = n + 1;
	      }
              if (ii == 12) {
                m = m + 2;
                n = n + 2;
	      }  
              if (jj == 12) {
                m = m + 2;
                n = n + 2;
	      }
              if (kk == 12) {
                m = m + 2;
                n = n + 2;
	      }
              if (ii == 13) {
                m = m + 1;
                n = n + 3;
	      }	       
              if (jj == 13) {
                m = m + 1;
                n = n + 3;
	      }
              if (kk == 13) {
                m = m + 1;
                n = n + 3;
	      }           
              if (ii == 14) n = n + 4;
              if (jj == 14) n = n + 4;
              if (kk == 14) n = n + 4;
              I[i][j][k] += J*C[0][ii]*C[1][jj]*C[2][kk]*cml[m][n];
            }
          }
        }
      }
    }
  }

  return(I);
}





void tableElements(void **elementsIn,int nElements)
{
  struct constantElement **elements;
  int i,j;

  /* Recast the sources to correct type: */
  elements = (struct constantElement**) elementsIn;

  fprintf(stdout,"WARNING: tableSources:\n"
	         "         Writing temporary output\n");

  fprintf(stdout," %d elements read \n",nElements);
  fprintf(stdout," center (x,y,z), radius \n");
  for (i=0;i<nElements;i++){
    fprintf(stdout,"Element %d:",i);
    if(elements[i]->crackedNeighbourIndex != i){
      fprintf(stdout," Cracked element. Other half is element %d.\n",
	     elements[i]->crackedNeighbourIndex);
    }
    else {
      fprintf(stdout,"\n");
    }
    for (j=0;j<elements[i]->shape;j++){
      fprintf(stdout,"   (%f,%f,%f)\n",
	      elements[i]->corners[j][0],
	      elements[i]->corners[j][1],
	      elements[i]->corners[j][2]);
    }
    fprintf(stdout,"   (%f,%f,%f) <- centroid\n",
	    elements[i]->centroid[0],
	    elements[i]->centroid[1],
	    elements[i]->centroid[2]);
   
  }
  return;
}

/*
* ==================================================================== 
* Calculate the LOCAL moments of a triangular flat panels.
* The moments are caluculated in a local curvilinear coordinate system 
* and taken with regard to the local origin (located at the first 
* corner of the panel.
* This routine was taken from the Stokes code and was originally 
* coded by Wenjing Ye. 
*
* This routine could be called for every panel. However, no geometry 
* information of the panel is used, so these values may as well be 
* calculated once and stored. Then they can be reused for every 
* panel. In either case, the difference in speed/memory will probably 
* not be a lot.
* 
* The moment M[i][j] is evaluated as the integral of eta1^i*eta2^j 
* over the triangular panel (i = 0,1,....m, j = 0,1,...n):
*
*            eta1
*              ^
*           1 -|
*              |\
*              | \ 
*              |  \
*              |   \
*              |    \
*              |     \
*              |      \
*           0 - --------> eta2
*       Origin |       |
*              0       1
*
* Parameters:
* m:        Input, maximum degree eta1^m
* n:        Input, maximum degree eta2^n
*           Note: Often n=m=2 will be used(?)
* moments:  Output. moments[i][j] contains the values of eta1^i*eta2^j
*           integrated over the panel.
*           moments must be of size at least (m+1)x(n+1): the last 
*           element accessed will be moments[m][n].
* ==================================================================== */ 
static void localElementMoments(int m, int n, double **moments)
{
  int i, j;
  /* Compute the first row and the first column of Mij. 
   * These are simply the integrals of eta1^i over the triangular 
   * panel, i.e. (using for short (x,y) for (eta1,eta2)): 
   * $\int_{x=0}^1 \int_{y=0}^{1-x} x^i dy dx $
   * This integral has the value 1/((i+1)(i+2)):
   */
  for (i=0; i<=m; i++) moments[i][0] = 1.0 / ((double) ((i+1)*(i+2)));
  for (j=0; j<=n; j++) moments[0][j] = 1.0 / ((double) ((j+1)*(j+2)));

  /* Compute the rest of matrix using a recursion relation. 
   * We want to compute 
   *  $I=\int_{x=0}^1 \int_{y=0}^{1-x} x^i y^j dy dx $
   * This recursion relation can be obtained by the following steps:
   * 1) Evaluate the INNER integral to get:
   *     $I=(1/(j+1)) \int_0^1 x^i (1-x)^{j+1} dx$
   * 2) Evaluate this integral using, e.g. Gradsteyn & Ryzhik
   *     (Eq. 2.151 in my version [5th edn.]).
   * The recursion relation follows after very few calculations.       */
  for (j=1; j<=n; j++) {
    for (i=1; i<=m; i++) {
      moments[i][j] = ((double) i) * moments[i-1][j] / ((double) (i+j+2));
    }
  }
  return;
}


