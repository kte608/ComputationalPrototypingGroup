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
/* Test routines. This should *not* be part of the 
 * End product code! */

void testINDEX(void);

/* To make life easy on the programmer, get information on what 
 * the grid and the elements are: */
# include "Grid.h" 
# include "Elements.h" 
# include "Sparse.h" 

void testINDEX(void)
{
  /* Test INDEX and inverse mapping INDEX_IJK */
  int i,j,k, s1,s2,s3, n, i1,j1,k1;
  
  s1=20;s2=30;s3=40;

  for (i=1; i<= s1; i++){
    for (j=1; j<= s2; j++){
      for (k=1; k<= s3; k++){
	n = INDEX(i,j,k, s1,s2,s3);
	INDEX_IJK_OLD(n, i1,j1,k1, s1,s2,s3);
	if (0) {  
	  /* This works fine. But it could spell trouble that floating
	   * point operations are used */
	  int left;
	  i1    = 1 + (int) floor( ((double)n)/( (double) s2*s3)); 
	  left = n - (i1-1)*s2*s3;                                
	  j1    = 1 + (int) floor( ((double)left)/( (double) s3)); 
	  k1    = left - (j1-1)*s3 +1;

	  printf("Shouldn't execute this!\n");                         
	}
	
	/* Test with integer operations only */
	if (0) { 
	  int left;
	  /* Integer division truncates any fractional part ! */
	  i1    = 1 + n/(s2*s3);
	  left =  n - (i1-1)*s2*s3;
	  j1    = 1+  left/s3;
	  k1    = left - (j1-1)*s3 +1;
	  assert(n==INDEX(i1,j1,k1,s1,s2,s3));

	  printf("Shouldn't execute this!\n");                         
	}
	if ( i1!=i || j1!=j || k1!=k ){
	  printf(" DISCREPANCY: %d, (%d,%d,%d), (%d,%d,%d)\n",
		 n, i,j,k, i1,j1,k1);
	}
      }
    }
  }
printf("TEST completed! n=0:%d\n",n);
  _EXIT_;
}


void testProject2grid(void *gridIn, void *projectMat, 
		      void *interpMat, int ncols)
{
  /* Test the projection onto the grid */
  struct grid *grid;
  fftGridData *gridData, *g;
  double *src, *eval, *data1;
  int igsrc, i,j,k, ip,jp,kp, ic,jc,kc, found, isrc, *isrcs,
    startRow, endRow, imax, index;
  gridPoint ****gp;
  void *srcp;
  sparseMat *pmat;
  
  /* Cast input data to correct type: */
  grid = (struct grid *) gridIn;
  /* Extract grid data grom the grid: */
  g = gridData = (fftGridData *) grid->griddata;
  /* Local shorthand: */
  gp = grid->points;


  printf(" grid info: nx=%d, ny=%d, nz=%d\n",grid->nx,grid->ny,grid->nz);
  printf(" data info: s1=%d, s2=%d, s3=%d\n",g->s1,g->s2,g->s3);

  /* Pick one source to examine */ 
  for (i=0,found=0;i<grid->nx;i++){
    for (j=0;j<grid->nx;j++){
      for (k=0;k<grid->nx;k++){
	if (gp[i][j][k]!=NULL && gp[i][j][k]->nAssociatedElements[0]>0){
	  ip=i;
	  jp=j;
	  kp=k;
	  /* First source on this point: */
	  igsrc = gp[ip][jp][kp]->associatedElements[0][0];
	  printf("Using grid-source %d associated with point  %d %d %d \n",
		 igsrc,ip,jp,kp);
	  printf("This element should have center (%g,%g,%g)\n"
		 " and radius %g\n",
		 grid->gridElements[0][igsrc]->center[0],
		 grid->gridElements[0][igsrc]->center[1],
		 grid->gridElements[0][igsrc]->center[2],
		 grid->gridElements[0][igsrc]->radius);
	  /* pointer to this source */
	  srcp = grid->gridElements[0][igsrc]->element;
	  /*nsrcs = numSourceUnknowns(srcp);*/
	  isrcs = sourceIndices(srcp);
	  isrc = isrcs[0]; /* First direct-index of this element */
	  printf("Direct-index of this source is: %d\n",isrc);
	  found = 1;
	  break;
	}
      }
      if (found) break;
    }
    if (found) break;
  }
  
  /* Allocate storage for sources to be projected */
  src = (double *) calloc(ncols, sizeof(double));
  /* Make the chosen source be a unit source */
  src[isrc]=1.0;

  /* Make sure all the grid data is zeroed out */
  for (i=1;i<=gridData->s1;i++){
    for (j=1;j<=gridData->s2;j++){
      for (k=1;k<=gridData->s3;k++){
	gridData->data[i][j][k] = 0.0;
      }
    }
  }
  /* Project onto grid */
  projectOntoGrid(grid, src, projectMat);

  /* Test for nonzero grid charges */
  data1 = &(gridData->data[1][1][1]);
  imax = gridData->s1*gridData->s2*gridData->s3;
  for (index=0;index<imax;index++){
    if (data1[index]!=0.0){
      INDEX_IJK(index, i,j,k,g->s1,g->s2,g->s3)
      printf("A: charge %d at [%d][%d][%d] has nonzero value: %16.8e\n",
	     index,i-1,j-1,k-1,data1[index]);
    }
  }

  /* Test for nonzero grid charges */
  for (i=0;i<grid->nx;i++){
    for (j=0;j<grid->ny;j++){
      for (k=0;k<grid->nz;k++){
	if (gridData->data[i+1][j+1][k+1]!=0.0){
	  ic=i;jc=j;kc=k;
	  printf("charge at [%d][%d][%d] has nonzero value: %16.8e\n",
		 i,j,k,gridData->data[i+1][j+1][k+1]);
	  /* if (gp[i][j][k]==NULL)
	     printf(" ...however [%d][%d][%d] has no associated elements\n",
	     i,j,k);*/
	}
      }
    }
  }

  /* Examine the projection matrix - because something was obviously wrong */
  pmat = (sparseMat *) projectMat;
  for(j = 0; j < pmat->n; j++) {
    startRow = pmat->colIndex[j];
    endRow   = pmat->colIndex[j+1]; 
    for(i = startRow;  i < endRow; i++) {
      printf("Projection matrix has nonzero value at (%d,%d)\n",
	     pmat->rowIndex[i],j);
    }
  }

  printf("the projection matrix now seems fine.\n"
	 "Test the interpolation matrix.\n");

  printf(" Interpolation test 1: Only eval #%d should be nonzero\n",isrc);
  /* Allocate storage for evals to be interpolated
   * (exploit that the number of sources and evals is the same)        */
  eval = (double *) calloc(ncols, sizeof(double));

  /* Interpolate */
  interpolateFromGrid(gridIn,eval,interpMat);

  /* Check non-zero evals: */
  for (i=0;i<ncols;i++){
    if (eval[i]!=0.0)
      printf("Eval #%d nonzero, value %16.8e\n",i,eval[i]);
  }

  printf(" Interpolation test 2: Every eval EXCEPT #%d should be nonzero\n",isrc);
  /* Zero out the evals */
  for (i=0;i<ncols;i++){
    eval[i]=0.0;
  }
  /* Make every grid potential unity *except* the one related to 
   * the eval in question */
  for (i=1;i<=gridData->s1;i++){
    for (j=1;j<=gridData->s2;j++){
      for (k=1;k<=gridData->s3;k++){
	gridData->data[i][j][k] = 1.0;
      }
    }
  }
  gridData->data[ic+1][jc+1][kc+1]=0.0;
  /* Interpolate */
  interpolateFromGrid(gridIn,eval,interpMat);

  /* Check non-zero evals: */
  for (i=0;i<ncols;i++){
    if (eval[i]!=0.0)
      printf("Eval #%d nonzero, value %16.8e\n",i,eval[i]);
  }

  /* Find a column of the grid-interaction matrix */
  /* Zero out sources and evals */
  for (i=0;i<ncols;i++)
    eval[i]=src[i]=0.0;
  src[0]=1.0;

  projectOntoGrid(gridIn, src, projectMat);
  grid2gridCalculate(gridIn);
  interpolateFromGrid(gridIn,eval,interpMat);
  
  for (i=0;i<ncols;i++)
    printf("grid col %d: %16.8e\n",0,eval[i]);


  return;
}

void testGrid2grid(void *gridIn)
{ 
/* Test the grid-to-grid interactions. 
 * This part seems OK */
 
  struct grid *grid;
  fftGridData *gridData;
  int i,j,k, ic,jc,kc;
  double delta, gval, errval;

  /* Cast input data to correct type */
  grid = (struct grid *) gridIn;

  delta = grid->dx;
  assert(delta==grid->dy);
  assert(delta==grid->dz);

  /* Extract grid data grom the grid: */
  gridData = (fftGridData *) grid->griddata;


  /* Make sure all the grid data is zeroed out */
  for (i=1;i<=gridData->s1;i++){
    for (j=1;j<=gridData->s2;j++){
      for (k=1;k<=gridData->s3;k++){
	gridData->data[i][j][k] = 0.0;
      }
    }
  }


  /* Include a single grid charge */
  ic=jc=kc= grid->nx;
  gridData->data[ic][jc][kc] = 1.0;
  
  /* Make grid-to-grid interactions */
  grid2gridCalculate(gridIn);

  /* Test that the grid charges now has the correct format 
   * (real part of the grid only, I guess) */
  errval=0.0;
  for (i=1;i<=grid->nx;i++){
    for (j=1;j<=grid->ny;j++){
      for (k=1;k<=grid->nz;k++){
	if (i!=ic || j!=jc || k!=kc){
	  gval = 1.0/sqrt(grid->dx*grid->dx*((double) (i-ic)*(i-ic))
			  +grid->dy*grid->dy*((double) (j-jc)*(j-jc))
			  +grid->dz*grid->dz*((double) (k-kc)*(k-kc)));
	  errval = MAX(abs(gval-gridData->data[i][j][k]),errval);
	}
      }
    }
  }
  printf("testGrid2grid: maximum error = %16.8e\n",errval);
  printf("testGrid2grid: Visual inspection, DX=%16.8e:\n",grid->dx);
  for (i=1;i<=grid->nx;i++){
    printf("%d %16.8e\n",i,gridData->data[i][jc][kc]);
  }


  return;
}
