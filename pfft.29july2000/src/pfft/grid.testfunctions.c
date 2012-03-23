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
# include "Calcp.h"
/* 
* ==================================================================== 
* The following part is an attempt to accellerate the precorrection
* step. Information is need for projection, interpolation AND direct 
* stencils. For each grid point in the direct stencil make a small 
* matrix with the grid-to-grid interactions between the projection 
* stencil centered at the central node and the interpolation stencil 
* centered at this particular direct-stencil grid point. 
* The precorrection consists of calculating INTERP*G2G*PROJ and 
* subtracting the sum (over projection and interpolation points) 
* from each matrix entry. Basically this may be written as
*
*        \sum_i \sum_j I_i G_{ij} P_j   = I*G*P
*
* or, equivalently as:
*
*        \sum_i (I_i ( \sum_j P_j G_{ij} ) ) 
* 
* For all elements, which are neighbours "through the grid" the small
* matrix G_{ij} is thus precomputed for the grid-to-grid interactions.
* Large elements, which have their own (linked) interaction lists 
* need special treatment for elements in those lists.
*
* NOTE: The matrix G2G here is (interpreted as) the transpose of what 
* is used in some other places of the program. This shouldn't matter 
* at all - but please DO be careful with verifying any modifications!
*
* The intent of the small matrices is to speed up the precorrection 
* calculations. Note that the approach ONLY works if the kernel function
* is on the form G(x,y) = G(x-y). If the kernel depends on the actual 
* position in space, then another precorrection approach should be 
* chosen (the entire pFFT process fails if the kernel is not on 
* convolution form).
* This routine can possibly be speeded up by calling the routine 
* supplyKernelDataValues after the grid is set up, but before the kernel 
* is convolved. In this proposed setting the kernel values will be 
* taken off the (precomputed) values on the grid, rather than being 
* computed locally here. Further, this approach will avoid any 
* discepancy that there may be between kernel data in this routine 
* and kernel data on the grid.
* ==================================================================== */ 
static 
void setupNearFieldGrid2grid(struct grid *grid, 
			     int refreshGridSizeOnly){
  int imat,irow,icol;
  int *iInterpPts, *jInterpPts, *kInterpPts;
  /* Shorthands */
  double dx,dy,dz, *thisrow;
  struct gridStencil *directStencil, *projStencil, *interpStencil;

  /* Setup shorthands */
  directStencil = grid->directStencil;
  projStencil   = grid->elementGridStencil;
  interpStencil = grid->elementGridStencil;
  dx = grid->dx;
  dy = grid->dy;
  dz = grid->dz;

  if (refreshGridSizeOnly){
    /* Everything should be allocated. Only new values are needed.     */
    assert(grid->numNearGrid2grid == directStencil->nStencilPoints);
    assert(grid->nNearGrid2grid == projStencil->nStencilPoints);
    assert(grid->mNearGrid2grid == interpStencil->nStencilPoints);
  }
  else {
    /* Set the size parameters and allocatememory for the matrices     */
    /* One matrix per point in the direct stencil                      */
    grid->numNearGrid2grid = grid->directStencil->nStencilPoints;
    /* One matrix column per point in the projection stencil           */
    grid->mNearGrid2grid = projStencil->nStencilPoints;
    /* One matrix row per point in the interpolation stencil        */
    grid->nNearGrid2grid = interpStencil->nStencilPoints;
    
    /* Allocate matrices                                               */
    grid->nearGrid2grid = (double ***) 
      calloc(grid->numNearGrid2grid, sizeof(double **));
    /* Allocate rows in each matrix and columns for each row.          */
    for (imat=0; imat<grid->numNearGrid2grid; imat++){
      grid->nearGrid2grid[imat] = (double **) 
	calloc(grid->nNearGrid2grid, sizeof(double *));
      /* For this matrix allocate the rows (# cols in each)            */
      for (irow=0; irow<grid->nNearGrid2grid; irow++){
	grid->nearGrid2grid[imat][irow] = (double *) 
	  calloc(grid->mNearGrid2grid, sizeof(double));
      }
    }    
  }

  /* Temporary memory for storing the interpolation points 
   * relative to the projection points.                                */
  iInterpPts = (int *) calloc(grid->nNearGrid2grid,sizeof(int));
  jInterpPts = (int *) calloc(grid->nNearGrid2grid,sizeof(int));
  kInterpPts = (int *) calloc(grid->nNearGrid2grid,sizeof(int));

  /* Setup the values of the G2G matrices                              */
  /* Go through the direct points making sure the ordering of the 
   * matrices matches the ordering of the points in the direct stencil 
   * (if this was not true, then more bookkeeping would be necessary.) */
  for (imat=0; imat<grid->numNearGrid2grid; imat++){
    /* Get the location of the interpolation points centered at this 
     * direct point:                                                   */
    for (irow=0; irow<grid->nNearGrid2grid; irow++){
      iInterpPts[irow] = grid->directStencil->ipNeighbours[imat]
	                +grid->elementGridStencil->ipNeighbours[irow];
      jInterpPts[irow] = grid->directStencil->jpNeighbours[imat]
	                +grid->elementGridStencil->jpNeighbours[irow];
      kInterpPts[irow] = grid->directStencil->kpNeighbours[imat]
	                +grid->elementGridStencil->kpNeighbours[irow];
    }
    /* Actually fill this matrix. 
     * Note that this computation may possibly be speeded up by
     * precomputing the values and just stuffing the matrix with   
     * the values.                                                     */
    for (irow=0; irow<grid->nNearGrid2grid; irow++){
      thisrow = grid->nearGrid2grid[imat][irow];
      for (icol=0; icol<grid->mNearGrid2grid; icol++){
	/* Calculate the kernel for this distance.                     */
	/* Zero kernel for zero distance:                              */
	if (iInterpPts[irow]-projStencil->ipNeighbours[icol] == 0 &&
	    jInterpPts[irow]-projStencil->jpNeighbours[icol] == 0 &&
	    kInterpPts[irow]-projStencil->kpNeighbours[icol] == 0)
	  thisrow[icol] = 0.0;
	/* Get kernel in all other cases:                              */
	else
	  thisrow[icol] = 
	    kernel(dx*(iInterpPts[irow]-projStencil->ipNeighbours[icol]),
		   dy*(jInterpPts[irow]-projStencil->jpNeighbours[icol]),
		   dz*(kInterpPts[irow]-projStencil->kpNeighbours[icol]));
      } /* Next entry in this row  */
    }   /* Next row in this matrix */
  }     /* Next matrix             */

  /* Free temporary memory */
  FREE(iInterpPts);
  FREE(jInterpPts);
  FREE(kInterpPts);

  return;
}
static 
void setupNearFieldGrid2gridUnion(struct grid *grid, 
				  int refreshGridSizeOnly){
  char fctName[] = "setupNearFieldGrid2gridUnion";
  int idir, iinterp, icol;
  int ithis,jthis,kthis;
  int foundThisPoint;
  /* Temporary storage for brute-force sorting */
  struct localGridIndexLink *compressedList, *thisLink, *nextLink, *newLink;
  int nDistinctPointsFound;
  /* Shorthands */
  double dx,dy,dz, *thisrow;
  struct gridStencil *directStencil, *projStencil, *interpStencil;
  /* Setup shorthands */
  directStencil = grid->directStencil;
  projStencil   = grid->elementGridStencil;
  interpStencil = grid->elementGridStencil;
  dx = grid->dx;
  dy = grid->dy;
  dz = grid->dz;

  if (refreshGridSizeOnly){
    /* Everything should be allocated. Only new values are needed.     */
    assert(grid->nRowsNearGrid2gridUnion >= directStencil->nStencilPoints);
    assert(grid->nColsNearGrid2gridUnion == projStencil->nStencilPoints);
    assert(grid->nRowsG2gUnionRowNumbers == directStencil->nStencilPoints);
    assert(grid->nColsG2gUnionRowNumbers == interpStencil->nStencilPoints);
    compressedList = grid->localGridIndexLink;
  }
  else {
    /* Needs to allocate stuff. Also, needs to compress the "most likely" 
     * interpolation points into one vector */

    /* Find unique interpolation points (use brute force - for now don't 
     * worry about searching a lot!)                                   */
    /* Allocate array to store "compressed" index for each 
     * interpolation point                                             */
    grid->g2gUnionRowNumbers 
      = (int**) calloc(directStencil->nStencilPoints,sizeof(int*));
    for (idir=0; idir<directStencil->nStencilPoints; idir++)
      grid->g2gUnionRowNumbers[idir] 
	= (int*) calloc(interpStencil->nStencilPoints,sizeof(int));

    /* Make a linked list of indices and compressed numbers */
    nDistinctPointsFound = 0; 
    compressedList = NULL;
    /* Compress (sort) all the "most likely" interpolation points:     */
    for (idir=0; idir<directStencil->nStencilPoints; idir++){
      for (iinterp=0; iinterp<interpStencil->nStencilPoints;iinterp++){
	/* For verification purposes, set index to negative one. 
	 * At the end of compression, the indices should be positive 
	 * (or zero).                                                  */
	grid->g2gUnionRowNumbers[idir][iinterp]=-1;
	/* Start at the beginning of the list */
	nextLink = compressedList;
	foundThisPoint=0;
	while (!foundThisPoint){
	/* Go through the list */
	  ithis = directStencil->ipNeighbours[idir]
	         +interpStencil->ipNeighbours[iinterp];
	  jthis = directStencil->jpNeighbours[idir]
	         +interpStencil->jpNeighbours[iinterp];
	  kthis = directStencil->kpNeighbours[idir]
	         +interpStencil->kpNeighbours[iinterp];
	  if (nextLink==NULL){
	    /* So we did not find this point in the list. Thus it is 
	     * a new unique point. Generate an entry for it!           */
	    newLink = (struct localGridIndexLink*) 
	      calloc(1,sizeof(struct localGridIndexLink));
	    /* Set the values */
	    newLink->ip = ithis;
	    newLink->jp = jthis;
	    newLink->kp = kthis;
	    newLink->compressedIndex = nDistinctPointsFound;
	    newLink->next = compressedList;
	    /* Update the number of distinct points found              */
	    nDistinctPointsFound++;
	    /* Let this be the newest (first) link in the list!        */
	    compressedList = newLink;
	    /* Store the compressed index value with this 
	     * interpolation point                                     */
	    grid->g2gUnionRowNumbers[idir][iinterp]
	      = newLink->compressedIndex;
	    /* So, we "found" the point ;)                             */
	    foundThisPoint=1;
	    /* The following is not really necessary: */
	    newLink = NULL;
	  }
	  /* So, this link exists! Check it out! */
	  else if ((ithis==nextLink->ip)&&
		   (jthis==nextLink->jp)&&
		   (kthis==nextLink->kp)){
	    /* This position (in the linked list) is the same as the 
	     * present position we are examining. So, just store the 
	     * compressed index with the present interpolation point   */
	    grid->g2gUnionRowNumbers[idir][iinterp]
	      = nextLink->compressedIndex;
	    /* And we actually found the point                         */
	    foundThisPoint=1;
	  }
	  else {
	    /* So, we haven't found it yet. No big deal - just try 
	     * the next link in the list.                              */
	    nextLink = nextLink->next;
	  }
	} /* Repeat until found this point                             */
	/* Make sure that the present index was given a value!         */
	assert(grid->g2gUnionRowNumbers[idir][iinterp]>=0);
      } /* Next point in the interpolation stencil                     */
    } /* Next point in the direct stencil                              */
    fprintf(stdout,
	    "%s: Found %d distinct most likely interpolation points\n",
	    fctName,nDistinctPointsFound);
    /* The points are now compressed. 
     * Store the linked list of compressed indices with the grid       */
    grid->localGridIndexLink=compressedList;

    /* Allocate memory for storing grid-to-grid interactions           */
    grid->nRowsNearGrid2gridUnion = nDistinctPointsFound;
    grid->nColsNearGrid2gridUnion = projStencil->nStencilPoints;
    grid->nearGrid2gridUnion = 
      (double **) calloc(grid->nRowsNearGrid2gridUnion,sizeof(double*));
    for(iinterp=0;iinterp<grid->nRowsNearGrid2gridUnion;iinterp++)
       grid->nearGrid2gridUnion[iinterp] = 
      (double *) calloc(grid->nColsNearGrid2gridUnion,sizeof(double));

    /* Setup work memory of the needed size */
    grid->nearGrid2gridUnionWork = 
      (double *) calloc(grid->nRowsNearGrid2gridUnion,sizeof(double));
    
  } /* End of memory setup and sorting */

  /* Now set the values of the grid-to-grid matrix for each (distinct) 
   * pair of projection and interpolation points                       */
  /* Go through the list of distinct interpolation points */
  thisLink =  grid->localGridIndexLink;
  nDistinctPointsFound=0;
  while (thisLink!=NULL){
    /* Step the number of found interpolation points (for verification) */
    nDistinctPointsFound++;
    iinterp = thisLink->compressedIndex;
    ithis   = thisLink->ip;
    jthis   = thisLink->jp;
    kthis   = thisLink->kp;
    /* Step through the projection points */
    for (icol=0;icol<projStencil->nStencilPoints;icol++){
      /* Get and store the kernel for this pair of interpolation 
       * and projection points.                                      */
      /* Zero kernel for zero distance:                              */
      if (ithis-projStencil->ipNeighbours[icol] == 0 &&
	  jthis-projStencil->jpNeighbours[icol] == 0 &&
	  kthis-projStencil->kpNeighbours[icol] == 0)
	grid->nearGrid2gridUnion[iinterp][icol] = 0.0;
      /* Get kernel in all other cases:                              */
      else
	grid->nearGrid2gridUnion[iinterp][icol] = 
	  kernel(dx*(ithis-projStencil->ipNeighbours[icol]),
		 dy*(jthis-projStencil->jpNeighbours[icol]),
		 dz*(kthis-projStencil->kpNeighbours[icol]));
    } /* Next column (projection point) */
    /* Examine next link in list */
    thisLink = thisLink->next;
  } /* Next interpolation point in linked list */
  /* Check that we have been through the right number of 
   * interpolation points                                            */
  assert(nDistinctPointsFound==grid->nRowsNearGrid2gridUnion);

  return;
} /* End of routine setupNearFieldGrid2gridUnion */


/* 
================================================================= 
Make a list of elements and their neighbours.
================================================================= */ 
static
void testTableElementNeighbours(int elementType,
				struct grid *grid,
				struct gridStencil *stencil)
{
  int iae;
  struct neighbourLinkedList 
    *thisLink;     /* Used for accessing the linked lists         */
  /* Show the obtained lists of neighbours */
  if (elementType==0){
    printf(" LIST OF NEIGHBOURS TO SOURCES \n");
    printf(" Source #   Neighbours \n");
  }
  else{
    printf(" LIST OF NEIGHBOURS TO EVALS   \n");
    printf(" Eval #     Neighbours \n");
  }
  for (iae=0; iae < grid->nGridElements[elementType]; iae++){
    if (stencil->gridElementNeighbours[elementType][iae]->neighbourLinkedList!=NULL){
      printf("     %4d  ",iae);
      thisLink = stencil
                   ->gridElementNeighbours[elementType][iae]
                     ->neighbourLinkedList;
      while(thisLink!=NULL){
	printf(" %d,",thisLink->iNeighbour);
	thisLink = thisLink->nextNeighbour;
      }
      printf("\n");
    }
  }
  return;
}

/* 
================================================================= 
Make a list of the gridElements - bounding spheres.
================================================================= */ 
static
void tableGridElements(int elementType,struct grid *grid)
{
  int i;

  /* Write table header */
  fprintf(stdout," gridElement bounding spheres\n");

  /* Write table. Only include non-zero entries. */
  for (i=0; i<grid->nGridElements[elementType]; i++){
    printf(" (%f, %f, %f) %f \n",
	   grid->gridElements[elementType][i]->center[0],
	   grid->gridElements[elementType][i]->center[1],
	   grid->gridElements[elementType][i]->center[2],
	   grid->gridElements[elementType][i]->radius);
  }

  return ; 
}


