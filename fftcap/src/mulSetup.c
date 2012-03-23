/*!\page LICENSE LICENSE
 
Copyright (C) 2003 by the Board of Trustees of Massachusetts Institute of Technology, hereafter designated as the Copyright Owners.
 
License to use, copy, modify, sell and/or distribute this software and
its documentation for any purpose is hereby granted without royalty,
subject to the following terms and conditions:
 
1.  The above copyright notice and this permission notice must
appear in all copies of the software and related documentation.
 
2.  The names of the Copyright Owners may not be used in advertising or
publicity pertaining to distribution of the software without the specific,
prior written permission of the Copyright Owners.
 
3.  THE SOFTWARE IS PROVIDED "AS-IS" AND THE COPYRIGHT OWNERS MAKE NO
REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED, BY WAY OF EXAMPLE, BUT NOT
LIMITATION.  THE COPYRIGHT OWNERS MAKE NO REPRESENTATIONS OR WARRANTIES OF
MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT THE USE OF THE
SOFTWARE WILL NOT INFRINGE ANY PATENTS, COPYRIGHTS TRADEMARKS OR OTHER
RIGHTS. THE COPYRIGHT OWNERS SHALL NOT BE LIABLE FOR ANY LIABILITY OR DAMAGES
WITH RESPECT TO ANY CLAIM BY LICENSEE OR ANY THIRD PARTY ON ACCOUNT OF, OR
ARISING FROM THE LICENSE, OR ANY SUBLICENSE OR USE OF THE SOFTWARE OR ANY
SERVICE OR SUPPORT.
 
LICENSEE shall indemnify, hold harmless and defend the Copyright Owners and
their trustees, officers, employees, students and agents against any and all
claims arising out of the exercise of any rights under this Agreement,
including, without limiting the generality of the foregoing, against any
damages, losses or liabilities whatsoever with respect to death or injury to
person or damage to property arising from or out of the possession, use, or
operation of Software or Licensed Program(s) by LICENSEE or its customers.
 
*/


#include "mulGlobal.h"

cube *cstack[1024];		/* Stack used in several routines. */

/*
  sets up the partitioning of space and room for charges and expansions
*/
ssystem *mulInit(autom, depth, order, charges, gsize)
int autom;			/* ON => choose depth automatically */
int depth, order;
charge *charges;
int gsize; 
{
  ssystem *sys;
  int qindex=1, cindex=1;

  CALLOC(sys, 1, ssystem, ON, AMSC);
  sys->depth = depth;		/* overwritten below if autom = ON */
  sys->order = order;
  sys->gridsize = gsize-1; 

  sys->depth = placeq(autom, sys, charges); /* create cubes, put in charges */


/*   getrelations(sys); */ 		/* Get all the prnts and kids for each cube. */

  setPosition(sys);		/* Figures out position of cube center. */
/*  indexkid(sys, sys->cubes[0][0][0][0], &qindex, &cindex);  */ 
  indexcubes(sys, sys->cubes, &qindex, &cindex);  
				/* Index chgs and cubes. */

#if (1==0)
#if ADAPT == ON
  setExact(sys, multerms(sys->order)); /* Note cubes to be done exactly and
					   determine the number of nonempty
					   kids for each cube. */
#else
  setExact(sys, 0);		/* Note cubes to be done exactly and
					   determine the number of nonempty
					   kids for each cube. */
#endif
#endif 

  getnbrs(sys);			/* Get all the nearest neighbors. At bot level
				   add as nearest nbrs cubes in exact block. */
  newlinkcubes(sys);		/* Make linked-lists of direct, multis, and
				   locals to do at each level. */
  set_vector_masks(sys);  	/* set up sys->is_dummy and sys->is_dielec */
/*  setMaxq(sys);  */                 /* Calculates the max # chgs in cubes treated
				   exactly, and over lowest level cubes. */
/*  getAllInter(sys);	*/	/* Get the interaction lists at all levels. */

  return(sys);
}

/*
  places the charges using best number of levels for adaptive algorithm
  - returns number of levels used
  - can be called with flag == OFF to set depth = sys->depth
  - uses entire panel list, dieletric and dummy panels are included
     (excluding BOTH or DIELELC type panels if NUMDPT == 2, ie 
     two-point flux-density-differences are being use ),
     so the constraining type of exactness is usually
     local expansion exactness---using this method the switch to multipole
     exactness is done automatically if there are no dielectric panels
  - when the method without dielectric panels is set up, these loops will
     need to ignore all dielectric panels
  - this routine is still called to set automatic levels if ADAPT is OFF,
     ie even when the calculation is not adaptive, so results can be compared
*/
static int placeq(flag, sys, charges)
int flag;			/* ON => set depth automatically */
ssystem *sys;
charge *charges;
{
  int i, j, k, l, side, totalq, isexact, multerms(), depth;
  int xindex, yindex, zindex, limit = multerms(sys->order), compflag;
  int exact_cubes_this_level, cubes_this_level;
  double length0, length, exact_ratio;
  double minx, maxx, miny, maxy, minz, maxz, tilelength(), maxTileLength;
  charge *nextq, *compq;
  cube *****cubes, *nextc;
  char *hack_path();
  int wflag = 0; 
  int dj, dk, dl, gridj, gridk, gridl, minj, mink, minl, maxj, maxk, maxl, gridsize=sys->gridsize; 
  int nself,  size1, size2, size3; 
  double cdirect; 
  double cost[MAXDEP], slog, mincost, cf, cd,memuse[MAXDEP]; 
  int optlevel, mindepth, npts; 
  long ****ncubes, *nextnc; 
  /* Figure out the length of lev 0 cube and total number of charges. */
  limit = 30; 

  limit = gridsize * gridsize * gridsize; 
  nextq = charges;
  minx = maxx = nextq->x;
  miny = maxy = nextq->y;
  minz = maxz = nextq->z;

  for(totalq = 1, nextq = nextq->next; nextq != NULL;
      totalq++, nextq = nextq->next) {
    maxx = MAX(nextq->x, maxx);
    minx = MIN(nextq->x, minx);
    maxy = MAX(nextq->y, maxy);
    miny = MIN(nextq->y, miny);
    maxz = MAX(nextq->z, maxz);
    minz = MIN(nextq->z, minz);
  }


  /* Make sure cube isn't smaller than a tile. */

  for(maxTileLength = 0.0, nextq = charges;
      nextq != NULL; nextq = nextq->next) {
    length = tilelength(nextq);
    maxTileLength = MAX(maxTileLength, length);
  }
  maxx += 0.5 * maxTileLength;
  minx -= 0.5 * maxTileLength;
  maxy += 0.5 * maxTileLength;
  miny -= 0.5 * maxTileLength;
  maxz += 0.5 * maxTileLength;
  minz -= 0.5 * maxTileLength;
#if GNDPLN == ON 
  if (minz < 0.0) {
    minz = 0.01 * maxTileLength; 
    printf("warning: minz was adjusted to be < 0; check panels\n"); 
  }
#endif 

  sys->minx = minx;
  sys->miny = miny;
  sys->minz = minz;

#if GNDPLN  == ON 
  /* presume ground-plane located at z=0 */ 
  GroundPlaneHeight = minz; 
#endif

  /* (see below for test for panel size vs cube size) */

  length0 = MAX((maxx - minx), (maxy - miny));
  length0 = MAX((maxz - minz), length0);


  /* Create the vectors for storing the charges and coefficients. */
  CALLOC(sys->q, totalq+1, double, ON, AMSC);
  CALLOC(sys->p, totalq+1, double, ON, AMSC);

  /* set up mask vector: is_dummy[i] = TRUE => panel i is a dummy */
  CALLOC(sys->is_dummy, totalq+1, int, ON, AMSC);

  /* set up mask vector: is_dielec[i] = TRUE => panel i is on DIELEC or BOTH */
  CALLOC(sys->is_dielec, totalq+1, int, ON, AMSC);

  if(flag == ON) {		/* set depth of partitions automatically */
    /* alloc spine of cube pntr array - leave enough room for depth = MAXDEP */
    FCALLOC(ncubes, MAXDEP+1, long***, ON, ACUB); 

    /* allocate for levels 0, 1, and 2 (always used) */
    /* in FFT code, allocate for level 0 only */ 

    side = 1; 
    i=0; 

    /* side /= 2; */

    /* for each level > 2: allocate for full cubes, count panels in each 
       - quit loop if all lowest level cubes are exact */
    mindepth = 1; 
    for(isexact = FALSE; isexact == FALSE;  i++) {

      npts = 1 << i; 


      side = (npts -1) / gridsize; 

      printf("NPTS  = %d SIDE = %d\n", npts, side); 

      if (side < 1)  {
        mindepth = i+1; 
        continue; 
      }
      
      if(i > MAXDEP) {
	fprintf(stderr, 
		"placeq: out of cube pntr space - increase MAXDEP == %d\n", 
		MAXDEP);
	exit(0);
      }
      length = 1.01 * length0 / (double) side; 

      FCALLOC(ncubes[i], side, long**, ON, ACUB);
      for(j=0; j < side; j++) {
	FCALLOC(ncubes[i][j], side, long*, ON, ACUB);
	for(k=0; k < side; k++) {
	  FCALLOC(ncubes[i][j][k], side, long, ON, ACUB);
          for (l=0; l< side; l++) {
            ncubes[i][j][k][l] = 0; 
          }
	}
      }
      /* Count the number of charges per cube and allocate if needed */
      for(nextq = charges; nextq != NULL; nextq = nextq->next) {
	xindex = (nextq->x - minx) / length;
	yindex = (nextq->y - miny) / length;
	zindex = (nextq->z - minz) / length;
        ncubes[i][xindex][yindex][zindex]++; 
      }
    
      /* if the current lowest level is not exact, loop back until it is */
      /*    check for exactness of this level, get cube statistics */
      isexact = TRUE;
      cubes_this_level = 0;
      exact_cubes_this_level = 0;
      for(j = 0; j < side; j++) {
	for(k = 0; k < side; k++) {
	  for(l = 0; l < side; l++) {
            if(ncubes[i][j][k][l] > limit) {
              isexact = FALSE;
            }
            else { 
              if (ncubes[i][j][k][l] > 0) 
                exact_cubes_this_level++;
            }
            if (ncubes[i][j][k][l] > 0) { 
              cubes_this_level++;
            }
	  }
	}
      }

      minj = side; maxj = 0; 
      mink = side; maxk = 0; 
      minl = side; maxl = 0; 
      /* cost estimates */ 
      cdirect = 0; 
      for(j = 0; j < side; j++) {
	for(k = 0; k < side; k++) {
	  for(l = 0; l < side; l++) {
	    if(ncubes[i][j][k][l] >0) {
              nself = ncubes[i][j][k][l]; 
              if (j < minj) minj = j; 
              if (k < mink) mink = k; 
              if (l < minl) minl = l; 
              
              if (j > maxj) maxj = j; 
              if (k > maxk) maxk = k; 
              if (l > maxl) maxl = l; 
              cdirect += (double) nself * nself; 

              if (j > 0) 
                  cdirect += (double) nself * (ncubes[i][j-1][k][l]); 

              if (k > 0) 
                  cdirect += (double)nself * (ncubes[i][j][k-1][l]); 

              if (l > 0) 
                  cdirect += (double)nself * (ncubes[i][j][k][l-1]); 

              if (j < side-1) 
                  cdirect += (double)nself * (ncubes[i][j+1][k][l]); 

              if (k < side-1) 
                  cdirect += (double)nself * (ncubes[i][j][k+1][l]); 

              if (l < side-1) 
                cdirect += (double)nself * (ncubes[i][j][k][l+1]);              
            }
          }
        }
      }

      for(j = 0; j < side; j++) {
	for(k = 0; k < side; k++) {
          free(ncubes[i][j][k]);            
        }
        free(ncubes[i][j]);   
      }
      free(ncubes[i]);   
      
                
      dj = (maxj - minj +1) * gridsize + 1;
      dk = (maxk - mink +1) * gridsize + 1;
      dl = (maxl - minl +1) * gridsize + 1;
      for(gridj=1; gridj < dj; gridj *= 2);
      for(gridk=1; gridk < dk; gridk *= 2);
      for(gridl=1; gridl < dl; gridl *= 2);

      size1 = 2 * gridj; size2 = 2 * gridk; size3 = 2 * gridl;

     
      slog = log((double)size1) + log((double)size2) + log((double)size3); 
      slog /= log(2.0); 

/* costs in arbitrary units */ 
#define COSTD 4.9e-5
#define COSTF (0.25*5.1e-5)
      
      /* total cost = COSTD * cdirect + COSTF * size1*size2*size3 * \sum (log_2 size)  */ 
      cd = COSTD * cdirect; cf =  COSTF * size1 * size2*size3*slog; 
      cost[i] = cf + cd; 
      memuse[i] = cdirect + size1*size2*size3; 
      memuse[i] =  8.0 * memuse[i]/(double)1048576.0; 
      
      printf("LEVEL %d SETUP COST EST.: \n", i); 
      printf("      DIRECT INTERACTIONS: %15g\n", cdirect); 
      printf("      FFT SIZES = %10d %10d %10d\n", size1, size2, size3); 
      printf("      GRID DJ,DK,DL = %d,%d,%d\n", dj, dk, dl); 
      printf("      CHGS IN (%d,%d), (%d,%d), (%d,%d)\n", minj, maxj, mink, maxk,minl, maxl); 
      printf("      COSTS: %10g %10g\n", cd, cf);
      printf("      MEMORY: %10g\n", memuse[i]); 
      printf("\n\n"); 


      /*    decide whether to go down another level by checking exact/ttl */
      exact_ratio = (double)exact_cubes_this_level/(double)cubes_this_level;
#define EXRTSH2 0.95
      if(exact_ratio > EXRTSH2) 
	  isexact = TRUE;      	/* set up to terminate level build loop */

/* if we increase by this much we stop */ 
#define FAC 2.0 
      if (i > 1 && i > mindepth && cost[i] > FAC * cost[i-1] && memuse[i] > FAC * memuse[i-1]) 
        isexact  = TRUE; 

    }

    mincost = 1e20; 
    printf("ESTIMATED COSTS:\n"); 
    for (j=mindepth; j<i; j++) {
      printf("     %2d : %10g %10g %10g\n", j, cost[j], memuse[j], cost[j]*memuse[j]); 
      if (cost[j] * memuse[j] < mincost) {
        mincost = cost[j] * memuse[j]; 
        optlevel = j; 
      }
    }
    printf("optimal level is %d, cost %g, mem %g, product %g\n", 
           optlevel, cost[optlevel], memuse[optlevel], cost[optlevel]*memuse[optlevel]); 
    depth = i - 1;		/* the automatically set depth */
    side /= 2;

    depth = optlevel; 
    npts = 1 << depth; 
    side = (npts - 1) / gridsize; 

    length = 1.01 * length0 / (double) side; 

    printf("auto set depth to %d, side to %d, length %g\n", depth, side, length); 
    sys->depth = depth; 
    sys->side = side; 
    free(ncubes);   
  }


  /* Allocate the cubes, note calloc used because zeros everything. */
  depth = sys->depth;
  npts = 1 << depth; 
  side = (npts - 1) / gridsize; 

  CALLOC(cubes, sys->depth+1, cube****, ON, ACUB);
  i = depth; 
  CALLOC(cubes[i], side, cube***, ON, ACUB);
    for(j=0; j < side; j++) {
      CALLOC(cubes[i][j], side, cube**, ON, ACUB);
      for(k=0; k < side; k++) {
        CALLOC(cubes[i][j][k], side, cube*, ON, ACUB);
      }
    }
  length = 1.01 * length0 / (double) side; 


  /* Count the number of charges per cube. */
  for(nextq = charges; nextq != NULL; nextq = nextq->next) {
    xindex = (nextq->x - minx) / length;
    yindex = (nextq->y - miny) / length;
    zindex = (nextq->z - minz) / length;
    nextc = cubes[depth][xindex][yindex][zindex];
    if(nextc == NULL) {
      CALLOC(nextc, 1, cube, ON, ACUB);
      cubes[depth][xindex][yindex][zindex] = nextc;
      nextc->upnumvects = 1;
      CALLOC(nextc->upnumeles, 1, int, ON, ACUB);
      nextc->upnumeles[0] = 1;
    }
    else {
      nextc->upnumeles[0]++;
    }
  }

  sys->length = length;
  sys->side = side;
  sys->cubes = cubes;
  printf("minx = %g maxx = %g length = %g\n", sys->minx, maxx, sys->length); 
#if GNDPLN  == ON 
  GroundPlaneHeight /= length; 
  printf("Ground plane at z = %g  = %g cell spacings\n",   GroundPlaneHeight*length,GroundPlaneHeight); 
  printf("minz = %g maxz = %g \n", minz, maxz); 
  if (GroundPlaneHeight < 0) {
/*    printf("Warning: some charges are very close to ground plane. \n");  */ 
  }
#endif
  /* Allocate space for the charges. */
  for(j=0; j < side; j++) {
    for(k=0; k < side; k++) {
      for(l=0; l < side; l++) {
        nextc = sys->cubes[depth][j][k][l];
        if(nextc != NULL) {  /* Only fill out nonempty cubes. */
          if (nextc->upnumeles[0] <= 0) {
            printf("bad # charges in cube (%d,%d,%d), depth %d, # %d\n", 
                   j,k,l,depth,nextc->upnumeles[0]); 
            exit(0); 
          }
          /* Allocate for the charge ptrs, and get q vector pointer. */
	  CALLOC(nextc->chgs, nextc->upnumeles[0], charge*, ON, ACUB);
	  CALLOC(nextc->upnumeles, 1, int, ON, ACUB);
	/* Zero the numchgs to use as index. */
	  nextc->upnumeles[0] = 0;
	}
      }
    }
  }

  printf("maxTileLength ,  cube length = %g %g\n", maxTileLength ,  length); 

  /* Put the charges in cubes; check to make sure they are not too big. */
  for(nextq = charges; nextq != NULL; nextq = nextq->next) {
    if(tilelength(nextq) > length && wflag == 0) {
      wflag = 1; 
      fprintf(stderr,
	      "\nplaceq: Warning, a panel is larger than the cube supposedly containing it\n");
      fprintf(stderr,"  cube length = %g panel length = %g\n", 
	      length, tilelength(nextq));
    }
    xindex = (nextq->x - minx) / length;
    yindex = (nextq->y - miny) / length;
    zindex = (nextq->z - minz) / length;
    nextc = cubes[depth][xindex][yindex][zindex];

    /* check if current charge is same as those already in the cube `nextc' */
    for(compflag = FALSE, i = (nextc->upnumeles[0] - 1); i >= 0; i--) {
      compq = nextc->chgs[i];
      if((compq->x == nextq->x) &&
         (compq->y == nextq->y) && (compq->z == nextq->z)) {
	fprintf(stderr, "placeq: Warning, removing identical");
	if(compq->shape == 3) fprintf(stderr, " triangular");
	else if(compq->shape == 4) fprintf(stderr, " quadrilateral");
	else fprintf(stderr, " illegal-shape");
	fprintf(stderr, " panel\n  rmved ctr = (%g %g %g) surf = `%s'", 
		compq->x, compq->y, compq->z, hack_path(compq->surf->name));
	fprintf(stderr, " trans = (%g %g %g)\n", compq->surf->trans[0],
		compq->surf->trans[1], compq->surf->trans[2]);
	fprintf(stderr, "  saved ctr = (%g %g %g) surf = `%s'", 
		nextq->x, nextq->y, nextq->z, hack_path(nextq->surf->name));
	fprintf(stderr, " trans = (%g %g %g)\n", nextq->surf->trans[0],
		nextq->surf->trans[1], nextq->surf->trans[2]);
        /* Remove charge from linked list. */
        for(compq = charges; compq->next != nextq; compq = compq->next) {};
        compq->next = nextq->next;
        nextq = compq;
        compflag = TRUE;
      }
    }
    if(compflag == FALSE) nextc->chgs[nextc->upnumeles[0]++] = nextq;

  }
  return(depth);
}

/*
GetRelations allocates parents links the children. 
*/
getrelations(sys)
ssystem *sys;
{
cube *nextc, *parent, *****cubes = sys->cubes;
int i, j, k, l, side;
  for(i = sys->depth, side = sys->side; i >= 0; i--, side /= 2) {
    for(j=0; j < side; j++) {
      for(k=0; k < side; k++) {
	for(l=0; l < side; l++) {
          nextc = cubes[i][j][k][l];
	  if(nextc != NULL) {
	/* Get the parents and children pointers of nonempty cubes. */
	    if(i < sys->depth) {
	      nextc->numkids = 8; /* all cubes, even empties, are counted */
	      CALLOC(nextc->kids, nextc->numkids, cube*, ON, ACUB);
	      nextc->kids[0] = cubes[i+1][2*j][2*k][2*l]; /* empties get */
	      nextc->kids[1] = cubes[i+1][2*j][2*k][2*l+1]; /* null pointers */
	      nextc->kids[2] = cubes[i+1][2*j][2*k+1][2*l];
	      nextc->kids[3] = cubes[i+1][2*j][2*k+1][2*l+1];
	      nextc->kids[4] = cubes[i+1][2*j+1][2*k][2*l];
	      nextc->kids[5] = cubes[i+1][2*j+1][2*k][2*l+1];
	      nextc->kids[6] = cubes[i+1][2*j+1][2*k+1][2*l];
	      nextc->kids[7] = cubes[i+1][2*j+1][2*k+1][2*l+1];
	    }
	    if(i > 0) {
	      parent = cubes[i-1][j/2][k/2][l/2];
	      if(parent == NULL) {
		CALLOC(parent, 1, cube, ON, ACUB);
		cubes[i-1][j/2][k/2][l/2] = parent;
	      }
	      nextc->parent = parent;
	    }
	  }
	}
      }
    }
  }
}

/*
Set the position coordinates of the cubes.
*/
setPosition(sys)
ssystem *sys;
{
  int i, j, k, l;
  int side = sys->side;
  double length = sys->length;
  cube *nextc;

/* Mark the position of the lowest level cubes. */

  i = sys->depth; 

/*  for(i=sys->depth; i >= 0; i--, side /= 2, length *= 2.0) { */ 
    for(j=0; j < side; j++) {
      for(k=0; k < side; k++) {
	for(l=0; l < side; l++) {
	  nextc = sys->cubes[i][j][k][l];
	  if(nextc != NULL) {
	    nextc->x = length * ((double) j + 0.5) + sys->minx;
	    nextc->y = length * ((double) k + 0.5) + sys->miny;
	    nextc->z = length * ((double) l + 0.5) + sys->minz;
/*
	    if (nextc->chgs != NULL) 
	      printf("cube (%d,%d,%d)  center (%g, %g, %g)\n", j,k,l,nextc->x,nextc->y,nextc->z); 
*/
	    nextc->level = i;
	    nextc->j = j;
	    nextc->k = k;
	    nextc->l = l;
	  }
	}
      }
    }
/*  } */ 
}

/*
Recursive routine to give indexes to the charges so that those in each 
cube are contiguous. In addition, insure that the charges in each parent 
at each level are numbered contiguously.  This is used to support a 
psuedo-adaptive scheme.  Also get the pointer to the appropriate section
of the charge and potential vector.  Uses the eval vector for the potential
coeffs at the lowest level.  Also index the lowest level cubes.
*/
static indexkid(sys, dad, pqindex, pcindex)
ssystem *sys;
cube *dad;
int *pqindex, *pcindex;
{
  int i;
  
  if(dad != NULL) {
    if((dad->numkids == 0) && (dad->upnumvects > 0)) {
      CALLOC(dad->upvects, 1, double*, ON, ACUB);
      CALLOC(dad->nbr_is_dummy, 1, int*, ON, ACUB);
      dad->upvects[0] = &(sys->q[*pqindex]);
      dad->eval = &(sys->p[*pqindex]); /* changed from local to eval 17Feb90 */
      dad->nbr_is_dummy[0] = &(sys->is_dummy[*pqindex]);
      dad->is_dielec = &(sys->is_dielec[*pqindex]);
      dad->index = (*pcindex)++;
      for(i=0; i < dad->upnumeles[0]; i++) {
	(dad->chgs[i])->index = (*pqindex)++;
      }
    }
    else {
      for(i=0; i < dad->numkids; i++) {
	indexkid(sys, dad->kids[i], pqindex, pcindex);
      }
    }
  }
}

static indexcubes(sys, cubes, pqindex, pcindex)
ssystem *sys; 
cube *****cubes;
int *pcindex, *pqindex; 
{
  int j,k,l, side=sys->side, i=sys->depth, n; 
  cube *cp; 

  for (j=0; j<side; j++) {
    for (k=0; k< side; k++) {
      for (l=0; l<side; l++) {
        cp = cubes[i][j][k][l]; 
        if (cp != NULL) {
          if (cp->upnumeles[0]> 0) {
            CALLOC(cp->upvects, 1, double *, ON, ACUB); 
            CALLOC(cp->nbr_is_dummy, 1, int*, ON, ACUB); 
            cp->upvects[0] = &(sys->q[*pqindex]); 
            cp->eval = &(sys->p[*pqindex]); 
            cp->nbr_is_dummy[0] = &(sys->is_dummy[*pqindex]); 
            cp->is_dielec = &(sys->is_dielec[*pqindex]); 
            cp->index = (*pcindex)++; 
            for (n=0; n<cp->upnumeles[0]; n++) {
              cp->chgs[n]->index = (*pqindex)++; 
            }
          }
        }
      }
    }
  }
}

/* 
SetExact marks as exact those cubes containing fewer than numterms
number of charges.  If the number of charges in the kids is less than
numterms, the box is marked as exact and the charges are copied up.
In addition, the vector of local expansion coeffs is set to the
potential vector.  Otherwise, the number of nonzero kids is counted
and put in upnumvects as usual.  
*/
/* added 30Mar91: provisions for loc_exact and mul_exact */
setExact(sys, numterms)
ssystem *sys;
int numterms;
{
int i, j, k, l, m, n;
int side = sys->side;
int depth = sys->depth;
int numchgs, num_eval_pnts, first, multerms();
cube *nc, *nkid, *****cubes = sys->cubes;
int all_mul_exact, all_loc_exact, p, num_real_panels;

  for(i=depth; i > 0; i--, side /= 2) {
    for(j=0; j < side; j++) {
      for(k=0; k < side; k++) {
	for(l=0; l < side; l++) {
          nc = cubes[i][j][k][l];
	  if(nc != NULL) {
	    if(i == depth) {
	      ASSERT(nc->upnumvects != 0);
	      /* count the number of true panels in this cube */
	      num_real_panels = 0;
	      for(p = 0; p < nc->upnumeles[0]; p++) {
		if(!nc->chgs[p]->dummy) num_real_panels++;
	      }
	      if(num_real_panels <= numterms) {
		nc->mul_exact = TRUE;
		nc->multisize = nc->upnumeles[0];
	      }
	      else {
		nc->mul_exact = FALSE;
		nc->multisize = multerms(sys->order);
	      }
	      if(nc->upnumeles[0] <= numterms) {
		nc->loc_exact = TRUE;
		nc->localsize = nc->upnumeles[0];
	      }
	      else {
		nc->loc_exact = FALSE;
		nc->localsize = multerms(sys->order); 
	      }
	    }
	    else {  
	      /* Count the number of charges and nonempty kids. */
	      all_loc_exact = all_mul_exact = TRUE;
	      num_eval_pnts = numchgs = nc->upnumvects = 0;
	      for(m = 0; m < nc->numkids; m++) {
		nkid = nc->kids[m];
		if(nkid != NULL) {
		  nc->upnumvects += 1;
		  if(nkid->mul_exact == FALSE) all_mul_exact = FALSE;
		  else {
		    num_eval_pnts += nkid->upnumeles[0];
		    for(p = 0; p < nkid->upnumeles[0]; p++) {
		      if(!nkid->chgs[p]->dummy) numchgs++;
		    }
		  }
		  if(nkid->loc_exact == FALSE) all_loc_exact = FALSE;
		}
	      }
	      /* If all nonempty kids exact, # chgs <= # terms, mark exact, 
		 copy chgs, and promote pointers to charge and potential.  
		 Note EXPLOITS special ordering of the pot and charge vectors.
		 */
	      if(!all_mul_exact || (numchgs > numterms)) { /* multi req'd */
		nc->mul_exact = FALSE;
		nc->multisize = multerms(sys->order);
	      }
	      else if(all_mul_exact && (numchgs <= numterms)) { 
		nc->mul_exact = TRUE;
		nc->upnumvects = 1;
		CALLOC(nc->upvects, 1, double*, ON, ACUB);
		CALLOC(nc->upnumeles, 1, int, ON, ACUB);
		nc->upnumeles[0] = num_eval_pnts; /* was numchgs 30Mar91 */
		nc->multisize = num_eval_pnts; /* was numchgs */
		CALLOC(nc->chgs, num_eval_pnts, charge*, ON, ACUB);
		num_eval_pnts = 0;
		for(m=0, first=TRUE; m < nc->numkids; m++) {
		  nkid = nc->kids[m]; 
		  if(nkid != NULL) {
		    if(first == TRUE) {
		      nc->upvects[0] = nkid->upvects[0];
		      if(nc->nbr_is_dummy == NULL)
			  CALLOC(nc->nbr_is_dummy, 1, int*, ON, ACUB);
		      nc->nbr_is_dummy[0] = nkid->nbr_is_dummy[0];
		      first = FALSE;
		    }
		    for(n=0; n < nkid->upnumeles[0]; n++) {
		      nc->chgs[num_eval_pnts++] = nkid->chgs[n];
		    }
		  }
		}
	      }

	      /* do the same for local expansion */
	      /* if local exact, must be multi exact => no promotion reqd */
	      if(!all_loc_exact || (num_eval_pnts > numterms)) { /* le req'd */
		nc->loc_exact = FALSE;
		nc->localsize = multerms(sys->order);
	      }
	      else if(all_loc_exact && (num_eval_pnts <= numterms)) { 
		nc->loc_exact = TRUE;
		nc->localsize = num_eval_pnts;
	      }
	    }
	  }
	}
      }
    }
  }
}


/*
Find all the nearest neighbors.
At the bottom level, get neighbors due to a parents being exact.
*/
static getnbrs(sys)
ssystem *sys;
{
cube *nc, *np, *****cubes = sys->cubes;
int depth = sys->depth;
int i, j, k, l, m, n, p, side, es;
int numnbrs;

/* Return if depth = 0, no neighbors. */
  if(depth == 0) return;

/*
At the every level, get the nearest nbrs combined with nbrs due to parents
being exact.
*/

side = sys->side; 
i = sys->depth; 

  /* exactness for local expansion is checked - nbrs used only in dwnwd pass */
/*  for(i = 1, side = 2; i <= depth; i++, side *= 2) {  */ 
    for(j=0; j < side; j++) {
      for(k=0; k < side; k++) {
	for(l=0; l < side; l++) {
	  nc = cubes[i][j][k][l];
	  if(nc != NULL) {
	    /* Find sidelength of exact cube. */
            /*
	    for(es=1, np=nc->parent; np->loc_exact==TRUE; 
                                     np = np->parent, es *= 2); 
                                     */ 
            es = 1; 
	  
	    /* Stack up the nearest nbrs plus nbrs in exact cube. */
	    numnbrs = 0;
	    for(m = MIN((j-NNBRS), es * (j/es));
		m < MAX((j+NNBRS+1), es * (1 + (j / es))); m++) {
	      for(n = MIN((k-NNBRS), es * (k/es));
		  n < MAX((k+NNBRS+1), es * (1 + (k/es))); n++) {
		for(p = MIN((l-NNBRS), es * (l/es));
		    p < MAX((l+NNBRS+1), es * (1+(l/es))); p++) {
		  if( (m >= 0) && (n >= 0) && (p >= 0)
		     && (m < side) && (n < side) && (p < side)
		     && ((m != j) || (n != k) || (p != l))
		     && (cubes[i][m][n][p] != NULL)) {
		    cstack[numnbrs++] = cubes[i][m][n][p];
		  }
		}
	      }
	    }
	    nc->numnbrs = numnbrs;
	    if(nc->numnbrs > 0) CALLOC(nc->nbrs, numnbrs, cube*, ON, ACUB);
	    for(m=numnbrs-1; m >= 0; m--) nc->nbrs[m] = cstack[m];
	  }
	}
      }
    }
/*  } */ 
}

/*
  returns the number of charges in the lowest level cubes contained in "cp"
*/
int cntDwnwdChg(cp, depth)
int depth;			/* number of lowest level */
cube *cp;
{
  int i;
  int cnt=0;
  cube *kidc;

  if(cp->level == depth) return(cp->upnumeles[0]);
  else for(i = 0; i < cp->numkids; i++) 
      cnt += cntDwnwdChg(cp->kids[i], depth);
  return(cnt);
}

/* 
Set up the links between cubes requiring multi work on each level, one
for the cubes requiring local expansion work, one for the cubes requiring
direct methods and one for cubes with potential evaluation points. 
Note, upnumvects and exact must be set!!!
*/
static linkcubes(sys)
ssystem *sys;
{
  cube *nc, **plnc, **pdnc, **pmnc, *****cubes = sys->cubes;
  int i, j, k, l, cnt = 0;
  int dindex, side, depth=sys->depth, numterms=multerms(sys->order);

  /* Allocate the vector of heads of cubelists. */
  CALLOC(sys->multilist, sys->depth+1, cube*, ON, ACUB);
  CALLOC(sys->locallist, sys->depth+1, cube*, ON, ACUB);

  pdnc = &(sys->directlist);
  for(dindex = 1, i=0, side = 1; i <= sys->depth; i++, side *= 2) {
    pmnc = &(sys->multilist[i]);
    plnc = &(sys->locallist[i]);
    for(j=0; j < side; j++) {
      for(k=0; k < side; k++) {
	for(l=0; l < side; l++) {
	  nc = cubes[i][j][k][l];
	  if(nc != NULL) {
	    /* Do the multi expansion if the cube is not treated exactly. */
	    if(i > 1) {		/* no multis over the root cube and its kids */
	      if(nc->mul_exact == FALSE) { /* exact -> mul_exact 1Apr91 */
		CALLOC(nc->multi, numterms, double, ON, ACUB);
		*pmnc = nc;
		pmnc = &(nc->mnext);
	      }
	    }
	    
	    /* Do the local expansion on a cube if it has chgs inside, and it's
	     not exact and not the root (lev 0) nor one of its kids (lev 1). */
	    if(i > 1) {		/* no locals with level 0 or 1 */
	      if(nc->loc_exact == FALSE) { /* exact -> loc_exact 1Apr91 */
		*plnc = nc;
		plnc = &(nc->lnext);
		CALLOC(nc->local, numterms, double, ON, ACUB);
	      }
	    }

	    /* Add to direct list if at bot level and not empty. */
	    if(i == depth) { 
	      *pdnc = nc;  /* For the direct piece, note an index. */
	      pdnc = &(nc->dnext);
	      nc->dindex = dindex++;
	    }
	  }
	}
      }
    }
  }
}


static newlinkcubes(sys)
ssystem *sys;
{
  cube *nc, **plnc, **pdnc, **pmnc, *****cubes = sys->cubes;
  int i, j, k, l, cnt = 0;
  int dindex, side, depth=sys->depth, numterms=multerms(sys->order);

  /* Allocate the vector of heads of cubelists. */
  CALLOC(sys->multilist, sys->depth+1, cube*, ON, ACUB);
  CALLOC(sys->locallist, sys->depth+1, cube*, ON, ACUB);

  pdnc = &(sys->directlist);
  i = sys->depth; side = sys->side; 
  dindex = 1; 
  for(j=0; j < side; j++) {
    for(k=0; k < side; k++) {
      for(l=0; l < side; l++) {
        nc = cubes[i][j][k][l];
	  if(nc != NULL) {
	    /* Add to direct list if at bot level and not empty. */
	      *pdnc = nc;  /* For the direct piece, note an index. */
	      pdnc = &(nc->dnext);
	      nc->dindex = dindex++;
          }
      }
    }
  }
}


/*
Determine maximum number of chgs contained in a single cube.
*/
static setMaxq(sys)
ssystem *sys;
{
  int i, j, k, l, side, p, kids_are_exact, all_null, depth = sys->depth;
  int mul_maxq, mul_maxlq, loc_maxq, loc_maxlq, num_chgs, real_panel_cnt;
  cube *nc, *****cubes = sys->cubes;

  mul_maxq = mul_maxlq = loc_maxq = loc_maxlq = 0;
  for(i = 1, side = 2; i <= depth; i++, side *= 2) {
    for(j=0; j < side; j++) {
      for(k=0; k < side; k++) {
	for(l=0; l < side; l++) {
	  nc = cubes[i][j][k][l];
	  if(nc != NULL) {
	    if(nc->mul_exact == TRUE) {
	      num_chgs = 0;
	      for(p = 0; p < nc->upnumeles[0]; p++) {
		if(!nc->nbr_is_dummy[0][p]) num_chgs++;
	      }
	      mul_maxq = MAX(mul_maxq, num_chgs);
	      if(i == depth) mul_maxlq = MAX(mul_maxlq, num_chgs); 
	    }
	    if(nc->loc_exact == TRUE) {
	      loc_maxq = MAX(loc_maxq, nc->upnumeles[0]);
	      if(i == depth) loc_maxlq = MAX(loc_maxlq,nc->upnumeles[0]); 
	    }
	  }
	}
      }
    }
  }
  sys->loc_maxq = loc_maxq;	/* max evaluation points, over all loc_exact */
  sys->loc_maxlq = loc_maxlq;	/* max evaluation pnts, over lowest level */
  sys->mul_maxq = mul_maxq;	/* max panels, over all mul_exact cubes */
  sys->mul_maxlq = mul_maxlq;	/* max panels, over lowest level cubes */

  /* find the maximum #panels in all non-exact cubes w/exact (or no) kids */
  sys->max_panel = 0;
  for(j = 2; j <= depth; j++) {
    for(nc = sys->multilist[j]; nc != NULL; nc = nc->mnext) {
      if(nc->level == depth) {
	real_panel_cnt = 0;
	for(i = 0; i < nc->upnumeles[0]; i++) {
	  if(!nc->nbr_is_dummy[0][i]) real_panel_cnt++;
	}
	sys->max_panel = MAX(sys->max_panel, real_panel_cnt);
      }
      else {
	kids_are_exact = all_null = TRUE;
	real_panel_cnt = 0;
	for(i = 0; i < nc->numkids && kids_are_exact; i++) {
	  if(nc->kids[i] == NULL) continue;
	  all_null = FALSE;
	  if(!nc->kids[i]->mul_exact) kids_are_exact = FALSE;
	  else {		/* count real panels */
	    for(l = 0; l < nc->kids[i]->upnumeles[0]; l++) {
	      if(!((nc->kids[i]->nbr_is_dummy[0])[l])) real_panel_cnt++;
	    }
	  }
	}
	if(kids_are_exact && !all_null) {
	  sys->max_panel = MAX(sys->max_panel, real_panel_cnt);
	}
      }
    }
  }

  /* find the maximum #eval points in all non-exact cubes w/exact children */
  sys->max_eval_pnt = 0;
  for(j = 2; j <= depth; j++) {
    for(nc = sys->locallist[j]; nc != NULL; nc = nc->lnext) {
      if(nc->level == depth) {
	sys->max_eval_pnt = MAX(sys->max_eval_pnt, nc->upnumeles[0]);
      }
      else {
	kids_are_exact = all_null = TRUE;
	real_panel_cnt = 0;
	for(i = 0; i < nc->numkids && kids_are_exact; i++) {
	  if(nc->kids[i] == NULL) continue;
	  all_null = FALSE;
	  if(!nc->kids[i]->loc_exact) kids_are_exact = FALSE;
	  else real_panel_cnt += nc->kids[i]->upnumeles[0];
	}
      }
      if(kids_are_exact && !all_null)
	  sys->max_eval_pnt = MAX(sys->max_eval_pnt, real_panel_cnt);
    }
  }
}

/* 
  markup sets the flag to "flag" in the child and its nearest nbrs
*/
static markUp(child, flag)
cube *child;
int flag;
{
  int i,j;
  cube *nc, *np;

  child->flag = flag;
  for(i = 0; i < child->numnbrs; i++) {
    child->nbrs[i]->flag = flag;
  }
}

/* 
  forms the true interaction list (see also comment at mulMatEval())
   for cube "child", excluding only empty cubes
  -interaction list pointer is saved in the interList cube struct field
*/
static getInter(child)
cube *child;
{
  int i, j, vects, usekids, lc, jc, kc, ln, jn, kn;
  int numnbr = (child->parent)->numnbrs; /* number of neighbors */
  cube **nbrc = (child->parent)->nbrs; /* list of neighbor pointers */
  cube *sib;			/* pointer to sibling (same level as child) */
  cube **pstack = &(cstack[0]); /* temporary storage pointer */

  /* mark the child cube and all its neighbors */
  markUp(child, TRUE);

  /* unmarked children of child's parent's neighbors become the ilist */
  for(i = 0; i < numnbr; i++) { /* loop on neighbors */
    /* Check nbr's kids for a marked kid. */
    for(usekids = FALSE, j = 0; j < nbrc[i]->numkids; j++) { 
      sib = (nbrc[i]->kids)[j];
      if((sib != NULL) && (sib->flag == TRUE)) { usekids = TRUE; break; };
    }
    /* Use nbr if no kids marked. */
    /* ...and it's really not a 1st nrst nbr of the parent 
       - this stops parent-sized cubes from getting into the ilist
         when they have empty child-sized cubes that are 2nd or 1st
	 nrst nbrs of the child cube 
       - should work with NNBRS = 1 (never allows parent-sized in list)
         and NNBRS > 2 (but cannot allow greater than parent-sized)
       (29May90) */
#if ON == ON
    lc = (child->parent)->l;
    jc = (child->parent)->j;
    kc = (child->parent)->k;
    ln = nbrc[i]->l;
    jn = nbrc[i]->j;
    kn = nbrc[i]->k;
    if((RADINTER == ON) && (usekids == FALSE) &&
       ((lc-1 != ln && lc+1 != ln && lc != ln)
       || (jc-1 != jn && jc+1 != jn && jc != jn)
       || (kc-1 != kn && kc+1 != kn && kc != kn))) {  
      *pstack = nbrc[i];
      pstack++;
    }
#else				/* USE THIS PART FOR TESTING ONLY */
    if(RADINTER && (usekids == FALSE)) { /* PRODUCES INCORRECT ILISTS!!! */
      *pstack = nbrc[i];
      pstack++;
    }
#endif
    else for(j = 0; j < nbrc[i]->numkids; j++) { /* use nbr's kids. */
      sib = (nbrc[i]->kids)[j];	/* get sib of child cube of interest */
      if((sib != NULL) && (sib->flag == FALSE)) { 
	*pstack = sib;
	pstack++;
      }
    }
  }

  /* clear all the flags */
  markUp(child, FALSE);

  /* allocate and save the interaction list */
  child->interSize = vects = pstack - &(cstack[0]);
  if(vects > 0) CALLOC(child->interList, vects, cube*, ON, ACUB);
  for(j = 0; j < vects; j++) child->interList[j] = cstack[j];

  return(vects);		/* return number of interaction elements */
}

/*
  generates explicit, true interaction lists for all non-empty cubes w/lev > 1
*/
static getAllInter(sys)
ssystem *sys;
{
  int i, j, k, l, side, depth = sys->depth;
  cube *nc, *****cubes = sys->cubes;
  for(i = 2, side = 4; i <= depth; i++, side *= 2) {
    for(j=0; j < side; j++) {	/* loop through all cubes at levels > 1 */
      for(k=0; k < side; k++) {
	for(l=0; l < side; l++) {
	  nc = cubes[i][j][k][l];
	  if(nc != NULL) getInter(nc);
	}
      }
    }
  }
}

/* 
  inits the dummy and dielec mask vectors; used to tell which entries to skip
  - mask vectors are redundant (could read flags in charge struct) 
  - done for speed in potential eval loop
*/
static set_vector_masks(sys)
ssystem *sys;
{
  int i;
  cube *cp;

  for(cp = sys->directlist; cp != NULL; cp = cp->dnext) {
    for(i = 0; i < cp->upnumeles[0]; i++) {
      if(!(cp->nbr_is_dummy[0][i] = cp->chgs[i]->dummy))
	  cp->is_dielec[i] = 
	      (cp->chgs[i]->surf->type == DIELEC 
	       || cp->chgs[i]->surf->type == BOTH);
    }
  }

}





