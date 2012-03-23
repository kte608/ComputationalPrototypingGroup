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
  First attempt to define structures for use in the new pFFT code. 
  It seems like structure must be declared before they can be used
  (even within the same file).
*/

/* Structure for storing information on ordering of panels 
   within the grid. The information is primarily of use when 
   building the sparse matrix. */
struct panelOnGrid {
  int ***nPanelsOnPoint; /*  # panels associated with each grid point 
			     Size: Number of grid points: nx*ny*nz */

  int ****iPanelsOnPoint; /* (Pointer to) an array with the numbers of 
			     the panels associated with each grid point.
			     Size: Number of grid points: nx*ny*nz.
			     The arrays are stored together in 
			     panelNumbers, see below. */

  int *panelNumbers;      /* List of panel numbers sorted by grid point
			     association. Size: Number of panels */
};

/* The grid is for the FFT part. It is a 3D grid, but since it
   is oriented in the cardinal directions of the coordinate system, 
   the coordinates of each grid point need not to be stored. */
struct grid {
  int nx,ny,nz;    /* Number of nodes in each cardinal direction */

  double dx,dy,dz; /* Mesh size in each cardinal direction */

  double xmin,xmax,ymin,ymax,zmin,zmax;  /* Extend of mesh */

  double *x,*y,*z; /* Vector of coordinates in each cardinal direction 
		      *x,*y and *z should point to arrays of size nx,  
		      ny and nz respectively.Note that hopefully: 
		      xmin = *x+1, xmax = *x+nx, etc. */

  struct panelOnGrid panels; /* Panel associations, see above. */
};




/* ==== WHATEVER IS BELOW THIS LINE IS NOT USED ==== */

/* For each grid point make a list of panels associated 
   with the point. */
struct panelList {
  int n;            /* Number of panels in the list */
  int *panelNumber; /* List (vector) of the panel numbers 0 ... n-1 */
};

/* We don't need (I guess) a list of neighbouring grid points.
   Rather, for each panel there will be a list of neighbouring 
   panels, and on each grid point there will be a list of panels
   associated with the grid point. */



/* For each grid point there needs to be a neighbour-list of
   adjacent grid points, so that the panels on these points can 
   be found later??. The easiest way to do this indexing may be \
   by using a "neighbour-shape" or a "default list" such that 
   the structure of the grid may be exploited:
   Thus, if, say, i1,j1,k1+1 is a neighbour of i1,j1,k1 
             then i2,j2,k2+1 is a neighbour of i2,j2,k2. */

