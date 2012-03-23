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
/* FILE: GridBundles.h
 * Structures that are used by grid functions when 
 * passing data around (possibly to other types of 
 * functions as well). */

/* Bundle used in sourcesAndNeighbourOperate */
struct gridBundle_SANO_1 {
  void *origData; /* Original data for the target routine */
  double ***nearFieldG2G; /* Small grid-to-grid matrices for near-field 
			   * computations (and precorrection!)         */
  int *g2gIndex;          /* Index for element association with grid 
			   * corresponding to which G2G matrix to use. */
  int numG2G;             /* Number of grid-to-grid matrices           */
  int G2Gn;               /* # rows in each small G2G matrix           */
  int G2Gm;               /* # columns in each small G2G matrix        */
  /* The following entries are for the use of "compressed precomputed"
   * interpolation points */
  double **nearFieldUnionG2G; /* One skinny grid-to-grid matrix 
			       * for all commonly encountered pairs 
			       * of projection and interpolation points. */
  int nRowsNearFieldUnion;    /* Numer of rows in nearFieldUnionG2G */
  int **nearFieldUnionG2Gindices; /* For each interpolation point 
				   * [idirect][iinterp] an index in the 
				   * compressed set of points */
  double *unionWork; /* Work vector of length nRowsNearFieldUnion */
};
