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
/* FILE: Grid.h
*
*  Structures for grid data.
*
*/

/* Only include this file if it hasn't been included before: */
#if ! defined GridDotHIsIncluded
#define GridDotHIsIncluded

/* 
====================================================================== 
  STRUCTURES
====================================================================== */

/* This should be a data structure with info for each grid point.
*  The grid will have an [i][j][k] instance of this to store all 
*  the points. */
struct gridPoint {
  int nAssociatedElements[2];  /* Number of associated grid elements
				* 0: Sources; 
				* 1: Evals;
				*/
  int *associatedElements[2];  /* Lists of associated grid elements. */
};
typedef struct gridPoint gridPoint;


/* When finding neighbouring elements the number of 
*  elements in each list is initially unknown, and it 
*  can't easily be counted, since the planned double-sweep
*  may introduce multiple entries. Thus, temporary storage 
*  can be employed to keep the information until it
*  can be stored permanently. */
struct neighbourLinkedList {
  /* Number on a neighbouring element: */
  int iNeighbour; 
  /* Next element in the list:         */
  struct neighbourLinkedList *nextNeighbour; 
};

/* This is an element for the grid. It should contain all 
*  the information the grid needs to know about an element 
*  as well as a pointer (type void) to wherever more 
*  information is located. The grid routines should never 
*  need to access this extra information directly. */
struct gridElement {
  point center;    /* Center of bounding sphere */

  double radius;   /* Radius of bounding sphere */


  /* Level of the projection / interpolation stencil when used 
   * for the present element. In most cases this should be one, 
   * corresponding to unity*DX as the basis of the stencil points.
   * If an element is large, then the projection / interpolation
   * may be made using a "larger" size (but not more points) 
   * stencil, using projectLevel*DX as the basis. 
   * Maybe this should be part of the stencil rather than the element? */
  unsigned short int projectLevel; 

  /*int nNeighbours;  Number of element neighbours 
		    * (of opposite type) 
		    * This integer is not REALLY required 
		    * if the list is made as a linked list, 
		    * since the list is terminated with a 
		    * NULL pointer. */

  /* Starting entry of list of element neighbours 
  struct neighbourLinkedList *neighbourLinkedList;  */ 

  void *element;   /* This is where the remaining info is stored, 
		    * i.e. all the info the grid doesn't need explicitly, 
		    * such as the specific geometry and basis function 
		    * descriptions. */
};

/* This is for a grid-element a list of those neighbours, 
*  which cannot be reached "through the grid" 
*  (i.e. using neighbouring grid points).                              */
struct gridElementNeighbours {
    int nNeighbours; /* Number of element neighbours 
		    * (of opposite type) 
		    * This integer is not REALLY required 
		    * if the list is made as a linked list, 
		    * since the list is terminated with a 
		    * NULL pointer. */

  /* Starting entry of list of element neighbours */
  struct neighbourLinkedList *neighbourLinkedList;
};

/* This is for grid "stencils", used for defining neighbouring
*  grid points. For each stencil a separate set of close-by neighbours
*  is needed for each element. This set is obtained using the 
*  separationDistanceFactor. Grid-elements that are separated by less 
*  than separationDistanceFactor*radius are considered "close", 
*  no matter the size of the grid stencil.                             
*  This data structure also serves to store the interpolation/projection.
*/
struct gridStencil {
  int nStencilPoints;   /* Number of neighbouring grid points for 
		         * points well inside grid (i.e. the length of 
			 * the vectors ipNeighbours, jpNeighbours and 
			 * kpNeighbours).                              */
  int *ipNeighbours,    /* Lists of indices for neighbouring points    */
    *jpNeighbours,      /* based on the i,j,k index for the present    */
    *kpNeighbours;      /* grid point                                  */

  /* The following part is for projection / interpolation stencils     */
  /* int terms;	*/	/* number of poly terms. */
  int nBasis;           /* Number of basis functions (poly terms)      */
  double **bValues;	/* The values of the interpolation/projection 
			 * basis functions in each of the stencil points.
			 * Matrix is nStencilPoints x nBasis.
			 * At the time of writing this matrix is not 
			 * used, but the storage needed is very little 
			 * and it could come in handy later...         */
  double **weights;	/* The interpolation weights. 
			 * Relates grid points to polynomial terms. 
			 * Matrix is nBasis x nStencilPoints and is the 
			 * (pseudo-)inverse of bValues.                */
  int *imonoOrder,      /* Orders the monomials (x^i * y^j * z^k)      */
      *jmonoOrder,      /* used in the interpolation / projection      */
      *kmonoOrder;      /* basis.                                      */
  /* There is a need for scaling of the basis functions (monomials). 
   * If the monomials are not scaled, then the results become inaccurate 
   * as dx becomes small (or large) due to poor scaling of the columns 
   * of the "values" matrix. Thus, the monomials are written as powers
   * of the scaled variables (xs,ys,zs)=polyScale*(x-x0,y-y0,z-z0). 
   * Here (x0,y0,z0) is the origin of the interpolation (projection) 
   * stencil. Typically, polyScale=1/dx will be a good choice. In this 
   * case "values" contains only zero and +/-unity if dx=dy=dz. Obviously, 
   * this is the optimal scaling of the matrix. 
   * Note, that the scale must be used every time the basis functions 
   * are evaluated and (equivalently) element moments are found.       */
  double polyScale;
  /* The order of the polynomial function (3D) and the maximum 
   * order in any one direction of the monomials it consistes of.
   * If a n'th order consistent polynomial is used then 
   *   polyOrder = maxMonoOrder1D = n
   * However, if a tensor product of 1D monomials is used, then
   * maxMonoOrder1D = n, but polyOrder = 3*n                           */
  int polyOrder, maxMonoOrder1D;
  /* Note: If a polynomial basis (rather than a basis of monomials) 
   * is needed at a later time, then a matrix can be introduced, 
   * which relates the monomials to the polynomial basis. Basically the 
   * entries of this matrix will be the coefficients of each polynomial,
   * ordered according to (imonoOrder,jmonoOrder,kmonoOrder). 
   * At the time of writing only monomial basis functions are 
   * considered, so the handling of this matrix is seen as a 
   * hindrance more than a benefit. As a consequence the matrix 
   * is left for future implementation                                 */

  /* The following part is for direct / precondition stencils 
   * (these stencils need grid-element neighbours to be found, and 
   * interaction ranges etc.)                                          */

  /* Number to multiply by radius of bounding sphere to get a 
   * condition of neighbouring elements, i.e. two elements are 
   * neighbours if they do not have the same type and their 
   * bounding spheres are separated by less than separationDistanceFactor
   * times the maximum of the radii of their bounding spheres.         */
  double separationDistanceFactor;

  /* Arrays of grid element neighbours. 
   * One array for the sources [0] and one for the evals [1]. 
   * Note, that if the sources are a subset of the evals (or vise versa) 
   * such that they are defined by the same pointer, but differ in 
   * number, then two distinct neighbourlists are still needed.        */
  struct gridElementNeighbours **gridElementNeighbours[2];
  /* Flag for showing when the grid element neighbours have been found. 
   * 0: Not found 
   * 1: Found (still there may be no neighbours!)                      */
  int gridElementNeighboursAreFound;
  /* Flag for showing when the number of grid element neighbours have 
   * been found. 
   * 0: Not found 
   * 1: Estimated 
   * 2: Found (still there may be no neighbours!)                      */
  int nGridElementNeighboursAreFound;
  /* Number of neighbours (counted or estimated)                       */
  int nGridElementNeighbours;

  /* Minimum distance between two grid points that are *not* neighbours,
   * i.e. the difference of their indices is not in the set of 
   * (ipNeighbours,jpNeighbours,kpNeighbours) triplets). 
   * This value may be set to zero, but a higher value will be more 
   * efficient (still it should be no larger than the smallest possible
   * distance between two grid points that are not neighbours).       */  
  double minNotNeighbourDist; 
};

/* This structure is used for compressing the "most likely" 
 * interpolation points to subsequently speed up the precorrection 
 * computations.                                                  */
struct localGridIndexLink{
  int compressedIndex; /* Row index in the "compressed" matrix 
			* (see nearGrid2gridUnion above)          */
  int ip,jp,kp;        /* Grid indices relative to the centre of
			* a projection stencil                    */
  struct localGridIndexLink *next;
};

/* The grid is for the FFT part. It is a 3D grid, but since it
*  is oriented in the cardinal directions of the coordinate system, 
*  the coordinates of each grid point need not to be stored. */
struct grid {
  int nx,ny,nz;    /* Number of nodes in each cardinal direction */

  int p2x,p2y,p2z; /* nx = pow(2,p2x) etc. [IF nx,ny,nz are powers of two] */

  double dx,dy,dz; /* Mesh size in each cardinal direction */

  /* Extent of mesh. The MAXes can be found from the MINs, DELTAs 
   * and n's, but it may come in handy to store them. 
   * Storage is very small, also.                                      */
  double xmin,xmax,ymin,ymax,zmin,zmax; 

  /* For efficiency reasons the grid should always be a "tight fit" of 
   * the geometry in one of the three cardinal directions (possibly 
   * adding some extra layers of points for interpolation/projection 
   * reasons). The direction (0,1,2) of a tight fit is stored in 
   * iDirTightFit.                                                     */
  int iDirTightFit;

  /* Vector of coordinates in each cardinal direction. *x,*y and *z 
   * should point to arrays of size nx, ny and nz respectively.
   * Note that hopefully: xmin = *x+1, xmax = *x+nx, etc. These values 
   * could be calculated on the fly, but they are very cheap to store 
   * O(nx+ny+nz)<<O(N), so calculate once and reuse when necessary.    */
  double *x,*y,*z; 

  /* Needed stencils, one for the direct part and one for the 
   * preconditioning                                                   */
  /* JKW: project/interp stencil computed from subset of the direct stencil. */
  /* BB: Eventually I want stencils for interpolation and for projection 
   * to be independent. I want the freedom to choose these methods 
   * independently to minimize the computational efforts. 
   * Thus, elementGridStencil will be replaced by interpolationStencil 
   * and projectionStencil. */
  struct gridStencil *directStencil, *precondStencil, *elementGridStencil;

  gridPoint ****points; /* Array (3D) of pointers to grid points       */

  /* Sometimes (well, maybe rather often) the sources 
   * and evals will be the same, i.e. be defined by the 
   * same pointer and be of equal number.
   * In this case memory and CPU time may be reduced by
   * exploiting this fact. 
   * This flag is true (1) if there is two separate arrays
   * (or different numbers), false (0) if not */
  int sourceAndEvalDiffer;

  /* Arrays of grid elements. 
   * One array for the sources [0] and one for the evals [1]. 
   * These arrays contain the information needed for the grid, 
   * i.e. information on bounding sphere and neighbours. 
   * Note, that if the sources are a subset of the evals 
   * (or vise versa) such that they are defined by the same 
   * pointer, but differ in number, then two distinct 
   * neighbourlists are still needed. */
  struct gridElement **gridElements[2]; 

  /* Number of grid elements of each type [0], [1]. */
  int nGridElements[2];

  /* List of grid elements (sources and evals) ordered by 
   * grid point. 
   * Each grid point will contain a pointer to these 
   * arrays for proper access to the list.  */
  int *associatedElements[2]; 

  /* It may be convenient to keep a little statistical information 
   * around regarding the sizes of the grid-elements. For instance, 
   * it can be useful when determining the grid size to use. Note that 
   * these measures may not be very accurate - a fast summation will 
   * be used (accuracy might requirte an initial sorting of the 
   * elements based on radius). Store the sum of the radii and the 
   * sum of the squares of the radii                                   */
  double elementRadiusSum[2], elementRadiusSquareSum[2];

  /* Storage for the "kernel" (Green function) data "on the grid". 
   * Basically, the values of G(x-y) for all possible combinations 
   * of x of y. Typically values of 1/r with r=|x-y|. 
   * Note that, to increase the flexibility of the code, "the grid 
   * does not know" what the structure of the kernel is. Thus, a 
   * generic type pointer is used. 
   * At the time of writing, the kernel is laid out for FFT use 
   * (in "nfft.c"), but that could be changed at a later point.        */
  void *kernel;

  /* Storage for the grid data. Starts as projected sources, ends as
   * to be potentials on the grid. The potentials should then be 
   * interpolated from the grid to wherever they are needed. Note that 
   * the grid does not know the actual structure of the grid data.     */
  void *griddata;  

/* The following part is an attempt to accellerate the precorrection
 * step. Information is need for projection, interpolation AND direct 
 * stencils. For each grid point in the direct stencil make a small 
 * matrix with the grid-to-grid interactions between the projection 
 * stencil centered at the central node and the interpolation stencil 
 * centered at this particular direct-stencil grid point. 
 * The precorrection consists of calculating PROJ*G2G*INTERP and 
 * subtracting the sum (over projection and interpolation points) 
 * from each matrix entry. Basically this may be written as
 *
 *        \sum_i \sum_j P_i G_{ij} I_j  = P*G*I
 *
 * or, equivalently as:
 *
 *        \sum_j (I_j ( \sum_i P_i G_{ij} ) ) 
 * 
 * For all elements, which are neighbours "through the grid" the small
 * matrix G_{ij} is thus precomputed for the grid-to-grid interactions.
 * Large elements, which have their own (linked) interaction lists 
 * need special treatment for elements in those lists.                 */
  double ***nearGrid2grid;
  int numNearGrid2grid, /* Number of these small matrices, should equal 
			 * nStencilPoints for the direct stencil).     */
      nNearGrid2grid,   /* Number of rows in each small matrix, should 
			 * equal nStencilPoints for the PROJ stencil   */
      mNearGrid2grid;   /* Number of rows in each small matrix, should 
			 * equal nStencilPoints for the INTERP stencil */
  /* As an alternative the union of the most likely interpolation 
   * points may be found, and one (larger) local grid2grid matrix be 
   * set up (one column for each point in the projection stencil). 
   * Then G*P is calculated for a larger (maybe too large) G and 
   * I*(G*P) is performed for each I by chosing the "right" entries in 
   * G*P. */
  double **nearGrid2gridUnion; /* Union g2g matrix */
  int nRowsNearGrid2gridUnion; /* # rows in nearGrid2gridUnion = Union # INTERP */
  int nColsNearGrid2gridUnion; /* # columns in nearGrid2gridUnion = #PROJ     */
  double *nearGrid2gridUnionWork; /* Work vector (length nRowsNearGrid2gridUnion) 
				   * for storing G*P */
  int **g2gUnionRowNumbers;    /* g2gUnionRowNumbers[i][j] is the row 
				* entry (in nearGrid2gridUnion) for the 
				* j'th interpolation point at the i'th 
				* direct point                         */
  int nRowsG2gUnionRowNumbers; /* # rows in g2gUnionRowNumbers = # DIRECT */
  int nColsG2gUnionRowNumbers; /* # cols in g2gUnionRowNumbers = # INTERP */
  struct localGridIndexLink *localGridIndexLink;
};

/* 
================================================================= 
  MACROS
================================================================= */ 
# define gridPointExists(i,j,k,grid) ( \
     (i)>=0 && (i)<(grid->nx) &&         \
     (j)>=0 && (j)<(grid->ny) &&         \
     (k)>=0 && (k)<(grid->nz)          )



/* This is the last statement to make sure that 
 * all of the above is only included once */
#endif
/* ====== DONT PUT ANYTHING BELOW THIS LINE ====== */
