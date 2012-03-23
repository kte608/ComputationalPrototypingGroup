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
/* FILE: Elements.h
*
*  This file contains macros, type definitions and structures 
*  needed for element routines.
*
*/

/* 
====================================================================== 
  STRUCTURES
====================================================================== */ 

/* All information about what a source/eval element is. 
*  This thing will probably be redefined soon enough. 
*  Right now a "source" just contains a center and a 
*  radius corresponding to a bounding sphere. 
*  Obviously, later this will be removed and other 
*  information will be added, from which the 
*  center and radius can be calculated. */
struct dummyElement {  
  point center;
  double radius;
  int directIndex;
};

struct constantElement {
  int shape;         /* 3: Triangular, 4: Quadrilateral.               */
  point corners[4];  /* Corner point coordinates: P1 P2 P3 [P4]        */
  point centroid;    /* Centroid of element. This point may be used as 
		      * the collocation point.                         */
  int directIndex;   /* Index in the direct matrix. This can also be 
		      * known as the "internal" index. (It may differ 
		      * from the index used in the input file.)        */
  int crackedNeighbourIndex; /* If a input quadrilateral element does not
			      * lie in a plane, it may be "cracked" into
			      * two triangular elements along its 
			      * shortest diagonal. For output purposes
			      * it may be advantageous to store the 
			      * index of each neighbour - such that an 
			      * "averaging" or so can be made.
			      * If a panel is not cracked this index 
			      * equals the "internal" index.          */
  short int bcType;   /* Type of boundary condition on the boundary 
		       * where the element is located.
		       * 0: Neumann boundary condition
		       * 1: Diriclet boundary condition
		       * other ??                                      */
  short int boundaryNumber; /* The piece of boundary where the 
			     * present element belongs.
			     * Equivalent to "conducter number" in 
			     * conductance computations.               */
};

/* == Maybe an element should be as the following structure.
 * == Then, based on the value of elementType, the generic 
 * == pointer could be recast to the proper element type?? 
 *  Probably not a good idea? Rather use a more general 
 *  structure, and give e.g. the order of the element. 
 *  The size of each element can still vary without having to 
 *  recast parts of it.                                           */
struct elementIdea {
  int elementType;   /* Determines the how to interpret the rest,
		      * i.e. to what form to recast "elementData" */
  void *elementData; /* Pointer to the data of the element. 
		      * Note that the data could take one of 
		      * many forms, which differ in the choice 
		      * of geometry description, basis function 
		      * type (if any/could be Dirac delta function?)... */

  /* Maybe it would be better to use one type (int + void*) for 
   * the geometry and one for the basis? Then it would be easier 
   * to separate these two things and they could be handled by
   * different, each knowing only little about the other...? */
};




