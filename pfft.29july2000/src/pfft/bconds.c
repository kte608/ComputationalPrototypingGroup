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
/* This file is for routines that set up boundary condition 
 * types and values. These routines may be dependent on the 
 * specific problem being solved. 
 * Hopefully very little about the rest of the code needs 
 * to be known by these routines! */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <float.h>

#include "Global.h"

/* 
* ====================================================================
* Definitions of boundary condition types (Neuman or Diriclet) based 
* on the coordinate of the boundary point and possibly the number 
* of the boundary part that the point lies on.
*
* This routine is a companion to "boundaryConditionValue".
* "conditionType" is described in detail in that routine.
*
* Return value is 
*  0 for a Neumann boundary
*  1 for a Diriclet boundary.
* ==================================================================== */ 
int boundaryConditionType(point coordinate, 
			  int iboundary,
			  char *conditionType,
			  int idata){
  char fctName[]="boundaryConditionType";
  int itype;
  double epsi=1.0e-6;

  if (strcmp(conditionType,"conductor i") == 0){
    /* Set potential to 1 on iboundary i (conductor i)
     * and zero everywhere else. i is given in idata.                  */
    /* Always Dirichlet boundary conditions */
    itype = 1;
  }
  else if (strcmp(conditionType,"unity") == 0){
    /* Give unity to everything */
    /* Let the condition TYPE (depend on value of idata): */
    if (idata==0) itype = 0;
    else          itype = 1;
  }
  else if (strcmp(conditionType,"uniform") == 0){
    /* Use uniform flow (phi = U dot X) in a specified direction U.  */
    /* Let the condition TYPE (depend on value of idata): */
    if (idata==0)      itype = 0;
    else if (idata==1) itype = 1;
    else {
    /* Let a positive (or zero) value of the z coordinate denote a 
     * Dirichlet boundary. (Mixed boundary condition system) */
      if (coordinate[2]>-epsi) itype = 1;
      else                     itype = 0;
    }
  }
  else if (strcmp(conditionType,"constant") == 0){
    /* Let the condition TYPE (depend on value of idata): */
    if (idata==0)      itype = 0;
    else if (idata==1) itype = 1;
    else {
    /* Let a positive (or zero) value of the z coordinate denote a 
     * Dirichlet boundary. (Mixed boundary condition system) */
      if (coordinate[2]>-epsi) itype = 1;
      else                     itype = 0;
    }
  }
  else if (strcmp(conditionType,"translating sphere") == 0){
    /* Sphere (unit radius, center at origin) translating in deep 
     * water. Basically this is an external Neumann problem, but it can 
     * be used to test also other problems (exterior Diriclet or mixed).
     * At the moment hard coded to Neumann:                            */
    itype = 0;
  }
  else if (strcmp(conditionType,"sphere charge") == 0){
    /* Nonuniform potential on a sphere of unit radius centered on 
     * the origin */
    /* This isa Diriclet problem.                                      */
    itype = 1;
  }
  else {
    /* This is unknown = bad. Complain and exit!                       */
    fprintf(stderr,"%s: Error. Unknown boundary condition type: \"%s\"\n",
	    fctName,conditionType);
    _EXIT_;
  }
  return itype;
} /* End of routine boundaryConditionType */

/* 
* ====================================================================
* Definitions of boundary condition values based on the boundary 
* condition type, the on the coordinate of the boundary point and 
* possibly the number of the boundary part that the point 
* lies on [to be defined]. For Neumann conditions the normal direction 
* of the point is also needed.
*
* iData and nData are here as a means of communicating other stuff to 
* the boundaryConditions.
*
* conditionType and idata are included as means to communicate the 
* boundary condition types. 
* Presently a few simple cases have been implemented.
*
*  conditionType:
*    "unity"    - value unity everywhere
*    "uniform"  - Uniform potential flow. Phi=U*x hard coded U.
*    "constant" - unity on Dirichlet boundaries 
*                 and zero on Neumann boundaries.
*    "conductor i" - Set potential to unity on conductor 
*                 (boundary part) #i and to zero everywhere else.
*                 This may be used to find mutual capacitances,
*                 i.e. entries in "the capacitance matrix".
*                 The full capacitance matrix may be obtained by 
*                 repeating this process for every conducter i, and 
*                 at ach step integrating the obtained solution 
*                 (the charge) over each conducter j.
*                 i = idata.
*
* Special cases for a spherical boundary. For these cases analytical 
* solutions exist, so the accuracy of the method can be controlled by 
* comparing to the analytical results.
*    "translating sphere" - Hydrodynamics: sphere translating in deep 
*                 water. Disturbance potential only. Sphere of unit 
*                 radius, centered on the origin with inwards pointing 
*                 normal (out of the fluid). Note that for this case 
*                 x = -n (n being the normal and x the position vector).
*                 PHI  = x dot U/2 = -n dot U/2
*                 PHIn = -n dot U  =  x dot U     (?)
*    "sphere charge" - Electrostatics: Sphere given an uneven potential 
*                 (or charge) distribution. See e.g. FastCap or FFTcap 
*                 papers or Phillips (1997) PhD thesis. The case is *not* 
*                 one that is "realizable", i.e. it is an artificially
*                 constructed problem/solution.
*                 PHI  =   cos(theta) /2
*                 PHIn = 3*cos(theta) /(8*PI)    (=charge)
* ==================================================================== */ 
double boundaryConditionValue(int bcType, 
			      point coordinate, double *normal, 
			      int iboundary,
			      char *conditionType,
			      int idata){
  char fctName[]="boundaryConditionValue";
  double value;
  double U[3], epsi2, epsi2b;

  /* For later testing: */
  assert(value = DBL_MAX > 0.0);

  /* Check that the boundary condition type is of a known type: */
  if (bcType!=0 && bcType!=1){
    fprintf(stderr,"%s: Error. Unknown boundary condition type: %d\n"
	    "   Expecting 0 or 1\n",
	    fctName,bcType);
    _EXIT_;
  }

  if (strcmp(conditionType,"conductor i") == 0){
    /* Set potential to 1 on iboundary i (conductor i)
     * and zero everywhere else. i is given in idata.                  */
    if (iboundary==idata) value = 1.0; /* i'th conductor               */
    else value = 0.0;                  /* Everything else              */
  }
  else if (strcmp(conditionType,"unity") == 0){
    /* Give unity to everything (on a particular boundary?) */
    value = 1.0;
  }
  else if (strcmp(conditionType,"uniform") == 0){
    /* Use uniform flow (phi = U dot X) in a specified direction U.  */
    U[0] = 1.0;
    U[1] = 0.0;
    U[2] = 0.0;
    
    if (bcType==1){ /* Diriclet condition */
      value = 
	U[0]*coordinate[0]+
	U[1]*coordinate[1]+
	U[2]*coordinate[2];
    }
    else { /* (bcType==0) - Neumann condition */
      value = 
	U[0]*normal[0]+
	U[1]*normal[1]+
	U[2]*normal[2];
    }
  }
  else if (strcmp(conditionType,"constant") == 0){
    /* Use constant potential and zero derivative on boundaries.       */
    if (bcType==1) value = 1.0; /* Dirichlet condition                 */
    else value = 0.0;           /* Neumann condition (bcType==0)       */
  }
  else if (strcmp(conditionType,"translating sphere") == 0){
    /* Sphere (unit radius, center at origin) translating in deep 
     * water. U=(1,0,0) */
    epsi2=1.e-1;
    epsi2b=1.e-6;
    /* Make sure that the sphere is of unit radius and is 
     * centered on the origin. The discretization may be fairly 
     * coarse, and the centroid may not be exactly at the sphere, 
     * so a low accuracy is called for when doing the test.
     * All coordinate vectors should be of length (roughly) unity      */
    /*assert(ABS( coordinate[0]*coordinate[0]
	       +coordinate[1]*coordinate[1]
	       +coordinate[2]*coordinate[2]-1.0) <= epsi2);*/
    /* Translation direction - change to see effects of 
     * discretization direction                                        */
    U[0] = 1.0;
    U[1] = 0.0;
    U[2] = 0.0;
    if (bcType==1) /* Diriclet condition */
      value = 0.5*(coordinate[0]*U[0]
		  +coordinate[1]*U[1]
		  +coordinate[2]*U[2]); /* =-0.5*(normal * U) */
    else { /* Neumann condition 
	    * Test inwards pointing normal - high accuracy here */
      assert(ABS(coordinate[0]+normal[0]) < epsi2b);
      assert(ABS(coordinate[1]+normal[1]) < epsi2b);
      assert(ABS(coordinate[2]+normal[2]) < epsi2b);
      value = -(normal[0]*U[0]+normal[1]*U[1]+normal[2]*U[2]);
    }
  }
  else if (strcmp(conditionType,"sphere charge") == 0){
    /* Nonuniform potential on a sphere of unit radius centered on 
     * the origin */
    /* POT    = +/- (1/2)*cos(theta)
     * CHARGE = +/- (3/8PI)*cos(theta) */
    double costheta, rad;
    int idir;
    epsi2=1.e-1;
    epsi2b=1.e-6;
    idir = 2; /* "polar" direction for potential distribution.
	       * Use idir \in \{0,1,2\}                                */
    /* Make sure that the sphere is of unit radius and is 
     * centered on the origin. 
     * All coordinate vectors should be of length roughly unity.       */
    /*assert(ABS( coordinate[0]*coordinate[0]
	       +coordinate[1]*coordinate[1]
	       +coordinate[2]*coordinate[2]-1.0) <= epsi2);*/
    /* Actual computations: 
    rad = VNORM(coordinate);
    costheta = coordinate[idir]/rad;*/
    costheta = coordinate[idir];
      rad = sqrt( coordinate[0]*coordinate[0]
		  +coordinate[1]*coordinate[1]
		  +coordinate[2]*coordinate[2]);
      costheta = costheta/rad;
    if (bcType==1) {
      /* Diriclet condition */
      value = -0.5*costheta; 
    }
    else /* Neumann condition */
      value = -3.0*costheta/(8.0*PI);
  }
  else {
    /* This is unknown = bad. Complain and exit!                       */
    fprintf(stderr,"%s: Error. Unknown boundary condition type: \"%s\"\n",
	    fctName,conditionType);
    _EXIT_;
  }
  /* Make sure that value has been set: */
  assert(value != DBL_MAX);
  return value;
} /* End of routine boundaryConditionValue */

