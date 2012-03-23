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
/* FILE: Direct.h
*
*  Structures for the computations of the direct matrix.
*
*  Note that this does *not* include the structure of the 
*  matrix itself. This is presently defined in "Sparse.h".
*  Consequently, Sparse.h will be included by default in 
*  this file.
*
*/
/* Only include this file if it hasn't been included before:           */
#if ! defined DirectDotHIsIncluded
#define DirectDotHIsIncluded

/* Make sure that the structure                                        */
#include "Sparse.h"

/* 
====================================================================== 
  STRUCTURES
====================================================================== */ 


/* This is information to the integration routine passed along 
by the routine that sets up the direct matrix. The idea of this 
structure is, that a relatively large/complex data set can be 
"passed through" an "ignorant" grid routine (which doesn't need 
to know what these data are anyway).                              */
/* struct integrationInformation{
   } */








#endif
/* ====== DONT PUT ANYTHING BELOW THIS LINE ====== */
