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
/* FILE: Nfft.h
 * Just the prototypes for the functions in nfft.c
 */

void *createFFT(int s1, int s2, int s3, char *Type);

void freeFFT(void *dataIn);

void fillKernel(void *griddataIn, double dx, double dy, double dz);

void convolve(void *data, void *kernel, int datatype);

int dataIndex(void *gridIn, int s1, int s2, int s3);

void spMatMultiply(void *passedMat, 
		   double *src, double *dest, int matNum);

void projectOntoGridData(void *gridDataIn,double *src,void *projectMat);

void interpolateFromGridData(void *gridDataIn,double *src,void *interpMat);
