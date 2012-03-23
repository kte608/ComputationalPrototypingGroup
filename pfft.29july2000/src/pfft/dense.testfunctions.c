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
/* ==== TEST FUNCTIONS (not really needed) ==== */ 

/* Test matrix invertion routine */
void testPinv(void)
{
  double **A,val;
  int nRows,nCols, i,j, n;
  nRows=3;
  nCols=4;

  n = MAX(nRows,nCols);
printf("n=%d\n",n);
  /* Allocate matrix */
  A = (double **) calloc(n,sizeof(double *));
  for (i=0;i<n;i++){
    A[i] = (double *) calloc(n,sizeof(double));
  }

  /* Set up test matrix values */
  val = 1.0;
  for (j=0;j<nCols;j++){
    if (j<nRows) A[j][j] = 2.0;
    for (i=0;i<nRows;i++){
      A[i][j]+=val;
      val += 1.0;
    }
  }

  /* Display matrix */
  for (i=0;i<nRows;i++){
    for (j=0;j<nCols;j++){
      printf(" %16.8f",A[i][j]);
    }
    printf("\n");
  }
  printf("\n");

  /* Invert matrix */
  pinv(A,nRows,nCols);

  /* Display matrix */
  for (i=0;i<nCols;i++){
    for (j=0;j<nRows;j++){
     printf(" %16.8f",A[i][j]);
    }
    printf("\n");
  }

  printf("testPinv halting program\n");
  _EXIT_;
}
