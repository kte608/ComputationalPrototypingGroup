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
/* FILE: stdmath.c
*
*  This file contains basic routines which may be used 
*  anywhere in the program .
*
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include "Global.h"

/* The following routines still need testing and should 
 * be used with caution: 
 *   setupCubatureRuleSquare
 *   setupGaussLegendreRule
 */



/* Internally used in this file: */
static
int intQuickSortPartition(int* Array, int nStart, int nEnd);

static
int intQuickSortPartition2(int* Array, int* Array2, int nStart, int nEnd);

static
int doubleQuickSortPartition2(double* Array, int* Array2, 
			      int nStart, int nEnd);

/*
* ====================================================================
* Calculate a positive (integer) power of a small integer.
* In test mode (compiled with debug flag) exits before overflow occurs.
* ==================================================================== */
int intpow (int base, int n)
{
  int i, p, maxval;
  maxval=(int) floor((double) INT_MAX / (double) base);
  p=1;
  for (i=1;i<=n; i++){
    assert(p<maxval);
    p = p*base;
  }
  return p;
}

/*
* ====================================================================
* Round of double to nearest integer
* This routine usd to be called "round", but apparently that messed 
* up a Linux-on-Alpha box. Thus the new name.
* ==================================================================== */
double dblRound (double x)
{
  return floor(0.5 + x);
}

/*
* ====================================================================
* Sort an array of n integers (Array[0] ... Array[n-1]) by calling
* a particular sorting algorithm. 
* Pt. this is just a "driver" routine for intQuickSort.
* Later, if other sorting routines are implemented, this routine 
* can be changed to call the fasted currently implemented routine, 
* possibly dependent on the number n.
* ==================================================================== */
void intSort(int* Array, int n)
{
  /* Use Quicksort (the only implemented at the moment); */
  intQuickSort(Array, 0, n-1);
  return;
} /* End of intSort                                                    */

/*
* ==================================================================== 
* Sort an array of integers using the Quicksort algorithm. 
* "Quicksort is a sorting algorithm whose worst-case running 
*  time is O(n^2) on an input array of n numbers. In spite of 
*  the slow worst-case running time, quicksort is often the 
*  best practical choise for sorting, because it is remarkably 
*  efficient on the average: its expected running time is O(n log n), 
*  and the constant factors hidden in the O(n log n) notation 
*  are quite small. It also has the advantage of sorting in place..."
*       [quote from Cormen, Leiserson and Rivest (1990)]
* ==================================================================== */
void intQuickSort(int* Array, int nStart, int nEnd)
     /* Sort elements nStart through nEnd of Array, i.e.
      * on return Array[nStart] ... Array[nEnd] will 
      * contain sorted values 
      */
{
  int divider;
  assert(nStart<=nEnd);

  if (nStart<nEnd){
    divider = intQuickSortPartition(Array,nStart,nEnd);
    intQuickSort(Array,nStart,divider);
    intQuickSort(Array,divider+1,nEnd);
  }
  return;
  /* Alternative sorting algorithms:
   * HeapSort, InsertionSort, MergeSort. */
} /* End of intQuickSort                                               */
/*
* ==================================================================== 
* Routine used by intQuickSort.
* ==================================================================== */
static
int intQuickSortPartition(int* Array,int nStart,int nEnd)
{
  int i,j, x, keep;

  x = Array[nStart];
  i = nStart - 1;
  j = nEnd   + 1;

  while(1){
    for(; Array[--j] > x ;);
    for(; Array[++i] < x ;);
    if (i<j){ /* Exchange i and j entries */
      keep     = Array[i];
      Array[i] = Array[j];
      Array[j] = keep;
    }
    else return j;
  }
} /* End of intQuickSortPartition                                      */


/*
* ====================================================================
* Sort an array of n integers (Array[0] ... Array[n-1]) by calling
* a particular sorting algorithm.
* Also sort a companying int array in the same order. This will allow
* fast sort of e.g. pointers and structures, since the companying array
* may give the position in the original structure. 
* Note: Array2 is sorted as per Array, but the entries of Array2 do 
* not influence the sorting of Array or Array2. Only the entries of 
* Array matters for the sorting.
*
* Pt. this is just a "driver" routine for intQuickSort.
* Later, if other sorting routines are implemented, this routine 
* can be changed to call the fasted currently implemented routine, 
* possibly dependent on the number n.
* ==================================================================== */
void intSort2(int* Array, int* Array2, int n)
{
  /* Use Quicksort (the only implemented at the moment); */
  intQuickSort2(Array, Array2, 0, n-1);
  return;
} /* End of intSort2                                                   */

/*
* ==================================================================== 
* Sort an array of integers using the Quicksort algorithm. 
* Also sort a companying int array in the same order. This will allow
* fast sort of e.g. pointers and structures, since the companying array
* may give the position in the original structure. 
* Note: Array2 is sorted as per Array, but the entries of Array2 do 
* not influence the sorting of Array or Array2. Only the entries of 
* Array matters for the sorting.
*
* "Quicksort is a sorting algorithm whose worst-case running 
*  time is O(n^2) on an input array of n numbers. In spite of 
*  the slow worst-case running time, quicksort is often the 
*  best practical choise for sorting, because it is remarkably 
*  efficient on the average: its expected running time is O(n log n), 
*  and the constant factors hidden in the O(n log n) notation 
*  are quite small. It also has the advantage of sorting in place..."
*       [quote from Cormen, Leiserson and Rivest (1990)]
* ==================================================================== */
void intQuickSort2(int* Array, int* Array2, int nStart, int nEnd)
     /* Sort elements nStart through nEnd of Array, i.e.
      * on return Array[nStart] ... Array[nEnd] will 
      * contain sorted values 
      */
{
  int divider;
  assert(nStart<=nEnd);

  if (nStart<nEnd){
    divider = intQuickSortPartition2(Array,Array2,nStart,nEnd);
    intQuickSort2(Array,Array2,nStart,divider);
    intQuickSort2(Array,Array2,divider+1,nEnd);
  }
  return;
  /* Alternative sorting algorithms:
   * HeapSort, InsertionSort, MergeSort. */
} /* End of intQuickSort2                                              */

/*
* ==================================================================== 
* Routine used by intQuickSort2.
* ==================================================================== */
static
int intQuickSortPartition2(int* Array,int* Array2,int nStart,int nEnd)
{
  int i,j, x, keep;

  x = Array[nStart];
  i = nStart - 1;
  j = nEnd   + 1;

  while(1){
    for(; Array[--j] > x ;);
    for(; Array[++i] < x ;);
    if (i<j){ /* Exchange i and j entries */
      keep     = Array[i];
      Array[i] = Array[j];
      Array[j] = keep;
      /* Do the same for the companying array                          */
      keep     = Array2[i];
      Array2[i] = Array2[j];
      Array2[j] = keep;
    }
    else return j;
  }
} /* End of intQuickSortPartition2                                     */

/*
* ====================================================================
* Sort an array of n double's (Array[0] ... Array[n-1]) by calling
* a particular sorting algorithm.
* Also sort a companying int array in the same order. 
* Note: Array2 is sorted as per Array, but the entries of Array2 do 
* not influence the sorting of Array or Array2. Only the entries of 
* Array matters for the sorting.
*
* Pt. this is just a "driver" routine for doubleQuickSort2.
* Later, if other sorting routines are implemented, this routine 
* can be changed to call the fastest currently implemented routine, 
* possibly depending on the number n.
* ==================================================================== */
void doubleSort2(double* Array, int* Array2, int n)
{
  /* Use Quicksort (the only implemented at the moment); */
  doubleQuickSort2(Array, Array2, 0, n-1);
  return;
} /* End of doubleSort2                                                   */

/*
* ==================================================================== 
* Sort an array of doubles using the Quicksort algorithm. 
* Also sort a companying int array in the same order. 
* Note: Array2 is sorted as per Array, but the entries of Array2 do 
* not influence the sorting of Array or Array2. Only the entries of 
* Array matters for the sorting.
* ==================================================================== */
void doubleQuickSort2(double* Array, int* Array2, int nStart, int nEnd)
     /* Sort elements nStart through nEnd of Array, i.e.
      * on return Array[nStart] ... Array[nEnd] will 
      * contain sorted values 
      */
{
  int divider;
  assert(nStart<=nEnd);

  if (nStart<nEnd){
    divider = doubleQuickSortPartition2(Array,Array2,nStart,nEnd);
    doubleQuickSort2(Array,Array2,nStart,divider);
    doubleQuickSort2(Array,Array2,divider+1,nEnd);
  }
  return;
  /* Alternative sorting algorithms:
   * HeapSort, InsertionSort, MergeSort. */
} /* End of doubleQuickSort2                                              */

/*
* ==================================================================== 
* Routine used by doubleQuickSort2.
* ==================================================================== */
static
int doubleQuickSortPartition2(double* Array,int* Array2,int nStart,int nEnd)
{
  int i,j, keepint;
  double x, keep;

  x = Array[nStart];
  i = nStart - 1;
  j = nEnd   + 1;

  while(1){
    for(; Array[--j] > x ;);
    for(; Array[++i] < x ;);
    if (i<j){ /* Exchange i and j entries */
      keep  = Array[i];
      Array[i] = Array[j];
      Array[j] = keep;
      /* Do the same for the companying array                          */
      keepint   = Array2[i];
      Array2[i] = Array2[j];
      Array2[j] = keepint;
    }
    else return j;
  }
} /* End of doubleQuickSortPartition2                                  */

/* 
* ==================================================================== 
* Utility to delete a leading char from a string.
*
*                 - Written by Tom Korsmeyer.
* ==================================================================== */
char *delcr(char *str)
{
  int i, j, k;
  for(k = 0; str[k] != '\0'; k++) if(str[k] == '\n') { str[k] = '\0'; break; }
  for(i = 0; str[i] == ' ' || str[i] == '\t'; i++); /* count leading spaces */
  if(i > 0) {
    for(j = 0; str[j+i] != '\0'; j++) str[j] = str[j+i];
    str[j] = '\0';
  }
  for(k--; str[k] == ' ' || str[k] == '\t'; k--) str[k] = '\0';
  return(str);
}

/* 
* ==================================================================== 
* Cubature formulae for integration over a square [-1:1]x[-1:1].
* Returns points (x,y) and weights for the chosen scheme.
* The schemes have been found from the on-line 
* "Encyclopaedia of Cubature Formulas" at
* URL: http://www.cs.kuleuven.ac.be/~nines/research/ecf/ecf.html
* Rules have been copied using "cut and paste" to reduce the risk of 
* typos. Any rule of order N integrates exactly all polynomials of 
* degree less than or equal to N.
*
* degree: Degree of the cubature rule to use.
* irule:  Number of rule of this degree to use (normally use irule=1).
*         Giving a negative value of irule returns the number of 
*         points in cubature rule #-irule of degree "degree"
*         without accessing x,y,w.
*
* x: X-coordinates of cubature points
* y: Y-coordinates of cubature points
* w: Weights of the cubature rule. w[i] is the weight to use 
*    at point (x[i],y[i]), i.e. the cubature rule is:
*
*     INT \approx SUM_{i=0}^{order-1} w[i]*f(x[i],y[i])
*
* Since the area of the domain of integration is four, the weights 
* of each scheme should add up to four.
*
* NOTE: This routine has not been tested yet. Test it and make it 
* globally available...
* ==================================================================== */
int setupCubatureRuleSquare(int degree, int irule, 
				   double *x, double *y, double *w)
{ 
  char fctName[]="setupCubatureRuleSquare";  /* Name of this routine   */
  int i, npoints;
  i=0;
  npoints=0;

  /* Shorthand notations used in this function. 
   * They will be undefined towards the end. */
  /* Fully symmetric, 1-point resultant: (coord1,coord2) must be (0,0) */
#define SQUARE_ORIGIN1(coord1,coord2,weight){\
    assert((coord1)==0.0);\
    assert((coord2)==0.0);\
    x[i]= (coord1); y[i]= (coord2); w[i]= weight; i++;\
}
  /* Fully symmetric, 4-point resultant: coord1==coord2 */
#define SQUARE_FULLY4XY(coord1,coord2,weight){\
    assert((coord1)==(coord2));\
    x[i]= (coord1); y[i]= (coord2); w[i]= (weight); i++;\
    x[i]= (coord1); y[i]=-(coord2); w[i]= (weight); i++;\
    x[i]=-(coord1); y[i]= (coord2); w[i]= (weight); i++;\
    x[i]=-(coord1); y[i]=-(coord2); w[i]= (weight); i++;\
}
  /* Fully symmetric, 4-point resultant: One coordinate (say, coord2) is zero */
#define SQUARE_FULLY4ZERO(coord1,coord2,weight){\
    assert((coord2)==0.0);\
    x[i]= (coord1); y[i]= (coord2); w[i]= (weight); i++;\
    x[i]=-(coord1); y[i]= (coord2); w[i]= (weight); i++;\
    x[i]= (coord2); y[i]= (coord1); w[i]= (weight); i++;\
    x[i]= (coord2); y[i]=-(coord1); w[i]= (weight); i++;\
}
  /* Partial symmetric, 1-point resultant: coord2 is zero */
#define SQUARE_PARTIAL1(coord1,coord2,weight){\
    assert((coord2==0.0));\
    x[i]= (coord1); y[i]= (coord2); w[i]= (weight); i++;\
}
  /* Partial symmetric, 2-point resultant */
#define SQUARE_PARTIAL2(coord1,coord2,weight){\
    assert((coord2)!=0.0);\
    x[i]= (coord1); y[i]= (coord2); w[i]= (weight); i++;\
    x[i]= (coord1); y[i]=-(coord2); w[i]= (weight); i++;\
}
  /* Central symmetric, 2-point resultant */
#define SQUARE_CENTRAL2(coord1,coord2,weight){\
    x[i]= (coord1); y[i]= (coord2); w[i]= (weight); i++;\
    x[i]=-(coord1); y[i]=-(coord2); w[i]= (weight); i++;\
}
  /* Rectangular symmetry */
#define SQUARE_RECTANGULAR4(coord1,coord2,weight){\
    x[i]= (coord1); y[i]= (coord2); w[i]= (weight); i++;\
    x[i]= (coord1); y[i]=-(coord2); w[i]= (weight); i++;\
    x[i]=-(coord1); y[i]= (coord2); w[i]= (weight); i++;\
    x[i]=-(coord1); y[i]=-(coord2); w[i]= (weight); i++;\
}
  /* Rectangular symmetry, coord2=0 */
#define SQUARE_RECTANGULAR2YZERO(coord1,coord2,weight){\
    assert((coord2)==0.0);\
    x[i]= (coord1); y[i]= (coord2); w[i]= (weight); i++;\
    x[i]=-(coord1); y[i]= (coord2); w[i]= (weight); i++;\
}
  /* Rectangular symmetry, coord1=0 */
#define SQUARE_RECTANGULAR2XZERO(coord1,coord2,weight){\
    assert((coord1)==0.0);\
    x[i]= (coord1); y[i]= (coord2); w[i]= (weight); i++;\
    x[i]= (coord1); y[i]=-(coord2); w[i]= (weight); i++;\
}
  /* Rotational invariant generator */
#define SQUARE_ROTATIONAL4(coord1,coord2,weight){\
    x[i]= (coord1); y[i]= (coord2); w[i]= (weight); i++;\
    x[i]=-(coord2); y[i]= (coord1); w[i]= (weight); i++;\
    x[i]=-(coord1); y[i]=-(coord2); w[i]= (weight); i++;\
    x[i]= (coord2); y[i]=-(coord1); w[i]= (weight); i++;\
}
  switch (100*degree+abs(irule)) {
  case  101: npoints= 1; break;
  case  301: npoints= 4; break;
  case  302: npoints= 4; break;
  case  401: npoints= 6; break;
  case  501: npoints= 7; break;
  case  502: npoints= 8; break;
  case  601: npoints=10; break;
  case  602: npoints=10; break;
  case  701: npoints=12; break;
  case  702: npoints=12; break;
  case  801: npoints=16; break;
  case  802: npoints=16; break;
  case  901: npoints=17; break;
  case  902: npoints=17; break;
  case  903: npoints=18; break;
  case 1001: npoints=24; break;
  default:
    printf("%s: ERROR: Unknown rule #%d of degree %d.\n",
	   fctName,irule,degree);
    printf("This halts the program!\n");
    exit(1);
  } /* End switch */

  if (irule<0) return npoints;

  switch (100*degree+irule) {

  case 101:
    /* Region: Cube 
     * Dimension: 2 
     * Degree: 1 
     * Points: 1 
     * Structure: Fully symmetric 
     * Rule struct: 1 0 0 0 
     * Centroid formula 
     * REF: [Str71] */
    /*npoints=1;*/
    /* Generator: [ Origin ] */
    SQUARE_ORIGIN1(0.00000000000000000,0.00000000000000000, 4.00000000000000000);
    break;

  case 301:
    /* Region: Cube 
     * Dimension: 2 
     * Degree: 3 
     * Points: 4 
     * Structure: Fully symmetric 
     * Rule struct: 0 1 0 0  
     * REF: [Str71] */
    npoints=4;
    /* Generator: [ Fully symmetric ] */
    SQUARE_FULLY4ZERO(0.81649658092772603,0.00000000000000000, 1.00000000000000000);
    break;

  case 302:
    /* Region: Cube 
     * Dimension: 2 
     * Degree: 3 
     * Points: 4 
     * Structure: Fully symmetric 
     * Rule struct: 0 0 1 0 
     * Product Gauss formula  
     * REF: [Str71] */
    npoints=4;
    /* Generator: [ Fully symmetric ] */
    SQUARE_FULLY4XY(0.57735026918962576,0.57735026918962576, 1.00000000000000000);
    break;

  case 401:
    /* Region: Cube 
     * Dimension: 2 
     * Degree: 4 
     * Points: 6 
     * Structure: Partial symmetric 
     * Rule struct: 1 1 0 2 
     * REF: [WB86] */
    npoints=6;
    /* Generator: [ Origin ] */
    SQUARE_ORIGIN1(  0.00000000000000000,0.00000000000000000, 1.1428571428571428);
    /* Generator: [ Partial symmetry ] (= single node in this case) */
    SQUARE_PARTIAL1( 0.96609178307929590,0.00000000000000000, 0.43956043956043956);
    /* Generator: [ Partial symmetry ] */
    SQUARE_PARTIAL2( 0.45560372783619284,0.85191465330460049, 0.56607220700753210);
    /* Generator: [ Partial symmetry ] */
    SQUARE_PARTIAL2(-0.73162995157313452,0.63091278897675402, 0.64271900178367668);
    break;

  case 501:
    /* Region: Cube 
     * Dimension: 2 
     * Degree: 5 
     * Points: 7 
     * Structure: Symmetric 
     * Rule struct: 1 1 0 1 
     * Points on origin or unit circle 
     * Ref: [Str71] */
    npoints=7;
    /* Generator: [ Origin ] */
    SQUARE_ORIGIN1(     0.00000000000000000,0.00000000000000000, 1.1428571428571428);
    /* Generator: [ Central symmetry ] */
    SQUARE_CENTRAL2(    0.96609178307929590,0.00000000000000000, 0.31746031746031746);
    /* Generator: [ Rectangular symmetry ] */
    SQUARE_RECTANGULAR4(0.57735026918962576,0.77459666924148337, 0.55555555555555555);
    break;

  case 502:
    /* Region: Cube 
     * Dimension: 2 
     * Degree: 5 
     * Points: 8 
     * Structure: Fully symmetric 
     * Rule struct: 0 1 1 0 
     * REF: [Str71] */
    npoints=8;
    /* Generator: [ Fully symmetric ] */
    SQUARE_FULLY4ZERO(0.68313005106397322,0.00000000000000000, 0.81632653061224489);
    /* Generator: [ Fully symmetric ] */
    SQUARE_FULLY4XY(  0.88191710368819686,0.88191710368819686, 0.18367346938775510);
    break;

  case 601:
    /* Region: Cube 
     * Dimension: 2 
     * Degree: 6 
     * Points: 10 
     * Structure: Partial symmetric 
     * Rule struct: 0 2 0 4 
     * REF: [WB86] */
    npoints=10;
    /* Generator: [ Partial symmetry ] (= single node in this case) */
    SQUARE_PARTIAL1( 0.83640563369762560,0.00000000000000000, 0.45534324571417343);
    /* Generator: [ Partial symmetry ] (= single node in this case) */
    SQUARE_PARTIAL1(-0.35746016539130718,0.00000000000000000, 0.82739597320296549);
    /* Generator: [ Partial symmetry ] */
    SQUARE_PARTIAL2( 0.87210153119313059,0.88876401465476453, 0.14400088459964538);
    /* Generator: [ Partial symmetry ] */
    SQUARE_PARTIAL2( 0.30598516215542666,0.60485763946468502, 0.66825910426266514);
    /* Generator: [ Partial symmetry ] */
    SQUARE_PARTIAL2(-0.41027089946665783,0.95544750664106374, 0.22547400489067935);
    /* Generator: [ Partial symmetry ] */
    SQUARE_PARTIAL2(-0.87286931115687933,0.56545999343875404, 0.32089639678844064);
    break;

  case 602:
    /* Region: Cube 
     * Dimension: 2 
     * Degree: 6 
     * Points: 10 
     * Structure: Partial symmetric 
     * Rule struct: 0 2 0 4 
     * REF: [WB86] */
    npoints=10;
    /* Generator: [ Partial symmetry ] (= single node in this case) */
    SQUARE_PARTIAL1( 0.86983337525005900,0.00000000000000000, 0.39275059096434794);
    /* Generator: [ Partial symmetry ] (= single node in this case) */
    SQUARE_PARTIAL1(-0.47940635161211124,0.00000000000000000, 0.75476288124261053);
    /* Generator: [ Partial symmetry ] */
    SQUARE_PARTIAL2( 0.80283751620765670,0.86374282634615388, 0.20616605058827902);
    /* Generator: [ Partial symmetry ] */
    SQUARE_PARTIAL2( 0.26214366550805818,0.51869052139258234, 0.68999213848986375);
    /* Generator: [ Partial symmetry ] */
    SQUARE_PARTIAL2(-0.36309658314806653,0.93397254497284950, 0.26051748873231697);
    /* Generator: [ Partial symmetry ] */
    SQUARE_PARTIAL2(-0.89660863276245265,0.60897753601635630, 0.26956758608606100);
    break;

  case 701:
    /* Region: Cube 
     * Dimension: 2 
     * Degree: 7 
     * Points: 12 
     * Structure: Fully symmetric 
     * Rule struct: 0 1 2 0 
     * REF: [Str71] */
    npoints=12;
    /* Generator: [ Fully symmetric ] */
    SQUARE_FULLY4ZERO(0.92582009977255146,0.00000000000000000, 0.24197530864197530);
    /* Generator: [ Fully symmetric ] */
    SQUARE_FULLY4XY(  0.38055443320831565,0.38055443320831565, 0.52059291666739445);
    /* Generator: [ Fully symmetric ] */
    SQUARE_FULLY4XY(  0.80597978291859874,0.80597978291859874, 0.23743177469063023);
    break;

  case 702:
    /* Region: Cube 
     * Dimension: 2 
     * Degree: 7 
     * Points: 12 
     * Structure: Symmetric 
     * Rule struct: 0 1+1 0 2 
     * REF: [HP77] */
    npoints=12;
    /* Generator: [ Rectangular symmetry ] */
    SQUARE_RECTANGULAR2YZERO(0.52942280204265532,0.00000000000000000, 0.63585388344327977);
    /* Generator: [ Rectangular symmetry ] */
    SQUARE_RECTANGULAR2XZERO(0.00000000000000000,0.62704137378039531, 0.59001271542103076);
    /* Generator: [ Rectangular symmetry ] */
    SQUARE_RECTANGULAR4(     0.91711782231277058,0.54793120682809232, 0.21305721162094912);
    /* Generator: [ Rectangular symmetry ] */
    SQUARE_RECTANGULAR4(     0.61126876646532841,0.93884325665885830, 0.17400948894689560);
    break;

  case 801:
    /* Region: Cube 
     * Dimension: 2 
     * Degree: 8 
     * Points: 16 
     * Structure: Partial symmetric 
     * Rule struct: 1 3 0 6 
     * REF: [wb86] */
    npoints=16;
    /* Generator: [ Origin ] */
    SQUARE_ORIGIN1(  0.00000000000000000,0.00000000000000000, 0.055364705621439772);
    /* Generator: [ Partial symmetry ] (= single node in this case) */
    SQUARE_PARTIAL1( 0.75762917766050544,0.00000000000000000, 0.40438936872607541);
    /* Generator: [ Partial symmetry ] (= single node in this case) */
    SQUARE_PARTIAL1(-0.23687184225570189,0.00000000000000000, 0.53354660495263445);
    /* Generator: [ Partial symmetry ] (= single node in this case) */
    SQUARE_PARTIAL1(-0.98971792904452668,0.00000000000000000, 0.11705418878673920);
    /* Generator: [ Partial symmetry ] */
    SQUARE_PARTIAL2( 0.95052095564566688,0.63909130490036965, 0.12561441761374680);
    /* Generator: [ Partial symmetry ] */
    SQUARE_PARTIAL2( 0.66388273688563317,0.93706907692499044, 0.13654458473358846);
    /* Generator: [ Partial symmetry ] */
    SQUARE_PARTIAL2( 0.30421068172410450,0.53708353054149367, 0.48340847921125695);
    /* Generator: [ Partial symmetry ] */
    SQUARE_PARTIAL2(-0.23649671853611948,0.88718850644962463, 0.25252850642954369);
    /* Generator: [ Partial symmetry ] */
    SQUARE_PARTIAL2(-0.69895347608656355,0.49469882067019633, 0.36126232388217256);
    /* Generator: [ Partial symmetry ] */
    SQUARE_PARTIAL2(-0.90039077421157973,0.89749581827976753, 0.085464254086247091);
    break;

  case 802:
    /* Region: Cube 
     * Dimension: 2 
     * Degree: 8 
     * Points: 16 
     * Structure: Partial symmetric 
     * Rule struct: 0 2 0 7 
     * REF: [wb86] */
    npoints=16;
    /* Generator: [ Partial symmetry ] (= single node in this case) */
    SQUARE_PARTIAL1( 0.65956013196034176, 0.00000000000000000, 0.45027677630559029);
    /* Generator: [ Partial symmetry ] (= single node in this case) */
    SQUARE_PARTIAL1(-0.94914292304312538, 0.00000000000000000, 0.16657042677781274);
    /* Generator: [ Partial symmetry ] */
    SQUARE_PARTIAL2( 0.76505181955768362, 0.95250946607156228, 0.098869459933431422);
    /* Generator: [ Partial symmetry ] */
    SQUARE_PARTIAL2( 0.93697598108841598, 0.53232745407420624, 0.15369674714081197);
    /* Generator: [ Partial symmetry ] */
    SQUARE_PARTIAL2( 0.33365671773574759, 0.68473629795173504, 0.39668697607290278);
    /* Generator: [ Partial symmetry ] */
    SQUARE_PARTIAL2(-0.079583272377396852,0.23314324080140552, 0.35201436794569501);
    /* Generator: [ Partial symmetry ] */
    SQUARE_PARTIAL2(-0.27224008061253425, 0.92768331930611748, 0.18958905457779799);
    /* Generator: [ Partial symmetry ] */
    SQUARE_PARTIAL2(-0.61373535339802760, 0.45312068740374942, 0.37510100114758727);
    /* Generator: [ Partial symmetry ] */
    SQUARE_PARTIAL2(-0.88847765053597136, 0.83750364042281223, 0.12561879164007201);
    break;

  case 901:
    /* Region: Cube 
     * Dimension: 2 
     * Degree: 9 
     * Points: 17 
     * Structure: Rotational invariant 
     * Rule struct: 1 0 0 4 
     * REF: [Mol76] */
    npoints=17;
    /* Rotational invariant: 
     * If (x,y) is a point with weight w, 
     * then (-x,-y), (-y,x), (y,-x) are points with weight w */
    /* Generator: [ Origin ] */
    SQUARE_ORIGIN1(    0.00000000000000000, 0.00000000000000000, 0.52674897119341563);
    /* Generator: [ Rotational invariant ] */
    SQUARE_ROTATIONAL4(0.96884996636197772, 0.63068011973166885, 0.088879378170198706);
    /* Generator: [ Rotational invariant ] */
    SQUARE_ROTATIONAL4(0.75027709997890053, 0.92796164595956966, 0.11209960212959648);
    /* Generator: [ Rotational invariant ] */
    SQUARE_ROTATIONAL4(0.52373582021442933, 0.45333982113564719, 0.39828243926207009);
    /* Generator: [ Rotational invariant ] */
    SQUARE_ROTATIONAL4(0.076208328192617173,0.85261572933366230, 0.26905133763978080);
    break;

  case 902:
    /* Region: Cube 
     * Dimension: 2 
     * Degree: 9 
     * Points: 18 
     * Structure: Symmetric 
     * Rule struct: 0 1+2 0 3 
     * REF: [PH75b] */
    npoints=17;
    /* Generator: [ Rectangular symmetry ] */
    SQUARE_RECTANGULAR2YZERO(0.57882826011929170,0.00000000000000000, 0.40927359555433144);
    /* Generator: [ Rectangular symmetry ] */
    SQUARE_RECTANGULAR2XZERO(0.00000000000000000,0.97700090158004246, 0.10648011781560231);
    /* Generator: [ Rectangular symmetry ] */
    SQUARE_RECTANGULAR2XZERO(0.00000000000000000,0.39364057271848893, 0.45321488105170985);
    /* Generator: [ Rectangular symmetry ] */
    SQUARE_RECTANGULAR4(     0.87980721399752853,0.92797961509268528, 0.068416522462309305);
    /* Generator: [ Rectangular symmetry ] */
    SQUARE_RECTANGULAR4(     0.50445910315479838,0.75347199103161505, 0.27903384209687301);
    /* Generator: [ Rectangular symmetry ] */
    SQUARE_RECTANGULAR4(     0.91531235408227324,0.42299357094876513, 0.16806533822999587);
    break;

  case 903:
    /* Region: Cube 
     * Dimension: 2 
     * Degree: 9 
     * Points: 18 
     * Structure: Symmetric 
     * Rule struct: 0 1+2 0 3 
     * REF: [PH75b] */
    npoints=18;
    /* Generator: [ Rectangular symmetry ] */
    SQUARE_RECTANGULAR2YZERO(0.49471787965159623,0.00000000000000000, 0.45212398131214854);
    /* Generator: [ Rectangular symmetry ] */
    SQUARE_RECTANGULAR2XZERO(0.00000000000000000,0.98085697194664054, 0.10243215270991495);
    /* Generator: [ Rectangular symmetry ] */
    SQUARE_RECTANGULAR2XZERO(0.00000000000000000,0.48311469619727965, 0.45601422352687001);
    /* Generator: [ Rectangular symmetry ] */
    SQUARE_RECTANGULAR4(     0.93742666622066710,0.94145119299928430, 0.042853317248897088);
    /* Generator: [ Rectangular symmetry ] */
    SQUARE_RECTANGULAR4(     0.57077001686857404,0.79214654516847247, 0.25788406360659644);
    /* Generator: [ Rectangular symmetry ] */
    SQUARE_RECTANGULAR4(     0.89774224179848572,0.40001733897633692, 0.19397744037003970);
    break;

  case 1001:
    /* Region: Cube 
     * Dimension: 2 
     * Degree: 11 
     * Points: 24 
     * Structure: Rotational invariant 
     * Rule struct: 0 0 0 6 
     * REF: [CH88a][VC94b] */
    npoints=24;
    /* Rotational invariant: 
     * If (x,y) is a point with weight w, 
     * then (-x,-y), (-y,x), (y,-x) are points with weight w */
    /* Generator: [ Rotational invariant ] */
    SQUARE_ROTATIONAL4(0.98263922354085547, 0.69807610454956756, 0.048020763350723814);
    /* Generator: [ Rotational invariant ] */
    SQUARE_ROTATIONAL4(0.82577583590296393, 0.93948638281673690, 0.066071329164550595);
    /* Generator: [ Rotational invariant ] */
    SQUARE_ROTATIONAL4(0.18858613871864195, 0.95353952820153201, 0.097386777358668164);
    /* Generator: [ Rotational invariant ] */
    SQUARE_ROTATIONAL4(0.81252054830481310, 0.31562343291525419, 0.21173634999894860);
    /* Generator: [ Rotational invariant ] */
    SQUARE_ROTATIONAL4(0.52532025036454776, 0.71200191307533630, 0.22562606172886338);
    /* Generator: [ Rotational invariant ] */
    SQUARE_ROTATIONAL4(0.041658071912022368,0.42484724884866925, 0.35115871839824543);
    break;
    
  default:
    printf("%s: ERROR: The flow should never reach this statement!.\n",
	   fctName);
    printf("This halts the program!\n");
    exit(1);
  } /* End switch */
  /* Make a small test to see i at least the right number of nodes and 
   * weigths have been obtained */
  assert(i==npoints);
  /* References (these also taken from the "Encyclopedia"): 
   *
   *  [CH88a]   R. Cools and A. Haegemans, 
   *            Another step forward in searching for cubature formulae with 
   *            a minimal number of knots for the square, 
   *            Computing 40 (1988), 139--146. 
   *
   *  [HP77]:   A. Haegemans and R. Piessens, 
   *            Construction of cubature formulas of degree seven and nine symmetric
   *            planar regions, using orthogonal polynomials, 
   *            SIAM J. Numer. Anal. 14 (1977), 492--508. 
   *
   *  [Mol76]   H.M. M{\"o}ller,
   *            Kubaturformeln mit minimaler Knotenzahl, 
   *            Numer. Math. 25 (1976), 185--200. 
   *                                                       extra " to help emacs' font selection
   *  [PH75b]   R. Piessens and A. Haegemans, 
   *            Cubature formulas of degree nine for symmetric planar regions, 
   *            Math. Comp. 29 (1975), 810--815. 
   * 
   *  [Str71]:  A.H. Stroud, 
   *            Approximate calculation of multiple integrals, 
   *            Prentice-Hall, Englewood Cliffs, N.J., 1971. 
   * 
   *  [VC94]:   P. Verlinden and R. Cools, 
   *            The algebraic construction of a minimal cubature formula 
   *            of degree 11 for the square, 
   *            Cubature Formulas and their Applications (Russian) (Krasnoyarsk) 
   *            (M.V. Noskov, ed.), 1994, pp. 13--23.
   *
   *  [WB86]:   J. W. Wissman and T. Becker, 
   *            Partially symmetric cubature formulas for even degrees of exactness,
   *            SIAM J. Numer. Anal. 23 (1986), 676--685. 
   */

#undef SQUARE_ORIGIN1
#undef SQUARE_FULLY4XY
#undef SQUARE_FULLY4ZERO
#undef SQUARE_PARTIAL1
#undef SQUARE_PARTIAL2
#undef SQUARE_CENTRAL2
#undef SQUARE_RECTANGULAR4
#undef SQUARE_RECTANGULAR2YZERO
#undef SQUARE_RECTANGULAR2XZERO
#undef SQUARE_ROTATIONAL4

  return npoints;
} 

/* 
* ==================================================================== 
* Cubature formulae for integration over a triangle with corners
* (0,0) (0,1) and (0,1). This is a two-dimensional simplex.
* Returns points (x,y) and weights for the chosen scheme.
* The schemes have been found from the on-line 
* "Encyclopaedia of Cubature Formulas" at
* URL: http://www.cs.kuleuven.ac.be/~nines/research/ecf/ecf.html
* Rules have been copied using "cut and paste" to reduce the risk of 
* typos. Any rule of order N integrates exactly all polynomials of 
* degree less than or equal to N.
*
* degree: Degree of the cubature rule to use.
* irule:  Number of rule of this degree to use (normally use irule=1).
*         Giving a negative value of irule returns the number of 
*         points in cubature rule #-irule of degree "degree"
*         without accessing a,b,c,w.
*
* (a,b,c) Barycentric coordinates of cubture points (see below).
* w: Weights of the cubature rule. w[i] is the weight to use 
*    at point (x[i],y[i]), i.e. the cubature rule is:
*
*     INT \approx SUM_{i=0}^{order-1} w[i]*f(x[i],y[i]) 
*
* Since the area of the domain of integration is one half, the weights 
* of each scheme should add up to one half.
*
* Barycentric coordinates (a,b,c) are used to represent points in the 
* plane of the triangle. If the corners of the triangle are denoted 
* by (X1,X2,X3), then the Cartesian coordinates of a point are 
* (a*X1,b*X2,c*X3).
* Note that a+b+c=1 (the Barycentric coordinates contain redundant 
* information). This is exploited in the present routine, so that 
* only two of the three coordinates are given at any time. The last 
* coordinate is found by subtracting from unity the values of the other 
* two coordinates. 
*
* NOTE: The higher-order cubature formulae have been tested by finding 
* moments of flat triangular panels. All tested formulae gave correct 
* results for all the polynomials used in the integration (polynomial 
* orders up to six and cubature degree from six to eleven was tested). 
* (At the present the convergence orders have not been tested.)
* ==================================================================== */
int setupCubatureRuleTriangle(int degree, int irule, 
			       double *a, double *b, double *c, 
			       double *w)
{ 
  char fctName[]="setupCubatureRuleTriangle";  /* Name of this routine */
  int i, npoints, one=1.00000000000000000;
  double coord3;
  i=0;
  npoints=0;

  /* Shorthand notations used in this function. 
   * They will be undefined towards the end. */
  /* Fully symmetric, 1-point resultant: must be (1/3,1/3,1/3) */
#define SIMPLEX2_FULLY1(coord1,coord2,weight){ \
    assert(coord1==0.33333333333333333); \
    assert(coord1==coord2); \
    a[i]= coord1; b[i]=coord2; c[i]=one-coord1-coord2; w[i]=weight; i++;\
}
  /* Fully symmetric, 3-point resultant: can be written so that coord1=coord2. */
#define SIMPLEX2_FULLY3(coord1,coord2,weight){ \
    assert(coord1==coord2); \
    coord3=one-coord1-coord2;\
    a[i]= coord1; b[i]=coord2; c[i]=coord3; w[i]=weight; i++;\
    b[i]= coord1; c[i]=coord2; a[i]=coord3; w[i]=weight; i++;\
    c[i]= coord1; a[i]=coord2; b[i]=coord3; w[i]=weight; i++;\
}
  /* Fully symmetric, 6-point resultant */
#define SIMPLEX2_FULLY6(coord1,coord2,weight){ \
    assert(coord1!=coord2); \
    coord3=one-coord1-coord2;\
    a[i]= coord1; b[i]=coord2; c[i]=coord3; w[i]=weight; i++;\
    a[i]= coord1; c[i]=coord2; b[i]=coord3; w[i]=weight; i++;\
    b[i]= coord1; a[i]=coord2; c[i]=coord3; w[i]=weight; i++;\
    b[i]= coord1; c[i]=coord2; a[i]=coord3; w[i]=weight; i++;\
    c[i]= coord1; a[i]=coord2; b[i]=coord3; w[i]=weight; i++;\
    c[i]= coord1; b[i]=coord2; a[i]=coord3; w[i]=weight; i++;\
}
  /* Rotational symmetric (Ro3 invariant) 3-point resultant */
#define SIMPLEX2_ROTATIONAL3(coord1,coord2,weight){ \
    coord3=one-coord1-coord2;\
    a[i]= coord1; b[i]=coord2; c[i]=coord3; w[i]=weight; i++;\
    b[i]= coord1; c[i]=coord2; a[i]=coord3; w[i]=weight; i++;\
    c[i]= coord1; a[i]=coord2; b[i]=coord3; w[i]=weight; i++;\
}

  switch (100*degree+abs(irule)) {
  case  101: npoints= 1; break;
  case  201: npoints= 3; break;
  case  301: npoints= 4; break; /* Contains negative weights */
  case  302: npoints= 6; break;
  case  401: npoints= 6; break;
  case  501: npoints= 7; break;
  case  601: npoints=12; break;
  case  602: npoints=12; break;
  case  701: npoints=12; break;
  case  801: npoints=16; break;
  case  901: npoints=19; break;
  case 1001: npoints=25; break;
  case 1002: npoints=25; break;
  case 1101: npoints=28; break;
  default:
    printf("%s: ERROR: Unknown rule #%d of degree %d.\n Switch value %d\n",
	   fctName,irule,degree,100*degree+abs(irule));
    printf("This halts the program!\n");
    exit(1);
  }
    
  if (irule<0) return npoints;

  switch (100*degree+irule) {

  case 101:
    /* Region: Simplex 
     * Dimension: 2 
     * Degree: 1 
     * Points: 1 
     * Structure: Fully symmetric 
     * Rule struct: 0 0 0 1 0 0 
     * REF: [Str71] */
    /*npoints=1;*/
    /* Generator: [ Fully symmetric ] (= single node in this case) */
    SIMPLEX2_FULLY1(0.33333333333333333,0.33333333333333333, 0.50000000000000000);
    break;

  case 201:
    /* Region: Simplex 
     * Dimension: 2 
     * Degree: 2 
     * Points: 3 
     * Structure: Fully symmetric 
     * Rule struct: 0 0 0 0 1 0 
     * REF: [Str71] */
    /*npoints=3;*/
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.16666666666666666,0.16666666666666666, 0.16666666666666666);
    break;

  case 301:
    /* Region: Simplex 
     * Dimension: 2 
     * Degree: 3 
     * Points: 4 
     * Structure: Fully symmetric 
     * Rule struct: 0 0 0 1 1 0 
     * REF: [Str71] 
     *      NOTE: This rule contains negative weights! */
    /*npoints=4;*/
    /* Generator: [ Fully symmetric ] (= single node in this case) */
    SIMPLEX2_FULLY1(0.33333333333333333,0.33333333333333333, -0.28125000000000000);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.20000000000000000,0.20000000000000000,  0.26041666666666667);
    break;

  case 302:
    /* Region: Simplex 
     * Dimension: 2 
     * Degree: 3 
     * Points: 6 
     * Structure: Fully symmetric 
     * Rule struct: 0 0 0 0 0 1 
     * REF: [Str71] */
    /*npoints=6;*/
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY6(0.10903900907287721,0.23193336855303057, 0.083333333333333333);
    break;

  case 401:
    /* Region: Simplex 
     * Dimension: 2 
     * Degree: 4 
     * Points: 6 
     * Structure: Fully symmetric 
     * Rule struct: 0 0 0 0 2 0 
     * REF: [cow73][dun85b][lg78][lj75][moa74][sf73] */
    /*npoints=6;*/
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.091576213509770743,0.091576213509770743, 0.054975871827660933);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.44594849091596488, 0.44594849091596488,  0.11169079483900573);
    break;

  case 501:
    /* Region: Simplex 
     * Dimension: 2 
     * Degree: 5 
     * Points: 7 
     * Structure: Fully symmetric 
     * Rule struct: 0 0 0 1 2 0 
     * REF: [str71] */
    /*npoints=7;*/
    /* Generator: [ Fully symmetric ]  (= single node in this case) */
    SIMPLEX2_FULLY1(0.33333333333333333,0.33333333333333333, 0.11250000000000000);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.10128650732345633,0.10128650732345633, 0.062969590272413576);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.47014206410511508,0.47014206410511508, 0.066197076394253090);
    break;

  case 601:
    /* Region: Simplex 
     * Dimension: 2 
     * Degree: 6 
     * Points: 12 
     * Structure: Fully symmetric 
     * Rule struct: 0 0 0 0 2 1 
     * REF: [cow73][dun85b][lg78][lj75][moa74][sf73] */
    /*npoints=12;*/
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.063089014491502228,0.063089014491502228, 0.025422453185103408);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.24928674517091042, 0.24928674517091042,  0.058393137863189683);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY6(0.053145049844816947,0.31035245103378440,  0.041425537809186787);
    break;

  case 602:
    /* Region: Simplex 
     * Dimension: 2 
     * Degree: 6 
     * Points: 12 
     * Structure: Fully symmetric 
     * Rule struct: 0 0 0 0 2 1 
     * REF: [bec87] */
    /*npoints=12;*/
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.21942998254978296,0.21942998254978296, 0.085666562076490515);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.48013796411221504,0.48013796411221504, 0.040365544796515489);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY6(0.83900925971479105,0.14161901592396815, 0.020317279896830331);
    break;

  case 701:
    /* Region: Simplex 
     * Dimension: 2 
     * Degree: 7 
     * Points: 12 
     * Structure: Ro3 invariant 
     * Rule struct: 0 0 0 0 0 4 
     * REF: [gat88] */
    /*npoints=12;*/
    /* Generator: [ Ro3 invariant ] (Cyclic rotations) */
    SIMPLEX2_ROTATIONAL3(0.062382265094402118,0.067517867073916085, 0.026517028157436251);
    /* Generator: [ Ro3 invariant ] */
    SIMPLEX2_ROTATIONAL3(0.055225456656926611,0.32150249385198182,  0.043881408714446055);
    /* Generator: [ Ro3 invariant ] */
    SIMPLEX2_ROTATIONAL3(0.034324302945097146,0.66094919618673565,  0.028775042784981585);
    /* Generator: [ Ro3 invariant ] */
    SIMPLEX2_ROTATIONAL3(0.51584233435359177,0.27771616697639178, 0.067493187009802774);
    break;

  case 801:
    /* Region: Simplex 
     * Dimension: 2 
     * Degree: 8 
     * Points: 16 
     * Structure: Fully symmetric 
     * Rule struct: 0 0 0 1 3 1 
     * REF: [lj75][dun85b][lg78] */
    /*npoints=16;*/
    /* Generator: [ Fully symmetric ] (= single node in this case) */
    SIMPLEX2_FULLY1(0.33333333333333333, 0.33333333333333333,  0.072157803838893584);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.17056930775176020, 0.17056930775176020,  0.051608685267359125);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.050547228317030975,0.050547228317030975, 0.016229248811599040);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.45929258829272315, 0.45929258829272315,  0.047545817133642312);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY6(0.72849239295540428, 0.26311282963463811,  0.013615157087217497);
    break;

  case 901:
    /* Region: Simplex 
     * Dimension: 2 
     * Degree: 9 
     * Points: 19 
     * Structure: Fully symmetric 
     * Rule struct: 0 0 0 1 4 1 
     * REF: [lj75][dun85b][lg78] */
    /*npoints=19;*/
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY1(0.33333333333333333, 0.33333333333333333,  0.048567898141399416);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.48968251919873762, 0.48968251919873762,  0.015667350113569535);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.43708959149293663, 0.43708959149293663,  0.038913770502387139);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.18820353561903273, 0.18820353561903273,  0.039823869463605126);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.044729513394452709,0.044729513394452709, 0.012788837829349015);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY6(0.74119859878449802, 0.036838412054736283, 0.021641769688644688);
    break;

  case 1001:
    /* Region: Simplex 
     * Dimension: 2 
     * Degree: 10 
     * Points: 25 
     * Structure: Fully symmetric 
     * Rule struct: 0 0 0 1 2 3 
     * REF [lg78] */
    /*npoints=25;*/
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY1(0.33333333333333333, 0.33333333333333333,  0.039947252370619853);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.42508621060209057, 0.42508621060209057,  0.035561901116188667);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.023308867510000190,0.023308867510000190, 4.1119093452320977e-3);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY6(0.62830740021349255, 0.22376697357697300,  0.022715296148085009);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY6(0.61131382618139764, 0.35874014186443146,  0.018679928117152638);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY6(0.82107206998562937, 0.14329537042686714,  0.015443328442281994);
    break;

  case 1002:
    /* Region: Simplex 
     * Dimension: 2 
     * Degree: 10 
     * Points: 25 
     * Structure: Fully symmetric 
     * Rule struct: 0 0 0 1 2 3 
     * REF: [lg78] */
    /*npoints=25;*/
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY1(0.33333333333333333, 0.33333333333333333,  0.040871664573142983);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.14216110105656438, 0.14216110105656438,  0.022978981802372364);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.032055373216943512,0.032055373216943512, 6.6764844065747831e-3);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY6(0.53005411892734402, 0.32181299528883542,  0.031952453198212022);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY6(0.60123332868345924, 0.36914678182781098,  0.017092324081479714);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY6(0.80793060092287906, 0.16370173373718249,  0.012648878853644192);
    break;

  case 1101:
    /* Region: Simplex 
     * Dimension: 2 
     * Degree: 11 
     * Points: 28 
     * Structure: Fully symmetric 
     * Rule struct: 0 0 1 1 5 1 
     * REF: [lj75] */
    /*npoints=28;*/
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY6(0.85887028128263670, 0.14112971871736329,  3.6811918916502771e-3);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY1(0.33333333333333333, 0.33333333333333333,  0.043988650581116119);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.025989140928287395,0.025989140928287395, 4.3721557768680115e-3);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.094287502647922495,0.094287502647922495, 0.019040785996967468);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.49463677501721381, 0.49463677501721381,  9.4277240280656460e-3); 
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.20734338261451133, 0.20734338261451133,  0.036079848772369763);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.43890780570049209, 0.43890780570049209,  0.034664569352767949);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY6(0.67793765488259040, 0.044841677589130443, 0.020528157714644283);
    break;

  default:
    printf("%s: ERROR: The flow should never reach this statement!.\n",
	   fctName);
    printf("This halts the program!\n");
    exit(1);
  } /* End switch */
  /* Make a small test to see that at least the right number of nodes 
   * and weigths have been obtained */
  assert(i==npoints);
  /* References (these also taken from the "Encyclopedia"): 
   *
   *  [Bec87]:  T. Becker, 
   *            Konstruktion von interpolatorischen Kubaturformeln mit 
   *            Anwendungen in der Finit-Element-Methode, 
   *            Ph.D. thesis, Technische Hochschule Darmstadt, 1987. 
   *
   *  [Cow73]:  G.R. Cowper, 
   *            Gaussian quadrature formulas for triangles, 
   *            Internat. J. Numer. Methods Engrg. 7 (1973), 405--408.
   *
   *  [Dun85b]: D.A. Dunavant, 
   *            High degree efficient symmetrical Gaussian quadrature rules for the triangle, 
   *            Internat. J. Numer. Methods Engrg. 21 (1985), 1129--1148. 
   *
   *  [Gat88]:  K. Gatermann, 
   *            The construction of symmetric cubature formulas for the square and the triangle,
   *            Computing 40 (1988), 229--240. 
   *
   *  [LG78]:   M.E. Laursen and M. Gellert, 
   *            Some criteria for numerically integrated matrices and quadrature formulas
   *            for triangles, Internat. J. Numer. Methods Engrg. 12 (1978), 67--76. 
   *
   *  [LJ75]:   J.N. Lyness and D. Jespersen, 
   *            Moderate degree symmetric quadrature rules for the triangle, 
   *            J. Inst. Math. Appl. 15 (1975), 19--32. 
   *
   *  [Moa74]:  T. Moan, 
   *            Experiences with orthogonal polynomials and ``best'' numerical 
   *            integration formulas on a triangle; with particular reference to 
   *            finite element approximations, 
   *            Z. Angew. Math. Mech. 54 (1974), 501--508. 
   *
   *  [SF73]:   G. Strang and G.J. Fix, 
   *            An analysis of the finite element method, 
   *            Prentice-Hall, London, 1973. 
   *
   *  [Str71]:  A.H. Stroud, 
   *            Approximate calculation of multiple integrals, 
   *            Prentice-Hall, Englewood Cliffs, N.J., 1971. 
   */

#undef SIMPLEX2_FULLY1
#undef SIMPLEX2_FULLY3
#undef SIMPLEX2_FULLY6
#undef SIMPLEX2_ROTATIONAL3

  return npoints;
}

/* 
* ==================================================================== 
* Gauss-Legendre quadrature points and weights for the 1D integration 
* over the region [-1,1]. Different orders of the quadrature formulae 
* are included.
* Some of these results were compiled from Abramowitz and Stegun (1964) 
* [abramowitz:1964:dover].
*
* The arrays "points" and "weights" must be (at least) of size order.
* order:   Input - the order of the Gauss-Legendre scheme.
*           At the time of writing acccepted values are 1 through 24.
* points:  Output - the coordinate of the quadrature points.
* weights: Output - the quadrature weights of each Gauss point.
*
* NOTE: This routine has not been tested yet.
* ==================================================================== */
void setupGaussLegendreRule(int order, double *points, double *weights)
{
  char fctName[]="setupGaussLegendreRule";  /* Name of this routine    */
  int maxorder=24;
  /* The following part (select construct including data) is 
   * computer generated using gauss-routines from NETLIB. 
   * These routines could well be implemented here, except that they 
   * are written in FORTRAN. If somebody want to reimplement these 
   * routines, or simply pass them through an f2c compiler, then by all 
   * means do so. Until then - the tables here will have to suffice.   */
  /* BEGIN COMPUTER GENERATED SOURCE */
  switch (order){
  case   1:
  /* Gauss-Legendre rule, order=  1   */
  /*  Gauss-points: */
    points[  0]  =    0.000000000000000;
  /*  Gauss-weights: */
    weights[  0] =    2.000000000000000;
    break;
  case   2:
  /* Gauss-Legendre rule, order=  2   */
  /*  Gauss-points: */
    points[  0]  =   -0.577350269189626;
    points[  1]  =    0.577350269189626;
  /*  Gauss-weights: */
    weights[  0] =    1.000000000000000;
    weights[  1] =    1.000000000000000;
    break;
  case   3:
  /* Gauss-Legendre rule, order=  3   */
  /*  Gauss-points: */
    points[  0]  =   -0.774596669241484;
    points[  1]  =    0.000000000000000;
    points[  2]  =    0.774596669241483;
  /*  Gauss-weights: */
    weights[  0] =    0.555555555555556;
    weights[  1] =    0.888888888888889;
    weights[  2] =    0.555555555555556;
    break;
  case   4:
  /* Gauss-Legendre rule, order=  4   */
  /*  Gauss-points: */
    points[  0]  =   -0.861136311594053;
    points[  1]  =   -0.339981043584856;
    points[  2]  =    0.339981043584856;
    points[  3]  =    0.861136311594053;
  /*  Gauss-weights: */
    weights[  0] =    0.347854845137454;
    weights[  1] =    0.652145154862546;
    weights[  2] =    0.652145154862547;
    weights[  3] =    0.347854845137454;
    break;
  case   5:
  /* Gauss-Legendre rule, order=  5   */
  /*  Gauss-points: */
    points[  0]  =   -0.906179845938664;
    points[  1]  =   -0.538469310105683;
    points[  2]  =    0.000000000000000;
    points[  3]  =    0.538469310105683;
    points[  4]  =    0.906179845938664;
  /*  Gauss-weights: */
    weights[  0] =    0.236926885056189;
    weights[  1] =    0.478628670499367;
    weights[  2] =    0.568888888888889;
    weights[  3] =    0.478628670499367;
    weights[  4] =    0.236926885056189;
    break;
  case   6:
  /* Gauss-Legendre rule, order=  6   */
  /*  Gauss-points: */
    points[  0]  =   -0.932469514203152;
    points[  1]  =   -0.661209386466264;
    points[  2]  =   -0.238619186083197;
    points[  3]  =    0.238619186083197;
    points[  4]  =    0.661209386466265;
    points[  5]  =    0.932469514203152;
  /*  Gauss-weights: */
    weights[  0] =    0.171324492379170;
    weights[  1] =    0.360761573048139;
    weights[  2] =    0.467913934572690;
    weights[  3] =    0.467913934572691;
    weights[  4] =    0.360761573048138;
    weights[  5] =    0.171324492379171;
    break;
  case   7:
  /* Gauss-Legendre rule, order=  7   */
  /*  Gauss-points: */
    points[  0]  =   -0.949107912342758;
    points[  1]  =   -0.741531185599394;
    points[  2]  =   -0.405845151377397;
    points[  3]  =    0.000000000000000;
    points[  4]  =    0.405845151377397;
    points[  5]  =    0.741531185599394;
    points[  6]  =    0.949107912342758;
  /*  Gauss-weights: */
    weights[  0] =    0.129484966168870;
    weights[  1] =    0.279705391489276;
    weights[  2] =    0.381830050505119;
    weights[  3] =    0.417959183673470;
    weights[  4] =    0.381830050505119;
    weights[  5] =    0.279705391489277;
    weights[  6] =    0.129484966168870;
    break;
  case   8:
  /* Gauss-Legendre rule, order=  8   */
  /*  Gauss-points: */
    points[  0]  =   -0.960289856497536;
    points[  1]  =   -0.796666477413627;
    points[  2]  =   -0.525532409916329;
    points[  3]  =   -0.183434642495650;
    points[  4]  =    0.183434642495650;
    points[  5]  =    0.525532409916329;
    points[  6]  =    0.796666477413627;
    points[  7]  =    0.960289856497536;
  /*  Gauss-weights: */
    weights[  0] =    0.101228536290376;
    weights[  1] =    0.222381034453374;
    weights[  2] =    0.313706645877887;
    weights[  3] =    0.362683783378362;
    weights[  4] =    0.362683783378362;
    weights[  5] =    0.313706645877887;
    weights[  6] =    0.222381034453374;
    weights[  7] =    0.101228536290376;
    break;
  case   9:
  /* Gauss-Legendre rule, order=  9   */
  /*  Gauss-points: */
    points[  0]  =   -0.968160239507626;
    points[  1]  =   -0.836031107326636;
    points[  2]  =   -0.613371432700590;
    points[  3]  =   -0.324253423403809;
    points[  4]  =    0.000000000000000;
    points[  5]  =    0.324253423403809;
    points[  6]  =    0.613371432700591;
    points[  7]  =    0.836031107326636;
    points[  8]  =    0.968160239507626;
  /*  Gauss-weights: */
    weights[  0] =    0.081274388361574;
    weights[  1] =    0.180648160694858;
    weights[  2] =    0.260610696402936;
    weights[  3] =    0.312347077040003;
    weights[  4] =    0.330239355001260;
    weights[  5] =    0.312347077040003;
    weights[  6] =    0.260610696402936;
    weights[  7] =    0.180648160694857;
    weights[  8] =    0.081274388361574;
    break;
  case  10:
  /* Gauss-Legendre rule, order= 10   */
  /*  Gauss-points: */
    points[  0]  =   -0.973906528517172;
    points[  1]  =   -0.865063366688985;
    points[  2]  =   -0.679409568299025;
    points[  3]  =   -0.433395394129247;
    points[  4]  =   -0.148874338981631;
    points[  5]  =    0.148874338981631;
    points[  6]  =    0.433395394129248;
    points[  7]  =    0.679409568299025;
    points[  8]  =    0.865063366688985;
    points[  9]  =    0.973906528517172;
  /*  Gauss-weights: */
    weights[  0] =    0.066671344308688;
    weights[  1] =    0.149451349150581;
    weights[  2] =    0.219086362515982;
    weights[  3] =    0.269266719309997;
    weights[  4] =    0.295524224714754;
    weights[  5] =    0.295524224714753;
    weights[  6] =    0.269266719309996;
    weights[  7] =    0.219086362515981;
    weights[  8] =    0.149451349150581;
    weights[  9] =    0.066671344308689;
    break;
  case  11:
  /* Gauss-Legendre rule, order= 11   */
  /*  Gauss-points: */
    points[  0]  =   -0.978228658146057;
    points[  1]  =   -0.887062599768095;
    points[  2]  =   -0.730152005574049;
    points[  3]  =   -0.519096129206812;
    points[  4]  =   -0.269543155952345;
    points[  5]  =    0.000000000000000;
    points[  6]  =    0.269543155952345;
    points[  7]  =    0.519096129206812;
    points[  8]  =    0.730152005574049;
    points[  9]  =    0.887062599768096;
    points[ 10]  =    0.978228658146057;
  /*  Gauss-weights: */
    weights[  0] =    0.055668567116174;
    weights[  1] =    0.125580369464905;
    weights[  2] =    0.186290210927734;
    weights[  3] =    0.233193764591991;
    weights[  4] =    0.262804544510247;
    weights[  5] =    0.272925086777900;
    weights[  6] =    0.262804544510246;
    weights[  7] =    0.233193764591990;
    weights[  8] =    0.186290210927735;
    weights[  9] =    0.125580369464905;
    weights[ 10] =    0.055668567116174;
    break;
  case  12:
  /* Gauss-Legendre rule, order= 12   */
  /*  Gauss-points: */
    points[  0]  =   -0.981560634246719;
    points[  1]  =   -0.904117256370475;
    points[  2]  =   -0.769902674194305;
    points[  3]  =   -0.587317954286618;
    points[  4]  =   -0.367831498998180;
    points[  5]  =   -0.125233408511469;
    points[  6]  =    0.125233408511469;
    points[  7]  =    0.367831498998180;
    points[  8]  =    0.587317954286617;
    points[  9]  =    0.769902674194305;
    points[ 10]  =    0.904117256370475;
    points[ 11]  =    0.981560634246719;
  /*  Gauss-weights: */
    weights[  0] =    0.047175336386512;
    weights[  1] =    0.106939325995318;
    weights[  2] =    0.160078328543346;
    weights[  3] =    0.203167426723066;
    weights[  4] =    0.233492536538355;
    weights[  5] =    0.249147045813402;
    weights[  6] =    0.249147045813403;
    weights[  7] =    0.233492536538354;
    weights[  8] =    0.203167426723065;
    weights[  9] =    0.160078328543347;
    weights[ 10] =    0.106939325995318;
    weights[ 11] =    0.047175336386512;
    break;
  case  13:
  /* Gauss-Legendre rule, order= 13   */
  /*  Gauss-points: */
    points[  0]  =   -0.984183054718588;
    points[  1]  =   -0.917598399222978;
    points[  2]  =   -0.801578090733310;
    points[  3]  =   -0.642349339440340;
    points[  4]  =   -0.448492751036447;
    points[  5]  =   -0.230458315955135;
    points[  6]  =    0.000000000000000;
    points[  7]  =    0.230458315955135;
    points[  8]  =    0.448492751036447;
    points[  9]  =    0.642349339440340;
    points[ 10]  =    0.801578090733310;
    points[ 11]  =    0.917598399222978;
    points[ 12]  =    0.984183054718588;
  /*  Gauss-weights: */
    weights[  0] =    0.040484004765316;
    weights[  1] =    0.092121499837728;
    weights[  2] =    0.138873510219786;
    weights[  3] =    0.178145980761952;
    weights[  4] =    0.207816047536885;
    weights[  5] =    0.226283180262897;
    weights[  6] =    0.232551553230874;
    weights[  7] =    0.226283180262896;
    weights[  8] =    0.207816047536889;
    weights[  9] =    0.178145980761946;
    weights[ 10] =    0.138873510219788;
    weights[ 11] =    0.092121499837728;
    weights[ 12] =    0.040484004765316;
    break;
  case  14:
  /* Gauss-Legendre rule, order= 14   */
  /*  Gauss-points: */
    points[  0]  =   -0.986283808696812;
    points[  1]  =   -0.928434883663573;
    points[  2]  =   -0.827201315069765;
    points[  3]  =   -0.687292904811685;
    points[  4]  =   -0.515248636358154;
    points[  5]  =   -0.319112368927890;
    points[  6]  =   -0.108054948707344;
    points[  7]  =    0.108054948707344;
    points[  8]  =    0.319112368927890;
    points[  9]  =    0.515248636358154;
    points[ 10]  =    0.687292904811685;
    points[ 11]  =    0.827201315069765;
    points[ 12]  =    0.928434883663574;
    points[ 13]  =    0.986283808696812;
  /*  Gauss-weights: */
    weights[  0] =    0.035119460331752;
    weights[  1] =    0.080158087159760;
    weights[  2] =    0.121518570687903;
    weights[  3] =    0.157203167158194;
    weights[  4] =    0.185538397477938;
    weights[  5] =    0.205198463721295;
    weights[  6] =    0.215263853463158;
    weights[  7] =    0.215263853463158;
    weights[  8] =    0.205198463721295;
    weights[  9] =    0.185538397477938;
    weights[ 10] =    0.157203167158194;
    weights[ 11] =    0.121518570687903;
    weights[ 12] =    0.080158087159760;
    weights[ 13] =    0.035119460331752;
    break;
  case  15:
  /* Gauss-Legendre rule, order= 15   */
  /*  Gauss-points: */
    points[  0]  =   -0.987992518020486;
    points[  1]  =   -0.937273392400706;
    points[  2]  =   -0.848206583410428;
    points[  3]  =   -0.724417731360170;
    points[  4]  =   -0.570972172608539;
    points[  5]  =   -0.394151347077563;
    points[  6]  =   -0.201194093997435;
    points[  7]  =    0.000000000000000;
    points[  8]  =    0.201194093997435;
    points[  9]  =    0.394151347077564;
    points[ 10]  =    0.570972172608539;
    points[ 11]  =    0.724417731360170;
    points[ 12]  =    0.848206583410427;
    points[ 13]  =    0.937273392400706;
    points[ 14]  =    0.987992518020486;
  /*  Gauss-weights: */
    weights[  0] =    0.030753241996117;
    weights[  1] =    0.070366047488108;
    weights[  2] =    0.107159220467172;
    weights[  3] =    0.139570677926154;
    weights[  4] =    0.166269205816994;
    weights[  5] =    0.186161000015563;
    weights[  6] =    0.198431485327111;
    weights[  7] =    0.202578241925561;
    weights[  8] =    0.198431485327112;
    weights[  9] =    0.186161000015562;
    weights[ 10] =    0.166269205816993;
    weights[ 11] =    0.139570677926153;
    weights[ 12] =    0.107159220467173;
    weights[ 13] =    0.070366047488108;
    weights[ 14] =    0.030753241996118;
    break;
  case  16:
  /* Gauss-Legendre rule, order= 16   */
  /*  Gauss-points: */
    points[  0]  =   -0.989400934991650;
    points[  1]  =   -0.944575023073232;
    points[  2]  =   -0.865631202387832;
    points[  3]  =   -0.755404408355003;
    points[  4]  =   -0.617876244402644;
    points[  5]  =   -0.458016777657227;
    points[  6]  =   -0.281603550779259;
    points[  7]  =   -0.095012509837637;
    points[  8]  =    0.095012509837637;
    points[  9]  =    0.281603550779259;
    points[ 10]  =    0.458016777657227;
    points[ 11]  =    0.617876244402644;
    points[ 12]  =    0.755404408355003;
    points[ 13]  =    0.865631202387832;
    points[ 14]  =    0.944575023073232;
    points[ 15]  =    0.989400934991650;
  /*  Gauss-weights: */
    weights[  0] =    0.027152459411754;
    weights[  1] =    0.062253523938648;
    weights[  2] =    0.095158511682492;
    weights[  3] =    0.124628971255534;
    weights[  4] =    0.149595988816577;
    weights[  5] =    0.169156519395003;
    weights[  6] =    0.182603415044923;
    weights[  7] =    0.189450610455068;
    weights[  8] =    0.189450610455069;
    weights[  9] =    0.182603415044923;
    weights[ 10] =    0.169156519395003;
    weights[ 11] =    0.149595988816577;
    weights[ 12] =    0.124628971255535;
    weights[ 13] =    0.095158511682493;
    weights[ 14] =    0.062253523938647;
    weights[ 15] =    0.027152459411755;
    break;
  case  17:
  /* Gauss-Legendre rule, order= 17   */
  /*  Gauss-points: */
    points[  0]  =   -0.990575475314417;
    points[  1]  =   -0.950675521768768;
    points[  2]  =   -0.880239153726986;
    points[  3]  =   -0.781514003896801;
    points[  4]  =   -0.657671159216691;
    points[  5]  =   -0.512690537086477;
    points[  6]  =   -0.351231763453876;
    points[  7]  =   -0.178484181495848;
    points[  8]  =    0.000000000000000;
    points[  9]  =    0.178484181495848;
    points[ 10]  =    0.351231763453876;
    points[ 11]  =    0.512690537086477;
    points[ 12]  =    0.657671159216691;
    points[ 13]  =    0.781514003896801;
    points[ 14]  =    0.880239153726986;
    points[ 15]  =    0.950675521768767;
    points[ 16]  =    0.990575475314417;
  /*  Gauss-weights: */
    weights[  0] =    0.024148302868548;
    weights[  1] =    0.055459529373987;
    weights[  2] =    0.085036148317178;
    weights[  3] =    0.111883847193404;
    weights[  4] =    0.135136368468526;
    weights[  5] =    0.154045761076809;
    weights[  6] =    0.168004102156450;
    weights[  7] =    0.176562705366993;
    weights[  8] =    0.179446470356207;
    weights[  9] =    0.176562705366993;
    weights[ 10] =    0.168004102156450;
    weights[ 11] =    0.154045761076810;
    weights[ 12] =    0.135136368468526;
    weights[ 13] =    0.111883847193404;
    weights[ 14] =    0.085036148317180;
    weights[ 15] =    0.055459529373988;
    weights[ 16] =    0.024148302868548;
    break;
  case  18:
  /* Gauss-Legendre rule, order= 18   */
  /*  Gauss-points: */
    points[  0]  =   -0.991565168420931;
    points[  1]  =   -0.955823949571398;
    points[  2]  =   -0.892602466497556;
    points[  3]  =   -0.803704958972523;
    points[  4]  =   -0.691687043060353;
    points[  5]  =   -0.559770831073948;
    points[  6]  =   -0.411751161462843;
    points[  7]  =   -0.251886225691506;
    points[  8]  =   -0.084775013041736;
    points[  9]  =    0.084775013041735;
    points[ 10]  =    0.251886225691505;
    points[ 11]  =    0.411751161462843;
    points[ 12]  =    0.559770831073948;
    points[ 13]  =    0.691687043060353;
    points[ 14]  =    0.803704958972523;
    points[ 15]  =    0.892602466497556;
    points[ 16]  =    0.955823949571398;
    points[ 17]  =    0.991565168420931;
  /*  Gauss-weights: */
    weights[  0] =    0.021616013526483;
    weights[  1] =    0.049714548894970;
    weights[  2] =    0.076425730254889;
    weights[  3] =    0.100942044106287;
    weights[  4] =    0.122555206711479;
    weights[  5] =    0.140642914670650;
    weights[  6] =    0.154684675126265;
    weights[  7] =    0.164276483745833;
    weights[  8] =    0.169142382963143;
    weights[  9] =    0.169142382963144;
    weights[ 10] =    0.164276483745834;
    weights[ 11] =    0.154684675126265;
    weights[ 12] =    0.140642914670650;
    weights[ 13] =    0.122555206711478;
    weights[ 14] =    0.100942044106288;
    weights[ 15] =    0.076425730254889;
    weights[ 16] =    0.049714548894970;
    weights[ 17] =    0.021616013526483;
    break;
  case  19:
  /* Gauss-Legendre rule, order= 19   */
  /*  Gauss-points: */
    points[  0]  =   -0.992406843843584;
    points[  1]  =   -0.960208152134830;
    points[  2]  =   -0.903155903614818;
    points[  3]  =   -0.822714656537143;
    points[  4]  =   -0.720966177335229;
    points[  5]  =   -0.600545304661681;
    points[  6]  =   -0.464570741375961;
    points[  7]  =   -0.316564099963630;
    points[  8]  =   -0.160358645640226;
    points[  9]  =    0.000000000000000;
    points[ 10]  =    0.160358645640225;
    points[ 11]  =    0.316564099963630;
    points[ 12]  =    0.464570741375961;
    points[ 13]  =    0.600545304661681;
    points[ 14]  =    0.720966177335229;
    points[ 15]  =    0.822714656537143;
    points[ 16]  =    0.903155903614818;
    points[ 17]  =    0.960208152134830;
    points[ 18]  =    0.992406843843585;
  /*  Gauss-weights: */
    weights[  0] =    0.019461788229726;
    weights[  1] =    0.044814226765699;
    weights[  2] =    0.069044542737642;
    weights[  3] =    0.091490021622450;
    weights[  4] =    0.111566645547334;
    weights[  5] =    0.128753962539337;
    weights[  6] =    0.142606702173606;
    weights[  7] =    0.152766042065860;
    weights[  8] =    0.158968843393954;
    weights[  9] =    0.161054449848784;
    weights[ 10] =    0.158968843393954;
    weights[ 11] =    0.152766042065859;
    weights[ 12] =    0.142606702173607;
    weights[ 13] =    0.128753962539336;
    weights[ 14] =    0.111566645547334;
    weights[ 15] =    0.091490021622450;
    weights[ 16] =    0.069044542737641;
    weights[ 17] =    0.044814226765699;
    weights[ 18] =    0.019461788229727;
    break;
  case  20:
  /* Gauss-Legendre rule, order= 20   */
  /*  Gauss-points: */
    points[  0]  =   -0.993128599185095;
    points[  1]  =   -0.963971927277914;
    points[  2]  =   -0.912234428251326;
    points[  3]  =   -0.839116971822219;
    points[  4]  =   -0.746331906460151;
    points[  5]  =   -0.636053680726515;
    points[  6]  =   -0.510867001950827;
    points[  7]  =   -0.373706088715419;
    points[  8]  =   -0.227785851141645;
    points[  9]  =   -0.076526521133497;
    points[ 10]  =    0.076526521133498;
    points[ 11]  =    0.227785851141645;
    points[ 12]  =    0.373706088715420;
    points[ 13]  =    0.510867001950827;
    points[ 14]  =    0.636053680726515;
    points[ 15]  =    0.746331906460151;
    points[ 16]  =    0.839116971822219;
    points[ 17]  =    0.912234428251326;
    points[ 18]  =    0.963971927277914;
    points[ 19]  =    0.993128599185095;
  /*  Gauss-weights: */
    weights[  0] =    0.017614007139152;
    weights[  1] =    0.040601429800388;
    weights[  2] =    0.062672048334108;
    weights[  3] =    0.083276741576705;
    weights[  4] =    0.101930119817240;
    weights[  5] =    0.118194531961518;
    weights[  6] =    0.131688638449177;
    weights[  7] =    0.142096109318382;
    weights[  8] =    0.149172986472604;
    weights[  9] =    0.152753387130726;
    weights[ 10] =    0.152753387130726;
    weights[ 11] =    0.149172986472605;
    weights[ 12] =    0.142096109318382;
    weights[ 13] =    0.131688638449177;
    weights[ 14] =    0.118194531961519;
    weights[ 15] =    0.101930119817241;
    weights[ 16] =    0.083276741576705;
    weights[ 17] =    0.062672048334109;
    weights[ 18] =    0.040601429800386;
    weights[ 19] =    0.017614007139153;
    break;
  case  21:
  /* Gauss-Legendre rule, order= 21   */
  /*  Gauss-points: */
    points[  0]  =   -0.993752170620389;
    points[  1]  =   -0.967226838566306;
    points[  2]  =   -0.920099334150400;
    points[  3]  =   -0.853363364583317;
    points[  4]  =   -0.768439963475678;
    points[  5]  =   -0.667138804197413;
    points[  6]  =   -0.551618835887220;
    points[  7]  =   -0.424342120207439;
    points[  8]  =   -0.288021316802401;
    points[  9]  =   -0.145561854160895;
    points[ 10]  =    0.000000000000000;
    points[ 11]  =    0.145561854160895;
    points[ 12]  =    0.288021316802401;
    points[ 13]  =    0.424342120207439;
    points[ 14]  =    0.551618835887220;
    points[ 15]  =    0.667138804197412;
    points[ 16]  =    0.768439963475678;
    points[ 17]  =    0.853363364583317;
    points[ 18]  =    0.920099334150401;
    points[ 19]  =    0.967226838566306;
    points[ 20]  =    0.993752170620389;
  /*  Gauss-weights: */
    weights[  0] =    0.016017228257775;
    weights[  1] =    0.036953789770852;
    weights[  2] =    0.057134425426858;
    weights[  3] =    0.076100113628379;
    weights[  4] =    0.093444423456028;
    weights[  5] =    0.108797299167148;
    weights[  6] =    0.121831416053729;
    weights[  7] =    0.132268938633345;
    weights[  8] =    0.139887394791072;
    weights[  9] =    0.144524403989970;
    weights[ 10] =    0.146081133649691;
    weights[ 11] =    0.144524403989970;
    weights[ 12] =    0.139887394791074;
    weights[ 13] =    0.132268938633337;
    weights[ 14] =    0.121831416053729;
    weights[ 15] =    0.108797299167148;
    weights[ 16] =    0.093444423456033;
    weights[ 17] =    0.076100113628379;
    weights[ 18] =    0.057134425426858;
    weights[ 19] =    0.036953789770853;
    weights[ 20] =    0.016017228257774;
    break;
  case  22:
  /* Gauss-Legendre rule, order= 22   */
  /*  Gauss-points: */
    points[  0]  =   -0.994294585482399;
    points[  1]  =   -0.970060497835428;
    points[  2]  =   -0.926956772187174;
    points[  3]  =   -0.865812577720300;
    points[  4]  =   -0.787816805979208;
    points[  5]  =   -0.694487263186683;
    points[  6]  =   -0.587640403506910;
    points[  7]  =   -0.469355837986757;
    points[  8]  =   -0.341935820892084;
    points[  9]  =   -0.207860426688221;
    points[ 10]  =   -0.069739273319723;
    points[ 11]  =    0.069739273319722;
    points[ 12]  =    0.207860426688221;
    points[ 13]  =    0.341935820892084;
    points[ 14]  =    0.469355837986757;
    points[ 15]  =    0.587640403506911;
    points[ 16]  =    0.694487263186682;
    points[ 17]  =    0.787816805979208;
    points[ 18]  =    0.865812577720300;
    points[ 19]  =    0.926956772187174;
    points[ 20]  =    0.970060497835428;
    points[ 21]  =    0.994294585482399;
  /*  Gauss-weights: */
    weights[  0] =    0.014627995298272;
    weights[  1] =    0.033774901584814;
    weights[  2] =    0.052293335152684;
    weights[  3] =    0.069796468424521;
    weights[  4] =    0.085941606217068;
    weights[  5] =    0.100414144442881;
    weights[  6] =    0.112932296080540;
    weights[  7] =    0.123252376810512;
    weights[  8] =    0.131173504787062;
    weights[  9] =    0.136541498346015;
    weights[ 10] =    0.139251872855632;
    weights[ 11] =    0.139251872855632;
    weights[ 12] =    0.136541498346015;
    weights[ 13] =    0.131173504787062;
    weights[ 14] =    0.123252376810513;
    weights[ 15] =    0.112932296080540;
    weights[ 16] =    0.100414144442881;
    weights[ 17] =    0.085941606217067;
    weights[ 18] =    0.069796468424520;
    weights[ 19] =    0.052293335152683;
    weights[ 20] =    0.033774901584814;
    weights[ 21] =    0.014627995298272;
    break;
  case  23:
  /* Gauss-Legendre rule, order= 23   */
  /*  Gauss-points: */
    points[  0]  =   -0.994769334997552;
    points[  1]  =   -0.972542471218115;
    points[  2]  =   -0.932971086826016;
    points[  3]  =   -0.876752358270442;
    points[  4]  =   -0.804888401618840;
    points[  5]  =   -0.718661363131951;
    points[  6]  =   -0.619609875763646;
    points[  7]  =   -0.509501477846008;
    points[  8]  =   -0.390301038030291;
    points[  9]  =   -0.264135680970345;
    points[ 10]  =   -0.133256824298466;
    points[ 11]  =    0.000000000000000;
    points[ 12]  =    0.133256824298466;
    points[ 13]  =    0.264135680970345;
    points[ 14]  =    0.390301038030291;
    points[ 15]  =    0.509501477846007;
    points[ 16]  =    0.619609875763646;
    points[ 17]  =    0.718661363131951;
    points[ 18]  =    0.804888401618840;
    points[ 19]  =    0.876752358270442;
    points[ 20]  =    0.932971086826016;
    points[ 21]  =    0.972542471218115;
    points[ 22]  =    0.994769334997552;
  /*  Gauss-weights: */
    weights[  0] =    0.013411859487142;
    weights[  1] =    0.030988005856979;
    weights[  2] =    0.048037671731085;
    weights[  3] =    0.064232421408526;
    weights[  4] =    0.079281411776719;
    weights[  5] =    0.092915766060035;
    weights[  6] =    0.104892091464541;
    weights[  7] =    0.114996640222411;
    weights[  8] =    0.123049084306730;
    weights[  9] =    0.128905722188082;
    weights[ 10] =    0.132462039404696;
    weights[ 11] =    0.133654572186107;
    weights[ 12] =    0.132462039404696;
    weights[ 13] =    0.128905722188082;
    weights[ 14] =    0.123049084306729;
    weights[ 15] =    0.114996640222411;
    weights[ 16] =    0.104892091464541;
    weights[ 17] =    0.092915766060034;
    weights[ 18] =    0.079281411776719;
    weights[ 19] =    0.064232421408526;
    weights[ 20] =    0.048037671731085;
    weights[ 21] =    0.030988005856980;
    weights[ 22] =    0.013411859487142;
    break;
  case  24:
  /* Gauss-Legendre rule, order= 24   */
  /*  Gauss-points: */
    points[  0]  =   -0.995187219997022;
    points[  1]  =   -0.974728555971309;
    points[  2]  =   -0.938274552002733;
    points[  3]  =   -0.886415527004401;
    points[  4]  =   -0.820001985973903;
    points[  5]  =   -0.740124191578554;
    points[  6]  =   -0.648093651936975;
    points[  7]  =   -0.545421471388839;
    points[  8]  =   -0.433793507626045;
    points[  9]  =   -0.315042679696163;
    points[ 10]  =   -0.191118867473617;
    points[ 11]  =   -0.064056892862605;
    points[ 12]  =    0.064056892862606;
    points[ 13]  =    0.191118867473616;
    points[ 14]  =    0.315042679696164;
    points[ 15]  =    0.433793507626045;
    points[ 16]  =    0.545421471388840;
    points[ 17]  =    0.648093651936976;
    points[ 18]  =    0.740124191578554;
    points[ 19]  =    0.820001985973903;
    points[ 20]  =    0.886415527004401;
    points[ 21]  =    0.938274552002733;
    points[ 22]  =    0.974728555971309;
    points[ 23]  =    0.995187219997021;
  /*  Gauss-weights: */
    weights[  0] =    0.012341229799987;
    weights[  1] =    0.028531388628934;
    weights[  2] =    0.044277438817419;
    weights[  3] =    0.059298584915437;
    weights[  4] =    0.073346481411081;
    weights[  5] =    0.086190161531953;
    weights[  6] =    0.097618652104114;
    weights[  7] =    0.107444270115966;
    weights[  8] =    0.115505668053725;
    weights[  9] =    0.121670472927803;
    weights[ 10] =    0.125837456346829;
    weights[ 11] =    0.127938195346752;
    weights[ 12] =    0.127938195346751;
    weights[ 13] =    0.125837456346829;
    weights[ 14] =    0.121670472927803;
    weights[ 15] =    0.115505668053726;
    weights[ 16] =    0.107444270115965;
    weights[ 17] =    0.097618652104115;
    weights[ 18] =    0.086190161531954;
    weights[ 19] =    0.073346481411081;
    weights[ 20] =    0.059298584915437;
    weights[ 21] =    0.044277438817421;
    weights[ 22] =    0.028531388628933;
    weights[ 23] =    0.012341229799987;
    break;
  /* END COMPUTER GENERATED SOURCE */
  default:
    printf("%s: Unknown order %d. Please use 1 through %d\n",
	   fctName,order,maxorder);
  } /* End switch */
    return;
}

/* 
* ==================================================================== 
* Utility to write a line of marks on screen to mark the progress
* of a part of the code.
* ==================================================================== */
void plotProgressMark(int iProgress, int iReset,
		      int iFirst,int iLast,
		      int iMarkLast,char mark){
  static int nMarks=62, iMark;
  int i;
  if (iReset>0){
    /* Set number of marks */
    /*nMarks = iReset;  Number of marks to set                         */
    iMark = 0;      /* Number of marks set already                     */
    if (iMarkLast){
      /* Print how many marks to expect                                */
      for (i=0;i<nMarks;i++) fprintf(stdout," ");
      fprintf(stdout,"v\n");
      fflush(stdout);
    }
  }
  else if (iReset<0){
    /* End of progress. Print new line. */
    printf("\n");
    
  }
  else {
    /* This is the standard. Plot a single mark to show the progress.  */
    if ((double)(iProgress-iFirst) /(double)(iLast-1) 
	>= (double)iMark/(double)nMarks){
      fprintf(stdout,"%c",mark);
      fflush(stdout);
      iMark++;
    }
  }
  return;
}

/* 
* ==================================================================== 
* Utility to write the memory consumption for a few chosen types to 
* screen.
* ==================================================================== */
void plotSizeOf(void){
  printf("plotSizeOf reporting memory needs for a few "
	 "commonly used types\n");
  printf(" INT                takes %2d bytes\n",sizeof(int));
  printf(" LONG INT           takes %2d bytes\n",sizeof(long int));
  printf(" SHORT INT          takes %2d bytes\n",sizeof(short int));
  printf(" UNSIGNED SHORT INT takes %2d bytes\n",sizeof(unsigned short int));
  printf(" FLOAT              takes %2d bytes\n",sizeof(float));
  printf(" DOUBLE             takes %2d bytes\n",sizeof(double));
  printf(" LONG DOUBLE        takes %2d bytes\n",sizeof(long double));
  printf(" VOID *             takes %2d bytes\n",sizeof(void *));
  return;
}
