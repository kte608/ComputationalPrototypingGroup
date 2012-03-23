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
#include <stdio.h>
#include <math.h>
void getEpsBnds(); 
int disRect(); 

#define DEFWID 1.0		/* default wire width */
#define DEFEFR 0.1		/* default edge-cell-width/inner-cell-width */
#define DEFNCL 3		/* default #cells on short side of faces */
#define DEFNWI 2		/* default problem is DEFNWI X DEFNWI bus xg */

#define CAPFRT 0		/* flags used in disSnakeLink() */
#define CAPEND 1
#define CAPALL 2
#define CAPNON 3

#define TRUE 1
#define FALSE 0

/*
  writes a quickif.c format dicretization of a bar given by 4 corners
  - corners must all be one edge away from corner1
  - returns the number of panels used
  - bar is set up to be part of a snakey bus that intertwines w/ straight bus
  - the short side with c1 on it is not paneled
  - nor is the section face in the c1,c2,c3 plane farthest from c1
  - this allows the bars to be butted together to make the snake conductors
*/
int disSnakeBar(fp, cond, edgefrac, ncells,
	     x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4)
int ncells;			/* number of cells on short side, > 2 */
int cond;			/* conductor number */
double edgefrac;       		/* edge cell widths =edgefrac*(inner widths) */
double x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4; /* 4 corners */
FILE *fp;
{
  int nsec = 4, sec, shortside, npanels = 0, i;
  double x12, y12, z12, x13, y13, z13, x14, y14, z14;
  double rx1, ry1, rz1, rx2, ry2, rz2, rx3, ry3, rz3, rx4, ry4, rz4;
  double ratio;

  /* each pass through the loop makes panels for two long faces opposite
     from each other - the ends are done after the loop */
  for(i = 0; i < 2; i++) {
    /* figure vectors to other corners relative to c1 */
    x12 = x2 - x1; y12 = y2 - y1; z12 = z2 - z1;
    x13 = x3 - x1; y13 = y3 - y1; z13 = z3 - z1;
    x14 = x4 - x1; y14 = y4 - y1; z14 = z4 - z1;

    /* figure which of the "bottom" sides is shorter */
    if(x12*x12 + y12*y12 + z12*z12 > x13*x13 + y13*y13 + z13*z13) shortside =3;
    else shortside = 2;

    /* discretize the two long faces - will be over/under 'wires' bars 
       so break into 'nsec' sections and discretize each */
    for(sec = 0; sec < nsec; sec++) {
      ratio = ((double)sec)/((double)nsec);
      if(shortside == 3) {
	rx1 = x1+ratio*x12; ry1 = y1+ratio*y12; rz1 = z1+ratio*z12;
	rx2 = rx1+x12/nsec; ry2 = ry1+y12/nsec; rz2 = rz1+z12/nsec;
	rx3 = rx1+x12/nsec+x13; ry3 = ry1+y12/nsec+y13; rz3 = rz1+z12/nsec+z13;
	rx4 = rx1+x13; ry4 = ry1+y13; rz4 = rz1+z13;
      }
      else {
	rx1 = x1+ratio*x13; ry1 = y1+ratio*y13; rz1 = z1+ratio*z13;
	rx2 = rx1+x12; ry2 = ry1+y12; rz2 = rz1+z12;
	rx3 = rx1+x13/nsec+x12; ry3 = ry1+y13/nsec+y12; rz3 =rz1+z13/nsec+z12;
	rx4 = rx1+x13/nsec; ry4 = ry1+y13/nsec; rz4 = rz1+z13/nsec;
      }
      /* face on reference side - don't panel last ref side section face */
      if(i == 1 || sec != nsec-1)
	  npanels += disRect(fp, cond, edgefrac, ncells,
			     rx1, ry1, rz1, rx2, ry2, rz2, 
			     rx3, ry3, rz3, rx4, ry4, rz4);
      /* opposite face */
      npanels += disRect(fp, cond, edgefrac, ncells,
			 rx1+x14, ry1+y14, rz1+z14, rx2+x14, ry2+y14, rz2+z14, 
			 rx3+x14, ry3+y14, rz3+z14, rx4+x14, ry4+y14, rz4+z14);
    }

    /* rotate reference coordinates to set up for last two long faces */
    if(shortside == 3) {
      x3 = x1; y3 = y1; z3 = z1;
      x1 = x4; y1 = y4; z1 = z4;
      x2 += x14; y2 += y14; z2 += z14;
      x4 += x13; y4 += y13; z4 += z13;
    }
    else {
      x2 = x1; y2 = y1; z2 = z1;
      x1 = x4; y1 = y4; z1 = z4;
      x3 += x14; y3 += y14; z3 += z14;
      x4 += x12; y4 += y12; z4 += z12;
    }
  }

  /* panel the end that doesnt have c1 on it */
  if(shortside == 3) {
    npanels += disRect(fp, cond, edgefrac, ncells,
		       x1+x12, y1+y12, z1+z12, x3+x12, y3+y12, z3+z12,
		       x4+x13+x12, y4+y13+y12, z4+z13+z12, 
		       x4+x12, y4+y12, z4+z12);
  }
  else {
    npanels += disRect(fp, cond, edgefrac, ncells,
		       x1+x13, y1+y13, z1+z13, x2+x13, y2+y13, z2+z13,
		       x4+x13+x12, y4+y13+y12, z4+z13+z12, 
		       x4+x13, y4+y13, z4+z13);
  }

  return(npanels);
}

/*
  writes a quickif.c format discretization of on link of a snakey conductor
  - c1, c2, c3, c4 are the corners of the first of 4 disSnakeBar's
  - the face determined by c1, c2, c3 is where rest of link is hooked on
  - ends are capped depending on 'capflag'
    CAPFRT => cap only front (c1, c4, c2 face - assuming |c1-c2| < |c1-c3|)
    CAPEND => cap only end (after U-shaped piece)
    CAPALL => cap both ends
    CAPNON => cap neither end
  - although disSnakeBar will accept either c1-c2 or c1-c3 as the short face,
    c1-c2 is assumed to be always shorter in this function
  - returns the number of panels produced
*/ 
int disSnakeLink(fp, cond, edgefrac, ncells, capflag,
		 x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4)
int ncells;			/* number of cells on short side, > 2 */
int cond;			/* conductor number */
/*int wires;			/* bar to be used in a wiresXwires xing */
int capflag;			/* tells if ends should be capped */
double edgefrac;       		/* edge cell widths =edgefrac*(inner widths) */
double x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4; /* 4 corners */
FILE *fp;
{
  int npanels = 0;
  double x12, y12, z12, x13, y13, z13, x14, y14, z14;
  double x1t, y1t, z1t, x2t, y2t, z2t, x3t, y3t, z3t, x4t, y4t, z4t;
  double xt, yt, zt;		/* translation vector */
  double norm, temp;

  /* figure vectors to other corners relative to c1 */
  x12 = x2 - x1; y12 = y2 - y1; z12 = z2 - z1;
  x13 = x3 - x1; y13 = y3 - y1; z13 = z3 - z1;
  x14 = x4 - x1; y14 = y4 - y1; z14 = z4 - z1;

  /* generate the four bars that make up the link */

  /* first "horizontal" bar assuming this shape: _|~| */
  npanels = disSnakeBar(fp, cond, edgefrac, ncells,
			x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);

  /* top "horizonal" bar - rotate first bar 180, then xlate over and up */
  /* rotate */
  x1t = x1+x12+x14; y1t = y1+y12+y14; z1t = z1+z12+z14;
  x3t = x3+x12+x14; y3t = y3+y12+y14; z3t = z3+z12+z14;
  x2t = x4; y2t = y4; z2t = z4;
  x4t = x2; y4t = y2; z4t = z2;
  /* translate "up" in negative c4-c1 direction, length |c3-c1|
     and "over" in positive c3-c1 direction, length |c3-c1| */
  norm = sqrt(x14*x14+y14*y14+z14*z14);	/* |c4-c1| */
  temp = sqrt(x13*x13+y13*y13+z13*z13);	/* |c3-c1| */
  xt = x13 - x14*temp/norm; 
  yt = y13 - y14*temp/norm; 
  zt = z13 - z14*temp/norm;
  x1t += xt; y1t += yt; z1t += zt;
  x2t += xt; y2t += yt; z2t += zt;
  x3t += xt; y3t += yt; z3t += zt;
  x4t += xt; y4t += yt; z4t += zt;
  npanels += disSnakeBar(fp, cond, edgefrac, ncells,
			 x1t, y1t, z1t, x2t, y2t, z2t, 
			 x3t, y3t, z3t, x4t, y4t, z4t);

  /* first "vertical" bar */
  x1t = x2+x13; y1t = y2+y13; z1t = z2+z13;
  x2t = x3; y2t = y3; z2t = z3;
  x3t = x4t; y3t = y4t; z3t = z4t;
  x4t = x1t-x13*norm/temp; y4t = y1t-y13*norm/temp; z4t = z1t-z13*norm/temp;
  npanels += disSnakeBar(fp, cond, edgefrac, ncells,
			 x1t, y1t, z1t, x2t, y2t, z2t, 
			 x3t, y3t, z3t, x4t, y4t, z4t);

  /* second "vertical" bar - first vert bar rotated and xlated down and over */
  /* rotate */
  x3t = x2t; y3t = y2t; z3t = z2t;
  x2t = x1t-x14*temp/norm; y2t = y1t-y14*temp/norm; z2t = z1t-z14*temp/norm;
  x1t = x3t-x14*temp/norm; y1t = y3t-y14*temp/norm; z1t = z3t-z14*temp/norm;
  x4t = x1t-x13*norm/temp; y4t = y1t-y13*norm/temp; z4t = z1t-z13*norm/temp;
  /* translate "over" in c3-c1 direction, length |c3-c1|
     and "down" in positive c4-c1 direction, length |c4-c1| */
  xt = x13 + x14;
  yt = y13 + y14;
  zt = z13 + z14;
  x1t += xt; y1t += yt; z1t += zt;
  x2t += xt; y2t += yt; z2t += zt;
  x3t += xt; y3t += yt; z3t += zt;
  x4t += xt; y4t += yt; z4t += zt;
  npanels += disSnakeBar(fp, cond, edgefrac, ncells,
			 x1t, y1t, z1t, x2t, y2t, z2t, 
			 x3t, y3t, z3t, x4t, y4t, z4t);

  /* cap front end, if required */
  if(capflag == CAPALL || capflag == CAPFRT) {
    npanels += disRect(fp, cond, edgefrac, ncells,
		       x1, y1, z1, x2, y2, z2, 
		       x1+x12+x14, y1+y12+x14, z1+z12+x14, x4, y4, z4);
  }

  /* cap end, if required */
  if(capflag == CAPALL || capflag == CAPEND) {
    npanels += disRect(fp, cond, edgefrac, ncells,
		       x3t-x14, y3t-y14, z3t-z14,
		       x3t-x14+x12, y3t-y14+y12, z3t-z14+z12,
		       x3t+x12, y3t+x12, z3t+x12,
		       x3t, y3t, z3t);
  }

  return(npanels);
}

/*
  writes a quickif.c format discretization of on link of a snakey conductor
  - does a backwards "L" shape, disSnakeLink does _|~| shape
  - c1, c2, c3, c4 are the corners of the first of 2 disSnakeBar's
  - the face determined by c1, c2, c3 is where rest of link is hooked on
  - ends are capped depending on 'capflag'
    CAPFRT => cap only front (c1, c4, c2 face - assuming |c1-c2| < |c1-c3|)
    CAPEND => cap only end 
    CAPALL => cap both ends
    CAPNON => cap neither end
  - although disSnakeBar will accept either c1-c2 or c1-c3 as the short face,
    c1-c2 is assumed to be always shorter in this function
  - returns the number of panels produced
*/ 
int disSnakeLLink(fp, cond, edgefrac, ncells, capflag,
		 x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4)
int ncells;			/* number of cells on short side, > 2 */
int cond;			/* conductor number */
/*int wires;			/* bar to be used in a wiresXwires xing */
int capflag;			/* tells if ends should be capped */
double edgefrac;       		/* edge cell widths =edgefrac*(inner widths) */
double x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4; /* 4 corners */
FILE *fp;
{
  int npanels = 0;
  double x12, y12, z12, x13, y13, z13, x14, y14, z14;
  double x1t, y1t, z1t, x2t, y2t, z2t, x3t, y3t, z3t, x4t, y4t, z4t;
  double xt, yt, zt;		/* translation vector */
  double norm, temp;

  /* figure vectors to other corners relative to c1 */
  x12 = x2 - x1; y12 = y2 - y1; z12 = z2 - z1;
  x13 = x3 - x1; y13 = y3 - y1; z13 = z3 - z1;
  x14 = x4 - x1; y14 = y4 - y1; z14 = z4 - z1;

  /* generate the four bars that make up the link */

  /* first "horizontal" bar assuming this shape: _| */
  npanels = disSnakeBar(fp, cond, edgefrac, ncells,
			x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);

  /* top "horizonal" bar - rotate first bar 180, then xlate over and up */
  /* THIS PART IS LEFT AS IS SO THAT t POINTS ARE AVAILABLE FOR BELOW */
  /* NOTE THAT disSnakeBar() CALL IS COMMENTED OUT */
  /* rotate */
  x1t = x1+x12+x14; y1t = y1+y12+y14; z1t = z1+z12+z14;
  x3t = x3+x12+x14; y3t = y3+y12+y14; z3t = z3+z12+z14;
  x2t = x4; y2t = y4; z2t = z4;
  x4t = x2; y4t = y2; z4t = z2;
  /* translate "up" in negative c4-c1 direction, length |c3-c1|
     and "over" in positive c3-c1 direction, length |c3-c1| */
  norm = sqrt(x14*x14+y14*y14+z14*z14);	/* |c4-c1| */
  temp = sqrt(x13*x13+y13*y13+z13*z13);	/* |c3-c1| */
  xt = x13 - x14*temp/norm; 
  yt = y13 - y14*temp/norm; 
  zt = z13 - z14*temp/norm;
  x1t += xt; y1t += yt; z1t += zt;
  x2t += xt; y2t += yt; z2t += zt;
  x3t += xt; y3t += yt; z3t += zt;
  x4t += xt; y4t += yt; z4t += zt;
  /*npanels += disSnakeBar(fp, cond, edgefrac, ncells,
			 x1t, y1t, z1t, x2t, y2t, z2t, 
			 x3t, y3t, z3t, x4t, y4t, z4t);*/

  /* first "vertical" bar */
  x1t = x2+x13; y1t = y2+y13; z1t = z2+z13;
  x2t = x3; y2t = y3; z2t = z3;
  x3t = x4t; y3t = y4t; z3t = z4t;
  x4t = x1t-x13*norm/temp; y4t = y1t-y13*norm/temp; z4t = z1t-z13*norm/temp;
  npanels += disSnakeBar(fp, cond, edgefrac, ncells,
			 x1t, y1t, z1t, x2t, y2t, z2t, 
			 x3t, y3t, z3t, x4t, y4t, z4t);

#if 1 == 0
  /* NOT USED HERE */
  /* second "vertical" bar - first vert bar rotated and xlated down and over */
  /* rotate */
  x3t = x2t; y3t = y2t; z3t = z2t;
  x2t = x1t-x14*temp/norm; y2t = y1t-y14*temp/norm; z2t = z1t-z14*temp/norm;
  x1t = x3t-x14*temp/norm; y1t = y3t-y14*temp/norm; z1t = z3t-z14*temp/norm;
  x4t = x1t-x13*norm/temp; y4t = y1t-y13*norm/temp; z4t = z1t-z13*norm/temp;
  /* translate "over" in c3-c1 direction, length |c3-c1|
     and "down" in positive c4-c1 direction, length |c4-c1| */
  xt = x13 + x14;
  yt = y13 + y14;
  zt = z13 + z14;
  x1t += xt; y1t += yt; z1t += zt;
  x2t += xt; y2t += yt; z2t += zt;
  x3t += xt; y3t += yt; z3t += zt;
  x4t += xt; y4t += yt; z4t += zt;
  npanels += disSnakeBar(fp, cond, edgefrac, ncells,
			 x1t, y1t, z1t, x2t, y2t, z2t, 
			 x3t, y3t, z3t, x4t, y4t, z4t);
#endif

  /* cap front end, if required */
  if(capflag == CAPALL || capflag == CAPFRT) {
    npanels += disRect(fp, cond, edgefrac, ncells,
		       x1, y1, z1, x2, y2, z2, 
		       x1+x12+x14, y1+y12+x14, z1+z12+x14, x4, y4, z4);
  }

  /* cap end, if required */
  if(capflag == CAPALL || capflag == CAPEND) {
    npanels += disRect(fp, cond, edgefrac, ncells,
		       x3t, y3t, z3t,
		       x3t-x12, y3t-y12, z3t-z12,
		       x3t+x14-x12, y3t+y14-y12, z3t+z14-z12,
		       x3t+x14, y3t+y14, z3t+z14);
  }

  return(npanels);
}

/* 
  does the rotation/translation to set up for the next snake backwards L link
*/
void doRotTrans(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4)
double *x1, *y1, *z1, *x2, *y2, *z2, *x3, *y3, *z3, *x4, *y4, *z4;
{
  double x12, y12, z12, x13, y13, z13, x14, y14, z14;
  double x1t, y1t, z1t, x2t, y2t, z2t, x3t, y3t, z3t, x4t, y4t, z4t;
  double xt, yt, zt;		/* translation applied between snake links */
  double temp, norm;

  x12 = *x2 - *x1; y12 = *y2 - *y1; z12 = *z2 - *z1;
  x13 = *x3 - *x1; y13 = *y3 - *y1; z13 = *z3 - *z1;
  x14 = *x4 - *x1; y14 = *y4 - *y1; z14 = *z4 - *z1;
  /* set up coordinates for next link - "upward" translation */
  /* LIFTED FROM disSnakeLink() */
  /* rotate */
  x1t = *x1+x12+x14; y1t = *y1+y12+y14; z1t = *z1+z12+z14;
  x3t = *x3+x12+x14; y3t = *y3+y12+y14; z3t = *z3+z12+z14;
  x2t = *x4; y2t = *y4; z2t = *z4;
  x4t = *x2; y4t = *y2; z4t = *z2;
  /* translate "up" in negative c4-c1 direction, length |c3-c1|
     and "over" in positive c3-c1 direction, length |c3-c1| */
  norm = sqrt(x14*x14+y14*y14+z14*z14);	/* |c4-c1| */
  temp = sqrt(x13*x13+y13*y13+z13*z13);	/* |c3-c1| */
  xt = x13 - x14*temp/norm; 
  yt = y13 - y14*temp/norm; 
  zt = z13 - z14*temp/norm;
  *x1 = x1t + xt; *y1 = y1t + yt; *z1 = z1t + zt;
  *x2 = x2t + xt; *y2 = y2t + yt; *z2 = z2t + zt;
  *x3 = x3t + xt; *y3 = y3t + yt; *z3 = z3t + zt;
  *x4 = x4t + xt; *y4 = y4t + yt; *z4 = z4t + zt;
}

/*
  writes a quickif.c format discretization of a snakey conductor
  - four corners specify the first disSnakeLink
  - wires = number of wires to cross, ie number of disSnakeLinks
  - |c2-c1| < |c3-c1| assumed
  - c1, c2, c3 plane has first link hook attached onto in
*/
int disSnake(fp, cond, edgefrac, ncells, wires,
	     x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4)
int ncells;			/* number of cells on short side, > 2 */
int cond;			/* conductor number */
int wires;			/* conductor setup for a wiresXwires xing */
/*int capflag;			/* tells if ends should be capped */
double edgefrac;       		/* edge cell widths =edgefrac*(inner widths) */
double x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4; /* 4 corners */
FILE *fp;
{
  int np = 0, i;
  
  /* make the snake-shaped conductor */
  if(wires == 1) np += disSnakeLLink(fp, cond, edgefrac, ncells, CAPALL,
				     x1, y1, z1, x2, y2, z2, 
				     x3, y3, z3, x4, y4, z4);
  else np += disSnakeLLink(fp, cond, edgefrac, ncells, CAPFRT,
			   x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);
  doRotTrans(&x1, &y1, &z1, &x2, &y2, &z2, &x3, &y3, &z3, &x4, &y4, &z4);
  for(i = 1; i < wires-1; i++) {
    np += disSnakeLLink(fp, cond, edgefrac, ncells, CAPNON,
			x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);
    doRotTrans(&x1, &y1, &z1, &x2, &y2, &z2, &x3, &y3, &z3, &x4, &y4, &z4);
  }
  if(wires > 1) np += disSnakeLLink(fp, cond, edgefrac, ncells, CAPEND,
				    x1, y1, z1, x2, y2, z2, 
				    x3, y3, z3, x4, y4, z4);

  return(np);
}

/*
  writes a quickif.c format discretization of a bar given by 4 corners
  - corners must all be one edge away from corner1
  - returns the number of panels used
*/
int disBar(fp, cond, edgefrac, ncells, wires,
	     x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4)
int ncells;			/* number of cells on short side, > 2 */
int cond;			/* conductor number */
int wires;			/* bar to be used in a wiresXwires xing */
double edgefrac;       		/* edge cell widths =edgefrac*(inner widths) */
double x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4; /* 4 corners */
FILE *fp;
{
  int nsec = 2*wires + 1, sec, shortside, npanels = 0, i;
  double x12, y12, z12, x13, y13, z13, x14, y14, z14;
  double rx1, ry1, rz1, rx2, ry2, rz2, rx3, ry3, rz3, rx4, ry4, rz4;
  double ratio;

  /* each pass through the loop makes panels for two long faces opposite
     from each other - the ends are done after the loop */
  for(i = 0; i < 2; i++) {
    /* figure vectors to other corners relative to c1 */
    x12 = x2 - x1; y12 = y2 - y1; z12 = z2 - z1;
    x13 = x3 - x1; y13 = y3 - y1; z13 = z3 - z1;
    x14 = x4 - x1; y14 = y4 - y1; z14 = z4 - z1;

    /* figure which of the "bottom" sides is shorter */
    if(x12*x12 + y12*y12 + z12*z12 > x13*x13 + y13*y13 + z13*z13) shortside =3;
    else shortside = 2;

    /* discretize the two long faces - will be over/under 'wires' bars 
       so break into 'wires' sections and discretize each */
    for(sec = 0; sec < nsec; sec++) {
      ratio = ((double)sec)/((double)nsec);
      if(shortside == 3) {
	rx1 = x1+ratio*x12; ry1 = y1+ratio*y12; rz1 = z1+ratio*z12;
	rx2 = rx1+x12/nsec; ry2 = ry1+y12/nsec; rz2 = rz1+z12/nsec;
	rx3 = rx1+x12/nsec+x13; ry3 = ry1+y12/nsec+y13; rz3 = rz1+z12/nsec+z13;
	rx4 = rx1+x13; ry4 = ry1+y13; rz4 = rz1+z13;
      }
      else {
	rx1 = x1+ratio*x13; ry1 = y1+ratio*y13; rz1 = z1+ratio*z13;
	rx2 = rx1+x12; ry2 = ry1+y12; rz2 = rz1+z12;
	rx3 = rx1+x13/nsec+x12; ry3 = ry1+y13/nsec+y12; rz3 =rz1+z13/nsec+z12;
	rx4 = rx1+x13/nsec; ry4 = ry1+y13/nsec; rz4 = rz1+z13/nsec;
      }
      /* face on reference side */
      npanels += disRect(fp, cond, edgefrac, ncells,
			 rx1, ry1, rz1, rx2, ry2, rz2, 
			 rx3, ry3, rz3, rx4, ry4, rz4);
      /* opposite face */
      npanels += disRect(fp, cond, edgefrac, ncells,
			 rx1+x14, ry1+y14, rz1+z14, rx2+x14, ry2+y14, rz2+z14, 
			 rx3+x14, ry3+y14, rz3+z14, rx4+x14, ry4+y14, rz4+z14);
    }

    /* rotate reference coordinates to set up for last two long faces */
    if(shortside == 3) {
      x3 = x1; y3 = y1; z3 = z1;
      x1 = x4; y1 = y4; z1 = z4;
      x2 += x14; y2 += y14; z2 += z14;
      x4 += x13; y4 += y13; z4 += z13;
    }
    else {
      x2 = x1; y2 = y1; z2 = z1;
      x1 = x4; y1 = y4; z1 = z4;
      x3 += x14; y3 += y14; z3 += z14;
      x4 += x12; y4 += y12; z4 += z12;
    }
  }

  /* panel the ends */
  if(shortside == 3) {
    npanels += disRect(fp, cond, edgefrac, ncells,
		       x1, y1, z1, x3, y3, z3,
		       x4+x13, y4+y13, z4+z13, x4, y4, z4);
    npanels += disRect(fp, cond, edgefrac, ncells,
		       x1+x12, y1+y12, z1+z12, x3+x12, y3+y12, z3+z12,
		       x4+x13+x12, y4+y13+y12, z4+z13+z12, 
		       x4+x12, y4+y12, z4+z12);
  }
  else {
    npanels += disRect(fp, cond, edgefrac, ncells,
		       x1, y1, z1, x2, y2, z2,
		       x4+x12, y4+y12, z4+z12, x4, y4, z4);
    npanels += disRect(fp, cond, edgefrac, ncells,
		       x1+x13, y1+y13, z1+z13, x2+x13, y2+y13, z2+z13,
		       x4+x13+x12, y4+y13+y12, z4+z13+z12, 
		       x4+x13, y4+y13, z4+z13);
  }

  return(npanels);
}

/*
  generates a wires crossing wires bus crossing example in quickif.c format
  - uses disRect() for all discretization of rectangular faces
*/
main(argc, argv)
int argc;
char *argv[];
{
  char **chkp, *chk;
  int wires, npanels = 0, ncells, cmderr = FALSE, i, cond;
  double edgefrac, pitch, strtod(), tem;
  double x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4; /* 4 corners */
  long strtol();
  FILE *fp, *fopen();

  /* load default parameters */
  pitch = 4*DEFWID;
  edgefrac = DEFEFR;
  ncells = DEFNCL;
  wires = DEFNWI;

  /* parse command line */
  chkp = &chk;			/* pointers for error checking */
  for(i = 1; i < argc && cmderr == FALSE; i++) {
    if(argv[i][0] != '-') {
      fprintf(stderr, "%s: illegal argument -- %s\n", argv[0], argv[i]);
      cmderr = TRUE;
      break;
    }
    else if(argv[i][1] == 'c') {
      wires = (int) strtol(&(argv[i][2]), chkp, 10);
      if(*chkp == &(argv[i][2]) || wires < 1) {
	fprintf(stderr, "%s: bad number of conductors/bus `%s'\n", 
		argv[0], &argv[i][2]);
	cmderr = TRUE;
	break;
      }
    }
    else if(argv[i][1] == 'n') {
      ncells = (int) strtol(&(argv[i][2]), chkp, 10);
      if(*chkp == &(argv[i][2]) || ncells < 1) {
	fprintf(stderr, 
		"%s: bad number of short side panels `%s'\n", 
		argv[0], &argv[i][2]);
	cmderr = TRUE;
	break;
      }
    }
    else if(argv[i][1] == 'w') {
      pitch = 2*strtod(&(argv[i][2]), chkp);
      if(*chkp == &(argv[i][2]) || pitch <= 0.0) {
	fprintf(stderr, "%s: bad wire width `%s'\n", 
		argv[0], &argv[i][2]);
	cmderr = TRUE;
	break;
      }
    }
    else if(argv[i][1] == 'e') {
      edgefrac = strtod(&(argv[i][2]), chkp);
      if(*chkp == &(argv[i][2]) || edgefrac < 0.0) {
	fprintf(stderr, "%s: bad edge panel fraction `%s'\n", 
		argv[0], &argv[i][2]);
	cmderr = TRUE;
	break;
      }
    }
    else {
      fprintf(stderr, "%s: illegal option -- %s\n", argv[0], &(argv[i][1]));
      cmderr = TRUE;
      break;
    }
  }

  if(cmderr == TRUE) {
    fprintf(stderr,
	    "Usage: %s [-c<conductors/bus(def=%d)>] [-w<wire width(def=%.3g)>]\n       [-n<num panels/wire width(def=%d)>] [-e<rel edge panel width(def=%.3g)>]\n", 
	    argv[0], DEFNWI, DEFWID, DEFNCL, DEFEFR);
    exit(0);
  }

  /* open output file */
  fp = stdout;

  /* write title */
  fprintf(fp, "0 %dX%d woven bus problem with %.3gm wires (n=%d e=%.3g)\n",
	  wires, wires, pitch/4.0, ncells, edgefrac);

  /* set up straight bars in bus crossing geometry */
  x1 = 2.0*pitch/4.0; y1 = 0.0; z1 = 3.0*pitch/4.0; /* setup first c1 */
  y2 = (pitch/4.0)*(2*wires+1); z2 = z1;	/* constant parts of c2 */
  y3 = 0.0; z3 = z1;	/* constant parts of c3 */
  y4 = 0.0; z4 = 2.0*pitch/4.0;		/* constant parts of c4 */
  for(cond = 1; cond <= wires; cond++) {
    x2 = x4 = x1; x3 = x1 - pitch/4.0; 
    npanels += disBar(fp, cond, edgefrac, ncells, wires,
		      x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);
    x1 += pitch;
  }

  /* setup snaked bars (step in y rather than x) */
  /* first "handle" has "lower" face on the x-y plane - gets rot'ed, xlated */
  x1 = 0.0; y1 = pitch/4.0; z1 = pitch/4.0; /* setup first c1 */
  x2 = 0.0; y2 = y1+y1; z2 = z1; /* setup first c2 */
  x3 = pitch; y3 = y1; z3 = z1;	/* setup first c3 */
  x4 = 0.0; y4 = y1; z4 = 0.0; /* setup first c4 */
  for(i = 1; cond <= 2*wires; cond++, i++) {
    npanels += disSnake(fp, cond, edgefrac, ncells, wires,
			x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);
    /* rotate */
    y1 = y2; z1 = z4;
    y3 = y2; z3 = z4;
    tem = y2; y2 = y4; y4 = tem;
    tem = z2; z2 = z4; z4 = tem;
    /* translate up or down (+-z) by pitch, over (+y) by pitch/2 */
    if(i % 2 != 0) {
      z1 += pitch; z2 += pitch; z3 += pitch; z4 += pitch;
    }
    else {
      z1 -= pitch; z2 -= pitch; z3 -= pitch; z4 -= pitch;
    }
    y1 += pitch/2.0; y2 += pitch/2.0; y3 += pitch/2.0; y4 += pitch/2.0;
  }

}


#define MAX(a,b) (((a)>(b))?(a):(b))


/*
  writes quad panel lines in quickif.c format to given file pointer
  - corners must be ordered around the rectangle, not across the diagonal
    (p1 must be opposite p3 and p2 must be opposite p4)
  - tries to panel specified rectangle with uniform inner and skinny edge cells
  - makes edge cell widths percentage of inner cell widths (10% best says A.R.)
  - ncells must be greater than zero, edge frac must be non-negative
  - ways to get uniform panels: edgefrac = 1.0, edgefrac = 0.0
    both give ncells panels on a side or ncells = 2 if want 2 on a side
  - ways to get one panel the same as the rectangle: anything with ncells = 1
  - should also work with parallelograms
  - returns the number of panels made
  - no_discr 
*/
int disRect(fp, cond, edgefrac, ncells, 
	     x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4)
int ncells;			/* number of cells on short side, > 2 */
int cond;			/* conductor number */
double edgefrac;       		/* edge cell widths =edgefrac*(inner widths) */
double x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4; /* 4 corners */
FILE *fp;
{
  int lflag, linernum, sinernum, npanels;
  double lside, sside, temp, sedgewid, sinerwid, linerwid;
  static double lepsilon, uepsilon;
  static int fstflag = 1;
  int no_discr = FALSE; 

  if(fp == NULL) {
    fprintf(stderr, "\ndisRect: bad output file pointer\n");
    exit(0);
  }

  if(no_discr) {
    fprintf(fp, "Q %d %.5e %.5e %.5e  %.5e %.5e %.5e",
	    cond, x1, y1, z1, x2, y2, z2);
    fprintf(fp, " %.5e %.5e %.5e  %.5e %.5e %.5e\n",
	    x3, y3, z3, x4, y4, z4);
    return(1);
  }

  /* setup bounds on machine precision on first call */
  if(fstflag == 1) {
    fstflag = 0;
    getEpsBnds(&uepsilon, &lepsilon);
  }

  /* find the sides */
  lflag = 2;			/* implies long side is p1-p2 side */
  lside = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
  sside = sqrt((x1-x4)*(x1-x4) + (y1-y4)*(y1-y4) + (z1-z4)*(z1-z4));
  if(lside < sside) {
    temp = lside;
    lside = sside;
    sside = temp;
    lflag = 4;
  }

  /* figure the short side cell widths */
  if(ncells == 1 || ncells == 2) edgefrac = 0.0;
  if(edgefrac > 0.0) {
    sinernum = ncells-2;	/* sinernum = number of inner cells */
    sinerwid = sside/(2.0*edgefrac + (double)sinernum);
    sedgewid = edgefrac*sinerwid;
  }
  else if(edgefrac == 0.0) {
    sinernum = ncells;	/* ncells = sinernum = #nonzero width cells */
    sinerwid = sside/((double)(ncells));
    sedgewid = 0.0;
  }
  else {
    fprintf(stderr, "\ndisRect: negative edge to inner panel ratio = %g\n",
	    edgefrac);
    exit(0);
  }

  /* figure the long side inner cell widths (edge cell widths are the same) */
  temp = (lside/sinerwid - 2.0*edgefrac);
  /* allow for cancellation error in temp (underestimates only)
     - otherwise possible for truncation to eliminate a wanted panel
     - if this screws up, the worst it'll do is produce an extra inner sec'n */
  /* add scaled epsilon to result to get close numbers over trunc. threshold */
  temp += (MAX(lside/sinerwid, 2.0*edgefrac)*uepsilon);
  linernum = temp;		/* truncate */
  linerwid = (lside-2.0*sedgewid)/((double)linernum); /* long side cell wdth */

  /* write out the quad cell lines */
  if(lflag == 4) npanels = wrQuadCells(fp, cond, sinernum, linernum, 
				       sinerwid, linerwid, sedgewid, 
				       x1, y1, z1, x2, y2, z2, x4, y4, z4);
  else npanels = wrQuadCells(fp, cond, sinernum, linernum, sinerwid, linerwid, 
		   sedgewid, x1, y1, z1, x4, y4, z4, x2, y2, z2);
  fprintf(fp, "*\n");

  return(npanels);
}

#define TOL 1e-10                /* tolerance on epsilon (machine precision) */
#define UPPEPS 1e-10              /* should be bigger than epsilon */
#define LOWEPS 1e-25             /* should be smaller than epsilon */

/*
  returns upper and lower bounds on machine precision (for doubles)
  - takes into account the limitations of the memory rep of a doubles
    by forcing the compiler to put all doubles in core
*/
void getEpsBnds(upper, lower)
double *upper, *lower;            /* current upper and lower epsilon bounds */
{
  double dif, tol, mid;
  double temp, one;

  double *difp = &dif;		/* to stop optimizer from putting */
  double *tolp = &tol;		/* variables in registers (not doing so */
  double *midp = &mid;		/* can lead to undully optomistic machine */
  double *tempp = &temp;	/* precision estimate that doesnt take */
  double *onep = &one;		/* memory storage rounding into account) */

  *upper = UPPEPS;
  *lower = LOWEPS;
  *onep = 1.0;

  *difp = *upper - *lower;
  *tolp = *difp/(*lower);
  while(*tolp > TOL) {
    *midp = (*upper + *lower)/2.0;
    *tempp = 1.0 + *midp;
    if(*tempp > *onep) *upper = *midp; 
    else *lower = *midp;
    *difp = *upper - *lower;
    *tolp = *difp/(*lower);
  }
}


/*
  write a set of quadralateral panels for a rectancle in quickif.c format
  - produces edge cells that are always sedgewid in width and either
    sinerwid or linerwid in length if on short or long side of rectangle
    - those on corners are sedgewidXsedgewid
  - produces inner cells that are all sinerwid long on the short side
    of the rectangle and linerwid long on the long side
  - there are sinernum inner panels along the short egde, linernum along long
  - zero with sedgewid ok but sinerwid and linerwid must always be nonzero
  - returns number of panels generated
*/

int wrQuadCells(fp, cond, sinernum, linernum, sinerwid, linerwid, sedgewid,
		 x1, y1, z1, x2, y2, z2, x4, y4, z4)
int cond;			/* conductor number */
int sinernum, linernum;
double sinerwid, linerwid, sedgewid;
double x1, y1, z1, x2, y2, z2, x4, y4, z4;
FILE *fp;
{
  int sncell, lncell, npanels = 0;
  double x12, y12, z12, x14, y14, z14; /* unit vector from p1 to p2, p4 */
  double norm2, norm4, x, y, z, temp;
  double xt, yt, zt;

  /* set up unit vectors - p1-p2 must be the short side */
  x12 = x2 - x1;
  y12 = y2 - y1;
  z12 = z2 - z1;
  norm2 = sqrt(x12*x12 + y12*y12 + z12*z12);
  x12 /= norm2;
  y12 /= norm2;
  z12 /= norm2;

  x14 = x4 - x1;
  y14 = y4 - y1;
  z14 = z4 - z1;
  norm4 = sqrt(x14*x14 + y14*y14 + z14*z14);
  x14 /= norm4;
  y14 /= norm4;
  z14 /= norm4;

  /* loop through long and short side points and make quads */
  for(sncell = 0; sncell < sinernum+2; sncell++) {
    for(lncell = 0; lncell < linernum+2; lncell++) {
      /* dump the quad line if cell has nonzero area
         - ie if not an edge cell or if an edge cell with nonzero width */
      if(sedgewid != 0.0 ||
	 (sncell != 0 && sncell != sinernum+1 &&
	  lncell != 0 && lncell != linernum+1 && sedgewid == 0.0)) {
	/* figure the current lower left point */
	x = x1;
	y = y1;
	z = z1;
	if(sncell > 0) {		/* add in short side edge cell width */
	  x += (sedgewid*x12);
	  y += (sedgewid*y12);
	  z += (sedgewid*z12);
	  if(sncell > 1) {	/* add in short side inner cell widths */
	    x += ((double)sncell-1)*(sinerwid*x12);
	    y += ((double)sncell-1)*(sinerwid*y12);
	    z += ((double)sncell-1)*(sinerwid*z12);
	  }
	}
	if(lncell > 0) {		/* add in long side edge cell width */
	  x += (sedgewid*x14);
	  y += (sedgewid*y14);
	  z += (sedgewid*z14);
	  if(lncell > 1) {	/* add in long side inner cell widths */
	    x += ((double)lncell-1)*(linerwid*x14);
	    y += ((double)lncell-1)*(linerwid*y14);
	    z += ((double)lncell-1)*(linerwid*z14);
	  }
	}

	fprintf(fp, "Q %d ", cond);
	/* dump point 1 */
	fprintf(fp, "%.5e %.5e %.5e ", x, y, z);
	/* dump point 2 */
	if(sncell == 0 || sncell == sinernum+1) { /* if short side edge cell */
	  xt = x+sedgewid*x12;	/* t coordinates are of point 2 */
	  yt = y+sedgewid*y12;
	  zt = z+sedgewid*z12;
	  fprintf(fp, "%.5e %.5e %.5e ", xt, yt, zt);
	}
	else {			/* if a short side inner cell */
	  xt = x+sinerwid*x12;
	  yt = y+sinerwid*y12;
	  zt = z+sinerwid*z12;
	  fprintf(fp, "%.5e %.5e %.5e ", xt, yt, zt);
	}
	/* dump point 3 (across from point 1 = (x,y,z)) and point 4 */
	if(lncell == 0 || lncell == linernum+1) {  /* if long side edge cell */
	  fprintf(fp, "%.5e %.5e %.5e ",
		  xt+sedgewid*x14, yt+sedgewid*y14, zt+sedgewid*z14); /* p3 */
	  fprintf(fp, "%.5e %.5e %.5e\n",
		  x+sedgewid*x14, y+sedgewid*y14, z+sedgewid*z14); /* p4 */
	}
	else {			/* if a long side inner cell */
	  fprintf(fp, "%.5e %.5e %.5e ",
		  xt+linerwid*x14, yt+linerwid*y14, zt+linerwid*z14);
	  fprintf(fp, "%.5e %.5e %.5e\n",
		  x+linerwid*x14, y+linerwid*y14, z+linerwid*z14);
	}
	npanels++;
      }
    }
  }
  return(npanels);
}
