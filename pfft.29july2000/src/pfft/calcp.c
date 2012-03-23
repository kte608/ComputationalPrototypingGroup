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
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <float.h>

#include "Calcp.h"

#define VOP(V,V1,op,V2) \
V[0] = V1[0] op V2[0];\
V[1] = V1[1] op V2[1];\
V[2] = V1[2] op V2[2]

#define VOPS(V, V1, op, S) \
V[0] = V1[0] op S;\
V[1] = V1[1] op S;\
V[2] = V1[2] op S

#define VOPSS(V, op, S) \
V[0] = V[0] op S;\
V[1] = V[1] op S;\
V[2] = V[2] op S

#define VASSIGN(V, V1) \
V[0] = V1[0]; \
V[1] = V1[1]; \
V[2] = V1[2]

#define VPRODUCT(V1,V2) (V1[0]*V2[0]+V1[1]*V2[1]+V1[2]*V2[2])

#define VNORM(V1) sqrt(V1[0]*V1[0]+V1[1]*V1[1]+V1[2]*V1[2])

#define PI 3.141592653589793

#define VCROSS(V, V1, V2) \
V[0] = V1[1]*V2[2] - V1[2]*V2[1]; \
V[1] = V1[2]*V2[0] - V1[0]*V2[2];\
V[2] = V1[0]*V2[1] - V1[1]*V2[0]

#define EQUIV_TOL 1.0e-9

#define MIN(a,b) ( ((a)>(b)) ? (b) : (a) )

/*
% C version of calcp, returns potential at evaluation point due
% to unit monopole and unit dipole uniformly distributed on a panel.
% Follows a left-hand rule (Clockwise ordered points has normal
% pointing up).
% panel -- vectors of panel vertices in rows of x,y,z (3 or 4 rows supported).
% evalpnts -- matrix of evaluation points, rows of x,y,z coordinates
% directions -- matrix of derivative directions, rows of x,y,z coordinates
% fss = the vector of potentials due to a monopole
% fds = the vector of potentials due to a panel normal dipole distribution
% area = panel area.
% centroid = panel centroid.
% Z = panel normal.
% fess = the derivative of the monopole potential at evalpnt along direction
% fess = the derivative of the dipole potential at evalpnt along direction
*/
void calcp(panel, verts, evalpnts, directions, numeval, fss, fds, fess, feds)
int verts, numeval;
double **panel, **evalpnts, **directions, *fss, *fds, *fess, *feds;
{
  double X[3], Y[3], Z[3], V0[3], V1[3], centroid[3], nrm[3];
  double edgeLength[4], ct[4], st[4], xmxv[4], ymyv[4], fe[4], r[4];
  double minDiagLength, minSideLength, dtol;
  double xri[4], yri[4];
  double npanel[4][3];
  double xc, yc, zc, v, fln, arg, xnorm, znorm, area;
  double y1, y3, x1, x3, zn, znabs, evalDistance, fs, fd;
  double c1, s1, c2, s2, s12, c12, val;
  double fsx, fsy, fdx, fdy, fdz;
  double fac, u1, u2, rr, fh1, fh2, fes, fed;
  int OK, i, evalindex, next;

  /* FIRST: The PANEL SETUP ************************* */

  /* Length of each side and the panel area. */
  for(minSideLength=DBL_MAX, i=0; i < (verts-1); i++) {
    VOP(V0, panel[i+1], -, panel[i]);
    edgeLength[i] = VNORM(V0);
    /* Minimum length of a panel side */
    minSideLength = MIN(minSideLength,edgeLength[i]);
  }
  VOP(V0, panel[0], -, panel[i]);
  edgeLength[i] = VNORM(V0);

  /* Calculate the panel coordinate system. */
  VOP(X, panel[2], -, panel[0]);
  xnorm =  VNORM(X);
  if(verts == 3) {
    VOP(Y, panel[1], -, panel[0]); 
  } else {
    VOP(Y, panel[1], -, panel[3]); 
  }
  /* Find minimum diagonal length (side length for triagonals) */
  if(verts==3) minDiagLength = minSideLength;
  else         minDiagLength = MIN(xnorm,VNORM(Y)); /*verts==4*/ 
  /* Set in-plane tolerance */
  dtol = EQUIV_TOL * minDiagLength;
  /* Calculate Z-axis is normal to two diags. */
  VCROSS(Z, X, Y);
  znorm = VNORM(Z);
  area = 0.5 * znorm;
  /* Normalize panel axes.  */
  VOPSS(Z, /, znorm);
  VOPSS(X, /, xnorm);
  VCROSS(Y, Z, X);

  /* Determine the centroid. */
  VOP(V0, panel[1], -, panel[0]);
  if(verts == 4) {
    VOP(V1, panel[3], -, panel[0]);
  } else  {
    VOP(V1, panel[2], -, panel[0]);
  }
  y1 = VPRODUCT(V0,Y);
  y3 = VPRODUCT(V1,Y);
  x1 = VPRODUCT(V0,X);
  x3 = VPRODUCT(V1,X);
  yc = (y1 + y3)/3.0;
  xc = (xnorm + ((x1 * y1 - x3 * y3)/(y1 - y3)))/3.0;
  VOPS(V0, X, *, xc);
  VOPS(V1, Y, *, yc);
  VOP(V0, V0, +, V1);
  VOP(centroid, panel[0], +, V0);

  /* Put panel vertices in newly defined coordinate system. */
  for (i=0; i < verts; i++) {
    VOP(V0, panel[i], -, centroid);
    npanel[i][0] = VPRODUCT(X, V0);
    npanel[i][1] = VPRODUCT(Y, V0);
    npanel[i][2] = VPRODUCT(Z, V0);
  }

  /* Check that panel is in the x-y plane.  */
  for(i=0; i < verts; i++) {
    if(fabs(npanel[i][2]) > (1.0e-8 * xnorm)) {
      printf("Coordinate transform failure!!");
      exit(0);
    }
  }

  /* Compute the contributions terms for each edge. */
  for(i=0; i < (verts-1); i++) {
    ct[i] = (npanel[i+1][0]-npanel[i][0])/edgeLength[i];
    st[i] = (npanel[i+1][1]-npanel[i][1])/edgeLength[i];
  }
  ct[i] = (npanel[0][0]-npanel[i][0])/edgeLength[i];
  st[i] = (npanel[0][1]-npanel[i][1])/edgeLength[i];


  /* Done with the PANEL SETUP!!!!!!!!!************************************ */


  /* Done with Setup, now loop through the evaluation points!! */
  for(evalindex = 0; evalindex < numeval; evalindex++) {

    /* First map the evaluation point into panel's coord system. */
    VOP(V1, evalpnts[evalindex], -, centroid);
    V0[0] = VPRODUCT(X,V1);
    V0[1] = VPRODUCT(Y,V1);
    zn = V0[2] = VPRODUCT(Z,V1);
    znabs=fabs(zn); 
    /* If we are in the plane of the panl (to some accuracy), then 
     * move out of the domain (with the normal) by a very small 
     * amount. This will give a 2pi term for the self-evaluation, 
     * which then is correctly taken care of.                          */
    if (znabs<dtol) {
      zn    =  0.5*dtol;
      znabs =  0.5*dtol;
    }

    evalDistance = VNORM(V0);

    /* Get derivative direction in coordinate system */
    if(directions != NULL) {
      VASSIGN(V1, directions[evalindex]);
      nrm[0] = VPRODUCT(X,V1);
      nrm[1] = VPRODUCT(Y,V1);
      nrm[2] = VPRODUCT(Z,V1);
    }
       
    /* Once per vertex computation. */
    OK = 1;
    for(i=0; i < verts; i++) {
      xc=V0[0]-npanel[i][0];
      yc=V0[1]-npanel[i][1];
      zc=V0[2]-npanel[i][2];  
      xmxv[i]=xc; 
      ymyv[i]=yc; 
      fe[i]=xc*xc+zc*zc;
      r[i]=sqrt(yc*yc+fe[i]);
      if(r[i] < (1.005*znabs)) {
        OK = 0; 
      }
      if(directions != NULL) {
        xri[i] = xmxv[i]/r[i];
        yri[i] = ymyv[i]/r[i];
      }
    }
    /* Potential and dipole formed by summing contributions from each edge */
    fs=0.0; 
    fd=0.0; 
    if(directions != NULL) {
      fsx = 0; fsy = 0;
      fdx= 0; fdy = 0; fdz = 0;
    }

    for(i=0; i <verts; i++) {
      if (i==(verts-1)) next=0;
      else next=i+1;

      /* v is the projection of the eval-i edge on the perpend to the side-i */
      /* Exploits the fact that corner points in panel coordinates.  */
      v=xmxv[i]*st[i] - ymyv[i]*ct[i];

      /* arg == zero if eval on next-i edge, but then v = 0. */
      arg=(r[i]+r[next]-edgeLength[i])/(r[i]+r[next]+edgeLength[i]);
      fln = -log(arg);
      if (arg>0.0) {
        fs = fs + v * fln;
        if (directions != NULL) {
          fac = (r[i]+r[next]-edgeLength[i]) * (r[i]+r[next]+edgeLength[i]);
          fac = v*(edgeLength[i]+ edgeLength[i])/fac;
          fsx = fsx + (fln*st[i] - fac*(xri[i] + xri[next]));
          fsy -= (fln*ct[i] + fac*(yri[i] + yri[next]));
          fdz -= (fac*( 1.0/r[i] + 1.0/r[next]));
        }
      }

      /* OK means eval not near a vertex normal, use Hess-Smith: */
      if (OK == 1) {
        s1=v*r[i];
        c1=znabs*(xmxv[i]*ct[i]+ymyv[i]*st[i]);
        s2=v*r[next];
        c2=znabs*(xmxv[next]*ct[i]+ymyv[next]*st[i]);
      } else { /* Near a vertex normal, use Newman  */
        s1=(fe[i]*st[i])-(xmxv[i]*ymyv[i]*ct[i]);
        c1=znabs*r[i]*ct[i];
        s2=(fe[next]*st[i])-(xmxv[next]*ymyv[next]*ct[i]);
        c2=znabs*r[next]*ct[i];
      }
  
      s12 = (s1*c2)-(s2*c1);
      c12 = (c1*c2)+(s1*s2);
      val = atan2(s12, c12);
      fd += val;
      if (directions != NULL) {
        u1 = xmxv[i]*ct[i] + ymyv[i]*st[i];
        u2 = xmxv[next]*ct[i]+ymyv[next]*st[i];
        if (OK == 0) {
          rr  = r[i]*r[i];
          fh1 = xmxv[i]*ymyv[i];
          fh2 = xmxv[next]*ymyv[next];
          fac = c1/((c1*c1+s1*s1)*rr );
          if(zn < 0.0) fac = -1.0 * fac;  /* Nick's 8/98 correction. */
          fdx = fdx + ((rr*v+fh1*u1)*fac);
          fdy = fdy - (fe[i]*u1*fac);
          rr  = r[next]*r[next];
          fac = c2/((c2*c2+s2*s2)*rr);
          if(zn < 0.0) fac = -1.0 * fac; /* Nick's 8/98 correction. */
          fdx = fdx - ((rr*v+fh2*u2)*fac);
          fdy = fdy + fe[next]*u2*fac;
        } else {
          fac = zn/(c1*c1+s1*s1);
          fdx = fdx + (u1*v*xri[i]+r[i]*ymyv[i])*fac;
          fdy = fdy + (u1*v*yri[i]-r[i]*xmxv[i])*fac;
          fac = zn/(c2*c2+s2*s2);
          fdx = fdx - ((u2*v*xri[next]+r[next]*ymyv[next])*fac);
          fdy = fdy - ((u2*v*yri[next]-r[next]*xmxv[next])*fac);
        }
      }
    }
    if (fd<0.0) fd += 2.0 * PI; 
    if (zn < 0) fd *= -1.0;
    fs -= zn*fd;

    if(directions != NULL ) {
      fsx = fsx - zn*fdx;
      fsy = fsy - zn*fdy;
      fes = nrm[0]*fsx + nrm[1]*fsy - nrm[2]*fd;
      fed = nrm[0]*fdx + nrm[1]*fdy + nrm[2]*fdz;
    }

    /* Save the computed potentials and derivatives. */
    if(fss != NULL) fss[evalindex] = fs;
    if(fds != NULL) fds[evalindex] = fd;
    if(fess != NULL) fess[evalindex] = fes;
    if(feds != NULL) feds[evalindex] = fed;
  }
} /* End of routine calcp */

double kernel(double x, double y, double z)
{
  double r;

  r = sqrt(x * x + y * y + z * z);
  if(r == 0.0) r = 0.0;
  else r = 1/r;
  return r;
}
