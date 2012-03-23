/***************************************************************************
 *   Copyright (C) 2005 by Xin Hu   *
 *   xinhu@mit.edu   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

//obselete program -- to be deleted 
 
#ifndef CALCPONEOVERR_H
#define CALCPONEOVERR_H

#include "panel.h"
#include "const_def.h"

#define MAX_REAL 1.0e+20;
#define PLANAR_TOL 1.0e-3
#define EQUIV_TOL 1.0e-9
#define TRUE true
#define FALSE false


class calcpOneOverR
{
public:
  calcpOneOverR(void);
  template<class T1, class T2>
  T2 calcp(const panel<T2>& panel_i, const panel<T1>& panel_j);
private:

  template<class T> struct charge
  {
    T corner[4][3];
    T length[4];
    T X[3],Y[3],Z[3];
    T moments[16];
    T x,y,z; //for centroid coord
    T max_diag, min_diag, area;
    int shape,order;
    T *multipoleM;
    T *multipoleD;
  };


  double max_my(double, double);
  double min_my(double, double);
  complex<double> max_my(complex<double>, complex<double>);
  complex<double> min_my(complex<double>, complex<double>);
  template<class T> void initcalcp(const panel<T>& p, charge<T>& c, const int order);
  template<class T>
  void Cross_Product(const T vector1[], const T vector2[], T result_vector[]);
  template<class T> T normalize(T vector[]);
  template<class T>  bool planarize(charge<T>& pq);
  template<class T> T Dot_Product(const T V1[],const T V2[]);
  template<class T> T DotP_Product(const T V1[], const T R, const T S, const T U);
  template<class T> void centroid(charge<T>& pq, T length);
  template<class T> void ComputeMoments(charge<T>& pp, const int order);
  template<class T>
  void calcspherical(charge<T>& pp, T** I, const int order);
  double** createBinom(int order);
  double** createPmn(int order);
  template<class T>
  T jacobid(int m,int mp, int n, T x,
            double **binom, int* fact);
  template <class T>
  void eulerAngles(const T X[], const T Y[], const T Z[],
                   T& pcosa, T& psina, T& pcosb, T& psignb, T& pcosg, T& psing);
  template<class T>
  void zrotate(const T cosa, const T sina, T** M, T** Mt, const int order);
  int* createFactorial(int size);
  bool less_than (complex<double> A, double B);
  bool less_than (complex<double> A, complex<double> B);
  bool less_than (double A, double B);
  bool less_equal_than (complex<double> A, double B);
  bool less_equal_than (double A, double B);
  double abs_my(double);
  complex<double> atan2_my(complex<double> val1,complex<double> val2 );
  double atan2_my(double val1,double val2 );
  complex<double> abs_my(complex<double>);
  int multerms(int order);
  int costerms(int);
  int sinterms(int);
  int sindex(int n, int m, int cterms);
  int index(int n, int m);
  int maxorderMom;
  int maxorderSphere;

};




template<class T1, class T2>
T2 calcpOneOverR::calcp(const panel<T2>& panel_i, const panel<T1>& panel_j)
{
  T2 x,y,z,xc,yc,zc,xn,yn,zn,zsq,xsq,ysq,rsq,diagsq,dtol;
  T2 fs, fd, *s, r2Inv, rInv, r3Inv, r5Inv, zr2Inv,ss1;
  T2 ss3, ss5, fdsum,s914, s813, s411,s512, s1215, r7Inv,r9Inv,ss7,ss9;
  T2 znabs, *corner, length, ct, st, v, arg;
  bool okay;
  int num2nd=0, num4th=0,numexact=0, i, next;
  T2 r[4], fe[4], xmxv[4], ymyv[4];
  T2 s1, c1, s2, c2, s12, c12, val;

  charge<T2> panel_src;
  initcalcp(panel_i,panel_src,0);

  charge<T1> panel_test;
  initcalcp(panel_j,panel_test,0);

  x=panel_test.x;
  y=panel_test.y;
  z=panel_test.z;

  xc = x - panel_src.x;
  yc = y - panel_src.y;
  zc = z - panel_src.z;

  xn = DotP_Product(panel_src.X, xc, yc, zc);
  yn = DotP_Product(panel_src.Y, xc, yc, zc);
  zn = DotP_Product(panel_src.Z, xc, yc, zc);

  zsq = zn * zn;
  xsq = xn * xn;
  ysq = yn * yn;
  rsq = zsq + xsq + ysq;
  diagsq = panel_src.max_diag * panel_src.max_diag;

  if(less_than((9.0 * diagsq),rsq))
  {
    fs = 0.0; fd = 0.0;
    s = panel_src.moments;
    /* First, second moments. */
    r2Inv = 1.0 / rsq;
    rInv = sqrt(r2Inv);
    r3Inv = r2Inv * rInv;
    r5Inv = r3Inv * r2Inv;
    zr2Inv = zn * r2Inv;
    ss1 = s[1] * rInv;
    ss3 = -(s[3] + s[10]) * r3Inv;
    ss5 = (xsq * s[10] + (xn * yn * s[7]) + ysq * s[3]) * r5Inv;
    fs = ss1 + (1./3.) * ss3 + ss5;
    fdsum = ss1 + ss3 + 5.0 * ss5;
    fd = zr2Inv * fdsum;


    if(less_than(rsq ,(36.0 * diagsq)))
    {
      /* Third and fourth moments added for diagsq/r2 between 40 and 150. */
      s914 = s[9] + s[14];
      s813 = s[8] + s[13];
      s411 = s[4] + s[11];
      s512 = s[5] + s[12];
      s1215 = s[12] + s[15];
      r7Inv = r5Inv * r2Inv;
      r9Inv = r7Inv * r2Inv;
      ss5 = (-xn * s813 - yn * s411 + 0.1 * (s512 + s1215)) * r5Inv;

      ss7 = ((5./3) *((xn * xsq * s[13] + yn * ysq * s[4])
                      + 3.0 * xn * yn * (xn * s[11]  +  yn * s[8]))
             - xsq * s1215 - ysq * s512 - xn * yn * s914) * r7Inv;

      ss9 = (7.0 * ((1./6) * (xsq * xsq * s[15] + ysq * ysq * s[5])
                    + xsq * ysq * s[12])
             + (7./3) * xn * yn * (xsq * s[14] + ysq * s[9])) * r9Inv;

      fs += ss5 + ss7 + ss9;
      fdsum = 5.0 * ss5 + 7.0 * ss7 + 9.0 * ss9;
      fd += zr2Inv * fdsum;
      /*ftk diagnostics
      		printf("fourth order source = %.8g dipole = %.8g\n", 
      		fs/panel->area, fd/panel->area); 
      */
      num4th++;
    }
    else num2nd++;
  }
  else
  {
    dtol = EQUIV_TOL * panel_src.min_diag;
    znabs = abs_my(zn);

    /* Always move the evaluation point a little bit off the panel. */
    if(less_than(znabs , dtol))
    {
      zn = 0.5 * dtol;  /* Half of dtol insures detection for zero dipole. */
      znabs = 0.5 * dtol;
    }

    /* Once per corner computations. */
    for(okay = TRUE, i=0; i < panel_src.shape; i++)
    {
      corner = panel_src.corner[i];
      xmxv[i] = xc = xn - corner[XX];
      ymyv[i] = yc = yn - corner[YY];
      zc = zn - corner[ZZ];
      fe[i] = xc * xc + zc * zc;
      r[i] = sqrt(yc * yc + fe[i]);
      if(less_than(r[i] ,(1.005 * znabs)))
      {  /* If r almost z, on vertex normal. */
        okay = FALSE;
      }
    }
    /* Once per edge computations. */
    fs = 0.0; fd = 0.0;
    for(i=0; i < panel_src.shape; i++)
    {
      if(i == (panel_src.shape - 1)) next = 0;
      else next = i + 1;

      /* Now calculate the edge contributions to a panel. */
      length = panel_src.length[i];
      ct = (panel_src.corner[next][XX] - panel_src.corner[i][XX]) / length;
      st = (panel_src.corner[next][YY] - panel_src.corner[i][YY]) / length;

      /* v is projection of eval-i edge onto perpend to next-i edge. */
      /* Exploits the fact that corner points in panel coordinates. */
      v = xmxv[i] * st - ymyv[i] * ct;

      /* arg == zero if eval on next-i edge, but then v = 0. */
      arg = (r[i] + r[next] - length)/(r[i] + r[next] + length);
	  if(less_than(0.0,arg)) fs -= v * log(arg);

      /* Okay means eval not near a vertex normal, Use Hess-Smith. */
      if(okay)
      {
        s1 = v * r[i];
        c1 = znabs * (xmxv[i] * ct + ymyv[i] * st);
        s2 = v * r[next];
        c2 = znabs * (xmxv[next] * ct + ymyv[next] * st);
      }
      /* Near a vertex normal, use Newman. */
      else
      {
        s1 = (fe[i] * st) - (xmxv[i] * ymyv[i] * ct);
        c1 = znabs * r[i] * ct;
        s2 = (fe[next] * st) - (xmxv[next] * ymyv[next] * ct);
        c2 = znabs * r[next] * ct;
      }

      s12 = (s1 * c2) - (s2 * c1);
      c12 = (c1 * c2) + (s1 * s2);
      val = atan2_my(s12, c12);
      fd += val;
    }

	/* Adjust the computed values. */

	if(less_than(fd , 0.0)) fd += 2*PI;
	if(less_than(zn , 0.0)) fd *= -1.0;

	/* If in the same plane as panel, fd = 0. */
	if(less_than(znabs , dtol)) {
		if(less_than(rsq , (dtol * dtol))) fd = -2*PI;
		else fd = 0.0;
	}
	fs -= zn * fd;
	numexact++;
  }
  
  /* Return values of the source and dipole, normalized by area. */
  fd *= -1.0;  /* Adjusts for opposite panel normal. */
  /* ftk diagnostics */

  fs = fs/panel_src.area;
  fd = fd/panel_src.area;
  
  //if(pfd != NULL) *pfd = fd;

//   if(less_than(fs , 0.0) ) {
// 	  fprintf(stderr, "FLE-clacp: fs (exact source) is less than zero \n");
// 	  fprintf(stderr, "FLE-calcp: Okay =Evaluation Point = \n");
// 	  fprintf(stderr, "FLE-calcp: Evaluation Point in local coords = \n");
// 	  fprintf(stderr, "FLE-calcp: Panel Description Follows\n");
// 	  exit(0);
//   }

  return (fs);
}


template<class T> void calcpOneOverR::initcalcp(const panel<T>& p, charge<T> & chg, const int order)
{
  int i,j,next;
  T sum, sum2, delta, length, maxlength, minlength;
  T length20, length31,  vtemp[3];


  maxlength = 0.0;
  minlength = MAX_REAL;
  for (i=0; i<4; i++)
  {
    chg.corner[i][0] = (p.get_apex(i+1))->x();
    chg.corner[i][1] = (p.get_apex(i+1))->y();
    chg.corner[i][2] = (p.get_apex(i+1))->z();
  }
  chg.shape=p.shape();

  //   for (i=0; i<4; i++)
  // 	  cout<<chg.corner[i][0] <<" "<< chg.corner[i][1] <<" "<<chg.corner[i][2]<<endl;

  for(i=0; i < p.shape(); i++)
  {
    if(i == (p.shape() -1)) next = 0;
    else next = i + 1;
    for(sum= 0., j = 0; j < 3; j++)
    {
      delta = chg.corner[next][j] - chg.corner[i][j];
      sum += delta * delta;
    }
    chg.length[i] = length = sqrt(sum);
    maxlength = max_my(maxlength, length);
    minlength = min_my(minlength, length);
  }
  //   for (i=0; i<4; i++) cout<<chg.length[i]<<endl;

  for(sum= 0.0, sum2 = 0.0, i = 0; i < 3; i++)
  {
    chg.X[i] = delta = chg.corner[2][i] - chg.corner[0][i];
    sum += delta * delta;
    if(p.shape() == 3) chg.Y[i] = chg.corner[1][i] - chg.corner[0][i];
    else
    {
      chg.Y[i] = delta = chg.corner[1][i] - chg.corner[3][i];
      sum2 += delta * delta;
    }
  }
  length20 = sqrt(sum);
  length31 = sqrt(sum2);

  //for (i=0; i<3; i++) cout<<chg.X[i]<<" "<<chg.Y[i]<<" "<<chg.Z[i]<<endl;


  /* Check on lengths for quad. */
  if(p.shape() == 3)
  {
    chg.max_diag = maxlength;
    chg.min_diag = minlength;
  }
  else
  {
    length = max_my(length20, length31);
    chg.max_diag = length;
    chg.min_diag = min_my(length20, length31);
  }
  //cout<<chg.max_diag<<" "<<chg.min_diag<<endl;

  /* Z-axis is normal to two diags. */
  Cross_Product(chg.X, chg.Y, chg.Z);
  chg.area = 0.5 * normalize(chg.Z);
  normalize(chg.X);
  //cout<<chg.area<<endl;

  /* Real Y-axis is normal to X and Z. */
  Cross_Product(chg.Z, chg.X, chg.Y);
  //   cout<<chg.X[0]<<chg.X[1]<<chg.X[2]<<endl;
  //   cout<<chg.Y[0]<<chg.Y[1]<<chg.Y[2]<<endl;
  //   cout<<chg.Z[0]<<chg.Z[1]<<chg.Z[2]<<endl;

  /* Project the corner points into the plane defined by edge midpoints. */
  if(planarize(chg) == FALSE)
    fprintf(stdout, "FLW-initcalcp: Panel skewed beyond the PLANAR_TOL tolerance.\n");

  /* Calculate the centroid. */
  centroid(chg, length20);
  //   cout<<chg.x<<" "<<chg.y<<" "<<chg.z<<endl;

  /* Put corners in the newly defined coord system. */
  for(i=0; i < chg.shape; i++)
  {
    chg.corner[i][XX] -= chg.x;
    chg.corner[i][YY] -= chg.y;
    chg.corner[i][ZZ] -= chg.z;
  }
  for(i=0; i < chg.shape; i++)
  {
    vtemp[XX] = Dot_Product(chg.corner[i], chg.X);
    vtemp[YY] = Dot_Product(chg.corner[i], chg.Y);
    vtemp[ZZ] = Dot_Product(chg.corner[i], chg.Z);
    chg.corner[i][0]=vtemp[XX];
    chg.corner[i][1]=vtemp[YY];
    chg.corner[i][2]=vtemp[ZZ];
    if(abs(chg.corner[i][ZZ]) > (EQUIV_TOL * abs(chg.min_diag)))
    {
      fprintf(stderr, "FLE-initcalcp: renormalized z=n" );
      exit(0);
    }
    chg.corner[i][ZZ] = 0.0;
  }
  //   for(i=0; i < chg.shape; i++)
  // 	  cout<<chg.corner[i][XX]<<" "<<chg.corner[i][YY]<<" "<<chg.corner[i][ZZ]<<endl;

  ComputeMoments(chg, order);

}

template<class T>
void calcpOneOverR::ComputeMoments(charge<T>& pp, const int order)
{
  int momOrder,i,j,nside,M1,M2,M,N, MN1,N1,MN2;
  T *XP[4], *YP[4], **I,*xp, *yp,*xpn, *ypn;
  T dx,dy,dydx,dxdy, SI;
  double CS[16] = { 0.0, 1.0, 1.0, 1.5, 1.5, 3.75, 1.0, 3.0,
                    1.5, 7.5, 1.5, 1.5, 3.75, 1.5, 7.5, 3.75 };

  momOrder = MAX(order,4);
  if(momOrder > maxorderMom)
  {
    for(i = 0; i < 4; i++)
    {
      XP[i]=(T*) malloc( (momOrder+3)*sizeof(T));
      YP[i]=(T*) malloc( (momOrder+3)*sizeof(T));
    }

    I=(T**)malloc((momOrder+1)*sizeof(T*));
    for(i=0; i <= momOrder; i++) I[i]=(T*) malloc((momOrder+1)*sizeof(T));
  }

  /* First zero out the Moments matrix. */
  for(i = 0; i <= momOrder; i++)
  {
    for(j = 0; j <= momOrder; j++) I[i][j] = (T) 0.0;
  }

  /* Compute powers of x and y at corner pts. */
  for(i = 0; i < pp.shape; i++)
  {
    xp = XP[i];
    yp = YP[i];
    xp[1] = pp.corner[i][XX];
    yp[1] = pp.corner[i][YY];
    for(j = 2; j <= momOrder+2; j++)
    {
      xp[j] = xp[j-1] * xp[1];
      yp[j] = yp[j-1] * yp[1];
    }
  }

  /* First moment, easy, just the panel area. */
  I[0][0] = pp.area;

  /* By using centroid, (1,0) and (0,1) are zero, so begin with (2,0). */
  for(nside = 0; nside < pp.shape; nside++)
  {
    xp = XP[nside];
    yp = YP[nside];
    if(nside == (pp.shape - 1))
    {
      xpn = XP[0];
      ypn = YP[0];
    }
    else
    {
      xpn = XP[nside + 1];
      ypn = YP[nside + 1];
    }


    dx = xpn[1] - xp[1];
    dy = ypn[1] - yp[1];

    if(abs(dx) >= abs(dy))
    {
      dydx = dy/dx;
      for(M = 2; M <= momOrder; M++)
      {
        M1 = M + 1;
        M2 = M + 2;

        SI = ((xpn[M1] * ypn[1]) - (xp[M1] * yp[1])) / (double) M1
             + dydx * (xp[M2] - xpn[M2]) / (double)(M1 * M2);
        I[M][0] += SI;

        for(N = 1; N <= M; N++)
        {
          N1 = N + 1;
          MN1 = M - N + 1;
          SI = (xpn[MN1] * ypn[N1] - xp[MN1] * yp[N1]) / (double)(MN1 * N1)
               - (dydx * (double)N * SI) / (double)MN1;
          I[M-N][N] += SI;
        }
      }
    }
    else
    {
      dxdy = dx/dy;
      for(M = 2; M <= momOrder; M++)
      {
        M1 = M + 1;
        M2 = M + 2;
        SI = (dxdy / (double)(M1 * M2)) * (ypn[M2] - yp[M2]);
        I[0][M] += SI;
        for(N = 1; N <= M; N++)
        {
          MN1 = M - N + 1;
          MN2 = MN1 + 1;
          SI = dxdy * ((xpn[N] * ypn[MN2] - xp[N] * yp[MN2]) /(double) (MN1 * MN2)
                       - ((double)N * SI /(double) MN1));
          I[N][M-N] += SI;
        }
      }
    }
  }
  /* Now Create the S vector for calcp. */
  for(i = 0, M = 0; M <= 4; M++)
  {
    for(N = 0; N <= (4 - M); N++)
    {
      i++;
      pp.moments[i] = I[M][N] * CS[i];
    }
  }
  //calcspherical(pp, I, order);

}

template<class T>
void calcpOneOverR::calcspherical(charge<T>& pp, T** I, const int order)
{
  T **Nmn, **Nrmn, **Ntmn, **Nrtmn;
  T **Mmn, **Mrmn, **Mtmn, **Mrtmn;
  T *multiM, *multiD;
  static double **Pmn, **Binom;
  double sqfac,sign;
  T sumc, msumc, msums, sums,cosa, sina, cosb, signb, cosg, sing;
  T dplus,dminus,rcoeff,icoeff, coeff;
  int*Fact;
  int i,j, n,m,r,mp,numterms, rterms, halfn,halfm,flrm, ceilm;


  //if(order > maxorderSphere)
  //{
  /* Allocate a temporary Multipole moments matrix, Mmn. */
  Mmn=(T**)malloc((order+1)*sizeof(T*));
  for(i = 0; i <= order; i++) Mmn[i]=(T*)malloc((order+1)*sizeof(T));
  Mtmn=(T**)malloc((order+1)*sizeof(T*));
  for(i = 0; i <= order; i++) Mtmn[i]=(T*)malloc((order+1)*sizeof(T));

  /* Allocate a rotated temporary Multipole matrix. */
  Mrmn=(T**)malloc((order+1)*sizeof(T*));
  for(i = 0; i <= order; i++) Mrmn[i]=(T*)malloc((order+1)*sizeof(T));
  Mrtmn=(T**)malloc((order+1)*sizeof(T*));
  for(i = 0; i <= order; i++) Mrtmn[i]=(T*)malloc((order+1)*sizeof(T));

  /* Allocate temporary spherical harmonic coefficient matrices for dipole
  terms, N*mn. */
  Nmn=(T**)malloc((order+1)*sizeof(T*));
  for(i = 0; i <= order; i++) Nmn[i]=(T*)malloc((order+1)*sizeof(T));
  Ntmn=(T**)malloc((order+1)*sizeof(T*));
  for(i = 0; i <= order; i++) Ntmn[i]=(T*)malloc((order+1)*sizeof(T));

  Nrmn=(T**)malloc((order+1)*sizeof(T*));
  for(i = 0; i <= order; i++) Nrmn[i]=(T*)malloc((order+1)*sizeof(T));
  Nrtmn=(T**)malloc((order+1)*sizeof(T*));
  for(i = 0; i <= order; i++) Nrtmn[i]=(T*)malloc((order+1)*sizeof(T));

  Binom = createBinom((order+1) * 2);
  Pmn = createPmn((order+1) * 2);
  Fact = createFactorial((order+1) * 2 + 1);

  maxorderSphere = order;
  //}

  pp.order = order;
  numterms = multerms(order);
  rterms = costerms(order);

  /* ftk allocate both source and dipole space */
  multiM=(T*)malloc(numterms*sizeof(T));
  multiD=(T*)malloc(numterms*sizeof(T));
  pp.multipoleM = multiM;
  pp.multipoleD = multiD;

  /* Just a thought. */
  for(n = 0; n <= order; n++)
  {
    for(m = 0; m <= n; m++)
    {
      Mrtmn[m][n] = Mtmn[m][n] = (T)0.0;
      Nrtmn[m][n] =  Ntmn[m][n] = (T)0.0;
      Mrmn[m][n] = Mmn[m][n] = (T)0.0;
      Nrmn[m][n] =  Nmn[m][n] = (T)0.0;
    }
  }

  /* Create the M^0_n multipole coefficients first. */
  Nmn[0][0] = 0.0;
  for(n = 0; n <= order; n++)
  {
    halfn = n / 2;
    if(2 * halfn == n)
    { /* n even. */
      for(sumc = 0.0, i=0; i <= halfn; i++)
      {
        sumc += Binom[halfn][halfn - i] * I[n - 2 * i][2 * i];
      }
      Mmn[0][n] = sumc * Pmn[0][n];
      /* ftk The dipole assignment is added. */
      if (n < order) Nmn[0][n+1] = (double)(n+1) * sumc * Pmn[0][n];
    }
    else
    {
      Mmn[0][n] = 0.0;
      if (n < order) Nmn[0][n+1] = 0.0;
    }
  }

  /* Create the rest. */
  for(n = 1; n <= order; n++)
  {
    for(m = 1; m <= n; m++)
    {
      halfn = (n - m) / 2;
      if((2 * halfn) == (n - m))
      { /* n-m even. */
        halfm = m / 2;
        if((2 * halfm) != m) { flrm = halfm; ceilm = halfm + 1;	}
        else { flrm = halfm; ceilm = halfm; }

        for(sumc = 0.0, sums = 0.0, i=0; i <= halfn; i++)
        {
          sign = -1.0;
          for(msumc = 0.0, r=0; r <= flrm; r++)
          {
            sign *= -1.0;
            msumc += sign * Binom[m][m - 2*r] * I[n - 2*(r+i)][2*(r+i)];
          }

          sign = 1.0;
          for(msums = 0.0, r=0; r <= (ceilm-1); r++)
          {
            sign *= -1.0;
            msums +=
              sign * Binom[m][m - (2*r+1)] * I[n - (2*(r+i)+1)][2*(r+i) + 1];
          }
          sumc = sumc+Binom[halfn][halfn - i] * msumc;
          sums = sums+Binom[halfn][halfn - i] * msums;
        }
        /* ftk The source assignment is unchanged. */
        Mmn[m][n]  = sumc *  2.0 * Pmn[m][n];
        Mtmn[m][n] = sums * -2.0 * Pmn[m][n];
        /* ftk The dipole assignment is added. */
        if (n < order)
        {
          sqfac = sqrt( (double) (((n+1) * (n+1)) + (m*m)));
          Nmn[m][n+1]  = sqfac * sumc *  2.0 * Pmn[m][n];
          Ntmn[m][n+1] = sqfac * sums * -2.0 * Pmn[m][n];
        }
      }
    }
  }
  eulerAngles(pp.X, pp.Y, pp.Z, cosa, sina, cosb, signb, cosg, sing);

  /* First rotation about cosa, sina (we are doing both source and dipole). */
  zrotate(cosa, sina, Mmn, Mtmn, order);
  zrotate(cosa, sina, Nmn, Ntmn, order);
  /* Rotate theta by 2pi if z prime is in the -xpp half space. */
  /*ftk This operates on all coefficients, so N's are just included. */

  if(less_than(signb, 0.))
  {
    cosb *= -1.0;
    for(n = 0; n <= order; n++)
    {
      if((n % 2) == 0) { rcoeff = 1.0; icoeff = -1.0; }
      else { rcoeff = -1.0; icoeff = 1.0; }

      for(m = 0; m <= n; m++)
      {
        Mmn[m][n] *= rcoeff;
        Mtmn[m][n] *= icoeff;
        Nmn[m][n] *= rcoeff;
        Ntmn[m][n] *= icoeff;
      }
    }
  }


  for(n = 0; n <= order; n++)
  {
    for(m = 0; m <= n; m++)
    {
      Mrmn[m][n] = 0.0;
      Mrtmn[m][n] = 0.0;
      Nrmn[m][n] = 0.0;
      Nrtmn[m][n] = 0.0;
      for(mp = 0; mp <= n; mp++)
      {
        coeff = ((double) Fact[n - mp]) / ((double) Fact[n - m]);
        dplus = jacobid(mp, m, n, cosb, Binom, Fact);
        dminus = jacobid(-mp, m, n, cosb, Binom, Fact);
        if(m == 0)
        { /* This test has to do with conjugation properties. */
          icoeff = 0.0;
          if((mp % 2) == 0) rcoeff = 0.5 * coeff * (dplus + dminus);
          else rcoeff = 0.5 * coeff * (dplus - dminus);
          Mrmn[m][n] += rcoeff * Mmn[mp][n];
          Nrmn[m][n] += rcoeff * Nmn[mp][n];
        }
        else
        {
          if((mp % 2) == 0)
          {
            rcoeff = coeff * (dplus + dminus);
            icoeff = coeff * (dplus - dminus);
          }
          else
          {
            rcoeff = coeff * (dplus - dminus);
            icoeff = coeff * (dplus + dminus);
          }
          Mrmn[m][n] += rcoeff * Mmn[mp][n];
          Mrtmn[m][n] += icoeff * Mtmn[mp][n];
          Nrmn[m][n] += rcoeff * Nmn[mp][n];
          Nrtmn[m][n] += icoeff * Ntmn[mp][n];
        }
      }
    }
  }

  /* Final rotation about z by cosg, sing. (we are doing both source
  and dipole). */
  zrotate(cosg, -sing, Mrmn, Mrtmn, order);
  zrotate(cosg, -sing, Nrmn, Nrtmn, order);

  /* Copy the Multipole matrix into the multi vector for this panel. */
  int temp, temp3, temp2;
  for(n = 0; n <= order; n++)
  {
    temp=index(n,0);
    multiM[temp] = Mrmn[0][n];
    multiD[temp] = Nrmn[0][n];
    for(m = 1; m <= n; m++)
    {
      temp2=index(n,m);
      temp3=sindex(n,m,rterms);
      multiM[temp2] = Mrmn[m][n];
      multiM[temp3] = Mrtmn[m][n];
      multiD[temp2] = Nrmn[m][n];
      multiD[temp3] = Nrtmn[m][n];
    }
  }

  /* Finally, scale to make M00 1, and flip dipole sign to adjust for panel
  normal. */
  for(n = 1; n < numterms; n++)
  {
    multiM[n] /= multiM[0];
    multiD[n] /= -multiM[0];
  }
  multiD[0] /= multiM[0]; /* This term is zero anyway.*/
  multiM[0] = 1.0;

}

template<class T>
T calcpOneOverR::jacobid(int m,int mp,
                         int n, T x,
                         double **binom, int* fact)
{
  double coeff = 1.0;
  T ax,val, part1, part2;
  int a, b, o, i;
  int temp, half, twosupm;

  /* First change m and mp to make sure m-mp > 0 and m+mp > 0. */
  if((m + mp) < 0)
  {
    m *= -1;
    mp *= -1;
    if(((mp - m) % 2) != 0) coeff *= -1.0;
  }
  if((m - mp) < 0)
  {  /* Must swap m and mp. */
    temp = m;
    m = mp;
    mp = temp;
    coeff = ((double) (fact[n-m] * fact[n+m]))
            / ((double) (fact[n-mp] * fact[n+mp]));
    if(((m - mp) % 2) != 0) coeff *= -1.0;
  }

  /* Compute 2^m. */
  for(twosupm = 1, i=1; i <= m; i++) twosupm *= 2;

  if((m - mp) == 0) part1 = 1.0;
  else part1 = pow(1.0 - x, 0.5 * (m - mp));
  if((m + mp) == 0) part2 = 1.0;
  else part2 = pow(1.0 + x, 0.5 * (m + mp));

  val = coeff * part1 * part2 / ((double) twosupm);

  /* The jacobi polynomial. */
  a = m - mp;
  b = m + mp;
  o = n - m;
  for(ax = 1.0, i=o; i > 0; i--)
  {
    ax = 1. - ax * (1.0 - x) * ((double) (o - i + 1) * (a + b + i + o))
         / ((double) (2 * i * (a + i)));
  }
  val *= binom[o + a][o] * ax;

  /*
  	printf("x=%g m=%d mp=%d, n=%d a=%d b=%d val=%g\n", x, m, mp, n, a, b, val);
  */
  return(val);
}

template<class T>
void calcpOneOverR:: zrotate(const T cosa, const T sina, T** M, T** Mt, const int order)
{
  T temp, cosm, sinm;
  int m, n;

  cosm = cosa;
  sinm = sina;
  for(m = 1; m <= order; m++)
  {
    for(n = m; n <= order; n++)
    {
      temp = M[m][n];
      M[m][n] = cosm * temp - sinm *  Mt[m][n];
      Mt[m][n] = cosm * Mt[m][n] + sinm * temp;
    }
    /* Compute  exp((m+1)a) = exp(a) * exp(ma). */
    temp = cosm;
    cosm = temp * cosa - sinm * sina;
    sinm = temp * sina + sinm * cosa;
  }
}



/* Calculate the Euler Angles between X, Y, Z and global coordinate system. */
template <class T>
void calcpOneOverR::eulerAngles(const T X[], const T Y[], const T Z[],
                                T& pcosa, T& psina, T& pcosb, T& psignb, T& pcosg, T& psing)
{
  T ypp[3], xpp[3], norm;

  /* First calculate ypp, the intersection of X-Y plane with global XY plane. */
  ypp[ZZ] = 0.0;
  if(Y[ZZ] != 0.0)
  {
    ypp[XX] = Y[ZZ] * X[XX] - X[ZZ] * Y[XX];
    ypp[YY] = Y[ZZ] * X[YY] - X[ZZ] * Y[YY];
  }
  else
  {
    ypp[XX] = Y[XX];
    ypp[YY] = Y[YY];
  }

  /* Force ypp to be in the positive global Y plane. */
  if(less_than(ypp[YY],0))
  {
    ypp[YY] *= -1.0;
    ypp[XX] *= -1.0;
  }
  //
  /* Normalize the ypp vector. Use fact that ypp[ZI] = 0. */
  norm = sqrt(ypp[XX] * ypp[XX] + ypp[YY] * ypp[YY]);
  ypp[XX] /= norm;
  ypp[YY] /= norm;

  Cross_Product(ypp, Z, xpp);


  /* Now calculate angles. */
  pcosa = ypp[XX] * Y[XX] + ypp[YY] * Y[YY];
  psina = ypp[XX] * X[XX] + ypp[YY] * X[YY];
  pcosb = Z[ZZ];
  if(less_equal_than(xpp[ZZ], 0.)) psignb = 1.0;
  else psignb = -1.0;
  pcosg = ypp[YY];
  psing = ypp[XX];
}



/* Calculates the cross product between two vectors. */
template<class T>
void calcpOneOverR::Cross_Product(const T vector1[], const T vector2[], T result_vector[])
{
  result_vector[XX] = vector1[YY]*vector2[ZZ] - vector1[ZZ]*vector2[YY];
  result_vector[YY] = vector1[ZZ]*vector2[XX] - vector1[XX]*vector2[ZZ];
  result_vector[ZZ] = vector1[XX]*vector2[YY] - vector1[YY]*vector2[XX];
}

template<class T>
T calcpOneOverR::normalize(T vector[])
{
  T length;
  int i;

  length = sqrt( vector[0]*vector[0]
                 + vector[1]*vector[1]
                 + vector[2]*vector[2]);

  for (i=0; i<3; i++) vector[i] = vector[i] / length;
  return length;
}

template<class T>
bool calcpOneOverR::planarize(charge<T>& pq)
{
  T origin[3], corner[3], delta[4][3], px, py, dx, dy, dz;
  int i, j, numcorners = pq.shape;
  T tolsq = PLANAR_TOL * pq.min_diag;

  tolsq *= tolsq;

  /* Triangular panels are just fine already. */
  if(numcorners != 4) return(TRUE);

  /* Pick edge midpoint as origin. */
  for(i=0; i < 3; i++) origin[i] = 0.5 * (pq.corner[1][i] + pq.corner[0][i]);
  //   for(i=0; i < 3; i++) cout<<origin[i]<<endl;

  for(i=0; i < numcorners; i++)
  {
    for(j=0; j < 3; j++) corner[j] = pq.corner[i][j] - origin[j];
    px = Dot_Product(corner, pq.X);
    py = Dot_Product(corner, pq.Y);

    //cout<<px<<" "<<py<<endl;

    dx = px * pq.X[XX] + py * pq.Y[XX] + origin[XX] - pq.corner[i][XX];
    dy = px * pq.X[YY] + py * pq.Y[YY] + origin[YY] - pq.corner[i][YY];
    dz = px * pq.X[ZZ] + py * pq.Y[ZZ] + origin[ZZ] - pq.corner[i][ZZ];

    // 	cout<<dx<<" "<<dy<<" "<<dz<<endl;

    delta[i][XX] = dx;
    delta[i][YY] = dy;
    delta[i][ZZ] = dz;
  }

  for(i=0; i < numcorners; i++)
    for(j=0; j < 3; j++)
      pq.corner[i][j] += delta[i][j];

  //   for(i=0; i < numcorners; i++)
  // 	  for(j=0; j < 3; j++)
  // 		  cout<<pq.corner[i][j];

  /* If moved beyond the tolerance, set flag. */
  if(abs(dx * dx + dy * dy + dz * dz) > abs(tolsq))
    return(FALSE);
  else
    return(TRUE);

}

template<class T>
T calcpOneOverR::Dot_Product(const T V1[],const T V2[])
{
  return V1[XX]*V2[XX]+V1[YY]*V2[YY]+V1[ZZ]*V2[ZZ];
}

template<class T>
T calcpOneOverR::DotP_Product(const T V1[], const T R, const T S, const T U)
{
  return V1[XX]*R+V1[YY]*S+V1[ZZ]*U;
}

template<class T>
void calcpOneOverR:: centroid(charge<T>& pp, T x2)
{
  T vertex1[3], vertex3[3];
  T sum, dl, x1, y1, x3, y3, xc, yc;
  int i;

  /* Use vertex 0 as the origin. */
  for(i=0; i< 3; i++)
  {
    vertex1[i] = pp.corner[1][i] - pp.corner[0][i];
    if(pp.shape == 4) vertex3[i] = pp.corner[3][i] - pp.corner[0][i];
    else vertex3[i] = pp.corner[2][i] - pp.corner[0][i];
  }

  /* Project into the panel axes. */
  y1 = Dot_Product(vertex1, pp.Y);
  y3 = Dot_Product(vertex3, pp.Y);
  x1 = Dot_Product(vertex1, pp.X);
  x3 = Dot_Product(vertex3, pp.X);

  yc = 1./3. * (y1 + y3);
  xc = 1./3. * (x2 + ((x1 * y1 - x3 * y3)/(y1 - y3)));

  pp.x = pp.corner[0][XX] + xc * pp.X[XX] + yc * pp.Y[XX];
  pp.y = pp.corner[0][YY] + xc * pp.X[YY] + yc * pp.Y[YY];
  pp.z = pp.corner[0][ZZ] + xc * pp.X[ZZ] + yc * pp.Y[ZZ];
}
#endif
