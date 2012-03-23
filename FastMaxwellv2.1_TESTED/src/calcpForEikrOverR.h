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
// compute slp=\int_s  d/rho^2 [1/(jk)^2(e^{-jkR}/R -1/R - e^{-jk|h|/|h|+1/|h|)] dr' 
// this is needed in the double volume computation to get A. 


// compute dlp=\int_s e^{-jkR}/R dr'  or   dlp=\int_s 1/R dr'.  
//This is needed for panel integration to compute potential P 
//double volume integration to get dAdx and dAdy calculation
                                                                  
 
 
#ifndef CALCPFOREIKROVERR_H_
#define CALCPFOREIKROVERR_H_
#include "panel.h"
#include "point3D.h"
#include <complex>
//#include <cmath.h>
//#include "complex_lib/cmath.h"
#include <iostream>
#include <cfloat>

template <class T1, class T2>
class calcpForEikrOverR
{
  static const size_t DEFAULT_NUM_GQ_PTS=24;
  enum CalcpEvalPntDirection {
    AGAINST_NORMAL_DIRECTION = +1,
    FOLLOW_NORMAL_DIRECTION = -1
  };

public:
  calcpForEikrOverR(void) {};
  calcpForEikrOverR( const complex<double> K_num, const bool needSLP, const bool needDLP,
                     const size_t n=DEFAULT_NUM_GQ_PTS);
  void init(const complex<double> K_num, const bool needSLP, const bool needDLP,
            const size_t n=DEFAULT_NUM_GQ_PTS);
  void operator() (panel<T1> src_p, point3D<T2> eval,
                   complex<double>& slp, complex<double>& dlp, complex<double>& zn);

private:
  void gaussQuadrature ( const size_t n, std::vector<double>& w, std::vector<double>& x);
  void setupPanelLocalCoordSystem(panel<T1>&, const point3D<T2> &,point3D<T1 > &);
  void getEvalPointDirection ( const point3D<T1 >& evalPnt);
  void setupCalcpTriangleInfo ( size_t i, const panel<T1>& src_panel);
  double getContributionSign (void);
  void compIntegralOnTriangle ( const point3D<T1 >& evalPnt,
                                std::complex<double>& slp,
                                complex<double> & grad);
  void compOuterIntegral ( const point3D<T1 >& evalPnt,
                           T1 distance,
                           T1 yStart,
                           std::complex<double>& slp,
                           complex<double> & dlp);
  complex<double> compSingleLayerInnerIntegral(const complex<double>& IK,
      const T1& Ra, const T1& R);
  complex<double> compDoubleLayerInnerIntegral(const complex<double>& IK,
      const T1& Ra, const T1& R);

  //for panel quadrature of eikr/r kernel
  complex<double> quad_panel(panel<T1>& src_p, const point3D<T1 >& eval, const int numQuadPts);




  bool greater_than(complex<double>, double);
  bool greater_than(double, double);
  bool greater_than( double A , complex<double> B);
  bool less_than(complex<double>, double);
  bool less_than(double, double);
  bool less_than(complex<double> A , complex<double> B);
  complex<double> abs_my(complex<double> A );
  double abs_my(double A );
  complex<double> atan_my(complex<double> A);
  double atan_my(double A);

  struct TriangleInfo_
  {
    point3D<T1 > PH_unit;
    point3D<T1 > AH;
    point3D<T1 > HB;
    point3D<T1 > X;
    point3D<T1 > Y;
    T1 PH_len, AB_len, AH_len, area;
  };
  TriangleInfo_ triangleInfo_;

  complex<double> K_;
  size_t NumGQpoint_;
  std::vector<double> GQweight_;
  std::vector<double> GQpoint_;

  bool needSingleLayerPotential_;
  bool needDoubleLayerPotential_;
  bool needNothing_;
  bool projectPointOnEdge_;
  CalcpEvalPntDirection evalPointDirection_;
};


template <class T1, class T2>
calcpForEikrOverR<T1, T2>::calcpForEikrOverR (
  const complex<double> K_num,
  const bool needSLP,
  const bool needDLP,
  const size_t n): K_(K_num), needSingleLayerPotential_(needSLP),
    needDoubleLayerPotential_(needDLP), projectPointOnEdge_(false)
{
  needNothing_ = (!needSLP) && (!needDLP);
  if (! needNothing_)
  {
    evalPointDirection_ = AGAINST_NORMAL_DIRECTION;

    if ((GQpoint_.size() != n) && (abs(K_)!=0))
    {
      // To avoid one of Gauss Quad point to be zero, please
      // use the even number of points;
      NumGQpoint_ = n;
      if(NumGQpoint_ % 2 != 0)
      {
        NumGQpoint_ += 1;
      }
      GQpoint_.resize(NumGQpoint_);
      GQweight_.resize(NumGQpoint_);
      gaussQuadrature(NumGQpoint_, GQweight_, GQpoint_);
    }
  }
}

template <class T1, class T2>
void calcpForEikrOverR<T1, T2>::init (
  const complex<double> K_num,
  const bool needSLP,
  const bool needDLP, const size_t n)
{
  K_=K_num;
  needSingleLayerPotential_=needSLP;
  needDoubleLayerPotential_=needDLP;
  projectPointOnEdge_=false;

  needNothing_ = (!needSLP) && (!needDLP);
  if (! needNothing_)
  {
    evalPointDirection_ = AGAINST_NORMAL_DIRECTION;
    if ((GQpoint_.size() != n) && (abs(K_)!=0))
    {
      // To avoid one of Gauss Quad point to be zero, please
      // use the even number of points;
      NumGQpoint_ = n;
      if(NumGQpoint_ % 2 != 0)
      {
        NumGQpoint_ += 1;
      }
      GQpoint_.resize(NumGQpoint_);
      GQweight_.resize(NumGQpoint_);
      gaussQuadrature(NumGQpoint_, GQweight_, GQpoint_);
    }
  }
}

template<class T1, class T2>
void calcpForEikrOverR<T1, T2>::
operator()(panel<T1> src_panel, point3D<T2> evalPt,
           complex<double>& slp, complex<double>& dlp, complex<double>& zn)
{
  double contributionSign;
  T1 totalArea;
  complex<double> slpLocal, dlpLocal;


  if (needNothing_) return;


 /// point3D<T1 > evalPt_new(evalPt.x(),evalPt.y(),evalPt.z() );
  point3D<T1 > evalPt_new(evalPt.x(),evalPt.y(),evalPt.z() );
  setupPanelLocalCoordSystem(src_panel, evalPt, evalPt_new);
  //src_panel.print_info();
  //cout<<evalPt_new.x()<<" "<<evalPt_new.y()<<" "<<evalPt_new.z()<<endl;


  //determine if simple quadrature is neceesary if distance is far enough
//     if(!needSingleLayerPotential_ && needDoubleLayerPotential_)
//     {
//       point3D<T1 > dist(evalPt.x()-src_panel.centroid().x(),
//                         evalPt.y()-src_panel.centroid().y(), evalPt.z()-src_panel.centroid().z());
//       T1 rc=length(dist);
//       T1 L1=length(src_panel.vertex(1)-src_panel.vertex(2));
//       T1 L2=length(src_panel.vertex(2)-src_panel.vertex(3));
//       T1 maxLength=max(abs(L2), abs(L1));
//   
//       if ( abs(rc) > abs(maxLength*8.))
//       {
//         dlp=quad_panel(src_panel, evalPt_new, 16);
//         slp=0.;
//         return;
//       }
//     }

  if( needDoubleLayerPotential_  )
    getEvalPointDirection(evalPt_new);

  slp = std::complex<double>(0., 0.);
  dlp = std::complex<double>(0., 0.);
  totalArea = 0.;
 // cout << "panel shape" << src_panel.shape() << endl;
  for(size_t ii=0; ii<src_panel.shape(); ++ii)
  {
    setupCalcpTriangleInfo(ii, src_panel);
    contributionSign = getContributionSign();
    //cout<<contributionSign<<endl;
    totalArea += triangleInfo_.area * contributionSign;
    //cout<<totalArea<<endl;

    compIntegralOnTriangle(evalPt_new, slpLocal, dlpLocal);

    if ( needSingleLayerPotential_ )
    {
      slp += slpLocal * contributionSign;
    }
    if ( needDoubleLayerPotential_ )
    {
      dlp += dlpLocal * contributionSign;
    }
  }

  zn=-1.*evalPt_new.z();

  if (less_than(totalArea, 0))
  {
    slp=-1.*slp;
    dlp=-1.*dlp;
    //zn=-1.*zn;
  }
  if (isnan(abs(slp)) || isnan(abs(dlp)))
    cout << totalArea << endl;
//cout<<zn<<" "<<slp<<endl;
  //cout<<dlp<<endl;

}


template<class T1, class T2>
void calcpForEikrOverR<T1, T2>::
setupPanelLocalCoordSystem (panel<T1>& src_panel, const point3D<T2>& evalPnt,
                            point3D<T1>& evalPnt_new)
{

  // 1) find the projection of evaluation point on the panel
  point3D<T1> evalPntProjection(evalPnt.x(), evalPnt.y(), evalPnt.z());
  evalPntProjection.transferGlobalToLocalCoord(src_panel.centroid(),
      src_panel.tangent1(),
      src_panel.tangent2(),
      src_panel.normal());

  T1 h = evalPntProjection.z();
  evalPntProjection = point3D<T1 >(evalPntProjection.x(), evalPntProjection.y(), 0.);

  // 2) shift back to the global coord system to be
  // compatible with all verteices which are in the globl coord.
  evalPntProjection.transferLocalToGlobalCoord(
    src_panel.centroid(),
    src_panel.tangent1(),
    src_panel.tangent2(),
    src_panel.normal() );

  // 3) use projection point as origin of the panel local coordinate system
  src_panel.shiftToLocalCoord(evalPntProjection);
  //evalPnt_new = point3D<complex<double> >(0., 0., h);
  evalPnt_new.transferGlobalToLocalCoord(evalPntProjection, src_panel.tangent1(),
                                         src_panel.tangent2(), src_panel.normal());
  //cout<<evalPnt_new.x()<<" "<<evalPnt_new.y()<<" "<<evalPnt_new.z()<<endl;
}


/**********************************************************************
  * getEvalPointDirection --
**********************************************************************/
template<class T1, class T2>
void calcpForEikrOverR<T1, T2>::getEvalPointDirection (
  const point3D<T1 >& evalPnt)
{
  if (abs(evalPnt.z()) == 0.)  //abs(evalPnt.z())
  {
    /*if ((evalPointDirection_ != FOLLOW_NORMAL_DIRECTION) &&
        (evalPointDirection_ != AGAINST_NORMAL_DIRECTION))
    {
      cout<<"setupPanelLocalCoordSystem, illegal evalPointDirection"<<endl;
    }*/
    evalPointDirection_ =AGAINST_NORMAL_DIRECTION;
  }
  else
  {
    // presumably, evalPointDirection has not been specified by the user,
    // it is specified here.
    if (greater_than(evalPnt.z() , 0.) )
    {
      evalPointDirection_ = AGAINST_NORMAL_DIRECTION;
    }
    else
    {
      evalPointDirection_ = FOLLOW_NORMAL_DIRECTION;
    }
  }
}


/**********************************************************************
  * setupCalcpTriangleInfo --
**********************************************************************/
template<class T1, class T2>
void calcpForEikrOverR<T1, T2>::setupCalcpTriangleInfo ( size_t i,
    const panel<T1>& src_panel)
{
  size_t j = (i+1) % src_panel.shape();
  point3D<T1> PA = src_panel.vertex(i);
  point3D<T1> PB = src_panel.vertex(j);
  //cout<<i<<endl;
  //PA.print_info();
  //PB.print_info();

  point3D<T1 > Y = PB - PA;
  //Y.print_info();
  triangleInfo_.AB_len = length(Y);
  if (abs(length(Y)) > 1)
    cout << length(Y) << endl;
  //cout<<(triangleInfo_.AB_len)<<endl;
  Y.normalize();
  //Y.print_info();

  point3D<T1 > Z(0., 0., 1.);
  triangleInfo_.X = crossProd(Y, Z);
  triangleInfo_.Y = Y;

  T1 project = - (Y * PA);
  triangleInfo_.AH = Y * project;
  triangleInfo_.AH_len = abs_my(project);
  //cout<<triangleInfo_.AH_len<<endl;

  triangleInfo_.PH_unit = PA + triangleInfo_.AH;
  triangleInfo_.PH_len = length(triangleInfo_.PH_unit);
  if( abs (triangleInfo_.PH_len) <
      (1e1*DBL_EPSILON * abs(triangleInfo_.AB_len)) )
  {
    triangleInfo_.PH_unit = point3D<T1 >(0., 0., 0.);
    triangleInfo_.PH_len = 0.;
  }
  else
  {
    triangleInfo_.PH_unit /= triangleInfo_.PH_len;
  }

  project = Y * PB;
  triangleInfo_.HB = Y * project;

  triangleInfo_.area = 0.5 * triangleInfo_.PH_len * triangleInfo_.AB_len;

  if (abs(triangleInfo_.area) == 0)
  {
    projectPointOnEdge_ = true;
  }
}

template<class T1, class T2>
double calcpForEikrOverR<T1, T2>::getContributionSign (void)
{
  double contributionSign = 1.;

  /* Assuming view point is AGAINST_NORMAL_VECTOR */
  if( abs(triangleInfo_.area) == 0. )
  {
    contributionSign = -1.;
  }
  else
  {
    complex<double> sign = triangleInfo_.X * triangleInfo_.PH_unit;
    if (greater_than(sign, 0.))
    {  /* counter-clock-wise */
      contributionSign = +1.;
    }
    else if (less_than(sign , 0.))
    {  /* clock-wise */
      contributionSign = -1.;
    }
  }

  /* flip the sign if the view point is not as assumed */
  if (evalPointDirection_ == FOLLOW_NORMAL_DIRECTION)
  {
    contributionSign *= -1.;
  }

  return contributionSign;
}

template<class T1, class T2>
void calcpForEikrOverR<T1, T2>::compIntegralOnTriangle (
  const point3D<T1 >& evalPnt,
  std::complex<double>& slp,
  complex<double> & dlp)
{
  // calculate slp, dI/dx & dI/dz. They are double integrals */
  T1 AH_len_withSign = triangleInfo_.AH * triangleInfo_.Y;
  if ( (greater_than(0., AH_len_withSign)) ||
       ( less_than(triangleInfo_.AB_len,triangleInfo_.AH_len))
       || (needDoubleLayerPotential_)&&(abs(K_)==0))
  {
    //cout<<"H not on AB"<<endl;
    // H is not on AB, use gauss quadrature on edge AB directly
    T1 distance = triangleInfo_.AB_len;
    T1 yStart = -AH_len_withSign;
    compOuterIntegral(evalPnt, distance, yStart, slp, dlp);
  }
  else
  {
    //cout<<"H on AB"<<endl;
    // H is on AB, use gauss quadrature on segment HA and HB separately
    // first on HA
    std::complex<double> slpLocal1;
    complex<double>  dlpLocal1;
    T1 distance = triangleInfo_.AH_len;
    T1 yStart = -distance;
    compOuterIntegral(evalPnt, distance, yStart,
                      slpLocal1, dlpLocal1);
    //cout<<length<<" "<<yStart<<endl;

    // then on HB
    std::complex<double> slpLocal2;
    complex<double> dlpLocal2;
    T1 HB_len = triangleInfo_.AB_len - triangleInfo_.AH_len;
    distance = HB_len;
    yStart = 0.;
    compOuterIntegral(evalPnt, distance, yStart,
                      slpLocal2, dlpLocal2);
    //cout<<length<<" "<<yStart<<endl;

    slp = slpLocal1 + slpLocal2;
    dlp = dlpLocal1 + dlpLocal2;
    //cout<<slp<<" "<<dlp<<endl;

  }
}

/**********************************************************************
 * compOuterIntegral --
**********************************************************************/
template<class T1, class T2>
void calcpForEikrOverR<T1,T2> ::compOuterIntegral (
  const point3D<T1>& evalPnt,
  T1 distance,
  T1 yStart,
  std::complex<double>& slp,
  complex<double> & dlp)
{
  T1 h = evalPnt.z();
  T1 d = triangleInfo_.PH_len;
  //cout<<h<<" "<<d<<endl;
  complex<double> IK = -1.*std::complex<double>(0., 1.) * K_;
  T1 Ra = abs_my(h);
  //   cout<<Ra<<endl;

  slp = std::complex<double>(0., 0.);
  dlp = complex<double> (0., 0.);

  if (abs(IK) !=0)
  {
    for (size_t jj=0; jj < NumGQpoint_; ++jj)
    {
      T1 w = distance * GQweight_[jj];
      T1 y = distance * GQpoint_[jj] + yStart;

      T1 rho = sqrt(d*d + y*y);
      T1 R = sqrt(rho*rho + h*h);
      //cout<<rho<<" "<<R<<endl;

      if ( abs(triangleInfo_.area) != 0 )
      {
	if (rho == 0.)
	  cout << rho << endl;
        if ( needSingleLayerPotential_ )
        {
          std::complex<double> slpLocal = compSingleLayerInnerIntegral(IK, Ra, R);
          slp += (slpLocal) * w * d / (rho*rho);
        }
        if ( needDoubleLayerPotential_ )
        {
          std::complex<double> dlpLocal = compDoubleLayerInnerIntegral(IK, Ra, R);
          dlp += (dlpLocal) * w * d / (rho*rho);
        }
      }
    }
  }
  else
  {
    if (needDoubleLayerPotential_)
    {
      //quasi-static approximation
      T1 l_minus=yStart;
      T1 l_plus=yStart+distance;
      //cout<<l_minus<<" "<<l_plus<<endl;
      T1 R_minus=sqrt(l_minus*l_minus+d*d+h*h);
      T1 R_plus=sqrt(l_plus*l_plus+d*d+h*h);
      T1 R0=sqrt(h*h+d*d);

      if (((R_minus+l_minus)==0.) || ((R0*R0+Ra*R_plus)==0.) || ((R0*R0+Ra*R_minus)==0.))
        dlp = dlp+0.;
      else
        dlp=dlp+d*log((R_plus+l_plus)/(R_minus+l_minus))
	  -Ra*( atan_my((d*l_plus)/(R0*R0+Ra*R_plus)) - atan_my((d*l_minus)/(R0*R0+Ra*R_minus)));

      if (isnan(abs(slp)) || isnan(abs(dlp)))
	cout << d << h << R0 << Ra << R_plus << l_plus << distance << R_minus << l_minus << yStart << endl;
    }
  }
}

template<class T1, class T2>
complex<double> calcpForEikrOverR<T1,T2> ::compSingleLayerInnerIntegral(const complex<double>& IK,
    const T1 & Ra, const T1& R)
{
  complex<double> slp;
  if ( abs(IK*(R-Ra)) > 1.e-3 )
  {
    if ((R == 0.) || (Ra == 0.))
      cout << R << Ra << endl; 
    slp=((exp(IK*R)-1.)/R-(exp(IK*Ra)-1.)/Ra)/(IK*IK);
  }
  else
  {
    slp=(R-Ra)/2.;
    complex<double> temp=1./2.;
    T1 temp1=R;
    T1 temp2=Ra;
    for (int m=3; m<=5; m++)
    {
      temp=(IK)*temp/ static_cast<double>(m);
      temp1=temp1*R;
      temp2=temp2*Ra;
      slp=slp+temp*(temp1-temp2);
    }
  }
  return slp;
}

template<class T1, class T2>
complex<double> calcpForEikrOverR<T1,T2>::compDoubleLayerInnerIntegral(const complex<double>& IK,
    const T1& Ra, const T1& R)
{
  complex<double> dlp;

  if ( abs(IK*(R-Ra)) > 1.e-3 )
  {
    dlp = (exp(IK*R) - exp(IK*Ra)) / IK;
    //cout<<dlp<<endl;
  }
  else
  {
    complex<double> tmp = R-Ra;
    dlp = tmp;
    for (int m=2; m<=4; m++)
    {
      tmp *= IK*(R-Ra) / static_cast<double>(m);
      dlp += tmp;
    }
    dlp *= exp(IK*Ra);
  }
  
  return dlp;
}

/**********************************************************************
  * panel quadrature
**********************************************************************/
template<class T1, class T2>
complex<double> calcpForEikrOverR<T1, T2> ::quad_panel(panel<T1>& src_p,
    const point3D<T1>& eval, const int numQuadPts)
{
  vector<point3D<T1> > quad_pts;
  vector<T1> quad_weights;
  T1 R;
  complex<double> val=0.;

  src_p.panel_quadrature(numQuadPts, quad_pts, quad_weights);


  for (int i=0; i<quad_pts.size(); i++)
  {
    R=length(eval-quad_pts[i]);
    val=val+exp(-1.*IMAG* K_*R)/R*quad_weights[i];
  }

  //cout<<"from quad "<<val<<endl;
  return val;

}

/**********************************************************************
  * gaussQuadrature --
**********************************************************************/
template<class T1, class T2>
void calcpForEikrOverR<T1, T2> ::gaussQuadrature (
  const size_t n, // number of points
  std::vector<double>& w,  // weights
  std::vector<double>& x)   // points
{
  double x1 = 0;
  double x2 = 1;

  int m = (n+1)/2;
  double xm = .5 * (x2 + x1);
  double xl = .5 * (x2 - x1);
  double pp;

  for (int i = 1; i <= m; i++)
  {
    // Initial value of the ith root
    double z = cos(PI * (i-.25) / (n+.5));

    // Get the ith root by Newton method
    bool converge = false;
    do
    {
      // Evaluate Leg polynomial. p1 is the value at z
      double p1 = 1.0;
      double p2 = 0.0;
      for (int j=1; j <= n; j++)
      {
        double p3 = p2;
        p2 = p1;
        p1 = ((2.0*j-1.0) * z * p2 - (j-1.0) * p3) / j;
      }

      // pp is the derivative
      pp = n * (z*p1 - p2) / (z*z - 1.0);
      double z1 = z;
      z = z1 - p1/pp;

      if (abs(z-z1) < 1e1*DBL_EPSILON )
      {
        converge = true;
      }

    }
    while( !converge);

    x[i-1] = xm - xl*z;
    x[n+1-i-1] = xm + xl*z;
    w[i-1] = 2. * xl / ((1 - z*z) * pp * pp);
    w[n+1-i-1] = w[i-1];
  }
}

template<class T1, class T2>
bool calcpForEikrOverR<T1, T2>:: greater_than(complex<double> A , double B)
{
  return real	(A)>B;
}
template<class T1, class T2>
bool calcpForEikrOverR<T1, T2>:: greater_than( double A , complex<double> B)
{
  return A>real(B);
}
template<class T1, class T2>
bool calcpForEikrOverR<T1, T2>:: greater_than(double A , double B)
{
  return A>B;
}
template<class T1, class T2>
bool calcpForEikrOverR<T1, T2>:: less_than(complex<double> A , double B)
{
  return real	(A)<B;
}
template<class T1, class T2>
bool calcpForEikrOverR<T1, T2>:: less_than(complex<double> A , complex<double> B)
{
  return real	(A)< real(B);
}
template<class T1, class T2>
bool calcpForEikrOverR<T1, T2>:: less_than(double A , double B)
{
  return A<B;
}
template<class T1, class T2>
double calcpForEikrOverR<T1, T2>:: abs_my(double A)
{
  return abs(A);
}

template<class T1, class T2>
complex<double> calcpForEikrOverR<T1, T2>:: abs_my(complex<double> A )
{
  //return complex<double>(abs(real(A)), -abs(imag(A)));
  if (real(A)<0)
    return -1.*A;
  else
    return A;
}

template<class T1, class T2>
complex<double> calcpForEikrOverR<T1, T2>::atan_my(complex<double> A)
{
  return IMAG/2.*log((IMAG+A)/(IMAG-A));
}

template<class T1, class T2>
double calcpForEikrOverR<T1, T2>::atan_my(double A)
{
  return atan(A);
}

#endif
