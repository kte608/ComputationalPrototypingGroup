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


//this program is the c version of \int_v 1/R dr 
// accept both real and complex inputs
// this program is a part of magnetic potential A computation

#ifndef CALCAONEOVERR_H
#define CALCAONEOVERR_H

#include "const_def.h"
#include "filament.h"
#include<complex>
#include<math.h>
#include<iostream>

#define LEN 4
#define WID 2
#define HEIGHT 1
#define DEG_TOL 1.1e-4
#define nearzero(x) (abs(x) < EPS)

class calcAOneOverR
{
  enum degen_type {brick = 0, flat = 1, skinny = 2, too_long = 3, too_short = 4,
                   short_flat = 5, short_skinny = 6, impossible = 7};
public:
  calcAOneOverR(void){};

  template <class T1, class T2>
  T2 mutual(const filament<T1>* f1, const filament<T2>* f2 );

  template <class T>
  T selfterm(const  filament<T>* fil);

  bool less_than (complex<double>, double);
  bool less_than (complex<double>, complex<double>);
  bool less_than(double, double);
  bool less_equal_than (complex<double>, double);
  bool less_equal_than (double, double);
  bool less_equal_than (complex<double>, complex<double>);
  bool greater_than (complex<double>, double);
  bool greater_than (double, double);
  bool greater_than (complex<double>, complex<double>);
  bool greater_equal_than (complex<double>, double);
  bool greater_equal_than (double, double);
  complex<double> abs_mydef(complex<double> val);
  double abs_mydef(double val);

private:
  template<class T> T vdotp(T v1[], T v2[]);
  template <class T1, class T2>
  T2 dist_betw_fils(const filament<T1>* f1, const filament<T2>* f2 , bool& paral);
  template <class T1, class T2>
  T2 min_endpt_sep(const filament<T1>* fil1, const filament<T2>* fil2);
  template<class T2, class T1>
  T2 dist_betw_pt_and_fil_cx(const filament<T2>* fil2, const T2 D[],
                             const T2 s[], const T2& DD, const filament<T1>* fil1,
                             const T2& t);
  template<class T2, class T1>
  T2 dist_betw_pt_and_fil(const filament<T1>* fil_line, const T2 D[],
                          const T2 s[], const T2& DD,
                          const filament<T2>* fil, const T2& t);
  template <class T1, class T2>
  T2 mutualfil(const filament<T1>* fil1, const filament<T2>* fil2);
  template<class T1, class T2>
  T2 magdiff(const filament<T1>* fil1, int node1,
             const filament<T2>* fil2, int node2);
  template<class T1, class T2>
  T2 dotprod(const filament<T1>* fil1,
             const filament<T2>* fil2);
  template<class T> T mut_rect(T len, T d);
  template <class T1, class T2>
  bool edges_parallel(const filament<T1>* fil_j, const filament<T2>* fil_m,
                      T2 wid1[], bool& whperp);
  template <class T1, class T2>
  T2 parallel_fils(const filament<T1>* fil_j, const filament<T2>* fil_m,
                   const bool whperp, T2 x_j[], T2 y_j[], const T2 dist);
  template <class T>
  int find_deg_dims(const filament<T>* fil);
  template <class T1, class T2>
  T2  exact_mutual(const filament<T1>* fil_j, const filament<T2>* fil_m,
                   const bool whperp, T2 x_j[], T2 y_j[],
                   const degen_type deg_j, const degen_type deg_m);
  template <class T>
  T brick_to_brick(const T& E, const T& a,const T& d, const T& P,
                   const T& b,const T& c,const T& l3,
                   const T& l1,const T& l2);
  template <class T>
  void fill_4(T vec[], const T& E, const T& a, const T& d);
  template<class T>
  T eval_eq(const T& x,const T& y,const T& z, const T& ref_len);
  template<class T>
  T log_term(const T& x, const T& xsq, const T& ysq, const T& zsq, const T& len);
  template<class T>
  T tan_term(const T& x, const T& xsq, const T& ysq, const T& zsq, const T& len);
  template <class T1, class T2>
  T2 fourfil(const filament<T1>* fil_j, const filament<T2>* fil_m);
  template<class T>
  void findfourfils(const filament<T>* fil, filament<T> subfils[]);


  complex<double> atanh_my(complex<double>);
  double atanh_my(double);
  complex<double> asinh_my(complex<double>);
  double asinh_my(double);
  complex<double> atan2_my(complex<double>,complex<double>);
  double atan2_my(double,double);
  complex<double> atan_my(complex<double>);
  double atan_my(double);
  complex<double> MAX_my(complex<double>, complex<double>);
  double MAX_my(double, double);
};

template <class T>
T calcAOneOverR::selfterm(const  filament<T>* fil)
{
  T w=fil->filWidth()/fil->filLength();
  T t=fil->filHeight()/fil->filLength();
  T r = sqrt(w*w+t*t);
  T aw = sqrt(w*w+1.0);
  T at = sqrt(t*t+1.0);
  T ar = sqrt(w*w+t*t+1.0);
  T z;


  z = 0.25 * ((1./w) * asinh(w/at) + (1./t) * asinh(t/aw) + asinh(1./r));
  z += (1./24.0) * ((t*t/w) * asinh(w/(t*at*(r+ar))) + (w*w/t) * asinh(t/(w*aw*(r+ar))) +
                   ((t*t)/(w*w)) * asinh(w*w/(t*r*(at+ar))) + ((w*w)/(t*t))*asinh(t*t/(w*r*(aw+ar))) +
                   (1.0/(w*t*t)) * asinh(w*t*t/(at*(aw+ar))) + (1.0/(t*w*w))*asinh(t*w*w/(aw*(at+ar))));
  z -= (1.0/6.0) * ((1.0/(w*t)) * atan(w*t/ar) + (t/w) * atan(w/(t*ar)) + (w/t) * atan(t/(w*ar)));
  z -= (1.0/60.0) * ( ((ar+r+t+at)*t*t)/((ar+r)*(r+t)*(t+at)*(at+ar))
                      + ((ar+r+w+aw)*(w*w)) / ((ar+r)*(r+w)*(w+aw)*(aw+ar))
                      + (ar+aw+1.+at)/((ar+aw)*(aw+1.)*(1.+at)*(at+ar)));
  z -= (1.0/20.0)*((1.0/(r+ar)) + (1.0/(aw+ar)) + (1.0/(at+ar)));

  z *= (2.0/PI);
  z *= fil->filLength();  /* this is inductance */

  return z*MU0;

}


template <class T1, class T2>
T2 calcAOneOverR::mutual(const filament<T1>* fil_j, const filament<T2>* fil_m )
{
  bool parallel=false;
  bool edge_par, whperp;
  T2 totalM,rj,rm;
  T2 widj[3], heightj[3];

  T2 dist = dist_betw_fils(fil_j, fil_m, parallel);
  rj = MAX_my(fil_j->filWidth(), fil_j->filHeight())/2.0;
  rm = MAX_my(fil_m->filWidth(), fil_m->filHeight())/2.0;

  if (greater_than(dist,MAX_my(rj,rm)*100.))
  {
    //      	  cout<<"here1"<<endl;
    totalM = mutualfil(fil_j, fil_m);
    return totalM;
  }
  else
  {
    if (parallel == true)
    {
      widj[0]=fil_j->w_dir().x();
      widj[1]=fil_j->w_dir().y();
      widj[2]=fil_j->w_dir().z();

      heightj[0]=fil_j->h_dir().x();
      heightj[1]=fil_j->h_dir().y();
      heightj[2]=fil_j->h_dir().z();

      // 	  cout<<widj[0]<<" "<<widj[1]<<" "<<widj[2]<<endl;
      // 	  cout<<heightj[0]<<" "<<heightj[1]<<" "<<heightj[2]<<endl;
    }
    else
    {
      T2 temp=fil_j->l_dir().x()*fil_m->l_dir().x()+
              fil_j->l_dir().y()*fil_m->l_dir().y()+
              fil_j->l_dir().z()*fil_m->l_dir().z();
      if ( less_than(abs(temp)/(fil_j->filLength()*fil_m->filLength()), EPS)  )
        /* fils are perpendicular */
        return 0.0;
    }
    edge_par = parallel
               == 1 && edges_parallel(fil_j,fil_m,widj,whperp);
    //cout<<edge_par<<endl;

    if (edge_par && (2*MAX(abs(rj),abs(rm))*10) > abs(dist))
    {
      //          			cout<<"here3"<<endl;
      totalM = parallel_fils(fil_j, fil_m, whperp, widj, heightj, dist);
      return totalM;
    }
    else
    {
      //      			cout<<"here4"<<endl;
      totalM = fourfil(fil_j, fil_m);
      return totalM;
    }
    return 0;
  }
}

template <class T1, class T2>
T2 calcAOneOverR:: fourfil(const filament<T1>* fil_j, const filament<T2>* fil_m)
{
  filament<T2>  subfilm[4];
  filament<T1> subfilj[4];
  T2 totalM;
  findfourfils(fil_j,subfilj);
  findfourfils(fil_m, subfilm);
  int i;

  //   for (i=0; i < 4; i++)
  //   {
  // 	  cout<<subfilj[i].loc0().x()<<" "<<subfilj[i].loc0().y()<<" "<<subfilj[i].loc0().z()<<endl;
  // 	  cout<<subfilj[i].loc1().x()<<" "<<subfilj[i].loc1().y()<<" "<<subfilj[i].loc1().z()<<endl;
  //   }

  totalM = 0.0;
  for(i = 0; i < 4; i++)
  {
    totalM += mutualfil(fil_j, &subfilm[i]);
    // 	cout<<totalM<<endl;
  }
  for(i = 0; i < 4; i++)
    totalM += mutualfil(&subfilj[i],fil_m);

  totalM += -2.0*mutualfil(fil_j, fil_m);

  totalM = totalM/6.0;
  return totalM;
}

template<class T>
void calcAOneOverR:: findfourfils(const filament<T>* fil, filament<T>  subfils[])
{
  int i;
  T wx,wy,wz,hx,hy,hz,mag;

  wx = fil->w_dir().x();
  wy = fil->w_dir().y();
  wz = fil->w_dir().z();
  hx = fil->h_dir().x();
  hy = fil->h_dir().y();
  hz = fil->h_dir().z();

  point3D<T> inc(fil->filWidth()*wx/2.,fil->filWidth()*wy/2.,fil->filWidth()*wz/2.);
  subfils[0].set_endpt0(fil->loc0()+inc);
  subfils[0].set_endpt1(fil->loc1()+inc);
  subfils[1].set_endpt0(fil->loc0()-inc);
  subfils[1].set_endpt1(fil->loc1()-inc);
  point3D<T> inc2(fil->filHeight()*hx/2.,fil->filHeight()*hy/2.,fil->filHeight()*hz/2.);
  subfils[2].set_endpt0(fil->loc0()+inc2);
  subfils[2].set_endpt1(fil->loc1()+inc2);
  subfils[3].set_endpt0(fil->loc0()-inc2);
  subfils[3].set_endpt1(fil->loc1()-inc2);

  for(i = 0; i < 4; i++)
    subfils[i].set_length(fil->filLength());


}



template <class T1, class T2>
T2 calcAOneOverR:: parallel_fils(const filament<T1>* fil_j, const filament<T2>* fil_m,
                                 const bool whperp, T2 x_j[], T2 y_j[], const T2 dist)
{
  enum degen_type deg_j, deg_m;

  /* find degenerate dimensions */
  deg_j =(degen_type) find_deg_dims(fil_j);
  deg_m =(degen_type) find_deg_dims(fil_m);
  //    cout<<deg_j<<" "<<deg_m<<endl;


  if (deg_j == brick && deg_m == brick)
  {
    /* no degenerate dimensions, both are bricks */
    return exact_mutual(fil_j, fil_m, whperp, x_j, y_j, deg_j, deg_m);
  }
  else
  {
    cout<<"need compute_for_degenerate routine"<<endl;
    //     return compute_for_degenerate(fil_j, fil_m, whperp, x_j, y_j,
    //							  deg_j, deg_m, dist);
  }
}


template <class T>
int calcAOneOverR::find_deg_dims(const filament<T>* fil)
{
  T max;
  max = MAX_my(fil->filLength(), fil->filWidth());
  max = MAX_my(max, fil->filHeight());

  return (less_than(fil->filLength()/max,DEG_TOL))*LEN +
         less_than(fil->filWidth()/max, DEG_TOL)*WID
         + less_than(fil->filHeight()/max, DEG_TOL)*HEIGHT;
}


template <class T1, class T2>
T2 calcAOneOverR:: exact_mutual(const filament<T1>* fil_j, const filament<T2>* fil_m,
                                const bool whperp, T2 x_j[], T2 y_j[],
                                const degen_type deg_j, const degen_type deg_m)
{
  T2 a,b,c,d,l1,l2;
  T2 z_j[3], ox, oy, oz, origin[3], endx, endy, endz,E,P;
  T2 l3_1, l3, totalM;
  int sign, a_deg,b_deg,c_deg,d_deg;

  z_j[XX] = fil_j->l_dir().x();
  z_j[YY] = fil_j->l_dir().y();
  z_j[ZZ] = fil_j->l_dir().z();

  //cout<<z_j[XX]<<" "<<z_j[YY]<<" "<<z_j[ZZ]<<endl;

  a = fil_j->filWidth();
  b = fil_j->filHeight();
  if (whperp == false)
  {
    c = fil_m->filHeight();
    d = fil_m->filWidth();
  }
  else
  {
    d = fil_m->filHeight();
    c = fil_m->filWidth();
  }
  //   cout<<a<<" "<<b<<" "<<c<<" "<<d<<endl;

  ox = origin[XX] = fil_j->loc0().x() - x_j[XX]*a/2.0 - y_j[XX]*b/2.0;
  oy = origin[YY] = fil_j->loc0().y() - x_j[YY]*a/2.0 - y_j[YY]*b/2.0;
  oz = origin[ZZ] = fil_j->loc0().z() - x_j[ZZ]*a/2.0 - y_j[ZZ]*b/2.0;
  //   cout<<ox<<" "<<oy<<" "<<oz<<" "<<endl;

  endx = fil_m->loc0().x() - ox;
  endy = fil_m->loc0().y() - oy;
  endz = fil_m->loc0().z()- oz;
  //   cout<<endx<<" "<<endy<<" "<<endz<<" "<<endl;


  E = x_j[XX]*endx+x_j[YY]*endy+x_j[ZZ]*endz - d/2.;
  P = y_j[XX]*endx+y_j[YY]*endy+y_j[ZZ]*endz - c/2.;
  l3 = z_j[XX]*endx+z_j[YY]*endy+z_j[ZZ]*endz;
  l3_1 = z_j[XX]*(fil_m->loc1().x() - ox)+z_j[YY]*(fil_m->loc1().y() - oy)+
         z_j[ZZ]*(fil_m->loc1().z() - oz);

  //   cout<<E<<" "<<P<<" "<<l3<<" "<<l3_1<<endl;

  l1 = fil_j->filLength();
  l2 = fil_m->filLength();

  if ( greater_than(abs(abs(l3 - l3_1) - l2)/l2 , EPS))
  {
    fprintf(stderr, "Huh?  filament not expected length\n");
    exit(1);
  }

  if (less_equal_than(l3,l3_1)) //might not be right
    sign = 1;
  else
  {
    sign = -1;
    l3 = l3_1;
  }
  //cout<<sign<<endl;
  a_deg = less_than(a/MAX_my(l1,b) , DEG_TOL);
  b_deg = less_than(b/MAX_my(l1,a) , DEG_TOL);
  c_deg = less_than(c/MAX_my(l2,d) , DEG_TOL);
  d_deg = less_than(d/MAX_my(l2,c) , DEG_TOL);

  //   cout<<a_deg<<" "<<b_deg<<" "<<c_deg<<" "<<d_deg<<endl;

  if (a_deg && b_deg && c_deg && d_deg)
  {
    /* two long filaments */
    cout<<"need fourfil"<<endl;
    totalM = fourfil(fil_j, fil_m)/(sign*MU0/(4*PI));
  }
  else if (!a_deg && b_deg && c_deg && d_deg)
  {
    /* one flat and one long */
    cout<<"need tape_to_fil"<<endl;
    //totalM = tape_to_fil(E + d/2, a, P - b/2 + c/2,l3,l1,l2) / a;
  }
  else if (!a_deg && b_deg && c_deg && !d_deg)
  {
    /* two flat */
    cout<<"need flat_to_flat"<<endl;
    //totalM = flat_to_flat_tape(E,a,d,P - b/2 + c/2 ,l3,l1,l2) / (a * d);
  }
  else if (!a_deg && b_deg && !c_deg && d_deg)
  {
    /* one flat and one skinny */
    cout<<"need flat_to_skinny_tape"<<endl;
    //totalM = flat_to_skinny_tape(E + d/2,a,P - b/2 ,c,l3,l1,l2) / (a * c);
  }
  else if (!a_deg && !b_deg && c_deg && d_deg)
  {
    cout<<"need brick_to_fil"<<endl;
    //totalM = brick_to_fil(E + d/2,a,P + c/2,b,l3,l1,l2) / (a * b);
  }
  else if (deg_j != brick || deg_m != brick)
  {
    fprintf(stderr,"Internal Error: Bad Degenerate filament a/l1<tol\n");
    fprintf(stderr,"Using fourfil() instead...\n");
    cout<<"need fourfil"<<endl;
    totalM = fourfil(fil_j, fil_m)/(sign*MU0/(4*PI));
  }
  else
  {
    /* all that are left are bricks */
    /* this is the full 3-D filament calculation, no degeneracies */
    //  cout<<brick_to_brick(E,a,d,P,b,c,l3,l1,l2)<<" "<<(a * b * c * d)<<endl;
    //  cout<<E<<" "<<a<<" "<<d<<" "<<P<<" "<<b<<" "<<c<<" "<<l3<<" "<<l1<<" "<<l2<<endl;
    totalM = brick_to_brick(E,a,d,P,b,c,l3,l1,l2)/ (a * b * c * d);
  }
  return sign*MU0/(4*PI)*totalM;

}

template <class T>
T calcAOneOverR:: brick_to_brick(const T& E, const T& a,const T& d, const T& P,
                                 const T& b,const T& c,const T& l3,
                                 const T& l1,const T& l2)
{
  T q[4], r[4], s[4], totalM;
  int i,j,k;
  double sign2;

  fill_4(q, E,a,d);
  //   cout<<q[0]<<" "<<q[1]<<" "<<q[2]<<" "<<q[3]<<endl;
  fill_4(r, P,b,c);
  //   cout<<r[0]<<" "<<r[1]<<" "<<r[2]<<" "<<r[3]<<endl;
  fill_4(s, l3,l1,l2);
  //   cout<<s[0]<<" "<<s[1]<<" "<<s[2]<<" "<<s[3]<<endl;

  totalM = 0;

  for(i = 0; i < 4; i++)
    for(j = 0; j < 4; j++)
      for(k = 0; k < 4; k++)
      {
        sign2 = ( (i+j+k)%2 == 0 ? 1 : -1);
        totalM += sign2*eval_eq(q[i],r[j],s[k], a);
      }

  return totalM;
}

template<class T>
T calcAOneOverR:: eval_eq(const T& x,const T& y,const T& z, const T& ref_len)
{
  static double one_60 = 1.0/60.0;
  static double one_6 = 1.0/6.0;
  static double one_24 = 1.0/24.0;
  T  len, xsq, ysq, zsq, retval;
  int num_nearzero;
  int num_nearzero_sq;
  T one_over_ref_len, one_over_ref_len_sq;

  one_over_ref_len = 1.0/MAX_my(ref_len, (abs_mydef(x) + abs_mydef(y) + abs_mydef(z)));
  one_over_ref_len_sq = (one_over_ref_len)*(one_over_ref_len);

  xsq = x*x;
  ysq = y*y;
  zsq = z*z;

  len = sqrt(xsq + ysq + zsq);
  retval = one_60*len
           *(xsq*(xsq - 3.*ysq) + ysq*(ysq - 3.*zsq) + zsq*(zsq - 3.*xsq));

  num_nearzero = nearzero(x*one_over_ref_len)
                 + nearzero(y*one_over_ref_len)
                 + nearzero(z*one_over_ref_len);

  num_nearzero_sq = nearzero(xsq*one_over_ref_len_sq)
                    + nearzero(ysq*one_over_ref_len_sq)
                    + nearzero(zsq*one_over_ref_len_sq);

  if (num_nearzero_sq < 2)
    retval += one_24*(log_term(x, xsq, ysq, zsq, len)
                      + log_term(y, ysq, xsq, zsq, len)
                      + log_term(z, zsq, xsq, ysq, len));

  if (num_nearzero < 1)
    retval -= one_6*(tan_term(x,y,z,zsq,len) + tan_term(x,z,y,ysq,len)
                     + tan_term(z,y,x,xsq,len));

  return retval;
}

template<class T>
T calcAOneOverR:: log_term(const T& x, const T& xsq, const T& ysq, const T& zsq, const T& len)
{
  return ((6.*zsq - ysq)*ysq - zsq*zsq)*x*log( (x + len)/sqrt(ysq + zsq));
}

template<class T>
T calcAOneOverR:: tan_term(const T& x, const T& y, const T& z, const T& zsq, const T& len)
{
  return  x*y*z*zsq*atan_my(x*y/(z*len));
}


template <class T>
void calcAOneOverR::fill_4(T vec[], const T& E, const T& a, const T& d)
{
  vec[0] = E - a;
  vec[1] = E + d - a;
  vec[2] = E + d;
  vec[3] = E;
}


template <class T1, class T2>
bool calcAOneOverR:: edges_parallel(const filament<T1>* fil_j, const filament<T2>* fil_m,
                                    T2 wid1[], bool& whperp)

{
  T2 wid_j[3], wid_m[3], prod,mj,mm ;

  wid_j[0]=fil_j->w_dir().x();
  wid_j[1]=fil_j->w_dir().y();
  wid_j[2]=fil_j->w_dir().z();
  wid_m[0]=fil_m->w_dir().x();
  wid_m[1]=fil_m->w_dir().y();
  wid_m[2]=fil_m->w_dir().z();


  mj = sqrt(wid_j[XX]*wid_j[XX]+wid_j[YY]*wid_j[YY]+wid_j[ZZ]*wid_j[ZZ]);
  mm = sqrt(wid_m[XX]*wid_m[XX]+wid_m[YY]*wid_m[YY]+wid_m[ZZ]*wid_m[ZZ]);
  prod = (wid_j[XX]*wid_m[XX]+wid_j[YY]*wid_m[YY]+wid_j[ZZ]*wid_m[ZZ])
         / (mm*mj);
  if (abs(prod) < EPS)
  {
    whperp = true;   /* width and height are perpend. to that on other fil*/
    return true;
  }
  else if (abs( abs(prod) - 1. ) < EPS)
  {
    whperp = false;
    return true;
  }
  else
    return false;

}

template <class T1, class T2>
T2 calcAOneOverR::dist_betw_fils(const filament<T1>* fil1, const filament<T2>* fil2 , bool& parallel)
{
  T2 s1[3], s2[3], s1ms2[3], D1[3], D2[3];
  T2 D1D1, D1D2, D2D2, D1s1s2, D2s1s2;
  T2 tmp1, tmp2, t1, t2;
  T2  x1,y1,z1,x2,y2,z2;

  s1[XX] = fil1->loc0().x();
  s1[YY] = fil1->loc0().y();
  s1[ZZ] = fil1->loc0().z();
  s2[XX] = fil2->loc1().x();
  s2[YY] = fil2->loc1().y();
  s2[ZZ] = fil2->loc1().z();

  s1ms2[XX] = s1[XX] - s2[XX];
  s1ms2[YY] = s1[YY] - s2[YY];
  s1ms2[ZZ] = s1[ZZ] - s2[ZZ];


  D1[XX] = fil1->loc1().x() - fil1->loc0().x();
  D1[YY] = fil1->loc1().y() - fil1->loc0().y();
  D1[ZZ] = fil1->loc1().z() - fil1->loc0().z();

  D2[XX] = fil2->loc1().x() - fil2->loc0().x();
  D2[YY] = fil2->loc1().y() - fil2->loc0().y();
  D2[ZZ] = fil2->loc1().z() - fil2->loc0().z();

  D1D1 = vdotp(D1,D1);
  D1D2 = vdotp(D1,D2);
  D2D2 = vdotp(D2,D2);
  D1s1s2 = vdotp(D1,s1ms2);
  D2s1s2 = vdotp(D2,s1ms2);

  tmp1 = D1D2*D1D2/D1D1 - D2D2;


  if (abs(tmp1/D2D2) < EPS)
  {
    parallel = true; /* fils are parallel */
    return min_endpt_sep(fil1,fil2);


  }
  else
    parallel = false;

  t2 = (D1D2*D1s1s2/D1D1 - D2s1s2)/tmp1;
  t1 = (t2*D1D2 - D1s1s2)/D1D1;

  if (less_equal_than(t1,1.) && greater_equal_than(t1,0.))
  {
    if (less_equal_than(t2,1.) && greater_equal_than(t2,0.))
    {
      x1=s1[0]+t1*D1[0];
      y1=s1[1]+t1*D1[1];
      z1=s1[2]+t1*D1[2];
      x2=s2[0]+t2*D2[0];
      y2=s2[1]+t2*D2[1];
      z2=s2[2]+t2*D2[2];
      return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
    }
    else
    {
      return dist_betw_pt_and_fil(fil1, D1, s1, D1D1, fil2,t2);
    }
  }
  else
  {
    if (less_equal_than(t2, 1.) && greater_equal_than(t2,0.))
    {
      /* nearest point along fil1 is outside line segment defining filament */
      return dist_betw_pt_and_fil_cx(fil2, D2, s2, D2D2, fil1, t1 );
    }
    else
      /* both point are out of range, just compare endpoints */
      return min_endpt_sep(fil1,fil2);
  }
}




template <class T2, class T1>
T2 calcAOneOverR::dist_betw_pt_and_fil_cx(const filament<T2>* fil_line, const T2 D[],
    const T2 s[], const T2& DD,
    const filament<T1>* fil, const T2& t)

{
  T2 e[3], sme[3], x,y,z;
  T2 tnew, Dsme;
  int idx;

  if (less_than(t,0.))
  {
    e[XX] = fil->loc0().x();
    e[YY] = fil->loc0().y();
    e[ZZ] = fil->loc0().z();
  }
  else if (greater_than(t,1.))
  {
    e[XX] = fil->loc1().x();
    e[YY] = fil->loc1().y();
    e[ZZ] = fil->loc1().z();
  }
  else
  {
    fprintf(stderr, "Internal err: dist_bet_pt_and_fil");
    exit(1);
  }
  sme[XX] = s[XX] - e[XX];
  sme[YY] = s[YY] - e[YY];
  sme[ZZ] = s[ZZ] - e[ZZ];
  Dsme = D[XX]*sme[XX]+D[YY]*sme[YY]+D[ZZ]*sme[ZZ];
  tnew = -Dsme/DD;
  if (less_equal_than(tnew, 1.) && greater_equal_than(tnew ,0.))
  {
    x=sme[0]+tnew*D[0];
    y=sme[1]+tnew*D[1];
    z=sme[2]+tnew*D[2];
    return sqrt(x*x + y*y + z*z);
  }
  else
  {
    if (less_than(tnew ,0.))
      return sqrt((e[XX]-fil_line->loc0().x())*(e[XX]-fil_line->loc0().x())+
                  (e[YY]- fil_line->loc0().y())*(e[YY]- fil_line->loc0().y())+
                  (e[ZZ]-fil_line->loc0().z())*(e[ZZ]-fil_line->loc0().z()));
    else
      return sqrt((e[XX]-fil_line->loc1().x())*(e[XX]-fil_line->loc1().x())+
                  (e[YY]- fil_line->loc1().y())*(e[YY]- fil_line->loc1().y())+
                  (e[ZZ]-fil_line->loc1().z())*(e[ZZ]-fil_line->loc1().z()));
  }
}

template <class T2, class T1>
T2 calcAOneOverR::dist_betw_pt_and_fil(const filament<T1>* fil_line, const T2 D[],
                                       const T2 s[], const T2& DD,
                                       const filament<T2>* fil, const T2& t)

{
  T2 e[3], sme[3], x,y,z;
  T2 tnew, Dsme;
  int idx;

  if (less_than(t,0.))
  {
    e[XX] = fil->loc0().x();
    e[YY] = fil->loc0().y();
    e[ZZ] = fil->loc0().z();
  }
  else if (greater_than(t,1.))
  {
    e[XX] = fil->loc1().x();
    e[YY] = fil->loc1().y();
    e[ZZ] = fil->loc1().z();
  }
  else
  {
    fprintf(stderr, "Internal err: dist_bet_pt_and_fil");
    exit(1);
  }
  sme[XX] = s[XX] - e[XX];
  sme[YY] = s[YY] - e[YY];
  sme[ZZ] = s[ZZ] - e[ZZ];
  Dsme = D[XX]*sme[XX]+D[YY]*sme[YY]+D[ZZ]*sme[ZZ];
  tnew = -Dsme/DD;
  if (less_equal_than(tnew, 1.) && greater_equal_than(tnew ,0.))
  {
    x=sme[0]+tnew*D[0];
    y=sme[1]+tnew*D[1];
    z=sme[2]+tnew*D[2];
    return sqrt(x*x + y*y + z*z);
  }
  else
  {
    if (less_than(tnew ,0.))
      return sqrt((e[XX]-fil_line->loc0().x())*(e[XX]-fil_line->loc0().x())+
                  (e[YY]- fil_line->loc0().y())*(e[YY]- fil_line->loc0().y())+
                  (e[ZZ]-fil_line->loc0().z())*(e[ZZ]-fil_line->loc0().z()));
    else
      return sqrt((e[XX]-fil_line->loc1().x())*(e[XX]-fil_line->loc1().x())+
                  (e[YY]- fil_line->loc1().y())*(e[YY]- fil_line->loc1().y())+
                  (e[ZZ]-fil_line->loc1().z())*(e[ZZ]-fil_line->loc1().z()));
  }
}


template <class T1, class T2>
T2 calcAOneOverR::min_endpt_sep(const filament<T1>* fil1, const filament<T2>* fil2)
{
  T2 min, tmp;

  min=magdiff(fil1, 0, fil2, 0);
  tmp=magdiff(fil1, 1, fil2, 0);
  if (abs(tmp) < abs(min)) min = tmp;
  tmp=magdiff(fil1, 0, fil2, 1);
  if (abs(tmp) < abs(min)) min = tmp;
  tmp=magdiff(fil1, 1, fil2, 1);
  if (abs(tmp) < abs(min)) min = tmp;
  return min;
}

template <class T1, class T2>
T2 calcAOneOverR::mutualfil(const filament<T1>* fil1, const filament<T2>* fil2)
{
  T2 R1, R2, R3, R4, maxR, minR, alpha, alpha2,R, M, blah;
  T2 R1sq, R2sq, R3sq, R4sq;
  T2 cose, realcos, sinsq,sine, tmp1, tmp2, tmp3, tmp4, tmp5;
  T2 omega;
  T2 Rx, Ry, Rz, ux, uy,uz,magu,dotp,dx,dy,dz,d,vx,vy,vz,x2_0,x2_1,vtemp;
  T2 u,v,u2,v2,l,m,x1_1,l2, m2, maxlength;
  double x1_0;
  double scaleEPS;
  int  signofM,realcos_error;
  double widj[3];

  R1=magdiff(fil1,1,fil2,1);
  R1sq =R1*R1;
  R2=magdiff(fil1,1,fil2,0);
  R2sq = R2*R2;
  R3 =magdiff(fil1,0,fil2,0);
  R3sq = R3*R3;
  R4=magdiff(fil1,0,fil2,1);
  R4sq = R4*R4;
  \

  maxR = minR = R1;
  if (less_than(maxR,R2)) maxR = R2;
  if (less_than(maxR,R3)) maxR = R3;
  if (less_than(maxR,R4)) maxR = R4;
  if (less_than(R2,minR)) minR = R2;
  if (less_than(R3,minR)) minR = R3;
  if (less_than(R4,minR)) minR = R4;
  //   cout<<maxR<<endl;
  //   cout<<minR<<endl;

  l = fil1->filLength();
  m = fil2->filLength();
  maxlength = MAX_my(l,m);

  //   cout<<maxlength<<endl;

  scaleEPS = abs(minR/maxlength)*10;
  if (scaleEPS < 1) scaleEPS = 1;
  if (scaleEPS > 100) scaleEPS = 100;

  alpha = R4*R4 - R3*R3 + R2*R2 - R1*R1;
  signofM = 1;
  //   cout<<alpha<<" "<<signofM<<" "<<scaleEPS<<endl;

  if ( (abs(R1) < EPS)||(abs(R2) < EPS)||(abs(R3) < EPS)||(abs(R4) < EPS) )
  {
    if (abs(R1) < EPS)  R = R3;
    else if (abs(R2) < EPS)  R = R4;
    else if (abs(R3) < EPS)  R = R1;
    else R = R2;

    M = MU0/(4*PI)*2*(dotprod(fil1,fil2)/(l*m))
        *(l*atanh_my(m/(l+R)) + m*atanh_my(l/(m+R)));

    return M;
  }

  cose = alpha/(2.*l*m);
  if (abs(cose) > 1) cose = (less_than(cose, 0.) ? -1.0 : 1.0);
  //blah = 1.0 - abs_mydef(cose);
  realcos = dotprod(fil1, fil2)/(l*m);
  if (abs(realcos) < EPS)
    return 0.0;
  if (greater_than(abs((realcos - cose)/cose), 0.1))
    if (realcos_error == 0)
    {
      fprintf(stderr, "Internal Warning: realcos =0");
      fprintf(stderr,"  This may be due to two filaments that are separated \n\
              by a distance 1e10 times their length\n");
      realcos_error = 1;
    }
  cose = realcos;

  //   cout<<cose<<endl;

  //tmp1 = abs( abs(cose) - 1);
  if ( abs( abs(cose) - 1) < EPS)
  {
    // 	  cout<<"here"<<endl;
    Rx = fil2->loc0().x() - fil1->loc0().x();  /* vector from fil1 to fil2 */
    Ry = fil2->loc0().y() - fil1->loc0().y();
    Rz = fil2->loc0().z() - fil1->loc0().z();
    ux = fil1->loc1().x() - fil1->loc0().x();
    uy = fil1->loc1().y() - fil1->loc0().y();
    uz = fil1->loc1().z() - fil1->loc0().z();
    //     	cout<<Rx<<" "<<Ry<<" "<<Rz<<endl;
    magu = sqrt(ux*ux + uy*uy + uz*uz);
    ux = ux/magu;    /* unit vector in direction of fil1 */
    uy = uy/magu;
    uz = uz/magu;
    // 	cout<<ux<<" "<<uy<<" "<<uz<<endl;

    dotp = ux*Rx + uy*Ry + uz*Rz;  /* component of R in direction of fil1 */
    //     	cout<<dotp<<endl;

    /* d vector is R vector without its component in the direction of fils */
    dx = Rx - dotp*ux;
    dy = Ry - dotp*uy;
    dz = Rz - dotp*uz;
    d = sqrt(dx*dx + dy*dy + dz*dz);
    //     cout<<dx<<" "<<dy<<" "<<dz<<" "<<d<<endl;


    /* let fil1 be the x axis, with node 0 being origin and u be */
    /* its positive direction */
    x1_0 = 0;
    x1_1 = l;
    /* (dotproduct just gives it correct sign) */
    vx =  (fil2->loc0().x() - ( fil1->loc0().x() + dx));
    vy =  (fil2->loc0().y() - ( fil1->loc0().y() + dy));
    vz =  (fil2->loc0().z() - ( fil1->loc0().z() + dz));
    x2_0 = vx*ux + vy*uy + vz*uz;
    vtemp = sqrt(vx*vx + vy*vy + vz*vz);

    /* same thing for x2_1 */
    vx =  (fil2->loc1().x() - ( fil1->loc0().x() + dx));
    vy =  (fil2->loc1().y() - ( fil1->loc0().y() + dy));
    vz =  (fil2->loc1().z() - ( fil1->loc0().z() + dz));
    x2_1 = vx*ux + vy*uy + vz*uz;

    if ( abs( (sqrt(vx*vx + vy*vy + vz*vz) - abs(x2_1))
              /(MAX(abs(abs(x2_0)+d),abs(abs(x2_1)+d)))) > EPS)
      printf("uh oh, segs don't seem parallel \n");
    if ( abs( (vtemp - abs(x2_0))/(MAX(abs(abs(x2_0)+d),abs(abs(x2_1)+d)))) > EPS)
      printf("uh oh, segs don't seem parallel\n");

    if ( abs(d) < EPS )
    { /* collinear! */
      M = (MU0/(4*PI))*(abs(x2_1 - x1_0)*log(abs(x2_1 - x1_0))
                        - abs(x2_1 - x1_1)*log(abs(x2_1 - x1_1))
                        - abs(x2_0 - x1_0)*log(abs(x2_0 - x1_0))
                        + abs(x2_0 - x1_1)*log(abs(x2_0 - x1_1)) );
      return M;
    }  /* end collinear */
    M = MU0/(4*PI)*(mut_rect(x2_1 - x1_1,d) - mut_rect(x2_1 - x1_0,d)
                    - mut_rect(x2_0 - x1_1,d) + mut_rect(x2_0 - x1_0,d) );

    return M;
  }

  l2 = l*l;
  m2 = m*m;
  alpha2 = alpha*alpha;

  u = l*(2.*m2*(R2sq -R3sq - l2) + alpha*(R4sq - R3sq - m2))
      / (4.*l2*m2 - alpha2);
  v = m* (2.*l2*(R4sq - R3sq - m2) + alpha*(R2sq - R3sq - l2))
      / (4.*l2*m2 - alpha2);

  u2 = u*u;
  v2 = v*v;

  d = (R3sq - u2 - v2 + 2.*u*v*cose);
  if (abs(d/(R3sq + u2 + v2 + 1.)*(maxlength*maxlength/(maxR*maxR))) < EPS)
    d = 0.0;
  d = sqrt(d);

  sinsq = 1.0 - cose*cose;
  if (abs(sinsq) < EPS) sinsq = 0.0;
  sine = sqrt(sinsq);
  tmp1 = d*d*cose;
  tmp2 = d*sine;
  tmp3 = sine*sine;

  if (abs(d) < EPS)
    omega = 0.0;   /* d is zero, so it doesn't matter */
  else
  {
    omega = atan2_my( (tmp1 + (u+l)*(v + m)*tmp3),(tmp2*R1))
            - atan2_my( (tmp1 + (u + l)*v*tmp3),(tmp2*R2))
            + atan2_my( (tmp1 + u*v*tmp3),(tmp2*R3))
            - atan2_my( (tmp1 + u*(v + m)*tmp3),(tmp2*R4) );
  }
  tmp4 =(  (u+l)*atanh_my( m/(R1 + R2))
           +(v+m)*atanh_my( l/(R1 + R4))
           -    u*atanh_my( m/(R3 + R4))
           -    v*atanh_my( l/(R2 + R3))  );

  if ( abs(sine) < 1e-150) tmp5 = 0.0;
  else tmp5 = omega*d/sine;

  M = MU0/(4*PI)*cose*(2.*tmp4 - tmp5 );

  return M;
}

template<class T>
T calcAOneOverR:: mut_rect(T len, T d)
{
  T temp,temp1;
  temp = sqrt(len*len + d*d);
  temp1 = len*asinh_my(len/d) ;
  return temp - temp1;
}


template<class T>
T calcAOneOverR::vdotp(T v1[], T v2[])
{
  return v1[XX]*v2[XX] + v1[YY]*v2[YY] + v1[ZZ]*v2[ZZ];
}



template<class T1, class T2>
T2 calcAOneOverR::dotprod(const filament<T1>* fil1,
                          const filament<T2>* fil2)
{
  return(  (fil1->loc1().x() - fil1->loc0().x())*(fil2->loc1().x() - fil2->loc0().x())
           + (fil1->loc1().y() - fil1->loc0().y())*(fil2->loc1().y() - fil2->loc0().y())
           + (fil1->loc1().z() - fil1->loc0().z())*(fil2->loc1().z() - fil2->loc0().z()) );
}

template<class T1, class T2>
T2 calcAOneOverR:: magdiff(const filament<T1>* fil1, int node1,
                           const filament<T2>* fil2, int node2)
{

  if ((node1==0) && (node2==0))
    return length(point3D<T2> ( (fil1->loc0().x()-fil2->loc0().x()),
                                (fil1->loc0().y()-fil2->loc0().y()),
                                fil1->loc0().z()-fil2->loc0().z()));
  if ((node1==1) && (node2==0))
    return length(point3D<T2> ( (fil1->loc1().x()-fil2->loc0().x()),
                                (fil1->loc1().y()-fil2->loc0().y()),
                                fil1->loc1().z()-fil2->loc0().z()));
  if ((node1==0) && (node2==1))
	  return length(point3D<T2> ( (fil1->loc0().x()-fil2->loc1().x()),
					(fil1->loc0().y()-fil2->loc1().y()),
					fil1->loc0().z()-fil2->loc1().z()));
  if ((node1==1) && (node2==1))
	  return length(point3D<T2> ( (fil1->loc1().x()-fil2->loc1().x()),
					(fil1->loc1().y()-fil2->loc1().y()),
					fil1->loc1().z()-fil2->loc1().z()));


}

#endif
