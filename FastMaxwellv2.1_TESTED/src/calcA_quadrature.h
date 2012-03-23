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
//calcAEikrOverR handles computation of slp= \int_v \int V e^{-jkr}/r dr'dr
//dIdx= \int_v \int V d/dx(e^{-jkR}/R) dr'dr
// dIdy= \int_v \int V d/dy(e^{-jkR}/R) dr'dr
// r can be real or complex
// this program is a part of magnetic potential A computation

#ifndef CALCAQUADRATURE_H_
#define CALCAQUADRATURE_H_

#include<iostream>
#include<math.h>

template <class T1, class T2>
class calcAQuadrature
{
public:
  calcAQuadrature(void){};
  void init(const complex<double> K_num, const bool needSLP, const bool needDLP);
  void operator() (filament<T1> src, filament<T2> eval,
                   complex<double>& slp, complex<double>& dIdx, complex<double>& dIdy);
  void selfterm (filament<T1> fil,complex<double>& slp);

private:
  void compute_fullQuadrature(filament<T1> src, filament<T2> eval,
                              complex<double>& slp, complex<double>& dIdx,  complex<double>& dIdy);
  complex<double> K_;
  bool needSingleLayerPotential_;
  bool needDoubleLayerPotential_;
  bool needNothing_;


};



template <class T1, class T2>
void calcAQuadrature<T1, T2>::init (
  const complex<double> K_num,
  const bool needSLP,
  const bool needDLP)
{
  K_=(K_num);
  needSingleLayerPotential_=(needSLP);
  needDoubleLayerPotential_=(needDLP);

  needNothing_ = (!needSLP) && (!needDLP);

}


template<class T1, class T2>
void calcAQuadrature<T1, T2>::
operator()(filament<T1> src, filament<T2> eval,
           complex<double>& slp, complex<double>& dIdx,  complex<double>& dIdy)
{
  if (needNothing_) return;

  slp=0., dIdx=0., dIdy=0.;


  compute_fullQuadrature(src, eval, slp, dIdx, dIdy);


  if ( needSingleLayerPotential_ )
  {
    slp=MU0*slp/(4*PI*src.filCrossA()*eval.filCrossA());
  }
  if ( needDoubleLayerPotential_ )
  {
    dIdx=MU0*dIdx/(4*PI*src.filCrossA()*eval.filCrossA());
    dIdy=MU0*dIdy/(4*PI*src.filCrossA()*eval.filCrossA());
  }

}


template<class T1, class T2>
void calcAQuadrature<T1, T2>::compute_fullQuadrature(filament<T1> src, filament<T2> eval,
    complex<double>& slp, complex<double>& dIdx, complex<double>& dIdy)
{
  slp=dIdx=dIdy=0.;
  //not equipped to handle singularity
  size_t numPtsW=3;
  size_t numPtsH=3;
  size_t numPtsL=300;

  vector<point3D<T1> > quad_pts_src;
  vector<T1> quad_weights_src;
  vector<point3D<T2> > quad_pts_eval;
  vector<T2> quad_weights_eval;
  T1 R;
  complex<double> dlpLocal=0.;

  //obtain volume quadrature pts for the eval. and src. filament
  quad_pts_eval.resize(numPtsW*numPtsH*numPtsL);
  quad_weights_eval.resize(numPtsW*numPtsH*numPtsL);

  quad_pts_src.resize(numPtsW*numPtsH*numPtsL);
  quad_weights_src.resize(numPtsW*numPtsH*numPtsL);

  eval.volume_quadrature(numPtsW, numPtsH, numPtsL, quad_pts_eval, quad_weights_eval);
  src.volume_quadrature(numPtsW, numPtsH, numPtsL, quad_pts_src, quad_weights_src);

  for (int i=0; i< numPtsW*numPtsH*numPtsL; i++)
  {
    for (int j=0; j<numPtsW*numPtsH*numPtsL; j++)
    {
      R=length(point3D<T1>( quad_pts_eval[i].x()-quad_pts_src[j].x(),
                            quad_pts_eval[i].y()-quad_pts_src[j].y(),
                            quad_pts_eval[i].z()-quad_pts_src[j].z()));

	  if (needSingleLayerPotential_)
      	slp+=exp(-IMAG*K_*R)/R*quad_weights_eval[i]*quad_weights_src[j];
	  if (needDoubleLayerPotential_){
	  dlpLocal=(-IMAG*K_/R-1./(R*R))*exp(-IMAG*K_*R);
      dIdx+=(quad_pts_eval[i].x()-quad_pts_src[j].x())/R*dlpLocal
            *quad_weights_eval[i]*quad_weights_src[j];

      dIdy+=(quad_pts_eval[i].y()-quad_pts_src[j].y())/R*dlpLocal
            *quad_weights_eval[i]*quad_weights_src[j];
	  }
    }
  }

}

//note that self-term is only encountered for real filaments and only slp is needed
template<class T1, class T2>
void calcAQuadrature<T1, T2>::selfterm(filament<T1> src,complex<double>& slp)
{

  slp=0.;

  size_t numPtsW=1;
  size_t numPtsH=1;
  size_t numPtsL=10;

  vector<point3D<T1> > quad_pts_src;
  vector<T1> quad_weights_src;
  complex<double> val,temp2;
  int temp;
  double R;
  calcAOneOverR calcA_static;


  //obtain volume quadrature pts for the eval. and src. filament
  quad_pts_src.resize(numPtsW*numPtsH*numPtsL);
  quad_weights_src.resize(numPtsW*numPtsH*numPtsL);
  src.volume_quadrature(numPtsW, numPtsH, numPtsL, quad_pts_src, quad_weights_src);


  if (abs(K_)>0.)
  {
    double threshold=0.0005*(1./(pow(abs(K_),1.25)));
    for (int i=0; i< numPtsW*numPtsH*numPtsL; i++)
    {
      for (int j=0; j<numPtsW*numPtsH*numPtsL; j++)
      {
        R=length(point3D<T1>( quad_pts_src[i].x()-quad_pts_src[j].x(),
                              quad_pts_src[i].y()-quad_pts_src[j].y(),
                              quad_pts_src[i].z()-quad_pts_src[j].z()));

        if (abs(R)>threshold)
          val=(exp(-IMAG*K_*R)-1.)/R;
        else
        {
          val=0.;
          for (int n=1; n<5; n++)
          {
            temp=1;
            temp2=1;
            for (int m=1; m<n+1; m++)
            {
              temp=temp*n;
              temp2=temp2*(-IMAG*K_);
            }
            val=val+1./temp*temp2*pow(R, n-1);
          }
          //cout<<val<<endl;
        }
        slp=slp+val*quad_weights_src[i]*quad_weights_src[j];
      }
    }
  }

  slp=MU0*slp/(4*PI*src.filCrossA()*src.filCrossA());
  temp2 = calcA_static.selfterm(&src);
  slp=slp + temp2;
}




#endif
