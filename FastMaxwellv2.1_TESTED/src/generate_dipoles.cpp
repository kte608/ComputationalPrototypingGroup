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
#include "generate_dipoles.h"
#include "const_def.h"
#include <math.h>
#include <stdlib.h>
//#include "/usr/local/include/itpp/itbase.h"
#include "../lib/itpp-4.0.1/itpp/itbase.h"

/*
#include "it++/mat.h"
#include "it++/svd.h"
#include "it++/matfunc.h"
#include "it++/eigen.h"
#include "it++/ls_solve.h"
*/

Dipoles::Dipoles(const LayerEnv& layers)
{
  MAX_CONDUCTORSYS_RANGE[0]=0.01;
  MAX_CONDUCTORSYS_RANGE[1]=0.01;
  MAX_CONDUCTORSYS_RANGE[2]=0.05;
  THRESHOLD=1e-4;

  if (layers.num_layers()==1)
  {
    dipoles_flag_=false;
  }
  else
  {
    substrateHeight=layers.get_subHeight();
    dipoles_flag_=true;
    k1=layers.list_layers[0].k();
    k2=layers.list_layers[1].k();
    n=k2/k1;
    //cout << n << endl;
	
//	cout<<"Approximate coes and exps for Gazz"<<endl;
    approx_integral_2levels(G_A_ZZ, coes_G_A_ZZ, exps_G_A_ZZ);

    for (int i=0; i<coes_G_A_ZZ.size(); i++) 
      coes_G_A_ZZ_MOD.push_back(-coes_G_A_ZZ[i]);
    
    //	cout<<"Approximate coes and exps for Gaxx"<<endl;
    approx_integral_2levels(G_A_XX, coes_G_A_XX, exps_G_A_XX);
    //	cout<<"Approximate coes and exps for Gaxz"<<endl;
    //	cout<<"Mistake in the approximation because "
    //    <<"reflection coefficient is not normalized"<<endl;
    approx_integral_2levels(G_A_XZ, coes_G_A_XZ, exps_G_A_XZ);
    //divide all a's by k0
    for (int i=0; i<coes_G_A_XZ.size(); i++) {
      coes_G_A_XZ[i]=(1./k1)*coes_G_A_XZ[i];
      //cout<<exps[i]<<endl;
    }
    
    //	cout<<"Approximate coes and exps for G_phi_h"<<endl;

    approx_integral_2levels(G_PHI_H, coes_G_PHI_H, exps_G_PHI_H);
    //for (int i=0; i<coes_G_PHI_H.size(); i++)
    //  cout<<exps_G_PHI_H[i]<<endl;
    //	cout<<"Approximate coes and exps for modified Gaxx"<<endl;
    approx_integral_2levels(G_A_XX_MOD, coes_G_A_XX_MOD, exps_G_A_XX_MOD);
    
    // 	cout<<coes_G_PHI_H.size()<<endl;
    //     for (int i=0; i<coes_G_A_XX.size(); i++)
    //       cout<<coes_G_PHI_H[i]<<" "<<exps_G_PHI_H[i]<<endl;
    
    //divide all a's by k0
    for (int i=0; i<coes_G_A_XX_MOD.size(); i++) {
      coes_G_A_XX_MOD[i]=(-TWO/k1/k1)*coes_G_A_XX_MOD[i];
      //cout<<exps[i]<<endl;
    }
    
  }
}

void Dipoles::setup(const LayerEnv& layers)
{
  MAX_CONDUCTORSYS_RANGE[0]=0.01;
  MAX_CONDUCTORSYS_RANGE[1]=0.01;
  MAX_CONDUCTORSYS_RANGE[2]=0.05;
  THRESHOLD=1e-4;

  if (layers.num_layers()==1)
  {
    dipoles_flag_=false;
  }
  else
  {
    substrateHeight=layers.get_subHeight();
    dipoles_flag_=true;
    k1=layers.list_layers[0].k();
    k2=layers.list_layers[1].k();
    n=k2/k1;

	cout<<"Approximate coes and exps for Gazz"<<endl;
    approx_integral_2levels(G_A_ZZ, coes_G_A_ZZ, exps_G_A_ZZ);

  for (int i=0; i<coes_G_A_ZZ.size(); i++)
	  coes_G_A_ZZ_MOD.push_back(-coes_G_A_ZZ[i]);
  
  cout<<"Approximate coes and exps for Gaxx"<<endl;
  approx_integral_2levels(G_A_XX, coes_G_A_XX, exps_G_A_XX);
  cout<<"Approximate coes and exps for Gaxz"<<endl;
  //cout<<"Mistake in the approximation because "
  //  <<"reflection coefficient is not normalized"<<endl;
  approx_integral_2levels(G_A_XZ, coes_G_A_XZ, exps_G_A_XZ);
  //divide all a's by k0
  for (int i=0; i<coes_G_A_XZ.size(); i++) 
    coes_G_A_XZ[i]=(1./k1)*coes_G_A_XZ[i];
  
  cout<<"Approximate coes and exps for G_phi_h"<<endl;
  approx_integral_2levels(G_PHI_H, coes_G_PHI_H, exps_G_PHI_H);
  cout<<"Approximate coes and exps for modified Gaxx"<<endl;
  approx_integral_2levels(G_A_XX_MOD, coes_G_A_XX_MOD, exps_G_A_XX_MOD);
  
  //divide all a's by k0
  for (int i=0; i<coes_G_A_XX_MOD.size(); i++) {
    coes_G_A_XX_MOD[i]=(-TWO/k1/k1)*coes_G_A_XX_MOD[i];
    //cout<<exps[i]<<endl;
  }

  }
}


//determine if the remainder of kernel-qs is small enough to use only
//quasistatic approximation
void Dipoles::approx_integral_2levels( dipoleType type,
                                       vector<Complex>& coes, vector<Complex>& exps)
{
  vector<Complex> a, b, a2, b2;	
  Complex T01, T02, period1, period2;
  int NSample1, NSample2, degRough, degFine;
  Complex quasiCoe;
  bool isQuasi, isQuasiEnough;
  int i;
  double err1, err2;
  
  
  determine_quasi(type, isQuasi, isQuasiEnough, quasiCoe);
  //cout<<real(quasiCoe)<<" "<<imag(quasiCoe)<<endl;
  
  if (isQuasiEnough==false)
    {
      initialize(type, T01, T02, period1, period2, NSample1, NSample2);
      //cout<<real(T01)<<" "<<imag(T01)<<endl;
      //cout<<real(T02)<<" "<<imag(T02)<<endl;
      //cout<<real(period1)<<" "<<imag(period1)<<endl;
      //cout<<real(period2)<<" "<<imag(period2)<<endl;
      //cout<<NSample1<<" "<<NSample2<<endl;
      
      sample_cap1(type, NSample1, period1, T02, a, b, degRough);
      //     for (int i=0; i<a.size(); i++)
      //    	cout<<a[i]<<" "<<b[i]<<endl;
      
      sample_cap2(type, NSample2, period2, T02, a, b, a2, b2,degFine);
      // 	for (int i=0; i<a2.size(); i++)
      //     	cout<<a2[i]<<" "<<b2[i]<<endl;
    }
  
  //concatenating the exps and coes from the 2-level process
  //cout << k1 << endl;
  for (i=0; i<a.size();i++)
    {
      if (abs(b[i]) < 4*PI) {
	//cout << b[i] << endl;
	coes.push_back(a[i]);
	exps.push_back(b[i]);
      }
    } 
  for (i=0; i<a2.size(); i++)
    {
      if (abs(b2[i]) < 4*PI) {
	//cout << b2[i] << endl;
	coes.push_back(a2[i]);
	exps.push_back(b2[i]);
      }
    }

  //select_relevant_dipoles(coes, exps);
  if (isQuasi==true)
    {
      coes.insert(coes.begin(), quasiCoe);
      exps.insert(exps.begin(), ZERO);
    }
  
  //divide all b's by k0
  for (int i=0; i<exps.size(); i++) {
    exps[i]=exps[i]/k1;
    //cout<<exps[i]<<endl;
  }
  
  
}

void Dipoles::select_relevant_dipoles(vector<Complex>& coes, vector<Complex>& exps)
{
  Complex val, maxDis;
  if (coes.size()!=0)
  {
    for (int i=0; i<coes.size(); i++)
    {
      maxDis=sqrt(MAX_CONDUCTORSYS_RANGE[0]*MAX_CONDUCTORSYS_RANGE[0]+
                  MAX_CONDUCTORSYS_RANGE[1]*MAX_CONDUCTORSYS_RANGE[1]+
                  (MAX_CONDUCTORSYS_RANGE[2]+exps[i])*(MAX_CONDUCTORSYS_RANGE[2]+exps[i]));
      val=coes[i]*exp(maxDis);
      if (abs(val) < THRESHOLD)
	{
	  coes.erase(coes.begin()+i);
	  exps.erase(exps.begin()+i);
	}
    }
  }
}

void Dipoles::determine_quasi(dipoleType type, bool& isQuasi,
                              bool& isQuasiEnough, Complex& qs)
{
  
  int numSamples=400;
  double lambdaMax= 5*max(real(n),1.);
  double lambdaDiv=lambdaMax/numSamples;
  double tMax=sqrt(lambdaMax*lambdaMax-1);
  double threshold=1e-4;
  Complex t, u0, u1, lambda, remainder;

  qs=obtain_qs_contribution(type);
  isQuasiEnough=true;

  if (qs==ZERO)
  {
    isQuasi=false;
    isQuasiEnough=false;
  }
  else
    isQuasi=true;

  if (isQuasi==true)
  {
    for (int i=0; i<numSamples; i++)
    {
      t=i*lambdaDiv;
      // u0=ONE-IMAG*t;
      u0=-IMAG*t+(ONE-t/tMax);
      lambda=sqrt(ONE-u0*u0);
      u1=sqrt(n*n-lambda*lambda);
      if (imag(u1) > 0) 
	{
	  u1 = -u1;
	  cout << imag(u1) << endl;
	}
      remainder=obtain_kernel_contribution(type, u0,u1)-qs;
      if (abs(remainder)>threshold)
      {
        isQuasiEnough=false;
        break;
      }
    }
  }
}


void  Dipoles::initialize(dipoleType type, Complex& T01,
                          Complex& T02, Complex& period1, Complex& period2,
                          int& NSample1, int& NSample2)
{
  double lambdaMax;

  lambdaMax=2*max(real(n),1.);
  T02=sqrt(lambdaMax*lambdaMax-1);
  
  lambdaMax=300*max(real(n),1.);
  T01=sqrt(lambdaMax*lambdaMax-1)-T02;

  NSample2=400;
  NSample1=400;

  period2=Complex(real(T02)/NSample2, imag(T02)/NSample2);
  period1=Complex(real(T01)/NSample1, imag(T01)/NSample1);


}


void Dipoles::sample_cap1(dipoleType type, int NSample, Complex period,
                          Complex T02, vector<Complex>& a,
                          vector<Complex>& b, int& deg)
{
  Complex t, u0, u1, lambda;
  vector<Complex> y;

  for (int i=0; i<NSample; i++)
  {
    t=Complex((i)*real(period), (i)*imag(period));
    u0=-IMAG*(t+T02);
    //u0=(ONE-IMAG*t);
    lambda=sqrt(ONE-u0*u0);
    u1=sqrt(n*n-lambda*lambda);
    if (imag(u1) > 0) 
      {
	u1 = -u1;
	cout << imag(u1) << endl;
      }
    y.push_back(obtain_kernel_contribution(type, u0,u1)
                -obtain_qs_contribution(type));
  }

  gpof(y,period, a,b,deg);


  for (int i=0; i<b.size(); i++)
  {
    b[i]=b[i]/IMAG;
    a[i]=a[i]*exp(-IMAG*b[i]*T02);
  }
}

/*double generate_dipoles::compute_err_cap1(dipoleType type,
  Complex period, Complex T01, Complex T02,
  vector<Complex>& a, vector<Complex>& b)
  {
  return 0.0;
  }*/



void Dipoles::sample_cap2(dipoleType type, int NSample2, Complex period2,
                          Complex T02, const vector<Complex>& a,
                          const vector<Complex>& b, vector<Complex>& a2,
                          vector<Complex>& b2,int& deg)
{
  Complex t, u0, u1, lambda, sum;
  vector<Complex> y;


  for (int i=0; i<NSample2; i++)
  {
    t=Complex((i)*real(period2), (i)*imag(period2));
    //u0=ONE+t;
    u0=-IMAG*t+(ONE-t/T02);
    lambda=sqrt(ONE-u0*u0);
    u1=sqrt(n*n-lambda*lambda);
    if (imag(u1) > 0) 
      {
	u1 = -u1;
	cout << imag(u1) << endl;
      }
    sum=ZERO;
    for (int i=0; i<b.size(); i++)
      sum=sum+a[i]*exp(-b[i]*u0);

    y.push_back(obtain_kernel_contribution(type, u0,u1)
                -obtain_qs_contribution(type)-sum);
  }

  gpof(y, period2, a2,b2,deg);

  for (int i=0; i<b2.size(); i++)
  {
    b2[i]=b2[i]*T02/(IMAG*T02+ONE);
    a2[i]=a2[i]*exp(b2[i]);
  }

}


complex<double> Dipoles::obtain_qs_contribution(dipoleType type)
{

  switch(type)
  {
  case G_A_XX:
    return ZERO;
  case G_A_XZ:
    return ZERO;
  case G_A_ZZ:
    return (n*n-ONE)/(n*n+ONE);
  case G_PHI_H:
    return (ONE-n*n)/(n*n+ONE);
  case G_A_XX_MOD:
    return ZERO;
  }
}


complex<double> Dipoles::obtain_kernel_contribution(dipoleType type,
    Complex u0, Complex u1)
{
  switch(type)
  {
  case G_A_XX:
    return (u0-u1)/(u0+u1);
  case G_A_XZ:
    return (TWO*IMAG*u0*(ONE-n*n))/((u0+u1)*(n*n*u0+u1));
  case G_A_ZZ:
    return (n*n*u0-u1)/(n*n*u0+u1);
  case G_PHI_H:
    return (u0-u1)/(u0+u1)+TWO*u0*u0*(ONE-n*n)/((u0+u1)*(n*n*u0+u1));
  case G_A_XX_MOD:
   // return TWO*(n*n-ONE)/((u0+u1)*(n*n*u0+u1));
    return (u0-u1)/(n*n*u0+u1);
  }
}


void Dipoles::gpof(const vector<Complex>& y, Complex period, vector<Complex>& a,
                   vector<Complex>&b, int& deg)
{
  int N=y.size();
  int L=static_cast<int> (floor(((double)N)/2));
  itpp::cvec temp(N-L), temp2(N-L);
  itpp::cmat Y1(N-L,L), Y2(N-L,L);
  itpp::cmat u,v;
  itpp::vec s;
  itpp::cvec eigZ;
  int i,j;


  for (int colIndex=0; colIndex<L; colIndex++)
  {
    for (j=colIndex; j<colIndex+N-L; j++)
      temp[j-colIndex]=y[j];
    for (j=colIndex+1; j<colIndex+1+N-L; j++)
      temp2[j-colIndex-1]=y[j];

    Y1.set_col(colIndex, temp);
    Y2.set_col(colIndex, temp2);
  }

  itpp::svd(Y1,u,s,v);

  for (i=0; i< s.size();i++)
  {
	if ((s[i]/s[0])< 1e-3)
    {
      deg=i;
      break;
    }
  }
  deg=min(deg, MAXDEG);
  //deg=deg+1;
  cout<<deg<<endl;

  itpp::vec s_reduced(deg);
  itpp::cmat u_reduced(N-L,deg), v_reduced(N-L,deg);
  for (i=0; i<deg; i++)
  {
    s_reduced[i]=s[i];
    u_reduced.set_col(i,u.get_col(i));
    v_reduced.set_col(i,v.get_col(i));
  }

  itpp::Mat<double> t=diag(1.0/s_reduced);
  itpp::cmat t2=hermitian_transpose(u_reduced)*Y2*v_reduced;

  itpp::cmat zMat(deg,deg);
  for (i=0; i<deg; i++)
  {
    for (j=0; j<deg; j++)
    {
      zMat(i,j)=t(i,i)*t2(i,j);
    }
  }
  eigZ=eig(zMat);


  itpp::cvec sPoles(deg);
  for (i=0; i<deg; i++)
    sPoles(i)=log(eigZ(i))/period;

  itpp::cmat A(N,deg);
  for (i=0; i<N; i++)
  {
    for (j=0; j<deg; j++)
    {
      A(i,j)=pow(eigZ(j),i);
    }
  }

  itpp::cvec Y(y.size());
  for (i=0; i<y.size(); i++)
    Y[i]=y[i];

  itpp::cvec coes(deg);
  coes=ls_solve_od(A, Y);

  //converting coes and sPoles to std::vectors
  for (i=0; i<coes.size(); i++)
  {
    a.push_back(coes[i]);
    b.push_back(sPoles[i]);
  }

}
