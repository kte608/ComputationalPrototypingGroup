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
#include "calcaoneoverr.h"


bool calcAOneOverR::less_equal_than (complex<double> A, double B)
{
  return A.real()<=B;
}
bool calcAOneOverR::less_than (complex<double> A, complex<double> B)
{
  return A.real()<B.real();
}
bool calcAOneOverR::less_equal_than (double A, double B)
{
  return A<=B;
}
bool calcAOneOverR::greater_equal_than (complex<double> A, double B)
{
  return A.real()>=B;
}
bool calcAOneOverR::greater_equal_than (double A, double B)
{
  return A>=B;
}
bool calcAOneOverR::less_than (complex<double> A, double B)
{
  return A.real()<B;
}
bool calcAOneOverR::less_than (double A, double B)
{
  return A<B;
}
bool calcAOneOverR::less_equal_than (complex<double> A, complex<double> B)
{
  return real(A)<=real(B);
}
bool calcAOneOverR::greater_than (complex<double> A, double B)
{
  return A.real()>B;
}
bool calcAOneOverR::greater_than (complex<double> A, complex<double> B)
{
  return real(A)>real(B);
}
bool calcAOneOverR::greater_than (double A, double B)
{
  return A>B;
}
complex<double> calcAOneOverR::atanh_my(complex<double> val)
{
  if ((ONE-val)==ZERO)
  {
    cout<<"can't perfrom atanh"<<endl;
    exit(1);
  }
  complex<double> var=(ONE+val)/(ONE-val);

  return (ONE/TWO)*log(var);
}
double calcAOneOverR::atanh_my(double val)
{
  return atanh(val);
}
complex<double> calcAOneOverR::asinh_my(complex<double> val)
{
  return log(val+sqrt(val*val+ONE));
}
double calcAOneOverR::asinh_my(double val)
{
  return asinh(val);
}
complex<double> calcAOneOverR::atan2_my(complex<double> val1,complex<double> val2 )
{
  return atan2(real(val1),real(val2));
}
double calcAOneOverR::atan2_my(double val1, double val2)
{
  return atan2(val1, val2);
}
complex<double> calcAOneOverR::atan_my(complex<double> val)
{
  return IMAG/2.*log((IMAG+val)/(IMAG-val));
}
double calcAOneOverR::atan_my(double val)
{
  return atan(val);
}
complex<double> calcAOneOverR::MAX_my(complex<double> A, complex<double> B)
{
  return( abs(A) > abs(B) ? (A) : (B) );
}

double calcAOneOverR:: MAX_my(double A, double B)
{
  return MAX(A,B);
}


complex<double> calcAOneOverR::abs_mydef(complex<double> val)
{
  if (real(val)<0.)
    return -val;
  else
    return val;
}


double calcAOneOverR::abs_mydef(double val)
{
  return abs(val);
}
