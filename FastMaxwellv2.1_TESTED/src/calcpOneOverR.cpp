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
#include "calcpOneOverR.h"


calcpOneOverR::calcpOneOverR(void)
{
  maxorderMom = -1;
  maxorderSphere = -1;
  printf("in initcalcp: Mom Spher %d %d\n",maxorderMom,maxorderSphere);
}

double calcpOneOverR::max_my(double A, double B)
{
  return MAX(A,B);
}
double calcpOneOverR::min_my(double A, double B)
{
  return ((A<B)? (A):(B));
}
complex<double> calcpOneOverR::max_my(complex<double> A, complex<double> B)
{
  return( abs(A) > abs(B) ? (A) : (B) );
}
complex<double> calcpOneOverR::min_my(complex<double> A, complex<double> B)
{
  return( abs(A) < abs(B) ? (A) : (B) );
}

bool calcpOneOverR::less_than (complex<double> A, double B)
{
  return A.real()<B;
}
bool calcpOneOverR::less_than (complex<double> A, complex<double> B)
{
  return A.real()<B.real();
}

double calcpOneOverR::abs_my(double A)
{
  return abs(A);
}
complex<double> calcpOneOverR::abs_my(complex<double> A)
{
  return complex<double>(abs(A.real()), abs(A.imag()));
	//return abs(A);
}

bool calcpOneOverR::less_than (double A, double B)
{
  return A<B;
}
bool calcpOneOverR::less_equal_than (complex<double> A, double B)
{
  return A.real()<=B;
}
bool calcpOneOverR::less_equal_than (double A, double B)
{
  return A<=B;
}
complex<double> calcpOneOverR::atan2_my(complex<double> val1,complex<double> val2 )
{	
	return atan2(real(val1),real(val2));
}
double calcpOneOverR::atan2_my(double val1, double val2)
{
	return atan2(val1, val2);
}

/* Create Binomial Array. */
double** calcpOneOverR::createBinom(int order)
{
  int m, n;
  double **Binom;

  Binom=(double**)malloc((order+1)*sizeof(double*));
  for(m=0; m <= order; m++)
    Binom[m]=(double*)malloc( (order+1)*sizeof(double));


  for(m=0; m <= order; m++)
  {
    Binom[m][0] = 1.0;
    for(n = 1; n <= order; n++)
    {
      Binom[m][n] = ((m - (n-1))  * Binom[m][n-1]) / n;
    }
  }
  return Binom;
}


double** calcpOneOverR::createPmn(int order)
{
  int n, m;
  double pmm;
  double **Pmn;

  Pmn=(double**)malloc( (order+1)*sizeof(double*));
  for(n=0; n <= order; n++)
    Pmn[n]=(double*)malloc( (order+2)*sizeof(double));

  pmm = 1.0;

  for(m=0; m <= order; m++)
  {
    Pmn[m][m] = pmm;
    Pmn[m][m+1] = 0.0;
    for(n = (m+2); n <= order; n++)
    {
      Pmn[m][n] = -((double) (m + n - 1)) / ((double) (n - m)) * Pmn[m][n-2];
    }
    pmm *= -(2.0 * m + 1.0);
  }
  return(Pmn);
}

/* Create Factorial Array. */
int* calcpOneOverR::createFactorial(int size)
{
  int i, *Fact;

  Fact=(int*)malloc(size*sizeof(int));

  Fact[0] = (int)1.0;
  for(i=1; i <= size; i++)
  {
    Fact[i] = i * Fact[i-1];
  }
  return(Fact);
}


int calcpOneOverR::multerms(int order)
{
  return(costerms(order) + sinterms(order));
}

/*
  returns number of cos(m*phi)-weighted terms in the real (not cmpx) multi exp
*/
int calcpOneOverR::costerms(int order)
{
  return(((order + 1) * (order + 2)) / 2);
}

/*
  returns number of sin(m*phi)-weighted terms in the real (not cmpx) multi exp
*/
int calcpOneOverR::sinterms(int order)
{
  return((((order + 1) * (order + 2)) / 2) - (order+1));
}

int calcpOneOverR::index(int n, int m)
{
  if(m > n)
  {
    fprintf(stderr, "FLE-index: m = %d > n = %d\n", m, n);
    exit(0);
  }
  if(n < 0 || m < 0)
  {
    fprintf(stderr, "FLE-index: n = %d or m = %d negative\n", n, m);
    exit(0);
  }
  return(m + (n*(n+1))/2);
}

/*
  gives the linear index into vector from n and m used by all routines dealing
  with moments (sine parts), e.g. Mn^m 
  assumes an array with all m = 0 (Mn^0) entries omitted to save space
  assumed entry order: (n,m) = (1,1) (2,1) (2,2) (3,1) (3,2) (3,3) (4,1)...
*/
int calcpOneOverR::sindex(int n, int m, int cterms)
{
  if(m > n)
  {
    fprintf(stderr, "FLE-sindex: m = %d > n = %d\n", m, n);
    exit(0);
  }
  if(n < 0 || m < 0)
  {
    fprintf(stderr, "FLE-sindex: n = %d or m = %d negative\n", n, m);
    exit(0);
  }
  if(m == 0)
  {
    fprintf(stderr, "FLE-sindex: attempt to index M%d^0\n", n);
    exit(0);
  }
  return(cterms + m + (n*(n+1))/2 - (n+1));
}
