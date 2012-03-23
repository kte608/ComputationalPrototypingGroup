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
#ifndef EXTRACTZ_H
#define EXTRACTZ_H
#include <complex>
#include <vector>
#include <iostream>

using namespace std;
/**
@author Xin Hu
*/
class extractZ
{
public:
  extractZ(void){};
  extractZ(const vector<complex<double> >& Ib,
           const vector<complex<double> >& Im,
           const complex<double>& adm,
           const complex<double>& Rsource,
           const double& freq);
  complex<double> Z(void)const {return Z_;};
  complex<double> Y(void)const {return Y_;};
  double R(void) const {return R_;};
  double X(void) const {return X_;};
  double C(void) const {return C_;};
  double L(void) const {return L_;};
  double Q(void) const {return Q_;};
private:
  double RSource_, ISource_;
  complex<double> Z_, Y_;
  double R_, X_,G_,B_, C_, L_, Q_;

};

#endif
