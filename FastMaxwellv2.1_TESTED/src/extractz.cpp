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
#include "extractz.h"
#include "const_def.h"



extractZ::extractZ(const std::vector<std::complex<double> >& Ib,
                   const std::vector<std::complex<double> >& Im,
                   const std::complex<double>& adm,
                   const std::complex<double>& Rsource,
                   const double& freq)
{
  ISource_=1.;
  RSource_=Rsource.real();
  double Vs=RSource_*ISource_;
  Z_=Vs/adm-RSource_;
  Y_=1./Z_;

  R_=real(Z_);
  X_=imag(Z_);
  L_=X_/2/PI/freq;

  G_=real(Y_);
  B_=imag(Y_);
  C_=B_/2/PI/freq;
  Q_=X_/R_;

}




