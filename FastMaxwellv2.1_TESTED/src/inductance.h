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
#ifndef INDUCTANCE_H
#define INDUCTANCE_H

//#include "spRowMat.h"
#include "../pfft/spRowMat.h"
#include "filament.h"
#include "generate_dipoles.h"
//#include "henry.h"
#include "calcaoneoverr.h"
#include "calcAEirkOverR.h"
#include "calcpForEikrOverR.h"

#include <complex>

using namespace std;

class Inductance
{
public:
  Inductance(void):isEMQS(true), isDipole(false), dipoles(NULL) {};
  void init (const bool& EMQS_flag, const bool& dipole_flag,const double epi);
  void fill_L(pfft::SpRowMat<complex<double> >& LMat,
              const int i, const int j,
              const filament<double> *f1,
              const filament<double> *f2);

  void Myfill_L(pfft::SpRowMat<complex<double> >& LMat,
              const int i, const int j,
              const filament<double> *f1,
              const filament<double> *f2);

  void Myfill_L_REDUCED(pfft::SpRowMat<complex<double> >& LMat,
              const int i, const int j,
              const filament<double> *f1,
              const filament<double> *f2);

  void set_dipoles(const Dipoles& dipoles_, const double freq);

private:
  complex<double> get_selfterm_dipoles( const filament<double> *f,
                                        const vector<complex<double> >& coes,
                                        const vector<complex<double> >& exps);
  complex<double> get_mutual_dipoles( const filament<double> *f1,
                                      const filament<double>* f2,
                                      const vector<complex<double> >& coes,
                                      const vector<complex<double> >& exps);
  void get_mutual_dipoles_asym(const filament<double> *f1,
                               const filament<double>* f2,
                               const vector<complex<double> >& coes,
                               const vector<complex<double> >& exps,
                               complex<double>& dAdx, complex<double>& dAdy);

  complex<double> Myget_selfterm_dipoles( const filament<double> *f,
                                        const vector<complex<double> >& coes,
                                        const vector<complex<double> >& exps,
                                        const vector<complex<double> >& coes2,
                                        const vector<complex<double> >& exps2);
  complex<double> Myget_mutual_dipoles( const filament<double> *f1,
                                      const filament<double>* f2,
                                      const vector<complex<double> >& coes,
                                      const vector<complex<double> >& exps,
                                      const vector<complex<double> >& coes2,
                                      const vector<complex<double> >& exps2);
  complex<double> Myget_mutual_dipoles_asym(const filament<double> *f1,
                               const filament<double>* f2,
                               const vector<complex<double> >& coes,
                               const vector<complex<double> >& exps);

  void get_cx_fil(const filament<double> *f,
                  const complex<double> cx_location,
                  filament<complex<double> >& cx_fil);
  void copy_fil(const filament<double> *f,
                  const complex<double> cx_location,
                  filament<complex<double> >& cx_fil);

  bool isEMQS;
  bool isDipole;
  const Dipoles *dipoles;
  double epi0, K_num;

  //static computations
  //for quasi-static computations
  calcAOneOverR calcA_static;
  //for cx filaments w/ asym computation (dAdx or dAdy) due to substrate
  calcAEikrOverR<complex<double>, complex<double> > calcA_X;

  //full wave computations
  //for real filaments
  calcAEikrOverR<double, double > calcA_fullWave;  
  //for cx filaments due to substrate
  //calcAEikrOverR<complex<double>, complex<double> > calcA_fullWave_cx;
  //calcAEikrOverR<complex<double>, complex<double> > calcA_fullWave_cx_XX;


};



#endif
