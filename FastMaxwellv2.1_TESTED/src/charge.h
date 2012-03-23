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
#ifndef CHARGE_H
#define CHARGE_H

#include "generate_dipoles.h"
//#include "spRowMat.h"
#include "../pfft/spRowMat.h"
//#include "henry.h"
#include "structure.h"
#include "calcpForEikrOverR.h"


class Charge
{
public:
  Charge(void):isEMQS(true), isDipole(false), dipoles(NULL), condStruc(NULL){};
  void init(const bool& EMQS_flag, const bool& dipole_flag,Structure& structure,
            const double epi);
  void fill_p(pfft::SpRowMat<complex<double> >& PMat,
              const int i,const int j,
              const panel<double>* p1, const panel<double>* p2);
  void Myfill_p(pfft::SpRowMat<complex<double> >& PMat,
              const int i,const int j,
              const panel<double>* p1, const panel<double>* p2);
  void set_dipoles(const Dipoles& dipoles_, const double freq);
private:
  //charge *build_charges_list();
  complex<double> get_phi_hor_dipoles(const panel<double>* p, const point3D<double>& eval,
					const vector<complex<double> >& coes,
					const vector<complex<double> >& exps,
					const int i, const int j);
  complex<double> get_phi_vert_dipoles(const panel<double>* p, const point3D<double>& eval,
						const vector<complex<double> >& coes,
						const vector<complex<double> >& exps,
						const int i, const int j);
  complex<double> Myget_phi_vert_dipoles(const panel<double>* p, const panel<double>* peval,
						const vector<complex<double> >& coes,
						const vector<complex<double> >& exps,
						const int i, const int j);

  void get_cx_panel(const panel<double>* p, 
							const complex<double> cx_location, 
							panel<complex<double> >& cx_panel);
  void get_r_panel(const panel<double>* p,
							const complex<double> cx_location,
							panel<double> & cx_panel);

  bool isEMQS;
  bool isDipole;
  const Dipoles *dipoles;
  Structure* condStruc;
  double epi0, K_num;
  
  calcpForEikrOverR<double, double > calcp_static;
  calcpForEikrOverR<complex<double>, double > calcp_static_cx;
  calcpForEikrOverR<double, double > calcp_fullWave;
  calcpForEikrOverR<complex<double>, double > calcp_fullWave_cx;
    
};

#endif
