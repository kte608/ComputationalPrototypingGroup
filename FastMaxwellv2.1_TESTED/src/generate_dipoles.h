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
#ifndef GENERATE_DIPOLES_H
#define GENERATE_DIPOLES_H

#include "layerenv.h"

#include <complex>
#include <vector>


using namespace std;


enum dipoleType {G_A_XX, G_A_XZ, G_A_YY, G_A_YZ, G_A_ZZ, G_PHI_H, G_PHI_V, G_A_XX_MOD};

class Dipoles {

friend class Inductance;
friend class Charge;
typedef complex<double>  Complex;
public:
	Dipoles(void){};
	Dipoles (const LayerEnv& layers);
	void setup(const LayerEnv& layers);
	bool dipoles_flag (void)const {return dipoles_flag_;};	

	vector<Complex> coes_G_A_XX;
	vector<Complex> exps_G_A_XX;

	vector<Complex> coes_G_A_XX_MOD;
	vector<Complex> exps_G_A_XX_MOD;


	vector<Complex> coes_G_A_XZ;
	vector<Complex> exps_G_A_XZ;

	//vector<Complex> coes_G_A_YZ;
	//vector<Complex> exps_G_A_YZ;

	vector<Complex> coes_G_A_ZZ;
	vector<Complex> coes_G_A_ZZ_MOD;
	vector<Complex> exps_G_A_ZZ;

	vector<Complex> coes_G_PHI_H;
	vector<Complex> exps_G_PHI_H;
	double substrateHeight;

private:
	
	void approx_integral_2levels(dipoleType type, 
								 vector<Complex>& coes,  
								 vector<Complex>& exps);
	
	void determine_quasi(dipoleType type, bool& isQuasi, 
						 bool& isQuasiEnough, Complex& qs);
	
	void initialize(dipoleType type, Complex& T01, Complex& T02, 
					Complex& period1 , Complex& period2,
					int& NSample1, int& NSample2);
	
	void sample_cap1(dipoleType type, int NSample, Complex period,
					 Complex T02, vector<Complex>& a,
					 vector<Complex>& b, int& deg);
	
	void sample_cap2(dipoleType type, int NSample2, Complex period2,
					 Complex T02, const vector<Complex>& a, 
					 const vector<Complex>& b, vector<Complex>& a2, 
					 vector<Complex>& b2,int& deg);
	
	//double compute_err_cap1(dipoleType, Complex, Complex, Complex, 
	//						vector<Complex>& ,vector<Complex>& );
	
	Complex obtain_qs_contribution(dipoleType type);
	
	Complex obtain_kernel_contribution(dipoleType type, Complex u0, Complex u1); 
	
	void select_relevant_dipoles(vector<Complex>& coes, vector<Complex>& exps);
	
	void gpof(const vector<Complex>& y, Complex period, vector<Complex>& a,
			  vector<Complex>&b, int& deg);
	

	Complex n, k1, k2;
	//Complex IMAG, ZERO, ONE, TWO;
	double MAX_CONDUCTORSYS_RANGE[3];
	double THRESHOLD;
	
	
	bool dipoles_flag_;
	
	
};

#endif
