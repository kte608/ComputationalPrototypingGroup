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
#ifndef DEFAULTS_H
#define DEFAULTS_H
#include <math.h>

/**
@author Xin Hu
*/

class Defaults
{
public:
  Defaults():isx_(0),isy_(0),isz_(0),ish_(0),isw_(0),issigma_(1),isnl_(1),
      isfnw_(1),isfnh_(1),ispnh_(1),ispnw_(1),isfrw_(1),isfrh_(1),ispfracw_(1),
	  ispfrach_(1),ispfracl_(1),isnep_(1),isepir_(0), isLayerSigma_(0), 
	  x_(0),y_(0),z_(0),h_(0),w_(0),sigma_(5.8e7),nl_(1),fnh_(1),fnw_(1),
	  pnh_(1),pnw_(1),frw_(2.),frh_(2.), pfracw_(0.1),pfrach_(0.1),
	  pfracl_(0.1), epir_(1), layerSigma_(0){};

  void set_isx(int val) {isx_=val;};
  void set_isy(int val) {isy_=val;};
  void set_isz(int val) {isz_=val;};
  void set_ish(int val) {ish_=val;};
  void set_isw(int val) {isw_=val;};
  void set_issigma(int val) {issigma_=val;};
  void set_isnl(int val) {isnl_=val;};
  void set_isfnh(int val) {isfnh_=val;};
  void set_isfnw(int val) {isfnw_=val;};
  void set_ispnh(int val) {ispnh_=val;};
  void set_ispnw(int val) {ispnw_=val;};
  void set_isfrw(int val) {isfrw_=val;};
  void set_isfrh(int val) {isfrh_=val;};
  void set_ispfracw(int val) {ispfracw_=val;};
  void set_ispfrach(int val) {ispfrach_=val;};
  void set_ispfracl(int val) {ispfracl_=val;};
  void set_isnep(int val) {isnep_=val;};
  void set_isepir(int val) {isepir_=val;};
  void set_isLayerSigma(int val) {isLayerSigma_=val;};
  
  void set_x(double val) {x_=val;};
  void set_y(double val) {y_=val;};
  void set_z(double val) {z_=val;};
  void set_h(double val) {h_=val;};
  void set_w(double val) {w_=val;};
  void set_sigma(double val) {sigma_=val;};
  void set_nl(int val) {nl_=val;};
  void set_fnh(int val) {fnh_=val;};
  void set_fnw(int val) {fnw_=val;};
  void set_pnh(int val) {pnh_=val;};
  void set_pnw(int val) {pnw_=val;};
  void set_frw(double val) {frw_=val;};
  void set_frh(double val) {frh_=val;};
  void set_pfracw(double val) {pfracw_=val;};
  void set_pfrach(double val) {pfrach_=val;};
  void set_pfracl(double val) {pfracl_=val;};
  void set_nep (int val) {nep_=val;};
  void set_epir(double val) {epir_=val;};
  void set_layerSigma(double val) {layerSigma_=val;};
  
  int isx(void) {return isx_;};
  int isy(void) {return isy_;};
  int isz(void) {return isz_;};
  int ish(void) {return ish_;};
  int isw(void) {return isw_;};
  int issigma(void) {return issigma_;};
  int isnl(void) {return isnl_;};
  int isfnh(void) {return isfnh_;};
  int isfnw(void) {return isfnw_;};
  int ispnh(void) {return ispnh_;};
  int ispnw(void) {return ispnw_;};
  int isfrh(void) {return isfrh_;};
  int isfrw(void) {return isfrw_;};
  int ispfrach(void) {return ispfrach_;};
  int ispfracw(void) {return ispfracw_;};
  int ispfracl(void) {return ispfracl_;};
  int isnep(void) {return isnep_;};
  int isepir(void) {return isepir_;};
  int isLayerSigma(void) {return isLayerSigma_;};

  
  double x(void) {return x_;};
  double y(void) {return y_;};
  double z(void) {return z_;};
  double h(void) {return h_;};
  double w(void) {return w_;};
  int nl(void) {return nl_;};
  int fnh(void) {return fnh_;};
  int fnw(void) {return fnw_;};
  int pnh(void) {return pnh_;};
  int pnw(void) {return pnw_;};
  double frw(void) {return frw_;};
  double frh(void) {return frh_;};
  double pfracw(void) {return pfracw_;};
  double pfrach(void) {return pfrach_;};
  double pfracl(void) {return pfracl_;};
  double sigma(void) {return sigma_;};
  int nep(void) {return nep_;};
  double epir(void) {return epir_;};
  double layerSigma(void){return layerSigma_;};


private:
  int isx_,
  isy_,
  isz_,
  ish_,
  isw_,
  issigma_,
  isnl_,
  isfnh_,
  isfnw_,
  ispnh_,
  ispnw_,
  isfrw_,           /* width ratio */
  isfrh_,           /* height ratio */
  ispfracw_,
  ispfrach_,
  ispfracl_,
  isnep_,
  isepir_,
  isLayerSigma_;


  double x_,y_,z_,h_,w_,sigma_;
  int nl_, fnh_, fnw_, pnh_, pnw_;
  double frw_, frh_;
  double pfracw_, pfrach_, pfracl_;
  int nep_;
  double epir_, layerSigma_;

};

#endif
