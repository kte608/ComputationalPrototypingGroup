/***************************************************************************
                          fullwaveVolumeGalerkin.h  -  description
                             -------------------
    begin                : Sat Feb 4 2006
    copyright            : (C) 2006 by root
    email                : tmoselhy@mit.edu
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef _FULLWAVE_GALERKIN_H_
#define _FULLWAVE_GALERKIN_H_

#include <vector>
#include <algorithm>
#include <stdexcept> //for exception handling
#include <utility> // for pair
#include "vector3D.h"
#include "kernelIntegration.h"
#include "../src/sysformation.h"
#include "../src/calcAEirkOverR.h"
#include "../src/calcaoneoverr.h"
#include "../src/filament.h"
#include "../src/generate_dipoles.h"



namespace pfft {

  template <class T1, class T2, class T3>

  class FullwaveVolumeGalerkin : public KernelIntegration {

  public:
    FullwaveVolumeGalerkin (void) {};
    FullwaveVolumeGalerkin (const std::vector<DifferentialOperator>&,
			    const std::vector<DifferentialOperator>&,
			    complex<double> KIn,
			    const bool& dipole_flag, const Dipoles& dipoles, 
			    const DyadicGreenFunction DyadicGreenFunc,
			    const SimulationType simuType=FULLWAVE, 
			    const int numGaussQuadraturePoint = 24);
    FullwaveVolumeGalerkin (const std::vector<DifferentialOperator>&,
			    const std::vector<DifferentialOperator>&);    

    const IntegrationScheme itegrationScheme(void) const { return GALERKIN; }
    
    void operator () (const T2& srcele, const T2& evalele);
    std::complex<double> result(const size_t integralIndex) const;

  private:

    
    bool isDipole_;
    
    DyadicGreenFunction DyadicGreenFunc_;
    SimulationType simuType_;
    const Dipoles *dipoles_;
    calcAEikrOverR<T1, T1> calcA;
    calcAEikrOverR<T3, T3> calcA_CX;
    calcAOneOverR calcA_static;
    complex<double> K_;
    enum ResultType { SINGLE_LAYER, D_DN_SINGLE_LAYER,
		      D_DX_SINGLE_LAYER, D_DY_SINGLE_LAYER, D_DZ_SINGLE_LAYER,
		      GRADIENT_SINGLE_LAYER,
		      DOUBLE_LAYER, D_DN_DOUBLE_LAYER};
    std::vector<ResultType> resultTypeList_;
    std::vector<ResultType> calcpResultTypeList_;

    std::complex<double> slp;
    std::complex<double> dlp;
    std::complex<double> dlp_dn;
    std::complex<double> slp_dn;
    pfft::vector3D<std::complex<double> > grad_slp;

    bool skipPanelIntegration_;

    void setupResultTypeList(void);
    ResultType mapDiffOperatorToResultType (const IntegralType integralType) const;

    bool IsWanted(ResultType rt) ;
    bool IsWantedFromCalcp(ResultType rt) ;
    
    void get_cx_fil(const filament<double> f,
		    const complex<double> cx_location,
		    filament<complex<double> >& cx_fil, bool isReal);
    
  };
  
  
  
  template <class T1, class T2, class T3>
    FullwaveVolumeGalerkin<T1, T2, T3>::
    FullwaveVolumeGalerkin (
			    const std::vector<DifferentialOperator>& outerOperatorList,
			    const std::vector<DifferentialOperator>& innerOperatorList,
			    complex<double> KIn,
			    const bool& dipole_flag, const Dipoles& dipoles, 
			    const DyadicGreenFunction DyadicGreenFunc,
			    const SimulationType simuType, const int numGaussQuadraturePoint)
    : KernelIntegration(outerOperatorList, innerOperatorList), K_(KIn),
    isDipole_(dipole_flag),  DyadicGreenFunc_(DyadicGreenFunc), simuType_(simuType)
    {
      if (isDipole_)
	dipoles_=&(dipoles);
      setupResultTypeList();
      calcA.init(K_, true, true);
      if (simuType_==EMQS) {
	calcA_CX.init(0., true, true);
      } else {
	calcA_CX.init(K_, true, true);
      }
      skipPanelIntegration_ = false;
      
    }

  template <class T1, class T2, class T3>
    FullwaveVolumeGalerkin<T1, T2, T3>::
    FullwaveVolumeGalerkin (
			    const std::vector<DifferentialOperator>& outerOperatorList,
			    const std::vector<DifferentialOperator>& innerOperatorList)
    : KernelIntegration(outerOperatorList, innerOperatorList)
    {
      
      setupResultTypeList();
      skipPanelIntegration_ = false;
      
    }
  
  
  /**********************************************************************
   * result --
   **********************************************************************/
  template <class T1, class T2, class T3>
  std::complex<double>
  FullwaveVolumeGalerkin<T1, T2, T3>::result (
					  const size_t integralIndex) const
  {
    switch (resultTypeList_[integralIndex]) {
    case SINGLE_LAYER:
      return slp;
      break;
    case D_DN_SINGLE_LAYER:
      return slp_dn;
      break;
    case D_DX_SINGLE_LAYER:
      return grad_slp.x();
      break;
    case D_DY_SINGLE_LAYER:
      return grad_slp.y();
      break;
    case D_DZ_SINGLE_LAYER:
      return grad_slp.z();
      break;
    case DOUBLE_LAYER:
      return dlp;
      break;
    case D_DN_DOUBLE_LAYER:
      return dlp_dn;
      break;
    default:
      throw std::domain_error("fullwaveCollocation.h: unknown resultType");
      break;
    }
  }

  /**********************************************************************
   * setupResultTypeList --
   **********************************************************************/
  template <class T1, class T2, class T3>
  void
  FullwaveVolumeGalerkin<T1, T2, T3>::setupResultTypeList (
						       void)
  {
    for (size_t i = 0; i < numIntegral(); i++) {
      resultTypeList_.push_back(mapDiffOperatorToResultType(integralType(i)));
    }

    calcpResultTypeList_ = resultTypeList_;
    for (size_t i = 0; i < numIntegral(); i++) {
      switch (calcpResultTypeList_[i]) {
      case D_DN_SINGLE_LAYER:
	calcpResultTypeList_[i] = GRADIENT_SINGLE_LAYER;
	break;
      case D_DX_SINGLE_LAYER:
	calcpResultTypeList_[i] = GRADIENT_SINGLE_LAYER;
	break;
      case D_DY_SINGLE_LAYER:
	calcpResultTypeList_[i] = GRADIENT_SINGLE_LAYER;
	break;
      case D_DZ_SINGLE_LAYER:
	calcpResultTypeList_[i] = GRADIENT_SINGLE_LAYER;
	break;
      }
    }
    sort(calcpResultTypeList_.begin(), calcpResultTypeList_.end());
    calcpResultTypeList_.erase(unique(calcpResultTypeList_.begin(),
				      calcpResultTypeList_.end()),
			       calcpResultTypeList_.end());
  }
  
  /**********************************************************************
   * mapDiffOperatorToResultType --
   **********************************************************************/
  template <class T1, class T2, class T3>
  typename FullwaveVolumeGalerkin<T1, T2, T3>::ResultType
  FullwaveVolumeGalerkin<T1, T2, T3>::
  mapDiffOperatorToResultType (
			       const IntegralType integralType) const
  {
    ResultType resultType;

    switch (integralType.first) {
    case NONE:
      switch(integralType.second) {
      case NONE:
	resultType = SINGLE_LAYER;
	break;
      case D_DN:
	resultType = DOUBLE_LAYER;
	break;
      case D_DX:
      case D_DY:
      case D_DZ:
	throw std::domain_error("fullwaveCollocation.h: can not handle this integral type");
      default:
	throw std::domain_error("fullwaveCollocation.h: unknown inner deferential operator");
      }
      break;
    case D_DN:
      switch(integralType.second) {
      case NONE:
	resultType = D_DN_SINGLE_LAYER;
	break;
      case D_DN:
	resultType = D_DN_DOUBLE_LAYER;
	break;
      case D_DX:
      case D_DY:
      case D_DZ:
	//	resultType = HYPER_SINGULAR;
	throw std::domain_error("fullwaveCollocation.h: can not handle this integral type");
      default:
	throw std::domain_error("fullwaveCollocation.h: unknown inner deferential operator");
      }
      break;
    case D_DX:
      switch(integralType.second) {
      case NONE:
	resultType = D_DX_SINGLE_LAYER;
	break;
      case D_DN:
      case D_DX:
	resultType = SINGLE_LAYER;
	break;
      case D_DY:
	resultType = SINGLE_LAYER;
	break;
      case D_DZ:
	throw std::domain_error("fullwaveCollocation.h: can not handle this integral type");
      default:
	throw std::domain_error("fullwaveCollocation.h: unknown inner deferential operator");
      }
      break;
    case D_DY:
      switch(integralType.second) {
      case NONE:
	resultType = D_DY_SINGLE_LAYER;
	break;
      case D_DN:
      case D_DX:
	resultType = SINGLE_LAYER;
	break;
      case D_DY:
	resultType = SINGLE_LAYER;
	break;
      case D_DZ:
	throw std::domain_error("fullwaveCollocation.h: can not handle this integral type");
      default:
	throw std::domain_error("fullwaveCollocation.h: unknown inner deferential operator");
      }
      break;
    case D_DZ:
      switch(integralType.second) {
      case NONE:
	resultType = D_DZ_SINGLE_LAYER;
	break;
      case D_DN:
      case D_DX:
      case D_DY:
      case D_DZ:
	throw std::domain_error("fullwaveCollocation.h: can not handle this integral type");
      default:
	throw std::domain_error("fullwaveCollocation.h: unknown inner deferential operator");
      }
      break;
    default:
      throw std::domain_error("fullwaveCollocation.h: unknown OUTER deferential operator");
      break;
    }
    return resultType;
  }

  /**********************************************************************
   * operator () --
   **********************************************************************/
  template <class T1, class T2, class T3>
  void
  FullwaveVolumeGalerkin<T1, T2, T3>::operator () (const T2& srcele, const T2& evalele)
    {
      complex<double> slptmp=0., dlptmp=0., dumy=0.;
      filament<double> srcfil, evalfil;
      srcfil = filament<double>(srcele.getfil().loc0(), srcele.getfil().loc1(), srcele.getfil().filWidth(),
				srcele.getfil().filHeight(), srcele.getfil().sigma());
      evalfil = filament<double>(evalele.getfil().loc0(), evalele.getfil().loc1(), 
				 evalele.getfil().filWidth(),
				 evalele.getfil().filHeight(), evalele.getfil().sigma());
      //calcp(srcPanel, evalPanel.centroid(), slp, dlp, grad_slp);
      //if (simuType_==EMQS) {
      slp = 0.;
      
      if ((srcele.getfil().l_dir().x()!=0) && (evalele.getfil().l_dir().x()!=0)) {
	filament<complex<double> > cx_fil, cx_fil2;
	if (length(srcele.getfil().filCentroid()-evalele.getfil().filCentroid())==0.)
	  slp += calcA_static.selfterm(&srcfil);
	else
	  slp += abs(calcA_static.mutual(&srcfil, &evalfil));
	if (isDipole_==true) {
	  const int numDipole=dipoles_->coes_G_A_XX.size();
	  get_cx_fil(evalele.getfil(), 0., cx_fil2, true);
	  for (int i=0; i<numDipole; i++) {
	      get_cx_fil(evalele.getfil(), dipoles_->exps_G_A_XX[i], cx_fil, false);
	      slptmp=0.;	
	      dlptmp=0.;
	      slptmp=calcA_static.mutual(&srcfil, &cx_fil);
	      if(real(slptmp)<0.) 
		slptmp *= -1.;
	      slp=slp+dipoles_->coes_G_A_XX[i]*slptmp;
	  }
	}
	slp *= (srcele.getfil().l_dir().x() * evalele.getfil().l_dir().x());
      }
      if ((srcele.getfil().l_dir().y()!=0) && (evalele.getfil().l_dir().y()!=0)) {
	filament<complex<double> > cx_fil, cx_fil2;
	if (length(srcele.getfil().filCentroid()-evalele.getfil().filCentroid())==0.)
	  slp += calcA_static.selfterm(&srcfil);
	else
	  slp += abs(calcA_static.mutual(&srcfil, &evalfil));
	if (isDipole_==true) {
	  const int numDipole=dipoles_->coes_G_A_XX.size();
	  get_cx_fil(evalele.getfil(), 0., cx_fil2, true);
	  for (int i=0; i<numDipole; i++) {
	    get_cx_fil(evalele.getfil(), dipoles_->exps_G_A_XX[i], cx_fil, false);
	    slptmp = 0.;	
	    dlptmp = 0.;
	    slptmp = calcA_static.mutual(&srcfil, &cx_fil);
	    if(real(slptmp)<0.) 
	      slptmp*=-1.;
	    slp = slp + dipoles_->coes_G_A_XX[i]*slptmp;
	  }
	}
	slp = slp * (srcele.getfil().l_dir().y() * evalele.getfil().l_dir().y());
      }
      if ((srcele.getfil().l_dir().x() !=0) && (evalele.getfil().l_dir().z()!=0)) {
	if (isDipole_==true) {
	  const int numDipole=dipoles_->coes_G_A_XZ.size();
	  filament<complex<double> > cx_fil, cx_fil2;
	  get_cx_fil(evalele.getfil(), 0., cx_fil2, true);
	  dlp = 0.;
	  for (int i=0; i<numDipole; i++)
	    {
	      get_cx_fil(srcele.getfil(), dipoles_->exps_G_A_XZ[i], cx_fil, false);
	      slptmp=0.;	dlptmp=0.;
	      calcA_CX(cx_fil2, cx_fil, slptmp, dlptmp, slptmp);
	      dlp=dlp+dipoles_->coes_G_A_XZ[i]*dlptmp;
	    }
	  }
	slp=slp+dlp * (srcele.getfil().l_dir().x() * evalele.getfil().l_dir().z());
      }
      if ((srcele.getfil().l_dir().z()!=0) && (evalele.getfil().l_dir().x()!=0)) {
	if (isDipole_==true) {
	  const int numDipole=dipoles_->coes_G_A_XZ.size();
	  filament<complex<double> > cx_fil, cx_fil2;
	  get_cx_fil(srcele.getfil(), 0., cx_fil2, true);
	  dlp = 0.;
	  for (int i=0; i<numDipole; i++)
	    {
	      get_cx_fil(evalele.getfil(), dipoles_->exps_G_A_XZ[i], cx_fil, false);
	      slptmp=0.;	dlptmp=0.;
	      calcA_CX(cx_fil2, cx_fil, slptmp, dlptmp, slptmp);
	      dlp=dlp+dipoles_->coes_G_A_XZ[i]*dlptmp;
	    }
	}
	slp=slp+dlp * (srcele.getfil().l_dir().z() * evalele.getfil().l_dir().x());
      } 
      if ((srcele.getfil().l_dir().y()!=0) && (evalele.getfil().l_dir().z()!=0)) {
	if (isDipole_==true) {
	  const int numDipole=dipoles_->coes_G_A_XZ.size();
	  filament<complex<double> > cx_fil, cx_fil2;
	  get_cx_fil(evalele.getfil(), 0., cx_fil2, true);
	  dlp = 0.;
	  for (int i=0; i<numDipole; i++)
	    {
	      get_cx_fil(srcele.getfil(), dipoles_->exps_G_A_XZ[i], cx_fil, false);
	      slptmp=0.;	dlptmp=0.;
	      calcA_CX(cx_fil2, cx_fil, slptmp, dlptmp, slptmp);
	      dlp=dlp + dipoles_->coes_G_A_XZ[i]*dlptmp;
	    }
	}
	slp=slp + dlp * (srcele.getfil().l_dir().y() * evalele.getfil().l_dir().z());
      }
      if ((srcele.getfil().l_dir().z()!=0) && (evalele.getfil().l_dir().y()!=0)) {
	if (isDipole_==true) {
	  const int numDipole=dipoles_->coes_G_A_XZ.size();
	  filament<complex<double> > cx_fil, cx_fil2;
	  get_cx_fil(srcele.getfil(), 0., cx_fil2, true);
	  dlp = 0.;
	  for (int i=0; i<numDipole; i++)
	    {
	      get_cx_fil(evalele.getfil(), dipoles_->exps_G_A_XZ[i], cx_fil, false);
	      slptmp=0.;	dlptmp=0.;
	      calcA_CX(cx_fil2, cx_fil, slptmp, slptmp, dlptmp);
	      dlp=dlp+dipoles_->coes_G_A_XZ[i]*dlptmp;
	    }
	}
	slp=slp+dlp * (srcele.getfil().l_dir().y() * evalele.getfil().l_dir().z());
      } 
      if ((srcele.getfil().l_dir().z()!=0) && (evalele.getfil().l_dir().z()!=0)) {
	complex<double> slpzz = 0.;
	//} else if (DyadicGreenFunc_==GZZ) {
	filament<complex<double> > cx_fil, cx_fil2;
	if (length(srcele.getfil().filCentroid()-evalele.getfil().filCentroid())==0.)
	  slpzz = calcA_static.selfterm(&srcfil);
	else
	  slpzz = calcA_static.mutual(&srcfil, &evalfil);
	if (isDipole_==true) {
	  get_cx_fil(evalele.getfil(), 0., cx_fil2, true);
	  const int numDipole=dipoles_->coes_G_A_ZZ.size();
	  for (int i=0; i<numDipole; i++) {
	    get_cx_fil(srcele.getfil(), dipoles_->exps_G_A_ZZ[i], cx_fil, false);
	    slptmp=0.;	dlptmp=0.;
	    slpzz = calcA_static.mutual(&srcfil, &cx_fil2);
	    slpzz = slpzz + dipoles_->coes_G_A_ZZ[i]*slptmp;
	  }
	}
	slp += slpzz * (srcele.getfil().l_dir().z() * evalele.getfil().l_dir().z());
      }
      slp /= 1e-7;
      dlp /= 1e-7;
      //cout << slp << " " << dlp << endl;
      //cout << srcele.getfil().l_dir().x() << srcele.getfil().l_dir().y() << srcele.getfil().l_dir().z() << evalele.getfil().l_dir().x() << evalele.getfil().l_dir().y() << evalele.getfil().l_dir().z() << endl;
    }
  
  
  /**********************************************************************
   * IsWanted --
   **********************************************************************/
  template <class T1, class T2, class T3>
  bool
  FullwaveVolumeGalerkin<T1, T2, T3>::IsWanted (
					    ResultType rt)
  {
    return (find(resultTypeList_.begin(), resultTypeList_.end(), rt)
	    != resultTypeList_.end());
  }

  /**********************************************************************
   * IsWantedFromCalcp --
   **********************************************************************/
  template <class T1, class T2, class T3>
  bool
  FullwaveVolumeGalerkin<T1, T2, T3>::IsWantedFromCalcp (
						     ResultType rt)
  {
    return (find(calcpResultTypeList_.begin(), calcpResultTypeList_.end(), rt)
	    != calcpResultTypeList_.end());
  }

/*************************************************************************/
/*************************************************************************/
  template <class T1, class T2, class T3>
    void
    FullwaveVolumeGalerkin<T1, T2, T3>::get_cx_fil(
						   const filament<double> f,
						   const complex<double> cx_location,
						   filament<complex<double> >& cx_fil, bool isReal)
    {
      if (isReal==false) {
	double subHeight=dipoles_->substrateHeight;	
	complex<double> x_val(f.loc0().x(),0);
	complex<double> y_val(f.loc0().y(),0);
	point_vec_operations::point3D<complex<double> > 
	  endPt0 (x_val, y_val,
		  (2*subHeight-f.loc0().z())+IMAG*cx_location );
	complex<double> x_val_pt2(f.loc1().x(),0);
	complex<double> y_val_pt2(f.loc1().y(),0);
	point_vec_operations::point3D<complex<double> >  
	  endPt1(x_val_pt2,y_val_pt2,
		 (2*subHeight-f.loc1().z())+IMAG*cx_location);
	double width= f.filWidth();
	double height=f.filHeight();
	double sigma=f.sigma();
	cx_fil.initialize(endPt0, endPt1, width, height, sigma);
      } else {
	complex<double> x_val(f.loc0().x(),0);
	complex<double> y_val(f.loc0().y(),0);
	complex<double> z_val(f.loc0().z(),0);
	point_vec_operations::point3D<complex<double> > endPt0 (x_val, y_val,z_val);
	complex<double> x_val_pt2(f.loc1().x(),0);
	complex<double> y_val_pt2(f.loc1().y(),0);
	complex<double> z_val_pt2(f.loc1().z(),0);
	point_vec_operations::point3D<complex<double> >  endPt1(x_val_pt2,y_val_pt2,z_val_pt2) ;
	double width= f.filWidth();
	double height=f.filHeight();
	double sigma=f.sigma();
	cx_fil.initialize(endPt0, endPt1, width, height, sigma);
      }
    }
  
  
} //namespace pfft

#endif
