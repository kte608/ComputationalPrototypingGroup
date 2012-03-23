/***************************************************************************
                          FullwaveSurfaceGalerkin.h  -  description
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

#ifndef _SufaceFW_GALERKIN_H_
#define _SufaceFW_GALERKIN_H_

#include <vector>
#include <algorithm>
#include <stdexcept> //for exception handling
#include <utility> // for pair
#include "vector3D.h"
#include "kernelIntegration.h"
#include "../src/calcpForEikrOverR.h"
#include "../src/panel.h"
#include "../src/generate_dipoles.h"
//#include "calcpForEikrOverR.h"

namespace pfft {


  template <class T1, class T2, class T3>

  class FullwaveSurfaceGalerkin : public KernelIntegration {

  public:
    FullwaveSurfaceGalerkin (void) {};
    
    FullwaveSurfaceGalerkin (const std::vector<DifferentialOperator>&,
			     const std::vector<DifferentialOperator>&,
			     complex<double> KIn,
			     const bool& dipole_flag, const Dipoles& dipoles, 
			     const DyadicGreenFunction DyadicGreenFunc,
			     const SimulationType simuType=EMQS, const int numGaussQuadraturePoint = 24);
    FullwaveSurfaceGalerkin (const std::vector<DifferentialOperator>&,
			     const std::vector<DifferentialOperator>&);

    const IntegrationScheme itegrationScheme(void) const { return GALERKIN; }    
    void operator () (const T2& srcele, const T2& evalele);
    std::complex<double> result(const size_t integralIndex) const;
    
  private:
    
    bool isDipole_;
    DyadicGreenFunction DyadicGreenFunc_;
    SimulationType simuType_;
    const Dipoles *dipoles_;
    calcpForEikrOverR<T1, T1> calcP;
    calcpForEikrOverR<T3, T1> calcP_CX;
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
    void get_cx_panel(const panel<double> p,
		      const complex<double> cx_location,
		      panel<complex<double> >& cx_panel, bool isReal);

  };


  template <class T1, class T2, class T3>
  FullwaveSurfaceGalerkin<T1, T2, T3>::
  FullwaveSurfaceGalerkin (
			   const std::vector<DifferentialOperator>& outerOperatorList,
			   const std::vector<DifferentialOperator>& innerOperatorList,
			   complex<double> KIn,
			   const bool& dipole_flag, const Dipoles& dipoles, 
			   const DyadicGreenFunction DyadicGreenFunc,
			   const SimulationType simuType, const int numGaussQuadraturePoint)
    : KernelIntegration(outerOperatorList, innerOperatorList),
    K_(KIn),
    isDipole_(dipole_flag), DyadicGreenFunc_(DyadicGreenFunc), simuType_(simuType)
    {
      if (isDipole_)
	dipoles_=&(dipoles);
      setupResultTypeList();
      
      calcP.init(0., false,true);
      calcP_CX.init(0., false,true);
      /*
      if (simuType_==EMQS) {
	calcP.init(0., false,true);
	calcP_CX.init(0., false,true);
      }
      else {
	calcP.init(K_, false,true);
	calcP_CX.init(K_, false,true);
      }
      */
      
      skipPanelIntegration_ = false;
      
    }
  


  template <class T1, class T2, class T3>
  FullwaveSurfaceGalerkin<T1, T2, T3>::
  FullwaveSurfaceGalerkin (
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
    FullwaveSurfaceGalerkin<T1, T2, T3>::result (
						 const size_t integralIndex) const
    {
      return slp;
    }
  

  /**********************************************************************
   * setupResultTypeList --
   **********************************************************************/
  template <class T1, class T2, class T3>
  void
  FullwaveSurfaceGalerkin<T1, T2, T3>::setupResultTypeList (
						       void)
  {
  }

  /**********************************************************************
   * mapDiffOperatorToResultType --
   **********************************************************************/
  template <class T1, class T2, class T3>
  typename FullwaveSurfaceGalerkin<T1, T2, T3>::ResultType
  FullwaveSurfaceGalerkin<T1, T2, T3>::
  mapDiffOperatorToResultType (
			       const IntegralType integralType) const
  {
  }

  /**********************************************************************
   * operator () --
   **********************************************************************/
  template <class T1, class T2, class T3>
    void
    FullwaveSurfaceGalerkin<T1, T2, T3>::operator () (const T2& srcele, const T2& evalele) {
    complex<double> slptmp, dlptmp;
    
    panel<complex<double> > cx_panel;
    size_t numPtsW=1;
    size_t numPtsH=1;
    vector<point_vec_operations::point3D<double > > OUTER_quad_pts;
    vector<double > OUTER_quad_weights;
    OUTER_quad_pts.resize(numPtsW*numPtsH);
    OUTER_quad_weights.resize(numPtsW*numPtsH);
    evalele.getpanel().panel_quadrature(numPtsW * numPtsH, OUTER_quad_pts, OUTER_quad_weights);
    complex<double> dlp_R=0., dlp_CX=0.;
    for (int j=0; j<OUTER_quad_pts.size(); j++) {
      dlp=0.;
      calcP(srcele.getpanel(), OUTER_quad_pts[j], slp, dlp, slp);
      dlp_R=dlp_R+dlp*OUTER_quad_weights[j];
    }
    dlp_R=dlp_R/srcele.getpanel().area()/evalele.getpanel().area();
    slp=dlp=dlp_R;
    if (isDipole_==true) {
      const int numDipole=dipoles_->coes_G_PHI_H.size();
      for (int i=0; i<numDipole; i++)
	{
	  get_cx_panel(srcele.getpanel(), dipoles_->exps_G_PHI_H[i], cx_panel, false);
	  dlp=0.;
	  for (int j=0; j<OUTER_quad_pts.size(); j++) {
	    dlptmp=0.; slptmp=0.;
	    calcP_CX(cx_panel, OUTER_quad_pts[j], slptmp, dlptmp, slptmp);
	    dlp=dlp+dlptmp*OUTER_quad_weights[j];
	  }
	  dlp_CX=dlp_CX+dipoles_->coes_G_PHI_H[i]*dlp/srcele.getpanel().area()/evalele.getpanel().area();
	}
      dlp=dlp_R+dlp_CX;
      slp=dlp;
      //cout << slp << dlp_R << dlp_CX << srcele.getpanel().area() << evalele.getpanel().area() << endl;
      //if (dlp_R==(0./0.,0.) | dlp_R == (1./0.,0.) | dlp_R == (-1./0.,0.))
      //   cout << dlp << " " << dlp_R <<  endl;
    }
}
  
  
  /**********************************************************************
   * IsWanted --
   **********************************************************************/
  template <class T1, class T2, class T3>
  bool
  FullwaveSurfaceGalerkin<T1, T2, T3>::IsWanted (
						 ResultType rt)
    {
    }
  
  /**********************************************************************
   * IsWantedFromCalcp --
   **********************************************************************/
  template <class T1, class T2, class T3>
  bool
  FullwaveSurfaceGalerkin<T1, T2, T3>::IsWantedFromCalcp (
						     ResultType rt)
  {
  }

/*************************************************************************/
/////////////////////////
  template <class T1, class T2, class T3>
  void
  FullwaveSurfaceGalerkin<T1, T2, T3>::get_cx_panel(const panel<double> p,
                          const complex<double> cx_location,
                          panel<complex<double> >& cx_panel, bool isReal)
{
  
  if (isReal==true) {
    
    vector<point_vec_operations::point3D<complex<double> > > cx_verts;
    cx_verts.resize(p.shape());
    for (int i=0; i<p.shape(); i++) {
      cx_verts[i]=point_vec_operations::point3D<complex<double> >
	(p.vertex(i).x(), p.vertex(i).y(), p.vertex(i).z());
    }
    cx_panel.initialize(cx_verts);
    
  } else {
    double subHeight=dipoles_->substrateHeight;
    vector<point_vec_operations::point3D<complex<double> > > cx_verts;
    cx_verts.resize(p.shape());
    for (int i=0; i<p.shape(); i++)
      {
	cx_verts[i]=point_vec_operations::point3D<complex<double> >
	  (p.vertex(i).x(), p.vertex(i).y(),
	   (2*subHeight-p.vertex(i).z())+IMAG*cx_location );
      }
    cx_panel.initialize(cx_verts);
  }
}
  
  
  
} //namespace pfft

#endif
