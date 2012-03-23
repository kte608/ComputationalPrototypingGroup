
/*-*-c++-*-*/
/******************************************************************************
 **                                                                          **
 ** Copyright (c) 2000 Massachusetts Institute of Technology, Cambridge, MA. **
 ** All rights reserved.                                                     **
 **                                                                          **
 ******************************************************************************
  ==========================================================================
  ============================= File Info ==================================

  Author: Zhenhai Zhu
  
  Description:

  Resources:

  See also:

  const static char cvsid[] = "$Id: eikrOverR.h,v 1.1.1.1 2002/03/12 20:08:56 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _EIKR_OVER_R_H_
#define _EIKR_OVER_R_H_

#include "vector3D.h"
#include <complex>
#include "kernelIntegration.h"


namespace pfft {

  template <class KT = std::complex<double> >

  class EikrOverRIMAGES {

  public:
    EikrOverRIMAGES(KT KIn) : K_(KIn) {}
    EikrOverRIMAGES(KT KIn, const Dipoles& dipoles,
		    const DyadicGreenFunction DyadicGreenFunc)
      : K_(KIn), dipoles_(dipoles), DyadicGreenFunc_(DyadicGreenFunc){}
    EikrOverRIMAGES(void) : K_(0) {}
    
    void set_shift(double dZ) { 
      AbsoluteShift_Z_ = dZ; 
    }
    void set_maxthreshold(double maxThreshold) {
      maxthreshold=maxThreshold;
    }
    
    std::complex<double> operator () (const double x, const double y, const double z) const
    {
      vector3D<double> r(x, y, z);
      if (length(r) <= maxthreshold) {
	//      if (length(r) == 0.) {
	return 0.;
      } else {
	return exp(std::complex<double>(0., 1.) * (-K_) * length(r)) / length(r);
      }
    }
    
    
    
    /**************************************
   ****************************************/
    std::complex<double> images (const double x, const double y, const double z) const
    {
      std::complex<double> value=std::complex<double>(0., 0.) ;
      std::complex<double> zTotal=std::complex<double>(0., 0.) ;
      vector3D<complex<double> > r(x, y, z);
      switch (DyadicGreenFunc_) {
      case GXXc:
	for (int i=0; i<dipoles_.exps_G_A_XX.size(); i++) {
	  zTotal=z-IMAG*dipoles_.exps_G_A_XX[i]+2.*AbsoluteShift_Z_ - 2.*dipoles_.substrateHeight;
	  vector3D<complex<double> > r(x, y, zTotal);
	  if (length(r) <= maxthreshold) {
	    //  if (length(r) == 0.) {
	    value = 0.;
	    return 0.;
	  } else {
	    value=value
	      +(dipoles_.coes_G_A_XX[i])*exp(std::complex<double>(0., 1.)*(-K_)*length(r))/length(r);
	  }
	}
	return   value;
	break;
      case GXZc:
	for (int i=0; i<dipoles_.exps_G_A_XZ.size(); i++) {
	  zTotal=z-IMAG*dipoles_.exps_G_A_XZ[i]+2.*AbsoluteShift_Z_
	    -2.*dipoles_.substrateHeight;
	  vector3D<complex<double> > r(x, y, zTotal);
	  if (length(r) == 0.) {
	    value=value+0.0;
	  } else {
	    value=value
	      +(dipoles_.coes_G_A_XZ[i])*exp(std::complex<double>(0., 1.)*(-K_)*length(r))/length(r);
	  }
	}
	return   value;
	break;
      case GZZ_phi:
	for (int i=0; i<dipoles_.exps_G_A_ZZ.size(); i++) {
	  zTotal=z-IMAG*dipoles_.exps_G_A_ZZ[i]+2.*AbsoluteShift_Z_
	    -2.*dipoles_.substrateHeight;
	  vector3D<complex<double> > r(x, y, zTotal);
	  if (length(r) == 0.) {
	    value=value+0.0;
	  } else {
	    value=value
	      +(dipoles_.coes_G_A_ZZ_MOD[i])*exp(std::complex<double>(0., 1.)*(-K_)*length(r))/length(r);
	  }
	}
	return   value;
	break;
      case GScalar:
	for (int i=0; i<dipoles_.exps_G_PHI_H.size(); i++) {
	  zTotal=z-IMAG*dipoles_.exps_G_PHI_H[i]+2.*AbsoluteShift_Z_
	    -2.*dipoles_.substrateHeight;
	  vector3D<complex<double> > r(x, y, zTotal);
	  if (length(r) == 0.) {
	    value=value+0.0;
	  } else {
	    value=value
	      +(dipoles_.coes_G_PHI_H[i])*exp(std::complex<double>(0., 1.)*(-K_)*length(r))/length(r);
	  }
	}
	return   value;
	break;
      case GZZ:
	for (int i=0; i<dipoles_.exps_G_A_ZZ.size(); i++) {
	  zTotal=z-IMAG*dipoles_.exps_G_A_ZZ[i]+2.*AbsoluteShift_Z_
	    -2.*dipoles_.substrateHeight;
	  vector3D<complex<double> > r(x, y, zTotal);
	  if (length(r) == 0.) {
	    value=value+0.0;
	  } else {
	    value=value
	      +(dipoles_.coes_G_A_ZZ[i])*exp(std::complex<double>(0., 1.)*(-K_)*length(r))/length(r);
	  }
	}
	return   value;
	break;
	
	
    default:
      errorMessage("eikrOverR_IMAGES.h::operator()",
		   "ERROR MUST CONTINUE REST OF CODE");
      }
    }
    /*******************************************
********************************************************************/
    

    std::complex<double> gxx (const double x, const double y, const double z) const
    {
      std::complex<double> value=std::complex<double>(0., 0.) ;
      std::complex<double> zTotal=std::complex<double>(0., 0.) ;
      vector3D<complex<double> > r(x, y, z);
      
      for (int i=0; i<dipoles_.exps_G_A_XX.size(); i++) {
	zTotal=z-IMAG*dipoles_.exps_G_A_XX[i]+2.*AbsoluteShift_Z_ - 2.*dipoles_.substrateHeight;
	vector3D<complex<double> > r(x, y, zTotal);
	if (length(r) <= maxthreshold) {
	  //  if (length(r) == 0.) {
	  value += 0.;
	  //return 0.;
	} else {
	  value=value
	    +(dipoles_.coes_G_A_XX[i])*exp(std::complex<double>(0., 1.)*(-K_)*length(r))/length(r);
	}
      }
      return   value;
    }
    
    
    std::complex<double> gxz (const double x, const double y, const double z) const
    {
      std::complex<double> value=std::complex<double>(0., 0.) ;
      std::complex<double> zTotal=std::complex<double>(0., 0.) ;
      vector3D<complex<double> > r(x, y, z);
      
      for (int i=0; i<dipoles_.exps_G_A_XZ.size(); i++) {
	zTotal=z-IMAG*dipoles_.exps_G_A_XZ[i]+2.*AbsoluteShift_Z_
	  -2.*dipoles_.substrateHeight;
	vector3D<complex<double> > r(x, y, zTotal);
	if (length(r) <= maxthreshold) {
	  value = value+0.0;
	  //return value;
	} else {
	  std::complex<double> ikR;
	  ikR = std::complex<double>(0., 1.)*(K_)*length(r);
	  value=value
	    +(dipoles_.coes_G_A_XZ[i])*exp(-ikR)*(-ikR-1.)*x/(length(r)*length(r)*length(r));
	}
      }
      return   value;
    }
    

    std::complex<double> gyz (const double x, const double y, const double z) const
    {
      std::complex<double> value=std::complex<double>(0., 0.) ;
      std::complex<double> zTotal=std::complex<double>(0., 0.) ;
      vector3D<complex<double> > r(x, y, z);
      
      for (int i=0; i<dipoles_.exps_G_A_XZ.size(); i++) {
	zTotal=z-IMAG*dipoles_.exps_G_A_XZ[i]+2.*AbsoluteShift_Z_
	  -2.*dipoles_.substrateHeight;
	vector3D<complex<double> > r(x, y, zTotal);
	if (length(r) <= maxthreshold) {
	  value = value+0.0;
	  //return value;
	} else {
	  std::complex<double > ikR; 
	  ikR = std::complex<double>(0., 1.)*(K_)*length(r);
	  value=value
	    +(dipoles_.coes_G_A_XZ[i])*exp(-ikR)*(-ikR-1.)*y/(length(r)*length(r)*length(r));
	}
      }
      return   value;
}
    
    std::complex<double> gphi_V (const double x, const double y, const double z) const
    {
      std::complex<double> value=std::complex<double>(0., 0.) ;
      std::complex<double> zTotal=std::complex<double>(0., 0.) ;
      vector3D<complex<double> > r(x, y, z);
      for (int i=0; i<dipoles_.exps_G_A_ZZ.size(); i++) {
	zTotal=z-IMAG*dipoles_.exps_G_A_ZZ[i]+2.*AbsoluteShift_Z_
	  -2.*dipoles_.substrateHeight;
	vector3D<complex<double> > r(x, y, zTotal);
	if (length(r) <= maxthreshold) {
	  value = value+0.0;
	  //return value;
	} else {
	  value=value
	    +(dipoles_.coes_G_A_ZZ_MOD[i])*exp(std::complex<double>(0., 1.)*(-K_)*length(r))/length(r);
	}
      }
      return   value;
    }

    std::complex<double> gphi_H (const double x, const double y, const double z) const
    {
      std::complex<double> value=std::complex<double>(0., 0.) ;
      std::complex<double> zTotal=std::complex<double>(0., 0.) ;
      vector3D<complex<double> > r(x, y, z);
      for (int i=0; i<dipoles_.exps_G_PHI_H.size(); i++) {
	zTotal=z-IMAG*dipoles_.exps_G_PHI_H[i]+2.*AbsoluteShift_Z_
	  -2.*dipoles_.substrateHeight;
	vector3D<complex<double> > r(x, y, zTotal);
	if (length(r) <= maxthreshold) {
	  value = value+0.0;
	  //return value;
	} else {
	  value=value
	    +(dipoles_.coes_G_PHI_H[i])*exp(std::complex<double>(0., 1.)*(-K_)*length(r))/length(r);
	}
      }
      return   value;
    }
    
    
    std::complex<double> gzz (const double x, const double y, const double z) const
    {
      std::complex<double> value=std::complex<double>(0., 0.) ;
      std::complex<double> zTotal=std::complex<double>(0., 0.) ;
      vector3D<complex<double> > r(x, y, z);
      for (int i=0; i<dipoles_.exps_G_A_ZZ.size(); i++) {
	zTotal=z-IMAG*dipoles_.exps_G_A_ZZ[i]+2.*AbsoluteShift_Z_
	  -2.*dipoles_.substrateHeight;
	vector3D<complex<double> > r(x, y, zTotal);
	if (length(r) <= maxthreshold) {
	  value = value+0.0;
	  //return value;
	} else {
	  value=value
	    +(dipoles_.coes_G_A_ZZ[i])*exp(std::complex<double>(0., 1.)*(-K_)*length(r))/length(r);
	}
      }
      return   value;
    }
    
    /*******************************************
********************************************************************/
    
    
    std::complex<double> all (const double x, const double y, const double zD,  const double zS) const
    {
      vector3D<complex<double> > r(x, y, zD);
      std::complex<double> value=std::complex<double>(0., 0.) ;
      std::complex<double> zTotal=std::complex<double>(0., 0.) ;
      switch (DyadicGreenFunc_) {
      case GXXc:
	if (length(r) <= maxthreshold) {
	  //  if (length(r) == 0.) {
	  value=0.;
	  return 0.;
	} else {
	  value=exp(std::complex<double>(0., 1.) * (-K_) * length(r)) / length(r);
	}
	for (int i=0; i<dipoles_.exps_G_A_XX.size(); i++) {
	  zTotal=zS-IMAG*dipoles_.exps_G_A_XX[i]+2.*AbsoluteShift_Z_
	    -2.*dipoles_.substrateHeight;
	  vector3D<complex<double> > r1(x, y, zTotal);
	  if (length(r1) == 0.) {
	    value=value+0.0;
	  } else {
	    value=value
	      +(dipoles_.coes_G_A_XX[i])*exp(std::complex<double>(0., 1.)*(-K_)*length(r1))/length(r1);
	  }
	}
	return   value;
	break;
      case GXZc:
	for (int i=0; i<dipoles_.exps_G_A_XZ.size(); i++) {
	  zTotal=zS-IMAG*dipoles_.exps_G_A_XZ[i]+2.*AbsoluteShift_Z_
	    -2.*dipoles_.substrateHeight;
	  //zTotal=zS-IMAG*dipoles_.exps_G_A_XZ[i];
	  vector3D<complex<double> > r(x, y, zTotal);
	  if (length(r) == 0) {
	    value=value+0.0;
	  } else {
	    value=value
	      +(dipoles_.coes_G_A_XZ[i])*exp(std::complex<double>(0., 1.)*(-K_)*length(r))/length(r);
	  }
	}
	return   value;
	break;
      case GXXa:
	//      if (length(r) == 0.) {
	//					 value=value+0.0;
	//      } else {
	//	value=exp(std::complex<double>(0., 1.) * (-K_) * length(r)) / length(r);
	//      }
	value=0.;
	for (int i=0; i<dipoles_.exps_G_A_XX_MOD.size(); i++) {
	  zTotal=zS-IMAG*dipoles_.exps_G_A_XX_MOD[i]+2.*AbsoluteShift_Z_
	    -2.*dipoles_.substrateHeight;
	  vector3D<complex<double> > r1(x, y, zTotal);
	  if (length(r1) == 0.) {
	    value=value+0.0;
	  } else {
	    value=value
	      +(dipoles_.coes_G_A_XX_MOD[i])*exp(std::complex<double>(0., 1.)*(-K_)*length(r1))/length(r1);
	  }
	}
	return   value;
	break;
      case GZZ_phi:
	if (length(r) == 0.) {
	  value=value+0.0;
	} else {
	  value=exp(std::complex<double>(0., 1.) * (-K_) * length(r)) / length(r);
	}
	for (int i=0; i<dipoles_.exps_G_A_ZZ.size(); i++) {
	  zTotal=zS-IMAG*dipoles_.exps_G_A_ZZ[i]+2.*AbsoluteShift_Z_
	    -2.*dipoles_.substrateHeight;
	  vector3D<complex<double> > r1(x, y, zTotal);
	  if (length(r1) == 0.) {
	    value=value+0.0;
	  } else {
	    value=value
	      +(dipoles_.coes_G_A_ZZ_MOD[i])*exp(std::complex<double>(0., 1.)*(-K_)*length(r1))/length(r1);
	  }
	}
	return   value;
	break;
      case GScalar:
	if (length(r) == 0.) {
	  value=value+0.0;
	} else {
	  value=exp(std::complex<double>(0., 1.) * (-K_) * length(r)) / length(r);
	}
	for (int i=0; i<dipoles_.exps_G_PHI_H.size(); i++) {
	  zTotal=zS-IMAG*dipoles_.exps_G_PHI_H[i]+2.*AbsoluteShift_Z_
	    -2.*dipoles_.substrateHeight;
	  vector3D<complex<double> > r1(x, y, zTotal);
	  if (length(r1) == 0.) {
	    value=value+0.0;
	  } else {
	    value=value
	      +(dipoles_.coes_G_PHI_H[i])*exp(std::complex<double>(0., 1.)*(-K_)*length(r1))/length(r1);
	  }
	}
	return   value;
	break;
      case GZZ:
	if (length(r) == 0.) {
	  value=value+0.0;
	} else {
	  value=exp(std::complex<double>(0., 1.) * (-K_) * length(r)) / length(r);
	}
	for (int i=0; i<dipoles_.exps_G_A_ZZ.size(); i++) {
	  zTotal=zS-IMAG*dipoles_.exps_G_A_ZZ[i]+2.*AbsoluteShift_Z_
	    -2.*dipoles_.substrateHeight;
	  vector3D<complex<double> > r1(x, y, zTotal);
	  if (length(r1) == 0.) {
	    value=value+0.0;
	  } else {
	    value=value
	      +(dipoles_.coes_G_A_ZZ[i])*exp(std::complex<double>(0., 1.)*(-K_)*length(r1))/length(r1);
	  }
	}
	return   value;
	break;
	
      default:
	errorMessage("eikrOverR_IMAGES.h::operator()",
		     "ERROR MUST CONTINUE REST OF CODE");
      }
    }
    /*******************************************
********************************************************************/
    
    
//
//    std::complex<double> operator () (const vector3D<double>& r) const {
//      if (length(r) == 0) {
//	return 0.0;
//      } else {
//	return  exp(std::complex<double>(0., 1.) * (-K_) * length(r)) / length(r);
//      }
//    }

    const KT KK(void) const { return K_; }
    
  private:
    KT K_; // k as in exp(ikR)/R
    Dipoles dipoles_;
    double AbsoluteShift_Z_;
    DyadicGreenFunction DyadicGreenFunc_;
    double maxthreshold;
  };
  
} //namespace pfft

#endif
