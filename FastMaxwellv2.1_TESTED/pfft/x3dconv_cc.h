/*-*-c++-*-*/
/******************************************************************************
 **                                                                          **
 ** Copyright (c) 2000 Massachusetts Institute of Technology, Cambridge, MA. **
 ** All rights reserved.                                                     **
 **                                                                          **
 ******************************************************************************
  ==========================================================================
  ============================= File Info ==================================

  Author: Ben Song
  
  Description:

  Resources:

  See also:

  const static char cvsid[] = "$Id: x3dconv_cc.h,v 1.7 2003/02/11 03:09:24 zhzhu Exp $";

  ==========================================================================
*/
#ifndef _X3DCONV_CC_H_
#define _X3DCONV_CC_H_

#include "x3dconv.h"
#include "TArray.h"
#include "rfft3d.h"

namespace pfft{
  
  template<class T3>
  class Fast3DConv<Cplx, Cplx, T3>{

  public:
    inline Fast3DConv(): isKernelSetup_(false), isKernelTransformed_(false),
			 isGridDataSetup_(false), isCalcKernelSetup_(false),
			 calcKernelPtr_(static_cast<T3*>(0)) {}

    inline void initialize(const size_t, const size_t, const size_t,
			   const double, const double, const double, const T3&);

    // test module:
    template<class Vector>
    void setupTestGridData(Vector&, bool randomOrNot = true);
  
    template<class VectorA, class VectorB>
    void getDirectPotential(const VectorA, 
			    VectorB&);

    template<class VectorA, class VectorB>
    bool compare(const VectorA&, const VectorB&);
    
    // pfft module:
    template<class KernelContainer>
    void setupKernel(KernelContainer&);

    template<class KernelContainer>
    void setupScalarImages(KernelContainer&);
    template<class KernelContainer>
    void setupKernelDiag(KernelContainer&, KernelContainer&);
    template<class KernelContainer>
    void setupKernelOffDiag(KernelContainer&, KernelContainer&);

    template<class KernelContainer>
    void fftOnKernel(KernelContainer&);

    template<class GridDataContainer>
  	void PermuteOnData(GridDataContainer&, const GridDataContainer&);


    template<class Vector, class GridDataContainer>
    void setupGridData(const Vector&, GridDataContainer&);
     
    template<class GridDataContainer, class Vector>
    void extractFastPotential(const GridDataContainer&, Vector&);

    template<class KernelContainer, class GridDataContainer>
    void operator() (const KernelContainer&,
		     GridDataContainer&, 
		     GridDataContainer&);

    template<class KernelContainer, class GridDataContainer>
    void operator() (const KernelContainer&, const KernelContainer&,
		     const KernelContainer&, const KernelContainer&,
		     GridDataContainer&, GridDataContainer&, 
		     GridDataContainer&, GridDataContainer&, 
		     GridDataContainer&, GridDataContainer&);

    template<class KernelContainer, class GridDataContainer>
    void operator() (const KernelContainer&,
		     GridDataContainer&, GridDataContainer&, 
		     GridDataContainer&, GridDataContainer&, 
		     GridDataContainer&, GridDataContainer&);
		       
  protected:
    // state variables
    size_t nx_, ny_, nz_;
    double dx_, dy_, dz_;
        
    bool isKernelSetup_, isKernelTransformed_;
    bool isGridDataSetup_, isCalcKernelSetup_; 

    const T3 *calcKernelPtr_;

    // working variables
    size_t nx2_, ny2_, nz2_; // 2 means double 
    size_t dim_[3], dim2_[3];
    double pace_[3];
    size_t length1D_, length1D8_, length1D4_;
  };

 
  /**********************************************************************
   * Fast3DConv<Cplx, Cplx, T3>::initialize --
   **********************************************************************/
  template<class T3>
  void Fast3DConv<Cplx, Cplx, T3>::initialize(const size_t n1,
					      const size_t n2,
					      const size_t n3,
					      const double d1,
					      const double d2,
					      const double d3,
					      const T3 &calcKernel)
  {
    nx_ = n1; ny_ = n2; nz_ = n3;
    dx_ = d1; dy_ = d2; dz_ = d3;
    
    pace_[0] = d1; pace_[1] = d2; pace_[2] = d3;

    nx2_ = nx_ * 2; ny2_ = ny_ * 2; nz2_ = nz_ * 2; 
    dim_[0] = nx_; dim_[1] = ny_; dim_[2] = nz_;
    dim2_[0] = nx2_; dim2_[1] = ny2_; dim2_[2] = nz2_;
    length1D_ = nx_ * ny_ * nz_;
    length1D8_ = nx2_ * ny2_ * nz2_;
    length1D4_ = length1D8_ / 2;

    calcKernelPtr_ = &calcKernel;
    isCalcKernelSetup_ = true;

  }
  
  
  /**********************************************************************
   * Fast3DConv<Cplx, Cplx, T3>::SetupKernel --
   **********************************************************************/
  template<class T3>
  template<class KernelContainer>
  void Fast3DConv<Cplx, Cplx, T3>
  ::setupKernel (KernelContainer& x)
  {
    assert(x.size() == length1D8_);
    pfft::setupKernelTopelitz(dim2_, pace_, x, *calcKernelPtr_);
    isKernelSetup_ = true;
  }

  /**********************************************************************
   * Fast3DConv<Cplx, Cplx, T3>::SetupKernelH --
   **********************************************************************/
  template<class T3>
  template<class KernelContainer>
  void Fast3DConv<Cplx, Cplx, T3>
  ::setupScalarImages (KernelContainer& x)
  {
    assert(x.size() == length1D8_);
    pfft::setupKernelHankel(dim2_, pace_, x, *calcKernelPtr_);
    isKernelSetup_ = true;
  }

  template<class T3>
  template<class KernelContainer>
  void Fast3DConv<Cplx, Cplx, T3>
  ::setupKernelDiag (KernelContainer& x, KernelContainer& z)
  {
    assert(x.size() == length1D8_);
    assert(z.size() == length1D8_);
    pfft::setupKerneldiag(dim2_, pace_, x, z, *calcKernelPtr_);
    isKernelSetup_ = true;
  }

  template<class T3>
  template<class KernelContainer>
  void Fast3DConv<Cplx, Cplx, T3>
  ::setupKernelOffDiag (KernelContainer& x, KernelContainer& z)
  {
    assert(x.size() == length1D8_);
    assert(z.size() == length1D8_);
    pfft::setupKerneloffdiag(dim2_, pace_, x, z, *calcKernelPtr_);
    isKernelSetup_ = true;
  }

  /**********************************************************************
   * Fast3DConv<Cplx, Cplx, T3>::Fast3DConv --
   **********************************************************************/
  template<class T3>
  template<class KernelContainer>
  void Fast3DConv<Cplx, Cplx, T3>
  ::fftOnKernel(KernelContainer& x)
  {
    assert(isKernelSetup_ == true);
   
    ::fftwnd_plan plan = ::fftw3d_create_plan(dim2_[0], dim2_[1], dim2_[2],
				FFTW_BACKWARD,
				FFTW_ESTIMATE|FFTW_IN_PLACE);

    // this always cause problem in our Itanium server. it return a null point for plan .
    // I don't know why. 
    // ::fftwnd_plan plan = ::fftwnd_create_plan(3, reinterpret_cast<int*>(dim2_), FFTW_BACKWARD, FFTW_ESTIMATE|FFTW_IN_PLACE);

    ::fftwnd_one(plan, x.p_, static_cast<fftw_complex*>(0));
    fftwnd_destroy_plan(plan);

    isKernelTransformed_ = true;
  } 

  /**********************************************************************
   * Fast3DConv<Cplx, Cplx, T3>::Fast3DConv --
   **********************************************************************/
  template<class T3>
  template<class Vector, class GridDataContainer>
  void Fast3DConv<Cplx, Cplx, T3>
  ::setupGridData(const Vector& gridData,
		  GridDataContainer& gridDataPadded)
  {
    assert(gridData.size() == length1D_);
    assert(gridDataPadded.size() == length1D4_);
    
    for(size_t ii=0; ii<nx_; ++ii){
      for(size_t jj=0; jj<ny2_; ++jj){
	for(size_t kk=0; kk<nz2_; ++kk){
	  int paddedIndex = kk + nz2_ * (jj + ny2_ * ii);
	  if(jj<ny_ && kk<nz_){
	    int index = kk + nz_ * (jj + ny_ * ii);
	    gridDataPadded[paddedIndex].re = gridData[index].real();
	    gridDataPadded[paddedIndex].im = gridData[index].imag();
	  } else {
	    gridDataPadded[paddedIndex].re = 0.;
	    gridDataPadded[paddedIndex].im = 0.;
	  }
	}
      }
    } // the grid data is zero padded, all these zeroes should
    // be ensured at its definition !!

    isGridDataSetup_ = true;
  }
 
  /**********************************************************************
   * Fast3DConv<Cplx, Cplx, T3>::Fast3DConv --
   **********************************************************************/
  template<class T3>
  template<class GridDataContainer, class Vector>
  void Fast3DConv<Cplx, Cplx, T3>
  ::extractFastPotential(const GridDataContainer& gridDataPadded,
			 Vector& gridData)
  {
    assert(gridData.size() == length1D_); 
    pfft::extractFastPotential(dim2_, gridDataPadded, gridData);
  }

  /**********************************************************************
   * Permutation --
   **********************************************************************/
  template<class T3>
  template<class GridDataContainer>
  void Fast3DConv<Cplx, Cplx, T3>
  ::PermuteOnData(GridDataContainer& P_gdata, const GridDataContainer& gdata)
  {
    
    /*   
    for(size_t kk=0; kk< gdata.size(); ++kk){
      P_gdata[kk] = 0.* gdata[kk];
      }*/
    
    for(size_t ii=0; ii<nx_; ++ii){
      for(size_t jj=0; jj<ny_; ++jj){
	for(size_t kk=0; kk<nz_; ++kk){
	  int index2 = kk + nz2_*(jj + ny2_*ii);
	  P_gdata[index2] = gdata[nz_-1-kk+nz2_*(jj + ny2_*ii)];
	}
      }
    }
    // copy(result.begin(), result.end(), Presult.begin());
  }
  
  /**********************************************************************
   * Fast3DConv<Cplx, Cplx, T3>::Fast3DConv --
   **********************************************************************/
  template<class T3>
  template<class KernelContainer, class GridDataContainer>
  void Fast3DConv<Cplx, Cplx, T3>
  ::operator()(const KernelContainer& kernel,
	       GridDataContainer& grid,
	       GridDataContainer& result)
  {
    assert(isKernelTransformed_ == true);
    assert(result.size() == length1D4_);
    
    int stride, dist, times;
    
    if(&result != &grid){
      copy(grid.begin(), grid.end(), result.begin());
    }

    fftw_plan plan;
 
    // transform along z-direction !
    times = ny_;
    stride = 1; dist = nz2_;
    for(size_t ii=0; ii<nx_; ++ii){
      plan = ::fftw_create_plan_specific(nz2_, FFTW_BACKWARD,
					 FFTW_ESTIMATE|FFTW_IN_PLACE,
					 result.p_ + ii * (ny2_ * nz2_),
					 stride, 
					 static_cast<fftw_complex*>(0), stride);

      ::fftw(plan, times, result.p_ + ii*(ny2_ * nz2_),
	   stride, dist, static_cast<fftw_complex*>(0), stride, dist);

      fftw_destroy_plan(plan);
    }
    
    times = nz2_;
    stride = nz2_; dist = 1;
    for(size_t ii=0; ii<nx_; ++ii){
      plan = ::fftw_create_plan_specific(ny2_, FFTW_BACKWARD,
					 FFTW_ESTIMATE|FFTW_IN_PLACE,
					 result.p_ + ii * (ny2_ * nz2_),
					 stride,
					 static_cast<fftw_complex*>(0),
					 stride);

      ::fftw(plan, times, result.p_ + ii*(ny2_ * nz2_),
	   stride, dist, static_cast<fftw_complex*>(0), stride, dist);

      fftw_destroy_plan(plan);
    }

    // 3-in-1
    // combined the last dimension transform(1d), multiplication, inverse transform,
    // then just put the half of the result back.
    fftw_complex *work
     = static_cast<fftw_complex*>(::operator new(nx2_ * sizeof(fftw_complex), nothrow));
    assert(work != 0);

    fftw_plan planPlus = fftw_create_plan(nx2_, FFTW_BACKWARD,
					  FFTW_ESTIMATE|FFTW_IN_PLACE);
    fftw_plan planMinus = fftw_create_plan(nx2_, FFTW_FORWARD,
					   FFTW_ESTIMATE|FFTW_IN_PLACE);
    for(size_t ii=0; ii<(ny2_ * nz2_); ++ii){
      for(size_t jj=0; jj<nx_; ++jj){
	work[jj] = result.p_[ii + jj*(ny2_*nz2_)];
      }
      for(size_t jj=0; jj<nx_; ++jj){
	work[jj + nx_].re = 0.0;
	work[jj + nx_].im = 0.0;
      }
      fftw_one(planPlus, work, 0);
      for(size_t jj=0; jj<nx2_; ++jj){
	work[jj] *= kernel.p_[ii + jj*(ny2_*nz2_)];
      }
      fftw_one(planMinus, work, 0);
      for(size_t jj=0; jj<nx_; ++jj){
	result.p_[ii + jj*(ny2_*nz2_)] = work[jj];
      }   
    }

    ::operator delete(work);

    fftw_destroy_plan(planPlus);
    fftw_destroy_plan(planMinus);

    times = nz2_;
    stride = nz2_; dist = 1;
    for(size_t ii=0; ii<nx_; ++ii){
      plan = ::fftw_create_plan_specific(ny2_, FFTW_FORWARD,
					 FFTW_ESTIMATE|FFTW_IN_PLACE,
					 result.p_ + ii * (ny2_ * nz2_),
					 stride,
					 static_cast<fftw_complex*>(0),
					 stride);
      
      ::fftw(plan, times, result.p_ + ii*(ny2_ * nz2_),
	     stride, dist, static_cast<fftw_complex*>(0), stride, dist);

      fftw_destroy_plan(plan);
    }

    times = ny_;
    stride = 1; dist = nz2_;
    for(size_t ii=0; ii<nx_; ++ii){
      plan = ::fftw_create_plan_specific(nz2_, FFTW_FORWARD,
				       FFTW_ESTIMATE|FFTW_IN_PLACE,
				       result.p_ + ii * (ny2_ * nz2_),
				       stride,
				       static_cast<fftw_complex*>(0),
				       stride);

      ::fftw(plan, times, result.p_ + ii*(ny2_ * nz2_),
	     stride, dist, static_cast<fftw_complex*>(0), stride, dist);

      fftw_destroy_plan(plan);
    }
    
    double scale_factor = 1./length1D8_;
    for(size_t ii=0; ii<length1D4_; ++ii){
      result[ii] *= scale_factor;
    }

  }




  /**********************************************************************
   * Fast3DConv<Cplx, Cplx, T3>::Fast3DConv --
   **********************************************************************/
  template<class T3>
  template<class KernelContainer, class GridDataContainer>
  void Fast3DConv<Cplx, Cplx, T3>
  ::operator()(const KernelContainer& kernelxx, const KernelContainer& kernelzz,
	       const KernelContainer& kernelxz, const KernelContainer& kernelyz,
	       GridDataContainer& gridx, GridDataContainer& gridy,GridDataContainer& gridz,
	       GridDataContainer& resultx, GridDataContainer& resulty,
	       GridDataContainer& resultz)
  {
    assert(isKernelTransformed_ == true);
    assert(resultx.size() == length1D4_);
    
    int stride, dist, times;
    
    if(&resultx != &gridx){
      copy(gridx.begin(), gridx.end(), resultx.begin());
    }
    if(&resulty != &gridy){
      copy(gridy.begin(), gridy.end(), resulty.begin());
    }
    if(&resultz != &gridz){
      copy(gridz.begin(), gridz.end(), resultz.begin());
    }

    fftw_plan planx;
 
    // transform along z-direction !
    times = ny_;
    stride = 1; dist = nz2_;
    for(size_t ii=0; ii<nx_; ++ii){
      planx = ::fftw_create_plan_specific(nz2_, FFTW_BACKWARD,
					 FFTW_ESTIMATE|FFTW_IN_PLACE,
					 resultx.p_ + ii * (ny2_ * nz2_),
					 stride, 
					 static_cast<fftw_complex*>(0), stride);

      ::fftw(planx, times, resultx.p_ + ii*(ny2_ * nz2_),
	   stride, dist, static_cast<fftw_complex*>(0), stride, dist);

      fftw_destroy_plan(planx);
    }
    
    times = nz2_;
    stride = nz2_; dist = 1;
    for(size_t ii=0; ii<nx_; ++ii){
      planx = ::fftw_create_plan_specific(ny2_, FFTW_BACKWARD,
					 FFTW_ESTIMATE|FFTW_IN_PLACE,
					 resultx.p_ + ii * (ny2_ * nz2_),
					 stride,
					 static_cast<fftw_complex*>(0),
					 stride);

      ::fftw(planx, times, resultx.p_ + ii*(ny2_ * nz2_),
	   stride, dist, static_cast<fftw_complex*>(0), stride, dist);

      fftw_destroy_plan(planx);
    }

    fftw_plan plany;
 
    // transform along z-direction !
    times = ny_;
    stride = 1; dist = nz2_;
    for(size_t ii=0; ii<nx_; ++ii){
      plany = ::fftw_create_plan_specific(nz2_, FFTW_BACKWARD,
					 FFTW_ESTIMATE|FFTW_IN_PLACE,
					 resulty.p_ + ii * (ny2_ * nz2_),
					 stride, 
					 static_cast<fftw_complex*>(0), stride);

      ::fftw(plany, times, resulty.p_ + ii*(ny2_ * nz2_),
	   stride, dist, static_cast<fftw_complex*>(0), stride, dist);

      fftw_destroy_plan(plany);
    }
    
    times = nz2_;
    stride = nz2_; dist = 1;
    for(size_t ii=0; ii<nx_; ++ii){
      plany = ::fftw_create_plan_specific(ny2_, FFTW_BACKWARD,
					 FFTW_ESTIMATE|FFTW_IN_PLACE,
					 resulty.p_ + ii * (ny2_ * nz2_),
					 stride,
					 static_cast<fftw_complex*>(0),
					 stride);

      ::fftw(plany, times, resulty.p_ + ii*(ny2_ * nz2_),
	   stride, dist, static_cast<fftw_complex*>(0), stride, dist);

      fftw_destroy_plan(plany);
    }

    fftw_plan planz;
 
    // transform along z-direction !
    times = ny_;
    stride = 1; dist = nz2_;
    for(size_t ii=0; ii<nx_; ++ii){
      planz = ::fftw_create_plan_specific(nz2_, FFTW_BACKWARD,
					 FFTW_ESTIMATE|FFTW_IN_PLACE,
					 resultz.p_ + ii * (ny2_ * nz2_),
					 stride, 
					 static_cast<fftw_complex*>(0), stride);

      ::fftw(planz, times, resultz.p_ + ii*(ny2_ * nz2_),
	   stride, dist, static_cast<fftw_complex*>(0), stride, dist);

      fftw_destroy_plan(planz);
    }
    
    times = nz2_;
    stride = nz2_; dist = 1;
    for(size_t ii=0; ii<nx_; ++ii){
      planz = ::fftw_create_plan_specific(ny2_, FFTW_BACKWARD,
					 FFTW_ESTIMATE|FFTW_IN_PLACE,
					 resultz.p_ + ii * (ny2_ * nz2_),
					 stride,
					 static_cast<fftw_complex*>(0),
					 stride);

      ::fftw(planz, times, resultz.p_ + ii*(ny2_ * nz2_),
	   stride, dist, static_cast<fftw_complex*>(0), stride, dist);

      fftw_destroy_plan(planz);
    }


    // 3-in-1
    // combined the last dimension transform(1d), multiplication, inverse transform,
    // then just put the half of the result back.
    fftw_complex *workx
     = static_cast<fftw_complex*>(::operator new(nx2_ * sizeof(fftw_complex), nothrow));
    assert(workx != 0);
    fftw_complex *worky
     = static_cast<fftw_complex*>(::operator new(nx2_ * sizeof(fftw_complex), nothrow));
    assert(worky != 0);
    fftw_complex *workz
     = static_cast<fftw_complex*>(::operator new(nx2_ * sizeof(fftw_complex), nothrow));
    assert(workz != 0);
    fftw_complex *work_x
     = static_cast<fftw_complex*>(::operator new(nx2_ * sizeof(fftw_complex), nothrow));
    assert(work_x != 0);
    fftw_complex *work_y
     = static_cast<fftw_complex*>(::operator new(nx2_ * sizeof(fftw_complex), nothrow));
    assert(work_y != 0);
    fftw_complex *work_z
     = static_cast<fftw_complex*>(::operator new(nx2_ * sizeof(fftw_complex), nothrow));
    assert(work_z != 0);

    fftw_plan planPlusx = fftw_create_plan(nx2_, FFTW_BACKWARD,
					  FFTW_ESTIMATE|FFTW_IN_PLACE);
    fftw_plan planMinusx = fftw_create_plan(nx2_, FFTW_FORWARD,
					   FFTW_ESTIMATE|FFTW_IN_PLACE);
    fftw_plan planPlusy = fftw_create_plan(nx2_, FFTW_BACKWARD,
					  FFTW_ESTIMATE|FFTW_IN_PLACE);
    fftw_plan planMinusy = fftw_create_plan(nx2_, FFTW_FORWARD,
					   FFTW_ESTIMATE|FFTW_IN_PLACE);
    fftw_plan planPlusz = fftw_create_plan(nx2_, FFTW_BACKWARD,
					  FFTW_ESTIMATE|FFTW_IN_PLACE);
    fftw_plan planMinusz = fftw_create_plan(nx2_, FFTW_FORWARD,
					   FFTW_ESTIMATE|FFTW_IN_PLACE);
    for(size_t ii=0; ii<(ny2_ * nz2_); ++ii){
      for(size_t jj=0; jj<nx_; ++jj){
	workx[jj] = resultx.p_[ii + jj*(ny2_*nz2_)];
	worky[jj] = resulty.p_[ii + jj*(ny2_*nz2_)];
	workz[jj] = resultz.p_[ii + jj*(ny2_*nz2_)];
      }
      for(size_t jj=0; jj<nx_; ++jj){
	workx[jj + nx_].re = 0.0;
	workx[jj + nx_].im = 0.0;
	worky[jj + nx_].re = 0.0;
	worky[jj + nx_].im = 0.0;
	workz[jj + nx_].re = 0.0;
	workz[jj + nx_].im = 0.0;
      }
      fftw_one(planPlusx, workx, 0);
      fftw_one(planPlusy, worky, 0);
      fftw_one(planPlusz, workz, 0);
      for(size_t jj=0; jj<nx2_; ++jj){
	work_x[jj] = workx[jj] * (kernelxx.p_[ii + jj*(ny2_*nz2_)]);
	work_x[jj] += workz[jj] * (kernelxz.p_[ii + jj*(ny2_*nz2_)]);
	work_y[jj] = worky[jj] * (kernelxx.p_[ii + jj*(ny2_*nz2_)]); 
	work_y[jj] += workz[jj] * (kernelxz.p_[ii + jj*(ny2_*nz2_)]);
	work_z[jj] = workx[jj] * (kernelxz.p_[ii + jj*(ny2_*nz2_)]);
	work_z[jj] += worky[jj] * (kernelyz.p_[ii + jj*(ny2_*nz2_)]);
	work_z[jj] += workz[jj] * (kernelzz.p_[ii + jj*(ny2_*nz2_)]);
      }
      fftw_one(planMinusx, work_x, 0);
      fftw_one(planMinusy, work_y, 0);
      fftw_one(planMinusz, work_z, 0);
      for(size_t jj=0; jj<nx_; ++jj){
	resultx.p_[ii + jj*(ny2_*nz2_)] = work_x[jj];
	resulty.p_[ii + jj*(ny2_*nz2_)] = work_y[jj];
	resultz.p_[ii + jj*(ny2_*nz2_)] = work_z[jj];
      }   
    }

    ::operator delete(workx);
    ::operator delete(worky);
    ::operator delete(workz);
    ::operator delete(work_x);
    ::operator delete(work_y);
    ::operator delete(work_z);

    fftw_destroy_plan(planPlusx);
    fftw_destroy_plan(planMinusx);
    fftw_destroy_plan(planPlusy);
    fftw_destroy_plan(planMinusy);
    fftw_destroy_plan(planPlusz);
    fftw_destroy_plan(planMinusz);

    times = nz2_;
    stride = nz2_; dist = 1;
    for(size_t ii=0; ii<nx_; ++ii){
      planx = ::fftw_create_plan_specific(ny2_, FFTW_FORWARD,
					 FFTW_ESTIMATE|FFTW_IN_PLACE,
					 resultx.p_ + ii * (ny2_ * nz2_),
					 stride,
					 static_cast<fftw_complex*>(0),
					 stride);
      
      ::fftw(planx, times, resultx.p_ + ii*(ny2_ * nz2_),
	     stride, dist, static_cast<fftw_complex*>(0), stride, dist);

      fftw_destroy_plan(planx);
    }

    times = ny_;
    stride = 1; dist = nz2_;
    for(size_t ii=0; ii<nx_; ++ii){
      planx = ::fftw_create_plan_specific(nz2_, FFTW_FORWARD,
				       FFTW_ESTIMATE|FFTW_IN_PLACE,
				       resultx.p_ + ii * (ny2_ * nz2_),
				       stride,
				       static_cast<fftw_complex*>(0),
				       stride);

      ::fftw(planx, times, resultx.p_ + ii*(ny2_ * nz2_),
	     stride, dist, static_cast<fftw_complex*>(0), stride, dist);

      fftw_destroy_plan(planx);
    }
    
   times = nz2_;
    stride = nz2_; dist = 1;
    for(size_t ii=0; ii<nx_; ++ii){
      plany = ::fftw_create_plan_specific(ny2_, FFTW_FORWARD,
					 FFTW_ESTIMATE|FFTW_IN_PLACE,
					 resulty.p_ + ii * (ny2_ * nz2_),
					 stride,
					 static_cast<fftw_complex*>(0),
					 stride);
      
      ::fftw(plany, times, resulty.p_ + ii*(ny2_ * nz2_),
	     stride, dist, static_cast<fftw_complex*>(0), stride, dist);

      fftw_destroy_plan(plany);
    }

    times = ny_;
    stride = 1; dist = nz2_;
    for(size_t ii=0; ii<nx_; ++ii){
      plany = ::fftw_create_plan_specific(nz2_, FFTW_FORWARD,
				       FFTW_ESTIMATE|FFTW_IN_PLACE,
				       resulty.p_ + ii * (ny2_ * nz2_),
				       stride,
				       static_cast<fftw_complex*>(0),
				       stride);

      ::fftw(plany, times, resulty.p_ + ii*(ny2_ * nz2_),
	     stride, dist, static_cast<fftw_complex*>(0), stride, dist);

      fftw_destroy_plan(plany);
    }

   times = nz2_;
    stride = nz2_; dist = 1;
    for(size_t ii=0; ii<nx_; ++ii){
      planz = ::fftw_create_plan_specific(ny2_, FFTW_FORWARD,
					 FFTW_ESTIMATE|FFTW_IN_PLACE,
					 resultz.p_ + ii * (ny2_ * nz2_),
					 stride,
					 static_cast<fftw_complex*>(0),
					 stride);
      
      ::fftw(planz, times, resultz.p_ + ii*(ny2_ * nz2_),
	     stride, dist, static_cast<fftw_complex*>(0), stride, dist);

      fftw_destroy_plan(planz);
    }

    times = ny_;
    stride = 1; dist = nz2_;
    for(size_t ii=0; ii<nx_; ++ii){
      planz = ::fftw_create_plan_specific(nz2_, FFTW_FORWARD,
				       FFTW_ESTIMATE|FFTW_IN_PLACE,
				       resultz.p_ + ii * (ny2_ * nz2_),
				       stride,
				       static_cast<fftw_complex*>(0),
				       stride);

      ::fftw(planz, times, resultz.p_ + ii*(ny2_ * nz2_),
	     stride, dist, static_cast<fftw_complex*>(0), stride, dist);

      fftw_destroy_plan(planz);
    }


    double scale_factor = 1./length1D8_;
    for(size_t ii=0; ii<length1D4_; ++ii){
      resultx[ii] *= scale_factor;
      resulty[ii] *= scale_factor;
      resultz[ii] *= scale_factor;
    }

  }


  /**********************************************************************
   * Fast3DConv<Cplx, Cplx, T3>::Fast3DConv --
   **********************************************************************/
  template<class T3>
  template<class KernelContainer, class GridDataContainer>
  void Fast3DConv<Cplx, Cplx, T3>
  ::operator()(const KernelContainer& kernel,
	       GridDataContainer& gridx, GridDataContainer& gridy,GridDataContainer& gridz,
	       GridDataContainer& resultx, GridDataContainer& resulty,
	       GridDataContainer& resultz)
  {
    assert(isKernelTransformed_ == true);
    assert(resultx.size() == length1D4_);
    
    int stride, dist, times;
    if(&resultx != &gridx){
      copy(gridx.begin(), gridx.end(), resultx.begin());
    }
    if(&resulty != &gridy){
      copy(gridy.begin(), gridy.end(), resulty.begin());
    }
    if(&resultz != &gridz){
      copy(gridz.begin(), gridz.end(), resultz.begin());
    }

    fftw_plan planx;
 
    // transform along z-direction !
    times = ny_;
    stride = 1; dist = nz2_;
    for(size_t ii=0; ii<nx_; ++ii){
      planx = ::fftw_create_plan_specific(nz2_, FFTW_BACKWARD,
					 FFTW_ESTIMATE|FFTW_IN_PLACE,
					 resultx.p_ + ii * (ny2_ * nz2_),
					 stride, 
					 static_cast<fftw_complex*>(0), stride);

      ::fftw(planx, times, resultx.p_ + ii*(ny2_ * nz2_),
	   stride, dist, static_cast<fftw_complex*>(0), stride, dist);

      fftw_destroy_plan(planx);
    }
    
    times = nz2_;
    stride = nz2_; dist = 1;
    for(size_t ii=0; ii<nx_; ++ii){
      planx = ::fftw_create_plan_specific(ny2_, FFTW_BACKWARD,
					 FFTW_ESTIMATE|FFTW_IN_PLACE,
					 resultx.p_ + ii * (ny2_ * nz2_),
					 stride,
					 static_cast<fftw_complex*>(0),
					 stride);

      ::fftw(planx, times, resultx.p_ + ii*(ny2_ * nz2_),
	   stride, dist, static_cast<fftw_complex*>(0), stride, dist);

      fftw_destroy_plan(planx);
    }

    fftw_plan plany;
 
    // transform along z-direction !
    times = ny_;
    stride = 1; dist = nz2_;
    for(size_t ii=0; ii<nx_; ++ii){
      plany = ::fftw_create_plan_specific(nz2_, FFTW_BACKWARD,
					 FFTW_ESTIMATE|FFTW_IN_PLACE,
					 resulty.p_ + ii * (ny2_ * nz2_),
					 stride, 
					 static_cast<fftw_complex*>(0), stride);

      ::fftw(plany, times, resulty.p_ + ii*(ny2_ * nz2_),
	   stride, dist, static_cast<fftw_complex*>(0), stride, dist);

      fftw_destroy_plan(plany);
    }
    
    times = nz2_;
    stride = nz2_; dist = 1;
    for(size_t ii=0; ii<nx_; ++ii){
      plany = ::fftw_create_plan_specific(ny2_, FFTW_BACKWARD,
					 FFTW_ESTIMATE|FFTW_IN_PLACE,
					 resulty.p_ + ii * (ny2_ * nz2_),
					 stride,
					 static_cast<fftw_complex*>(0),
					 stride);

      ::fftw(plany, times, resulty.p_ + ii*(ny2_ * nz2_),
	   stride, dist, static_cast<fftw_complex*>(0), stride, dist);

      fftw_destroy_plan(plany);
    }

    fftw_plan planz;
 
    // transform along z-direction !
    times = ny_;
    stride = 1; dist = nz2_;
    for(size_t ii=0; ii<nx_; ++ii){
      planz = ::fftw_create_plan_specific(nz2_, FFTW_BACKWARD,
					 FFTW_ESTIMATE|FFTW_IN_PLACE,
					 resultz.p_ + ii * (ny2_ * nz2_),
					 stride, 
					 static_cast<fftw_complex*>(0), stride);

      ::fftw(planz, times, resultz.p_ + ii*(ny2_ * nz2_),
	   stride, dist, static_cast<fftw_complex*>(0), stride, dist);

      fftw_destroy_plan(planz);
    }
    
    times = nz2_;
    stride = nz2_; dist = 1;
    for(size_t ii=0; ii<nx_; ++ii){
      planz = ::fftw_create_plan_specific(ny2_, FFTW_BACKWARD,
					 FFTW_ESTIMATE|FFTW_IN_PLACE,
					 resultz.p_ + ii * (ny2_ * nz2_),
					 stride,
					 static_cast<fftw_complex*>(0),
					 stride);

      ::fftw(planz, times, resultz.p_ + ii*(ny2_ * nz2_),
	   stride, dist, static_cast<fftw_complex*>(0), stride, dist);

      fftw_destroy_plan(planz);
    }


    // 3-in-1
    // combined the last dimension transform(1d), multiplication, inverse transform,
    // then just put the half of the result back.
    fftw_complex *workx
     = static_cast<fftw_complex*>(::operator new(nx2_ * sizeof(fftw_complex), nothrow));
    assert(workx != 0);
    fftw_complex *worky
     = static_cast<fftw_complex*>(::operator new(nx2_ * sizeof(fftw_complex), nothrow));
    assert(worky != 0);
    fftw_complex *workz
     = static_cast<fftw_complex*>(::operator new(nx2_ * sizeof(fftw_complex), nothrow));
    assert(workz != 0);
    fftw_complex *work_x
     = static_cast<fftw_complex*>(::operator new(nx2_ * sizeof(fftw_complex), nothrow));
    assert(work_x != 0);
    fftw_complex *work_y
     = static_cast<fftw_complex*>(::operator new(nx2_ * sizeof(fftw_complex), nothrow));
    assert(work_y != 0);
    fftw_complex *work_z
     = static_cast<fftw_complex*>(::operator new(nx2_ * sizeof(fftw_complex), nothrow));
    assert(work_z != 0);

    fftw_plan planPlusx = fftw_create_plan(nx2_, FFTW_BACKWARD,
					  FFTW_ESTIMATE|FFTW_IN_PLACE);
    fftw_plan planMinusx = fftw_create_plan(nx2_, FFTW_FORWARD,
					   FFTW_ESTIMATE|FFTW_IN_PLACE);
    fftw_plan planPlusy = fftw_create_plan(nx2_, FFTW_BACKWARD,
					  FFTW_ESTIMATE|FFTW_IN_PLACE);
    fftw_plan planMinusy = fftw_create_plan(nx2_, FFTW_FORWARD,
					   FFTW_ESTIMATE|FFTW_IN_PLACE);
    fftw_plan planPlusz = fftw_create_plan(nx2_, FFTW_BACKWARD,
					  FFTW_ESTIMATE|FFTW_IN_PLACE);
    fftw_plan planMinusz = fftw_create_plan(nx2_, FFTW_FORWARD,
					   FFTW_ESTIMATE|FFTW_IN_PLACE);
    for(size_t ii=0; ii<(ny2_ * nz2_); ++ii){
      for(size_t jj=0; jj<nx_; ++jj){
	workx[jj] = resultx.p_[ii + jj*(ny2_*nz2_)];
	worky[jj] = resulty.p_[ii + jj*(ny2_*nz2_)];
	workz[jj] = resultz.p_[ii + jj*(ny2_*nz2_)];
      }
      for(size_t jj=0; jj<nx_; ++jj){
	workx[jj + nx_].re = 0.0;
	workx[jj + nx_].im = 0.0;
	worky[jj + nx_].re = 0.0;
	worky[jj + nx_].im = 0.0;
	workz[jj + nx_].re = 0.0;
	workz[jj + nx_].im = 0.0;
      }
      fftw_one(planPlusx, workx, 0);
      fftw_one(planPlusy, worky, 0);
      fftw_one(planPlusz, workz, 0);
      for(size_t jj=0; jj<nx2_; ++jj){
	work_x[jj] = workx[jj] * (kernel.p_[ii + jj*(ny2_*nz2_)]);
	work_y[jj] = worky[jj] * (kernel.p_[ii + jj*(ny2_*nz2_)]); 
	work_z[jj] = workx[jj] * (kernel.p_[ii + jj*(ny2_*nz2_)]);
      }
      fftw_one(planMinusx, work_x, 0);
      fftw_one(planMinusy, work_y, 0);
      fftw_one(planMinusz, work_z, 0);
      for(size_t jj=0; jj<nx_; ++jj){
	resultx.p_[ii + jj*(ny2_*nz2_)] = work_x[jj];
	resulty.p_[ii + jj*(ny2_*nz2_)] = work_y[jj];
	resultz.p_[ii + jj*(ny2_*nz2_)] = work_z[jj];
      }   
    }

    ::operator delete(workx);
    ::operator delete(worky);
    ::operator delete(workz);
    ::operator delete(work_x);
    ::operator delete(work_y);
    ::operator delete(work_z);

    fftw_destroy_plan(planPlusx);
    fftw_destroy_plan(planMinusx);
    fftw_destroy_plan(planPlusy);
    fftw_destroy_plan(planMinusy);
    fftw_destroy_plan(planPlusz);
    fftw_destroy_plan(planMinusz);

    times = nz2_;
    stride = nz2_; dist = 1;
    for(size_t ii=0; ii<nx_; ++ii){
      planx = ::fftw_create_plan_specific(ny2_, FFTW_FORWARD,
					 FFTW_ESTIMATE|FFTW_IN_PLACE,
					 resultx.p_ + ii * (ny2_ * nz2_),
					 stride,
					 static_cast<fftw_complex*>(0),
					 stride);
      
      ::fftw(planx, times, resultx.p_ + ii*(ny2_ * nz2_),
	     stride, dist, static_cast<fftw_complex*>(0), stride, dist);

      fftw_destroy_plan(planx);
    }

    times = ny_;
    stride = 1; dist = nz2_;
    for(size_t ii=0; ii<nx_; ++ii){
      planx = ::fftw_create_plan_specific(nz2_, FFTW_FORWARD,
				       FFTW_ESTIMATE|FFTW_IN_PLACE,
				       resultx.p_ + ii * (ny2_ * nz2_),
				       stride,
				       static_cast<fftw_complex*>(0),
				       stride);

      ::fftw(planx, times, resultx.p_ + ii*(ny2_ * nz2_),
	     stride, dist, static_cast<fftw_complex*>(0), stride, dist);

      fftw_destroy_plan(planx);
    }
    
   times = nz2_;
    stride = nz2_; dist = 1;
    for(size_t ii=0; ii<nx_; ++ii){
      plany = ::fftw_create_plan_specific(ny2_, FFTW_FORWARD,
					 FFTW_ESTIMATE|FFTW_IN_PLACE,
					 resulty.p_ + ii * (ny2_ * nz2_),
					 stride,
					 static_cast<fftw_complex*>(0),
					 stride);
      
      ::fftw(plany, times, resulty.p_ + ii*(ny2_ * nz2_),
	     stride, dist, static_cast<fftw_complex*>(0), stride, dist);

      fftw_destroy_plan(plany);
    }

    times = ny_;
    stride = 1; dist = nz2_;
    for(size_t ii=0; ii<nx_; ++ii){
      plany = ::fftw_create_plan_specific(nz2_, FFTW_FORWARD,
				       FFTW_ESTIMATE|FFTW_IN_PLACE,
				       resulty.p_ + ii * (ny2_ * nz2_),
				       stride,
				       static_cast<fftw_complex*>(0),
				       stride);

      ::fftw(plany, times, resulty.p_ + ii*(ny2_ * nz2_),
	     stride, dist, static_cast<fftw_complex*>(0), stride, dist);

      fftw_destroy_plan(plany);
    }

   times = nz2_;
    stride = nz2_; dist = 1;
    for(size_t ii=0; ii<nx_; ++ii){
      planz = ::fftw_create_plan_specific(ny2_, FFTW_FORWARD,
					 FFTW_ESTIMATE|FFTW_IN_PLACE,
					 resultz.p_ + ii * (ny2_ * nz2_),
					 stride,
					 static_cast<fftw_complex*>(0),
					 stride);
      
      ::fftw(planz, times, resultz.p_ + ii*(ny2_ * nz2_),
	     stride, dist, static_cast<fftw_complex*>(0), stride, dist);

      fftw_destroy_plan(planz);
    }

    times = ny_;
    stride = 1; dist = nz2_;
    for(size_t ii=0; ii<nx_; ++ii){
      planz = ::fftw_create_plan_specific(nz2_, FFTW_FORWARD,
				       FFTW_ESTIMATE|FFTW_IN_PLACE,
				       resultz.p_ + ii * (ny2_ * nz2_),
				       stride,
				       static_cast<fftw_complex*>(0),
				       stride);

      ::fftw(planz, times, resultz.p_ + ii*(ny2_ * nz2_),
	     stride, dist, static_cast<fftw_complex*>(0), stride, dist);

      fftw_destroy_plan(planz);
    }


    double scale_factor = 1./length1D8_;
    for(size_t ii=0; ii<length1D4_; ++ii){
      resultx[ii] *= scale_factor;
      resulty[ii] *= scale_factor;
      resultz[ii] *= scale_factor;
    }

  }

  /**********************************************************************
   * Fast3DConv<Cplx, Cplx, T3>
   * ::SetupTestGridData<Vector> --
   **********************************************************************/
  template<class T3>
  template<class Vector>
  void Fast3DConv<Cplx, Cplx, T3>
  ::setupTestGridData (Vector& gridDataVec, bool randomOrNot)
  {
    assert(gridDataVec.size() == length1D_);
    for(size_t ii=0; ii<length1D_; ++ii){
      if(randomOrNot){
	gridDataVec[ii] = Cplx(mydrand(), mydrand());
      } else {
	gridDataVec[ii] = Cplx(1., 0.);
      }
    } 
  }

  /**********************************************************************
   * Fast3DConv<Cplx, Cplx, T3>::GetDirectPotential --
   **********************************************************************/
  template<class T3>
  template<class VectorA, class VectorB>
  void Fast3DConv<Cplx, Cplx, T3> 
  ::getDirectPotential (const VectorA gridData,
			VectorB& directPotential)
  {
    assert(directPotential.size() == length1D4_);
    assert(isCalcKernelSetup_);

    //    assert(&gridData != &directPotential);
      
    int tIndex, sIndex;
    for(int tx=0; tx<nx_; ++tx){
      for(int ty=0; ty<ny_; ++ty){
	for(int tz=0; tz<nz_; ++tz){
	  tIndex = tz + nz_ * (ty + ny_ * tx);
	  directPotential[tIndex] = 0.;
	  for(int sx=0; sx<nx_; ++sx){
	    for(int sy=0; sy<ny_; ++sy){
	      for(int sz=0; sz<nz_; ++sz){
		sIndex = sz + nz_ * (sy + ny_ * sx);
		directPotential[tIndex] 
		  += (*calcKernelPtr_)((tx-sx)*pace_[0], (ty-sy)*pace_[1], (tz-sz)*pace_[2]) * gridData[sIndex];
	      }
	    }
	  }
	}
      }
    }
  }

  
  /**********************************************************************
   * compare --
   **********************************************************************/
  template<class T3>
  template<class VectorA, class VectorB>
  bool Fast3DConv<Cplx, Cplx, T3>
  ::compare (
	     const VectorA& directPotential,
	     const VectorB& fastPotential)
	     
  {
    assert(directPotential.size() == fastPotential.size());

    double TOLERANCE = 1.e-6;
    double reErr;
    
    typename VectorA::const_iterator iterA = directPotential.begin();
    typename VectorA::const_iterator iterB = fastPotential.begin();
    while(iterA != directPotential.end()){
      reErr = abs(*iterA - *iterB) / abs(*iterA);
      if(reErr > TOLERANCE){
	return false;
      }
      ++iterA; ++iterB;
    }
      
    return true;
  }


} // namespace pfft

#endif



