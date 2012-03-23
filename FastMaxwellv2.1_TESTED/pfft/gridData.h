/*-*-c++-*-*/
/******************************************************************************
 **                                                                          **
 ** Copyright (c) 2000 Massachusetts Institute of Technology, Cambridge, MA. **
 ** All rights reserved.                                                     **
 **                                                                          **
 ******************************************************************************
  ==========================================================================
  ============================= File Info ==================================

  Author: Zhenhai Zhu, Ben Song
  
  Description:

  Resources:

  See also:

  const static char cvsid[] = "$Id: gridData.h,v 1.4 2002/10/01 21:17:37 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _GRID_DATA_H_
#define _GRID_DATA_H_

#include "TArray.h"
#include "x3dconv_rr.h"
#include "x3dconv_rc.h"
#include "x3dconv_cc.h"
#include "spRowMat.h"
#include "spColMat.h"
#include "../src/generate_dipoles.h"

namespace pfft {
  
  template <class KernelType, class DataType, class GreenFunc>
  class GridData : public Fast3DConv<KernelType, DataType, GreenFunc>  {
    
  public:
    GridData(void) {}
    GridData(const GreenFunc& greenFuncIn) : greenFunc_(greenFuncIn) {}
    GridData(const GreenFunc& greenFuncIn, const Dipoles& dipoles) : greenFunc_(greenFuncIn), dipoles_(dipoles) {}
    
    size_t memoryUsage (size_t, size_t, size_t) const;
    void fillKernel(void) { 
      setupKernel(gridKernel_); 
    }
    void fillKernelScalarImages(void) { 
      setupScalarImages (gridKernelImages_); 
    }
    void fillKerneldiag(void) { 
      setupKernelDiag(gridKernelxx_,gridKernelzz_); 
    }
    void fillKerneloffdiag(void) { 
      setupKernelOffDiag(gridKernelxz_,gridKernelyz_); 
    }
    void fftKernel(void) { 
      fftOnKernel(gridKernel_); 
    }
    void fftKernelScalarImages(void) { 
      fftOnKernel(gridKernelImages_); 
    }
    
    void fftKerneldyadic(void) { 
      fftOnKernel(gridKernelxx_); 
      fftOnKernel(gridKernelzz_); 
      fftOnKernel(gridKernelxz_); 
      fftOnKernel(gridKernelyz_); 
    }
    
    void PermutationImages(void)
    {

      PermuteOnData(PgridData, gridDataImages_);
      for (int i=0; i < PgridData.size(); i++)
	gridDataImages_[i]=PgridData[i];
    }
    
    
    
    void PermutationDyadic(void)
    {

      PermuteOnData(PgridData, gridDatax_);
      for (int i=0; i < PgridData.size(); i++) 
	gridDatax_[i]=PgridData[i];

      PermuteOnData(PgridData, gridDatay_);
      for (int i=0; i < PgridData.size(); i++)
	gridDatay_[i]=PgridData[i];

      PermuteOnData(PgridData, gridDataz_);
      for (int i=0; i < PgridData.size(); i++)
	gridDataz_[i]=PgridData[i];
    }
    
    void conv3D(void) { 
      operator() (gridKernel_, gridData_, gridData_); 
    }

    void conv3DScalarImages (void) { 
      operator() (gridKernelImages_, gridDataImages_, gridDataImages_); 
    }    
    
    void conv3Ddyadic(void) { 
      operator() (gridKernelxx_, gridKernelzz_, 
		  gridKernelxz_, gridKernelyz_, gridDatax_, gridDatay_, 
		  gridDataz_, gridDatax_, gridDatay_, gridDataz_); 
    }

    void conv3Dtoeplitz(void) { 
      operator() (gridKernel_, gridDatax_, gridDatay_, 
		  gridDataz_, gridDatax_, gridDatay_, gridDataz_); 
    }
    
    
    void allocate(const size_t nx, const size_t ny, const size_t nz, 
		  const double gridStep, int flag)
    {
      PgridData.resize(nx, ny, nz);
      gridKernel_.resize(nx, ny, nz);
      gridData_.resize(nx, ny, nz);
      //if (flag==1)
	gridKernelxx_.resize(nx, ny, nz);
	//else if (flag == 2) {
	//gridKernelxx_.resize(nx, ny, nz);
	gridKernelzz_.resize(nx, ny, nz);
	gridKernelxz_.resize(nx, ny, nz);
	gridKernelyz_.resize(nx, ny, nz);
	//}
      gridDatax_.resize(nx, ny, nz);
      gridDatay_.resize(nx, ny, nz);
      gridDataz_.resize(nx, ny, nz);
      initialize(nx, ny, nz, gridStep, gridStep, gridStep, greenFunc_);
    }

    void allocateScalar(const size_t nx, const size_t ny, const size_t nz, 
		  const double gridStep, int flag)
    {
      PgridData.resize(nx, ny, nz);
      gridKernel_.resize(nx, ny, nz);
      gridData_.resize(nx, ny, nz);
      // if (flag == 1) {
	gridKernelImages_.resize(nx, ny, nz);
	gridDataImages_.resize(nx, ny, nz);
	//}
      initialize(nx, ny, nz, gridStep, gridStep, gridStep, greenFunc_);
    }

    void outputGridData(std::ostream&) const;
    void outputGridKernel(std::ostream&) const;
    
    template <class ProjectMat, class SrcVec>
    void projectOnGrid(const ProjectMat& projectMat, const SrcVec& src) {
      matMultVec(gridData_, projectMat, src);
    }
    
    template <class InterpMat, class DesVec>
    void interpFromGrid(const InterpMat& interpMat, DesVec& des) const {
      matMultVec(des, interpMat, gridData_);
    }

    template <class ProjectMat, class SrcVec>
    void projectOnGridScalarImages(const ProjectMat& projectMat, const SrcVec& src) {
      matMultVec(gridDataImages_, projectMat, src);
    }
    
    template <class InterpMat, class DesVec>
    void interpFromGridScalarImages(const InterpMat& interpMat, DesVec& des) const {
      matMultVec(des, interpMat, gridDataImages_);
    }


    template <class ProjectMat, class SrcVec>
    void projectOnGridDyadic(const ProjectMat& projectMat, const SrcVec& srcx, const SrcVec& srcy, const SrcVec& srcz)
    {
      matMultVec(gridDatax_, projectMat, srcx);
      matMultVec(gridDatay_, projectMat, srcy);
      matMultVec(gridDataz_, projectMat, srcz);
    }

    template <class InterpMat, class DesVec>
    void interpFromGridDyadic(const InterpMat& interpMat, 
			      DesVec& desx, DesVec& desy, DesVec& desz) const 
    {
      matMultVec(desx, interpMat, gridDatax_);
      matMultVec(desy, interpMat, gridDatay_);
      matMultVec(desz, interpMat, gridDataz_);
    }
    
    
    const typename TypeSelect<DataType>::TYPE& operator[] (size_t i) const {
      return gridData_[i];
    }

  private:
    GreenFunc greenFunc_;
    Dipoles dipoles_;
    TArray<typename TypeSelect<KernelType>::TYPE, KernelTag> gridKernel_, 
      gridKernelxx_, gridKernelzz_, gridKernelxz_, gridKernelyz_, gridKernelImages_;
    TArray<typename TypeSelect<DataType>::TYPE, GridDataTag> gridData_, 
      PgridData, gridDatax_, gridDatay_, gridDataz_, gridDataImages_; 
    
  };
  
  /**********************************************************************
   * memoryEstimate --
   **********************************************************************/
  template <class KernelType, class DataType, class GreenFunc>
  size_t
  GridData<KernelType, DataType, GreenFunc>::
  memoryUsage (
	       size_t numPointX, 
	       size_t numPointY, 
	       size_t numPointZ) const
  {
    size_t nx = 2 * numPointX;
    size_t ny = 2 * numPointY;
    size_t nz = 2 * numPointZ;
    size_t gridDataUsage = sizeof(KernelType) * (nx * ny * numPointZ);
    size_t kernelUsage = sizeof(DataType) * (nx * ny * nz) * 1;
      
    return gridDataUsage + kernelUsage;
  }
  
  /**********************************************************************
   * outputGridData --
   **********************************************************************/
  template <class KernelType, class DataType, class GreenFunc>
  void
  GridData<KernelType, DataType, GreenFunc>::
  outputGridData (std::ostream& out) const
  {
    out<<" \n Grid Data \n"; 
    for (size_t ii=0; ii < gridData_.size(); ++ii) {
      out << gridData_[ii] << std::endl;
    }
    out << std::endl;
  }


  /**********************************************************************
   * outputGridKernel --
   **********************************************************************/
  template <class KernelType, class DataType, class GreenFunc>
  void
  GridData<KernelType, DataType, GreenFunc>::
  outputGridKernel (std::ostream& out) const
  {
    out<<" \n Grid Data \n";
    for (size_t ii=0; ii < gridData_.size(); ++ii) {
      out << gridKernel_[ii] << std::endl;
    }
    out << std::endl;
  }


} //namespace pfft

#endif








