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

  This is main header file for pfft module

  Resources:

  See also:

  const static char cvsid[] = "$Id: pfft.h,v 1.24 2003/07/10 19:14:26 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _PFFT_H_

#define _PFFT_H_

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "element.h"
#include "gridElement.h"
#include "grid.h"
#include "gridData.h"
#include "directMat.h"
#include "projectMat.h"
#include "interpMat.h"
#include "stpwatch.h"
#include "utils.h"

namespace pfft {

#ifndef ULLONG_MAX
#define ULLONG_MAX	18446744073709551615ULL
#endif

  //  const std::string PFFT_VERSION = "pfft++ v1.1";

  enum KernelSize { XY, XY_dipoles, XYZ, XYZ_dipoles, Scalar, Scalar_dipoles};


  template <class Kernel, class Data, class GreenFunc, class Integration>
  class Pfft {

    typedef Kernel KernelValueType;
    typedef Data DataType;
    typedef GreenFunc GreenFuncType;
    typedef Integration CalcpType;

  public:

    Pfft(void) : constructCompleted_(false), gridElementDeleted_(false) {}


   Pfft(
	 std::vector<element>& sourceElementList,
	 std::vector<element>& targetElementList,
	 const GreenFunc& greenFunc,
	 const Integration& integration,
	 const int directStencilSize,
	 const int projectStencilSize,
	 const int interpStencilSize, KernelSize kernelsize);


   Pfft(
	 std::vector<element>& sourceElementList,
	 std::vector<element>& targetElementList,
	 const Integration& integration,
	 const int directStencilSize = 3,
	 const int projectStencilSize = 1,
	 const int interpStencilSize = 1, KernelSize kernelsize = XYZ_dipoles);

    ~Pfft(void) {
      if (! gridElementDeleted_)
	deleteGridElement();
    }

    void constructKernelIndependent(
				    std::vector<element>& sourceElementList,
				    std::vector<element>& targetElementList);
    
    void constructKernelDependent(
				  const GreenFunc& greenFunc,
				  const Integration& integration);
    
    //pfft::SpColMat<complex<double> > getPreCondMat(void) { return directMat.getPreCondMat(); }
    std::vector<complex<double> > getPreCondMat(void) { 
      return directMat.getPreCondMat(); 
    }
    //SpColMat<Kernel>& getPreCondMat(void) { 
    //  return directMat[0];
    //}
    size_t memoryEstimate (void) const { return minMemoryUsage; }
    
    void deleteGridElement(void) {
      if (constructCompleted_) {
	delete srcGridElementPtr;      
	if (srcAndEvalDiff)
	  delete evalGridElementPtr;
      }
      gridElementDeleted_ = true;
    }

    template<class VecA, class VecB>
    void IEoperator(
		    VecA& des,
		    const DifferentialOperator outerOperator,
		    const DifferentialOperator innerOperator,
		    const VecB& src);
    
    template <class VecA, class VecB>
    void IEoperatorDyadic (
			   VecA& des, VecA& desx,  VecA& desy,  VecA& desz,
			   const DifferentialOperator outerOperator,
			   const DifferentialOperator innerOperator,
			   const VecB& src, const VecB& srcx, const VecB& srcy, const VecB& srcz);    
    
    
    template<class VecA, class VecB>
    void IEoperatorWithoutFastConvolution(
					  VecA& des,
					  const DifferentialOperator outerOperator,
					  const DifferentialOperator innerOperator,
					  const VecB& src);
    
  private:
    // state variables
    GreenFunc greenFunc;
    Integration integration;
    int directStencilSize;
    int projectStencilSize;
    int interpStencilSize;
    bool constructCompleted_;
    KernelSize kernelSize_;
    // bool GridCompleted_;
    
    // These are four most important matrices, results of pfft construction
    DirectMat<GreenFunc, Integration, Kernel> directMat;
    ProjectMat projectMat;
    InterpMat interpMat;
    GridData<Kernel, Data, GreenFunc> gridData;

    // working variables
    bool srcAndEvalDiff;
    size_t minMemoryUsage, currentMemoryUsage;
    // this must be trnasformed into array
    // please check the delete statement  
    GridElement* srcGridElementPtr;
    GridElement* evalGridElementPtr;
    Grid grid;
    Stencil projectStencil;
    Stencil interpStencil;
    Stencil directStencil;
    bool reportMatVecMultTime_;
    bool gridElementDeleted_;
    
    void checkDirectStencilSize (void);
    void initialize (std::vector<element>&, std::vector<element>&);
    void findOptimalNumGridPoint (void);
    void setupGridAndGridElement(void);


    size_t estimateMemoryByDirectPart (void);
    enum MemoryStatus { IMPROVING, WORSENING, WORSENING_BUT_KEEP_GOING };
    MemoryStatus checkMemoryUsage (void);

    void memoryReport (const size_t, const size_t) const;
    void memoryReport (const std::string matrixName, const double memory) const;
    void gridAndElementReport (void) const;
    void progressReport (const std::string message) const;
    void IEOperatorTimeReport (double time);
  };

  /**********************************************************************
   * Pfft --
   **********************************************************************/
  template <class Kernel, class Data, class GreenFunc, class Integration>
  Pfft<Kernel, Data, GreenFunc, Integration>::
  Pfft (
		 std::vector<element>& sourceElementList,
		 std::vector<element>& targetElementList,
		 const GreenFunc& greenFuncIn,
		 const Integration& integrationIn,
		 const int directStencilSizeIn,
		 const int projectStencilSizeIn,
		 const int interpStencilSizeIn,
		 KernelSize kernelSize)
  {
    greenFunc=(greenFuncIn); integration=(integrationIn);
    directStencilSize=(directStencilSizeIn);
    projectStencilSize=(projectStencilSizeIn);
    interpStencilSize=(interpStencilSizeIn);
    kernelSize_=(kernelSize);
    checkDirectStencilSize();    
    constructKernelIndependent(sourceElementList, targetElementList);
    
    constructKernelDependent(greenFuncIn, integrationIn);
    progressReport("pfft setup is done!");
    
    constructCompleted_ = true;
    reportMatVecMultTime_ = true;
    
    // gridElement is not needed after this point. 
    // This will save a lot of memory. But it is not possible if
    // the user calls the constructKernelDependent multiple times.
    // This constructor is for one frequency point simulation. So
    // it makes sense to delete gridElement here.
    deleteGridElement();
  }



  /**********************************************************************
   * Pfft --
   * The integration object could be a dummy one except that it 
   * should have right inner and outer operator list. This list
   * will be used to construct interpMat and projectMat. And after
   * the construction, both matrices are not dependent on the list
   * any more. So integration object could be overide later by a real one. 
   * So this function is still a kernel independent one.
   *
   * Note:
   * This constructor should be followed by the calling of 
   * constructKernelDependent() to complete the construction of a 
   * pfft object.
   **********************************************************************/
  template <class Kernel, class Data, class GreenFunc, class Integration>
  Pfft<Kernel, Data, GreenFunc, Integration>::
  Pfft (
	std::vector<element>& sourceElementList,
	std::vector<element>& targetElementList,
	const Integration& integrationIn,
	const int directStencilSizeIn,
	const int projectStencilSizeIn,
	const int interpStencilSizeIn,
	KernelSize kernelSize)
    : integration(integrationIn)
      , directStencilSize(directStencilSizeIn)
      , projectStencilSize(projectStencilSizeIn)
      , interpStencilSize(interpStencilSizeIn)
      , kernelSize_(kernelSize)
  {
    // GridCompleted_ = false;
    checkDirectStencilSize();
    constructKernelIndependent(sourceElementList, targetElementList);
    constructCompleted_ = false;
    // GridCompleted_ = true;
    reportMatVecMultTime_ = true;
    gridElementDeleted_ = false;
  }
  
  /**********************************************************************
   * checkDirectStencilSize --
   * Make sure directStencilSize >= InterpStencilSize + projectStencilSize
   **********************************************************************/
  template <class Kernel, class Data, class GreenFunc, class Integration>
  void
  Pfft<Kernel, Data, GreenFunc, Integration>::
  checkDirectStencilSize (
			  void)
  {
    int minDirectStencilSize = projectStencilSize + interpStencilSize;
    if (directStencilSize < minDirectStencilSize) {
      warningMessage("pfft.h::checkDirectStencilSize()",
		     "illegal direct stencil size, the new size is ", minDirectStencilSize);;
      directStencilSize = minDirectStencilSize;
    }    
  }

  /**********************************************************************
   * constructKernelIndependent --
   **********************************************************************/
  template <class Kernel, class Data, class GreenFunc, class Integration>
  void
  Pfft<Kernel, Data, GreenFunc, Integration>::
  constructKernelIndependent (
			      std::vector<element>& sourceElementList,
			      std::vector<element>& targetElementList)
  {
    progressReport("Start to setup kernel independent part of pfft ......");
    
    initialize(sourceElementList, targetElementList);
    setupGridAndGridElement();
    gridAndElementReport();
    
    interpStencil.setup(grid.gridStep());
    projectStencil.setup(grid.gridStep());
    directStencil.setup(grid.gridStep());
    TNT::stopwatch time;  time.start();
    interpMat.setup();
    progressReport("interpMat has been generated");
    timeReport("Time used for filling the interpolation matrix := ", time.read());
    
    time.reset();  time.start();
    projectMat.setup();
    progressReport("projectMat has been generated");
    timeReport("Time used for filling the projection matrix := ", time.read());
    
    progressReport("setup of pfft kernel independent part has been done!");
  }
  
  
  /**********************************************************************
   * constructKernelDependent --
   **********************************************************************/
  template <class Kernel, class Data, class GreenFunc, class Integration>
  void
  Pfft<Kernel, Data, GreenFunc, Integration>::
  constructKernelDependent (
			    const GreenFunc& greenFuncIn,
			    const Integration& integrationIn)
  {
    progressReport("Start to setup kernel dependent part of pfft ......");
    
    greenFunc = greenFuncIn;
    integration = integrationIn;
    
    greenFunc.set_maxthreshold(max(srcGridElementPtr->get_maxThreshold(),
				   evalGridElementPtr->get_maxThreshold()));
    directMat =
      DirectMat<GreenFunc, Integration, Kernel>(srcGridElementPtr,
						evalGridElementPtr,
						&integration, &greenFunc,
						&grid, &directStencil,
						&interpStencil,
						&projectStencil,
						&projectMat, &interpMat);
    
    TNT::stopwatch time;  time.start();

    if ( kernelSize_  == XY_dipoles ) {
      gridData = GridData<Kernel, Data, GreenFunc>(greenFunc);
      gridData.allocate(grid.numPointX(), grid.numPointY(), grid.numPointZ(),
			    grid.gridStep(), 1);
      gridData.fillKernel();
      gridData.fftKernel();

      gridData.fillKerneldiag();
      gridData.fillKerneloffdiag();
      gridData.fftKerneldyadic();
      
    } else if ( kernelSize_  == XYZ_dipoles ) {       

      gridData = GridData<Kernel, Data, GreenFunc>(greenFunc);
      gridData.allocate(grid.numPointX(), grid.numPointY(), grid.numPointZ(),
			    grid.gridStep(), 2);
      gridData.fillKernel();
      gridData.fftKernel();

      gridData.fillKerneldiag();
      gridData.fillKerneloffdiag();
      gridData.fftKerneldyadic();
      
    } else if ((kernelSize_  == XY) || (kernelSize_  == XYZ)) {
      gridData = GridData<Kernel, Data, GreenFunc>(greenFunc);
      gridData.allocate(grid.numPointX(), grid.numPointY(), grid.numPointZ(),
			grid.gridStep(), 0);
      gridData.fillKernel();
      gridData.fftKernel();

    } else if ( kernelSize_  == Scalar_dipoles ) {

      gridData = GridData<Kernel, Data, GreenFunc>(greenFunc);
      gridData.allocateScalar(grid.numPointX(), grid.numPointY(), grid.numPointZ(),
			      grid.gridStep(), 1);
      gridData.fillKernel();
      gridData.fftKernel();
      gridData.fillKernelScalarImages();
      gridData.fftKernelScalarImages();

    } else if  ( kernelSize_  == Scalar ) {

      gridData = GridData<Kernel, Data, GreenFunc>(greenFunc);
      gridData.allocateScalar(grid.numPointX(), grid.numPointY(), grid.numPointZ(),
			      grid.gridStep(), 0);
      gridData.fillKernel();
      gridData.fftKernel();
    }




    progressReport("gridKernel has been filled");
    progressReport("gridKernel has been FFTed");
    timeReport("Time used for filling the convolution matrix := ", time.read());
    timeReport("Time used for FFT the convolution matrix := ", time.read());
    
    time.reset();  time.start();
    directMat.setup();
    progressReport("directMat has been generated and pre-corrected");
    timeReport("Time used for forming and pre-correcting the direct matrix := ", time.read());
    
    constructCompleted_ = true;
    progressReport("setup of pfft kernel dependent part has been done!");
  }
  


  /**********************************************************************
   * IEoperator --
   **********************************************************************/
  template <class Kernel, class Data, class GreenFunc, class Integration>
  template <class VecA, class VecB>
  void
  Pfft<Kernel, Data, GreenFunc, Integration>::
  IEoperator (
	      VecA& des,
	      const DifferentialOperator outerOperator,
	      const DifferentialOperator innerOperator,
	      const VecB& src)
  {
    VecA des1(des.size());

    VecB desx(des.size()), desy(des.size()), desz(des.size());
    VecB srcx(src.size()), srcy(src.size()), srcz(src.size());
    
    if (constructCompleted_ == true) {
      size_t innerOperatorIndex = integration.findInnerOperatorIndex(innerOperator);
      size_t outerOperatorIndex = integration.findOuterOperatorIndex(outerOperator);
      size_t integralIndex = integration.findIntegralIndex(outerOperator, 
							   innerOperator);
      
      TNT::stopwatch time;  time.start();
      
      if (kernelSize_ == Scalar) {
	
	gridData.projectOnGrid(projectMat[innerOperatorIndex], src);
	gridData.conv3D();
	gridData.interpFromGrid(interpMat[outerOperatorIndex], des);
	matMultVecAdd(des, directMat[integralIndex], src);

      	IEOperatorTimeReport(time.read());

      } else if (kernelSize_ == Scalar_dipoles) {

	gridData.projectOnGrid(projectMat[innerOperatorIndex], src);
	gridData.conv3D();
	gridData.interpFromGrid(interpMat[outerOperatorIndex], des);
	matMultVecAdd(des, directMat[integralIndex], src);

	gridData.projectOnGridScalarImages(projectMat[innerOperatorIndex], src);
	gridData.conv3DScalarImages();
	gridData.PermutationImages();
	gridData.interpFromGridScalarImages(interpMat[outerOperatorIndex], des1);
	des += des1;
      	IEOperatorTimeReport(time.read());
      } 

    } else {
      pfft::errorMessage("pfft.h: IEoperator()",
			 "pfft object has not been fully constructed yet! Please call function construct() first.");
    }
  }
  
  
  template <class Kernel, class Data, class GreenFunc, class Integration>
  template <class VecA, class VecB>
  void
  Pfft<Kernel, Data, GreenFunc, Integration>::
  IEoperatorDyadic (
		    VecA& des, VecA& desx,  VecA& desy,  VecA& desz,
		    const DifferentialOperator outerOperator,
		    const DifferentialOperator innerOperator,
		    const VecB& src, const VecB& srcx, const VecB& srcy, const VecB& srcz)
  {
    VecB desx1(desx.size()), desy1(desy.size()),desz1(desz.size());    
    if (constructCompleted_ == true) {
      size_t innerOperatorIndex = integration.findInnerOperatorIndex(innerOperator);
      size_t outerOperatorIndex = integration.findOuterOperatorIndex(outerOperator);
      size_t integralIndex = integration.findIntegralIndex(outerOperator, 
							   innerOperator);      
      TNT::stopwatch time;  time.start();
      if ((kernelSize_ == XYZ_dipoles)||(kernelSize_ == XY_dipoles)) {

	gridData.projectOnGridDyadic(projectMat[innerOperatorIndex], srcx, srcy, srcz);
      	gridData.conv3Dtoeplitz();
      	gridData.interpFromGridDyadic(interpMat[outerOperatorIndex], desx, desy, desz);

	gridData.projectOnGridDyadic(projectMat[innerOperatorIndex], srcx, srcy, srcz);
	gridData.conv3Ddyadic();
	gridData.PermutationDyadic();
      	gridData.interpFromGridDyadic(interpMat[outerOperatorIndex], desx1, desy1, desz1);

	desx += desx1;
	desy += desy1;
	desz += desz1;	

      	matMultVec(des, directMat[integralIndex], src);
      	IEOperatorTimeReport(time.read());

      } /* else if (kernelSize_ == XY_dipoles) {

	gridData.projectOnGridDyadic(projectMat[innerOperatorIndex], srcx, srcy, srcz);
      	gridData.conv3Dtoeplitz();
      	gridData.interpFromGridDyadic(interpMat[outerOperatorIndex], desx, desy, desz);

	gridData.projectOnGridDyadic(projectMat[innerOperatorIndex], srcx, srcy, srcz);
	gridData.conv3Ddyadic();
	gridData.PermutationDyadic();
      	gridData.interpFromGridDyadic(interpMat[outerOperatorIndex], desx1, desy1, desz1);

	desx += desx1;
	desy += desy1;
	desz += desz1;	

      	matMultVec(des, directMat[integralIndex], src);
      	IEOperatorTimeReport(time.read());
	} */ else if ((kernelSize_ == XYZ) || (kernelSize_ == XY)) {

	gridData.projectOnGridDyadic(projectMat[innerOperatorIndex], srcx, srcy, desz);
      	gridData.conv3Dtoeplitz();
      	gridData.interpFromGridDyadic(interpMat[outerOperatorIndex], desx, desy, desz);

      	matMultVec(des, directMat[integralIndex], src);
      	IEOperatorTimeReport(time.read());

      } /* else if (kernelSize_ == XY) {

	gridData.projectOnGridDyadic(projectMat[innerOperatorIndex], srcx, srcy, desz);
      	gridData.conv3Dtoeplitz();
      	gridData.interpFromGridDyadic(interpMat[outerOperatorIndex], desx, desy, desz);

      	matMultVec(des, directMat[integralIndex], src);
      	IEOperatorTimeReport(time.read());
	} */

    } else {
      pfft::errorMessage("pfft.h: IEoperator()",
			 "pfft object has not been fully constructed yet! Please call function construct() first.");
    }
  }
  
  /**********************************************************************
   * IEoperatorWithoutFastConvolution -- for test only
   **********************************************************************/
  template <class Kernel, class Data, class GreenFunc, class Integration>
  template <class VecA, class VecB>
  void
  Pfft<Kernel, Data, GreenFunc, Integration>::
  IEoperatorWithoutFastConvolution (
				    VecA& des,
				    const DifferentialOperator outerOperator,
				    const DifferentialOperator innerOperator,
				    const VecB& src)
  {
    if (constructCompleted_ == true) {
      size_t innerOperatorIndex = integration.findInnerOperatorIndex(innerOperator);
      size_t outerOperatorIndex = integration.findOuterOperatorIndex(outerOperator);
      size_t integralIndex = integration.findIntegralIndex(outerOperator, 
							   innerOperator);

      TNT::stopwatch time;  time.start();
      gridData.projectOnGrid(projectMat[innerOperatorIndex], src);
      gridData.directConvolution();
      gridData.interpFromGrid(interpMat[outerOperatorIndex], des);
      matMultVecAdd(des, directMat[integralIndex], src);
      IEOperatorTimeReport(time.read());

    } else {
      errorMessage("pfft.h: IEoperatorWithoutFastConvolution()",
		   "pfft object has not been fully constructed yet! Please call function construct() first.");
    }
  }

  /**********************************************************************
   * initialize --
   * In find optimal grid size, the member function memoryEstimation
   * of all important objects, interpMat, projectMat, directMat, grid,
   * gridElement, and etc. are to be used. These objects are 
   * instantiated here. And since they depend upon each other, 
   * this function also let them "know each other". They
   * will remember other's name(pointer) and use them in their own setup.
   **********************************************************************/
  template <class Kernel, class Data, class GreenFunc, class Integration>
  void
  Pfft<Kernel, Data, GreenFunc, Integration>::
  initialize (
	      std::vector<element>& sourceElementList,
	      std::vector<element>& targetElementList)
  {
    srcAndEvalDiff = (sourceElementList != targetElementList);

    interpStencil = Stencil(INTERP_STENCIL, interpStencilSize);
    projectStencil = Stencil(PROJECT_STENCIL, projectStencilSize);
    directStencil = Stencil(DIRECT_STENCIL, directStencilSize);

    srcGridElementPtr = new GridElement(&sourceElementList);
    if (srcAndEvalDiff) {
      evalGridElementPtr = new GridElement(&targetElementList);
    } else {
      evalGridElementPtr = srcGridElementPtr;
    }
    //	if (GridCompleted_==false)
    {
      grid = Grid(srcGridElementPtr, evalGridElementPtr, &projectStencil, 
		  &interpStencil, &directStencil); 
    }
    interpMat = InterpMat(&interpStencil, evalGridElementPtr, &integration, 
			  &grid);
    projectMat = ProjectMat(&projectStencil, srcGridElementPtr, &integration, 
			    &grid);
    directMat = 
      DirectMat<GreenFunc, Integration, Kernel>(srcGridElementPtr, 
						evalGridElementPtr, 
						&integration, &greenFunc,
						&grid, &directStencil, 
						&interpStencil, 
						&projectStencil, 
						&projectMat, &interpMat);
    gridData = GridData<Kernel, Data, GreenFunc>(greenFunc);
  }
 

  /**********************************************************************
   * setupGridAndGridElement --
   **********************************************************************/

  template <class Kernel, class Data, class GreenFunc, class Integration>
  void
  Pfft<Kernel, Data, GreenFunc, Integration>::setupGridAndGridElement (
								       void)
  {

    //	if (GridCompleted_==false)
    {
      //  findOptimalNumGridPoint();
      progressReport("optimal grid size has been found");
      grid.setupGridSize();
    }
    
    srcGridElementPtr->mapElementToGridPoint(grid);
    if (srcAndEvalDiff) {
      evalGridElementPtr->mapElementToGridPoint(grid);  
    }
    srcGridElementPtr->setupElementNeighbor(*evalGridElementPtr, grid, 
					    directStencil);
    progressReport("grid and gridElement have been set up");

    srcGridElementPtr->countNumNeighbor(*evalGridElementPtr, grid, directStencil);
  }


  /**********************************************************************
   * findOptimalNumGridPoint --
   * Try to minimize the memory consumption by ajusting grid size or 
   * numbers of grid point along three directions.
   * Essentially this is a try-and-measure process.
   **********************************************************************/
  template <class Kernel, class Data, class GreenFunc, class Integration>
  void
  Pfft<Kernel, Data, GreenFunc, Integration>::findOptimalNumGridPoint (
								       void)
  {
    bool optimalAchieved = false;
    minMemoryUsage = ULLONG_MAX;
    Grid::GridStatus gridStatus = Grid::TOO_COARSE;
    MemoryStatus memoryStatus = IMPROVING;
    while ( (gridStatus != Grid::TOO_FINE) && 
	    (memoryStatus != WORSENING) ) {
      grid.setupGridSize();
      srcGridElementPtr->mapElementToGridPoint(grid);
      if (srcAndEvalDiff) evalGridElementPtr->mapElementToGridPoint(grid);
      gridStatus = grid.checkGridSize();
      if (gridStatus == Grid::SIZE_OK) {
	memoryStatus = checkMemoryUsage();
	if (memoryStatus == IMPROVING) {
	  grid.backupNumPoint();
	  grid.doubleNumPoint();
	  optimalAchieved = true;
	} else if (memoryStatus == WORSENING) {
	  grid.recoverNumPoint();
	} else if (memoryStatus == WORSENING_BUT_KEEP_GOING) {
	  grid.doubleNumPoint();
	}
      } else if (gridStatus == Grid::TOO_COARSE) {
	grid.doubleNumPoint();
      } else if (gridStatus == Grid::TOO_FINE) {
	grid.recoverNumPoint();
      }
    }
    if (! optimalAchieved) {
      warningMessage("pfft.h::findOptimalNumGridPoint()", 
		     "Fail to find the optimal grid! Use the default setting. This is almost equivalent to direct matrix vector product");
    }

  }
  
  /**********************************************************************
   * checkMemoryUsage --
   * 1) check if more memory is consumed
   * 2) check if the memory used by direct, projection and interpolation
   *    matrices (called direct part) is roughly equal to that used by 
   *    the convolution matrix and the working vector that stores the
   *    unknwons on the grid. 
   **********************************************************************/
  template <class Kernel, class Data, class GreenFunc, class Integration>
  typename Pfft<Kernel, Data, GreenFunc, Integration>::MemoryStatus
  Pfft<Kernel, Data, GreenFunc, Integration>::checkMemoryUsage (
								void)
  {
    // hard-coded parameter here
    const double threshold = 3.;

    // If gridData alone consumes more memory than previous betting, then no need to 
    // estimate other objects' memory usage
    size_t memoryUsageByGridData = gridData.memoryUsage(grid.numPointX(), 
							grid.numPointY(), 
							grid.numPointZ());
    if (memoryUsageByGridData > minMemoryUsage) {
      memoryReport(0, memoryUsageByGridData);
      currentMemoryUsage = memoryUsageByGridData;
      return WORSENING;
    }
    
    srcGridElementPtr->countNumNeighbor(*evalGridElementPtr, grid, directStencil);
    size_t memoryUsageByDirectPart = estimateMemoryByDirectPart();    
    memoryReport(memoryUsageByDirectPart, memoryUsageByGridData);
    currentMemoryUsage = memoryUsageByDirectPart + memoryUsageByGridData;
    double ratio = currentMemoryUsage / minMemoryUsage;
    if (ratio < 1.) {    
      minMemoryUsage = currentMemoryUsage;
      return IMPROVING;
    } else if (ratio >= threshold) {
      return WORSENING;
    } else {
      if (memoryUsageByDirectPart > memoryUsageByGridData) {
	return WORSENING_BUT_KEEP_GOING;
      } else {
	return WORSENING;
      }
    } 
  }

  /**********************************************************************
   * estimateMemoryByDirectPart --
   **********************************************************************/
  template <class Kernel, class Data, class GreenFunc, class Integration>
  size_t
  Pfft<Kernel, Data, GreenFunc, Integration>::estimateMemoryByDirectPart (
									  void)
  {
    size_t memByDirect = directMat.memoryEstimation();
    size_t memByProject = projectMat.memoryEstimation();
    size_t memByInterp = interpMat.memoryEstimation();
    size_t memByGridElement = srcGridElementPtr->memoryUsage();
    if (srcAndEvalDiff) {
      bool forSrcGridElement = false;
      memByGridElement += evalGridElementPtr->memoryUsage(forSrcGridElement);
    }
    // stencil
    // vanishingly small

#ifdef MEMORY_REPORT
    double Mb = 1024*1024;
    memoryReport("direct matrix", memByDirect/Mb);
    memoryReport("projection matrix", memByProject/Mb);
    memoryReport("interpolation matrix", memByInterp/Mb);
    memoryReport("gridElement", memByGridElement/Mb);
#endif

    return memByDirect + memByProject + memByInterp + memByGridElement;
  }

  /**********************************************************************
   * memoryReport --
   **********************************************************************/
  template <class Kernel, class Data, class GreenFunc, class Integration>
  void 
  Pfft<Kernel, Data, GreenFunc, Integration>::
  memoryReport (
		const std::string matrixName, 
		const double memory) const
  {
#ifdef MEMORY_REPORT
    std::streamsize prec = std::cout.precision();
    std::cout << setprecision(3)
	      << "memory used by " << matrixName << "(MB) := " << memory
	      << setprecision(prec) << std::endl;
#endif
  }

  /**********************************************************************
   * gridAndElementReport --
   **********************************************************************/
  template <class Kernel, class Data, class GreenFunc, class Integration>
  void
  Pfft<Kernel, Data, GreenFunc, Integration>::
  gridAndElementReport (
			void) const
  {
#ifdef PRINT_ELEMENT
    std::cout << "\tNumber of src elements := " 
	      << srcGridElementPtr->numElement() << std::endl;
    if (srcAndEvalDiff)
      std::cout << "\tNumber of eval elements := " 
		<< evalGridElementPtr->numElement() << std::endl;
    std::cout <<"\tMin src element size := " 
	      << srcGridElementPtr->minElementSize() << std::endl;
    std::cout <<"\tMax src element size := " 
	      << srcGridElementPtr->maxElementSize() << std::endl;
    std::cout <<"\tAverage src element size := "
	      << srcGridElementPtr->averageElementSize()
	      << std::endl;
    std::cout << "\tNumber of grids := (" << grid.numPointX() << "," 
	      << grid.numPointY() << "," << grid.numPointZ() << ")"
	      << "  grid size := " << grid.gridStep() << std::endl;
    std::cout << "\tMax number of elements mapped to one grid point := "
	      << srcGridElementPtr->maxNumElementMappedToOneGrid()
	      << std::endl;
    std::cout << "\tMin number of neighbors to one element := "
	      << srcGridElementPtr->minNumNeighbor()
	      << std::endl;
    std::cout << "\tMax number of neighbors to one element := "
	      << srcGridElementPtr->maxNumNeighbor()
	      << std::endl;
#endif
  }

  /**********************************************************************
   * progressReport --
   **********************************************************************/
  template <class Kernel, class Data, class GreenFunc, class Integration>
  void
  Pfft<Kernel, Data, GreenFunc, Integration>::
  progressReport (
		  const std::string message) const
  {
#ifdef DEBUG_PROGRESS    
    std::cout << std::endl << "\t" << message << std::endl;
#endif
  }

  /**********************************************************************
   * IEOperatorTimeReport --
   **********************************************************************/
  template <class Kernel, class Data, class GreenFunc, class Integration>
  void
  Pfft<Kernel, Data, GreenFunc, Integration>::
  IEOperatorTimeReport (
			double time)
  {
#ifdef PRINT_RUN_TIME
    if (reportMatVecMultTime_) {
      if (time < 0) time = 0.;
      std::cout << "\tTime used for IE operator on a vector := " 
		<< time << " (seconds)" << std::endl;
      reportMatVecMultTime_ = false;
    }
#endif
  }

  /**********************************************************************
   * memoryReport --
   **********************************************************************/
  template <class Kernel, class Data, class GreenFunc, class Integration>
  void 
  Pfft<Kernel, Data, GreenFunc, Integration>::
  memoryReport (
		const size_t memoryUsageByDirectPart, 
		const size_t memoryUsageByGridData) const
  {
#ifdef MEMORY_REPORT
    double md = static_cast<double>(memoryUsageByDirectPart)/1024/1024;
    double mg = static_cast<double>(memoryUsageByGridData)/1024/1024;
    std::streamsize prec = std::cout.precision();
    std::cout << setprecision(3)
	      << "memory usage (MB): DirectPart = " << md
	      << "   GridData = " << mg 
	      << "   total = " << md+mg 
	      << setprecision(prec) << std::endl;
#endif
  }


} //namespace pfft

#endif
