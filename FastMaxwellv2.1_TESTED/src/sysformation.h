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
#ifndef SYSFORMATION_H
#define SYSFORMATION_H


#include "structure.h"
#include "generate_dipoles.h"
//#include "spRowMat.h"
//#include "spColMat.h"
#include "layerenv.h"
//#include "henry.h"
#include "resistance.h"
//#include "superlu.h"
#include "inductance.h"
#include "charge.h"
#include "../pfft/element.h"
#include "../pfft/eikrOverR_IMAGES.h"
#include "../pfft/eikrOverR.h"
#include "../pfft/fullwaveSurfaceGalerkin.h"
#include "../pfft/fullwaveVolumeGalerkin.h"
#include "../pfft/pfft.h"
#include "../pfft/spRowMat.h"
//#include "../pfft/vector3D.h"


class sysFormation
{
public:

  sysFormation(Structure& strcture,
               const LayerEnv& layerenv,
               const SimulationType simuType, const bool debug);
  
  void setup_system(const Dipoles& dipoles,
                    const double freq);
  
  void solve_system_direct(vector<complex<double> >& Im, 
			   vector<complex<double> >& Ib,
			   complex<double>& adm,
			   complex<double>& Rsource);
  
  void Solver(const Dipoles& dipoles,
	      const double freq, const SimulationType simuType, 
	      const FormulationType formulationType,
	      const SolverType solverTypeIn);
  
  
  
  TNT::Vector<complex<double> > MatrixVectorProduct_C(TNT::Vector<complex<double> > RHS);
  
  //formulation
  
 private:
  void assemble_M_mat(void);
  void fill_RMat(void);
  void fill_LMat(void);
  void fill_PMat(void);
  
  void setup_freq_info(double freq);
  void setup_loads(vector<complex<double> >& zLoads);
  void setup_ZEM(const vector<complex<double> >& zLoads);
  void M_mult_ZEM_mult_MT (pfft::SpColMat<std::complex<double> >& mZEMmt,
                           const pfft::SpColMat<std::complex<double> >& ZEM,
                           const pfft::SpRowMat<double>& M);
  
  // for pfft:
  void setup_RHS(const vector<complex<double> >& zLoads);
  
  void pfftSolver(const Dipoles& dipoles,
		    const double freq,
		    const SimulationType simuType);
  
  void PreConditionerSetup(void);
  void setup_M_Elements(void);
  double get_shift(vector<pfft::element> ElementList);
  
  //void setup_pfft_C(const Dipoles& dipoles,
  //	    const double freq,
  //	    const SimulationType simuType);
  void setup_Pfft(const pfft::KernelSize, const pfft::KernelSize, 
		  const SimulationType simuType);
  
  TNT::Vector<complex<double> > setup_RHS_InitialGuess();
  
  
  template<class T>  
    void insert_subMatrix_spColMat(pfft::SpColMat<T>& systemMat,
				   const int startRow,
				   const int startCol,
				   const pfft::SpRowMat<T>&  mat);
  
  template<class T>  
    void insert_subMatrix_spColMat(pfft::SpColMat<T>& systemMat,
				   const int startRow,
				   const int startCol,
				   const pfft::SpColMat<T>&  mat);
  
  
  bool EMQS_flag, dipole_flag, debug_flag;
  
  //pfft
  typedef pfft::EikrOverRIMAGES<complex<double> > DynamicGreenFunc;
  typedef pfft::FullwaveSurfaceGalerkin<double, pfft::element, complex<double> > CalcPfullwave3D;
  typedef pfft::FullwaveVolumeGalerkin<double, pfft::element, complex<double> > CalcAfullwave3D;
  typedef pfft::Pfft<std::complex<double>, std::complex<double>, DynamicGreenFunc, CalcPfullwave3D> PfftPanel;
  typedef pfft::Pfft<std::complex<double>, std::complex<double>, DynamicGreenFunc, CalcAfullwave3D> PfftA;
  
  vector<pfft::element> ElementListX, ElementListY, ElementListZ, ElementList;
  vector<pfft::element> ElementListP, ElementListPH, ElementListPV;
  vector<pfft::DifferentialOperator> innerOperatorList;
  vector<pfft::DifferentialOperator> innerOperatorListX;
  vector<pfft::DifferentialOperator> innerOperatorListY;
  vector<pfft::DifferentialOperator> outerOperatorList;
  vector<pfft::DifferentialOperator> outerOperatorListX;
  vector<pfft::DifferentialOperator> outerOperatorListY;
  // pfft::DyadicGreenFunction DyadicGreenFunc;
  
  DynamicGreenFunc kernel_XX, kernel_YY, kernel_XX_MOD, kernel_YY_MOD, 
    kernel_XY, kernel_YX, kernel_ZZ, kernel_XZ, kernel_YZ, kernel_P;
  
  PfftA  pfftDyadic;
  PfftPanel pfftScalar;

  CalcAfullwave3D calcA, calcADyadic;
  CalcPfullwave3D calcP, calcPScalar;
  
  pfft::SpRowMat<double> Mpfft;
  pfft::SpRowMat<double> MTpfft;
  pfft::SpColMat<std::complex<double> > MpreCondMatMT;
  
  
  bool VerticalFilFlag;
  
  //freq-independent matrices
  pfft::SpRowMat<std::complex<double> >  LMat;
  pfft::SpRowMat<std::complex<double> > PMat;
  TNT::Vector<double> RVec;
  pfft::SpRowMat<double> RMat;
  pfft::SpRowMat<double> M;
  pfft::SpVec<complex<double> > Vs;
  vector<complex<double> > zLoads;
  
  //freq-dependent matrices
  pfft::SpColMat<std::complex<double> > ZEM;
  pfft::SpColMat<std::complex<double> > mZEMmt;
  vector<std::complex<double> >  Vm;
  TNT::Vector<complex<double> > RList;
  vector<point3D<double> > DirectionMatrix;
  
  Structure& condStruc;
  const LayerEnv& layers;
  //SuperLU superLU;
  //SuperLU superLU_PreCond;
  Inductance inductance;
  Charge chg;
  const SimulationType simuType_;
  int numPanels, numFilaments;
  int numLayers_;
  int port_location_;
  double omega;
  complex<double> k0;
  
};



template<class T>
void sysFormation::insert_subMatrix_spColMat(pfft::SpColMat<T>& systemMat,
					     const int startRow,
					     const int startCol,
					     const pfft::SpRowMat<T>&  mat)
{
  T value;
  for (size_t rowIndex = 0; rowIndex < mat.numRow(); rowIndex++)
    {
      for (size_t colIndex = 0; colIndex < mat[rowIndex].size(); colIndex++)
	{
	  value = mat[rowIndex].value(colIndex);
	  systemMat.insertElement(startRow+rowIndex,
				  startCol+mat[rowIndex].index(colIndex), value);
	}
    }
}


template<class T>
void sysFormation::insert_subMatrix_spColMat(pfft::SpColMat<T>& systemMat,
					     const int startRow,
					     const int startCol,
					     const pfft::SpColMat<T>&  mat)
{
  T value;
  for (size_t colIndex = 0; colIndex < mat.numCol(); colIndex++)
    {
      for (size_t rowIndex = 0; rowIndex < mat[colIndex].size(); rowIndex++)
	{
	  value = mat[colIndex].value(rowIndex);
	  systemMat.insertElement(startRow+mat[colIndex].index(rowIndex),
				  startCol+colIndex, value);
	}
    }
}


/**********************************************************************
 * fillDenseMatrix
 **********************************************************************/
template <class CalcpXX, class T>
void fillDenseMatrix (
		      CalcpXX& calcpXX,
		      const std::vector<pfft::element>& srcElementList,
		      const std::vector<pfft::element>& evalElementList,
		      TNT::Matrix<T>& mat)
{
  const size_t nRow = evalElementList.size();
  const size_t nCol = srcElementList.size();
  for (size_t rowIdx = 0; rowIdx < nRow; rowIdx ++) {
    for (size_t colIdx = 0; colIdx < nCol; colIdx ++) {
      calcpXX(srcElementList[colIdx], evalElementList[rowIdx]);
      mat[rowIdx][colIdx] = calcpXX.result(0);
    }
  }
}

/**********************************************************************
 * randomSourceGenerator --
 **********************************************************************/
template <class T>
void
randomSourceGenerator (
		       TNT::Vector<T>& v)
{
  for (size_t ii = 0; ii < v.size(); ii ++) {
    v[ii] = static_cast<T>(rand() / double(RAND_MAX));
  }
}

#endif


