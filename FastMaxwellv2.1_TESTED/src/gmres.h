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

  const static char cvsid[] = "$Id: gmres.h,v 1.10 2003/02/05 21:14:34 zhzhu Exp $";

  ==========================================================================
*/


#ifndef _GMRES_H_
#define _GMRES_H_

#include <complex>
#include "../pfft/cmat.h"
#include "../pfft/vec.h"
#include "sysformation.h"
#include "../lib/SuperLU_3.0/SRC/slu_zdefs.h"
#include "../lib/SuperLU_3.0/SRC/slu_util.h"
//#include "superlu.h"
/**********************************************************************
 * outputResid --
 **********************************************************************/


static void
outputResid (
	     const int outIt,
	     const int inIt,
	     const double resid)
{
#ifdef PRINT_GMRES_RESID
  if (outIt == 0) {
    std::cout << std::endl;
  }
  std::cout << inIt << "\t" 
	    << outIt << "th iter, residual = " << resid << std::endl;

#else
  const int printPoint = 10; 
  if (outIt == 0) {
    std::cout << std::endl;
  } else if ((outIt != 0) && (inIt == 0)) {    
    std::cout << "\touter resid = " << resid << std::endl;
  } else {    
    //    if (inIt == 1)  std::cout << std::endl;
    std::cout << outIt << " ";
    if (outIt % printPoint == 0) {
      std::cout << " inner resid = " << resid << std::endl;
    }
    flush(std::cout);
  }
#endif
}

/**********************************************************************
 * givens --
 **********************************************************************/
template <class T>
static void
givens (
	const T a, 
	const T b, 
	T& c, 
	T& s)
{
  if (b == 0.) {
    c = 1.;
    s = 0.;
  } else {
    double r1 = abs(a);
    double r2 = abs(b);
    double length = sqrt(r1*r1 + r2*r2);
    c = a / length;
    s = -b / length; 
  }
}

/**********************************************************************
 * givensRotate --
 **********************************************************************/
template <class T>
static void
givensRotate (
	      const T c, 
	      const T s, 
	      T& a, 
	      T& b)
{
  T a1 = a;
  T b1 = b;
  a = a1 * conj(c) - b1 * conj(s);
  b = c * b1 + s * a1;
}

/**********************************************************************
 * givensRotate --
 **********************************************************************/
static void
givensRotate (
	      const double c, 
	      const double s, 
	      double& a, 
	      double& b)
{
  double a1 = a;
  double b1 = b;
  a = a1 * c - b1 * s;
  b = c * b1 + s * a1;
}

/**********************************************************************
 * gmres --
 * 
 * On return flag: 
 * true -- converged, false -- diverged
 **********************************************************************/
//template <class Matrix, class Preconditioner, class T>

template <class MatrixVectorProd, class T, class Preconditioner>
bool
gmres (
       MatrixVectorProd &MZMtPfft,
       TNT::Vector<T> &x0, // initial guess and solution.
       const TNT::Vector<T> &b,
       Preconditioner &M,
       int maxiters = 300,
       int restart = 20,
       double tol = 1e-6)
{
  TNT::Vector<T> c(restart + 1);
  TNT::Vector<T> s(restart + 1);
  TNT::Vector<T> g(restart + 1);
  TNT::Vector<T> y(restart + 1);
  TNT::Matrix<T> H(restart+1, restart+1);

  int vectorLength = b.size();

  TNT::Vector<T> AP(vectorLength);
  TNT::Vector<T> P(vectorLength);
  std::vector<TNT::Vector<T> > bv(restart+1, TNT::Vector<T>(vectorLength));
  TNT::Vector<T> xLast(x0);

  int numRow;
  int numCol;
  int numNonZero;

  SuperMatrix A, L, U, X, B;
  doublecomplex *mat;
  int *colIndex, *colStart;
  int *perm_c, *perm_r, *etree;
  doublecomplex *superLUrhs, *superLUsol;
  double         rpg, rcond;
  
  SuperLUStat_t stat;
  
  double         *R, *C;
  double         *ferr, *berr;
  mem_usage_t    mem_usage;
  char           equed[1];
  
  
  //vector<T> AP_(vectorLength), M_AP_(vectorLength);
  //SuperLU superLU;
  
  //SuperLUsetup(&superLU,M);
  
  
  AP=MZMtPfft->MatrixVectorProduct_C(x0);
  //cout << AP << endl;
  //cout << "i am in step 3" << endl;
  //	matMultVec(AP,MZMtPfft,x0);
  AP=b-AP;
  
  numRow = M.numRow();
  numCol = M.numCol();
  numNonZero = M.numNonZero();
  
  mat = (doublecomplex *) doublecomplexMalloc(numNonZero);
  colIndex = (int *) intMalloc(numNonZero );
  colStart = (int *) intMalloc(numCol+1);
  
  int elementCount = 0;
  for (size_t col = 0; col < numCol; col++)
    {
      colStart[col] = elementCount;
      for (size_t row = 0; row < M[col].size(); row++)
	{
	  mat[elementCount].r = M[col].value(row).real();
	  mat[elementCount].i = M[col].value(row).imag();
	  colIndex[elementCount] = M[col].index(row);
	  elementCount++;
	}
    }
  colStart[M.numCol()] = elementCount;
  
  //create matrix A in the format expected by SuperLU
  zCreate_CompCol_Matrix(&A, numRow, numCol, numNonZero, mat, colIndex,
                         colStart, SLU_NC, SLU_Z, SLU_GE);
  

  int info;
  int lwork=0;
  void *work;
  double drop_tol = 0.;
  int relax      = sp_ienv(2);
  int panel_size = sp_ienv(1);
  
  //  if(!memoryAllocated) allocate();
  superLUrhs = (doublecomplex*) doublecomplexMalloc(numRow); 
  if ( !(perm_c = intMalloc(numCol)) ) ABORT("Malloc fails for perm_c[].");
  if ( !(perm_r = intMalloc(numRow)) ) ABORT("Malloc fails for perm_r[].");
  if (! (etree=intMalloc(numCol))) ABORT("Malloc fails for etree[].");
  if ( !(R = (double *) SUPERLU_MALLOC(A.nrow * sizeof(double))) )
    ABORT("SUPERLU_MALLOC fails for R[].");
  if ( !(C = (double *) SUPERLU_MALLOC(A.ncol * sizeof(double))) )
    ABORT("SUPERLU_MALLOC fails for C[].");
  if ( !(ferr = (double *) SUPERLU_MALLOC(sizeof(double))) )
    ABORT("SUPERLU_MALLOC fails for ferr[].");
  if ( !(berr = (double *) SUPERLU_MALLOC(sizeof(double))) )
    ABORT("SUPERLU_MALLOC fails for berr[].");
  //  memoryAllocated = true;
  /* Initialize the statistics variables. */
  StatInit(&stat);
  
  superlu_options_t options;
  set_default_options(&options);
  options.DiagPivotThresh=1.;
  
  options.PivotGrowth = YES;    /* Compute reciprocal pivot growth */
  options.ConditionNumber = YES;/* Compute reciprocal condition number */
  options.IterRefine = DOUBLE;  /* Perform single-precision refinement */
  
  int permc_spec = 2;
  get_perm_c(permc_spec, &A, perm_c);
  
  SuperMatrix AC;
  sp_preorder(&options, &A, perm_c, etree, &AC);
  
  if ( !(superLUsol = doublecomplexMalloc(numRow)) ) 
    ABORT("Malloc fails for superLUsol[].");
  if ( !(superLUrhs = doublecomplexMalloc(numRow)) ) 
    ABORT("Malloc fails for superLUrhs[].");
  for (size_t i = 0; i < numRow; i++){
    superLUrhs[i].r = AP[i].real();
    superLUrhs[i].i = AP[i].imag();
  }
  zCreate_Dense_Matrix(&B, numRow, 1, superLUrhs, numRow, SLU_DN, SLU_Z, SLU_GE);
  zCreate_Dense_Matrix(&X, numRow, 1, superLUsol, numRow, SLU_DN, SLU_Z, SLU_GE);
  
  /* ONLY PERFORM THE LU DECOMPOSITION */
  B.ncol = 0;  /* Indicate not to solve the system */
  
  zgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C,
	 &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
	 &mem_usage, &stat, &info);
  
  
  /* ------------------------------------------------------------
     NOW WE SOLVE THE LINEAR SYSTEM USING THE FACTORED FORM OF A.
     ------------------------------------------------------------*/
  options.Fact = FACTORED; /* Indicate the factored form of A is supplied. */
  B.ncol = 1;  /* Set the number of right-hand side */
  
  /* Initialize the statistics variables. */
  StatInit(&stat);
  
  zgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C,
	 &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
	 &mem_usage, &stat, &info);
  
  superLUsol = (doublecomplex*) ((DNformat*) X.Store)->nzval;
  
  for (size_t i = 0; i < numRow; i++)
    AP[i] = std::complex<double>(superLUsol[i].r, superLUsol[i].i);
  
  //cout << AP << endl;


  double outer_rnorm = two_norm(AP); 
  double bnorm = outer_rnorm; // assuming the initial guess is zero
  //cout << bnorm << endl;  
  int outIt = 0;
  outputResid(outIt, 0, outer_rnorm/bnorm);
  
  while (outIt < maxiters && outer_rnorm/bnorm > tol) {
    
    g[0] = outer_rnorm;
    AP /= outer_rnorm;
    bv[0] = AP;
    
    int inIt = 0;
    double inner_rnorm;
    do {
      
      //matMultVec(AP,MZMtPfft,x0);
      // if (MZMtPfft->Solver.solverTypeIn) matMultVec(AP,MZMtPfft,x0);
      //	if (MZMtPfft->Solver.formulationType==FORMULATION_A)
      //		AP=MZMtPfft->MatrixVectorProduct_A(bv[inIt]);
      //	else if (MZMtPfft->Solver.formulationType==FORMULATION_C)
      AP=MZMtPfft->MatrixVectorProduct_C(bv[inIt]);
      // tarek for beta release
      //for (int i=0; i<AP.size(); i++) {
      //	AP_[i]=AP[i];
      //}
      //  SuperLUsolve(M, AP, AP); // AP = precondition(P);
      
      if ( !(superLUsol = doublecomplexMalloc(numRow)) ) ABORT("Malloc fails for superLUsol[].");
      if ( !(superLUrhs = doublecomplexMalloc(numRow)) ) ABORT("Malloc fails for superLUrhs[].");
      for (size_t i = 0; i < numRow; i++){
	superLUrhs[i].r = AP[i].real();
	superLUrhs[i].i = AP[i].imag();
      }
      zCreate_Dense_Matrix(&B, numRow, 1, superLUrhs, numRow, SLU_DN, SLU_Z, SLU_GE);
      zCreate_Dense_Matrix(&X, numRow, 1, superLUsol, numRow, SLU_DN, SLU_Z, SLU_GE);
      options.Fact = FACTORED; /* Indicate the factored form of A is supplied. */
      B.ncol = 1;  /* Set the number of right-hand side */
      
      /* Initialize the statistics variables. */
      StatInit(&stat);
      
      zgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C,
	     &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
	     &mem_usage, &stat, &info);
      
      superLUsol = (doublecomplex*) ((DNformat*) X.Store)->nzval;
      
      for (size_t i = 0; i < numRow; i++)
	AP[i] = std::complex<double>(superLUsol[i].r, superLUsol[i].i);
      
      //    SuperLUsolve(&superLU, AP, AP); // AP = precondition(P);
      //for (int i=0; i<AP.size(); i++) {
      //	AP[i]=M_AP_[i];
      //}
      
      
      //	cout << AP.size() << endl;
      //	for (int i=0; i<AP.size(); i++) {
      //			cout << AP[i] <<  endl; }
      for (int j = 0; j <= inIt; j++) {
	H[inIt][j] = inner_prod(AP, bv[j]);
	AP -= bv[j] * H[inIt][j];
	// cout << H[inIt][j] << endl;
      }
      H[inIt][inIt+1] = two_norm(AP);   
      bv[inIt+1] = AP / H[inIt][inIt+1];
      
      for (int i = 0; i < inIt; i++) {
	givensRotate(c[i], s[i], H[inIt][i], H[inIt][i+1]);
      }
      givens(H[inIt][inIt], H[inIt][inIt+1], c[inIt], s[inIt]);
      givensRotate(c[inIt], s[inIt], H[inIt][inIt], H[inIt][inIt+1]);
      g[inIt+1] = 0.; // this should be zero
      givensRotate(c[inIt], s[inIt], g[inIt], g[inIt+1]);
      inner_rnorm = abs(g[inIt+1]);
      
      inIt++;
      outIt++;
      
      outputResid(outIt, inIt, inner_rnorm/bnorm);
      
    } while (inIt < restart && inner_rnorm/bnorm > tol);
    
    // it is incremented by one before exit. Adjustment is needed here
    inIt --;
    
    // Compute solution, note, H is H[col][row]. 
    // backsolve for y: H*y = g
    for (int k = 0; k <= inIt; ++k) {
      y[k] = g[k];
    }
    for(int i = inIt; i >= 0; i--) {
      //cout << H[i][i] << endl;
      y[i] /= H[i][i];
      for(int j = i-1; j >= 0; j--) {
	y[j] -= H[i][j] * y[i]; 
      }
    }
    
    // x = x0 + bv * y
    for(int j=0; j <= inIt; j++) {
      x0 +=  bv[j] * y[j];
    }
    
    // A.matMultVec(AP, x0); // AP = A * x0;
    
    
    //matMultVec(AP,MZMtPfft,x0);
    //	if (MZMtPfft->Solver.formulationType==FORMULATION_A)
    //		AP=MZMtPfft->MatrixVectorProduct_A(x0);
    //	else if (MZMtPfft->Solver.formulationType==FORMULATION_C)
    AP=MZMtPfft->MatrixVectorProduct_C(x0);
    
    AP=b-AP;
    // tarek for beta release
    
    //for (int i=0; i<AP.size(); i++) {
    //	AP_[i]=AP[i];
    //}
    //    SuperLUsolve(M, AP, AP); // AP = precondition(P);
    if ( !(superLUsol = doublecomplexMalloc(numRow)) ) ABORT("Malloc fails for superLUsol[].");
    if ( !(superLUrhs = doublecomplexMalloc(numRow)) ) ABORT("Malloc fails for superLUrhs[].");
    for (size_t i = 0; i < numRow; i++){
      superLUrhs[i].r = AP[i].real();
      superLUrhs[i].i = AP[i].imag();
    }
    zCreate_Dense_Matrix(&B, numRow, 1, superLUrhs, numRow, SLU_DN, SLU_Z, SLU_GE);
    zCreate_Dense_Matrix(&X, numRow, 1, superLUsol, numRow, SLU_DN, SLU_Z, SLU_GE);
    options.Fact = FACTORED; /* Indicate the factored form of A is supplied. */
    B.ncol = 1;  /* Set the number of right-hand side */
    
    /* Initialize the statistics variables. */
    StatInit(&stat);
    
    zgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C,
           &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
           &mem_usage, &stat, &info);
    
    superLUsol = (doublecomplex*) ((DNformat*) X.Store)->nzval;
    
    for (size_t i = 0; i < numRow; i++)
      AP[i] = std::complex<double>(superLUsol[i].r, superLUsol[i].i);
    
    
    //  SuperLUsolve(&superLU, AP, AP); // AP = precondition(P);
    //for (int i=0; i<AP.size(); i++) {
    //	AP[i]=M_AP_[i];
    //}
    
    outer_rnorm = two_norm(AP); 
    outputResid(outIt, 0, outer_rnorm/bnorm);
  }
  
  std::cout << "\n\t Number of Gmres iterations := " << outIt << std::endl;
  std::cout << "\t Final residual := " << outer_rnorm/bnorm << std::endl;
  StatFree(&stat);
  SUPERLU_FREE (etree);
  SUPERLU_FREE (perm_r);
  SUPERLU_FREE (perm_c);
  SUPERLU_FREE (R);
  SUPERLU_FREE (C);
  SUPERLU_FREE (ferr);
  SUPERLU_FREE (berr);
  // SUPERLU_FREE (superLUrhs);
  // SUPERLU_FREE (superLUsol);
  // SUPERLU_FREE (mat);
  Destroy_CompCol_Matrix(&A);
  Destroy_SuperNode_Matrix(&L);
  Destroy_CompCol_Matrix(&U);
  Destroy_SuperMatrix_Store(&B);
  Destroy_SuperMatrix_Store(&X);
  Destroy_CompCol_Permuted(&AC);
  if (outer_rnorm/bnorm < tol) {
    return true;
  } else {
    return false;
  }
}



#endif

    
