/***************************************************************************
*   Copyright (C) 2005 by Xin Hu   * *   xinhu@mit.edu   * * * * This
program is free software; you can redistribute it and/or modify * *
it under the terms of the GNU General Public License as published by
* *   the Free Software Foundation; either version 2 of the License,
or     * *   (at your option) any later version.  * * * * This program
is distributed in the hope that it will be useful, * * but WITHOUT ANY
WARRANTY; without even the implied warranty of * * MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the * * GNU General Public
License for more details.  * * * *   You should have received a copy
of the GNU General Public License     * *   along with this program;
if not, write to the * *   Free Software Foundation, Inc., * *   59
Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*
***************************************************************************/
#include "sysformation.h" 
#include "filament.h" 
#include "const_def.h"
#include "gmres.h" 
#include "../pfft/cmat.h" 
#include "../pfft/vec.h"


sysFormation::sysFormation(Structure& structure, const LayerEnv&
                           layerenv, const SimulationType simuType,
                           const bool debug)
  : condStruc(structure),
    layers(layerenv), simuType_(simuType),
    EMQS_flag(false), dipole_flag(false),
    debug_flag(debug) {
  
  pfft::KernelSize DyadSize, ScalarSize;

  numPanels=condStruc.num_panels();
  numFilaments=condStruc.num_filaments();
  numLayers_=layers.num_layers();
  
  VerticalFilFlag = false;
  
  if ((simuType_==EMQS) && (numLayers_==1)) { 
    EMQS_flag=true;
    dipole_flag=false; 
  } 
  else if ((simuType_==EMQS) && (numLayers_>1)) 
    {
      EMQS_flag=true; 
      dipole_flag=true; 
    } 
  else if ((simuType_==FULLWAVE) && (numLayers_==1)) 
    { 
      EMQS_flag=false; 
      dipole_flag=false; 
    } 
  else if ((simuType_==FULLWAVE) && (numLayers_>1)) 
    { 
      EMQS_flag=false;
      dipole_flag=true; 
    } 
  else 
    cout<<"unknown options "<<endl;
  
  //initialized inductance and charge calculations
  inductance.init(EMQS_flag,dipole_flag,layers.get_epir_topmostlayer()*EPI0);
  chg.init(EMQS_flag, dipole_flag, condStruc, layers.get_epir_topmostlayer()*EPI0);

  
  
  //if EMQS mode is used and there is no substrate, precompute the L
  //and P matrices which are freq-independent.  
  //if ((EMQS_flag==true) && (dipole_flag==false)) { 
  // TAREK MOSELHY 
  //   fill_LMat(); 
  //   fill_PMat(); 
  //}
 
  //R and M matrices should only be computed as well
  //
  fill_RMat();
  assemble_M_mat();

  // Tarek Moselhy November 09, 2006
  setup_M_Elements();
  if (numLayers_ == 1) {
    ScalarSize = pfft::Scalar; 
    if (VerticalFilFlag == false) {
      DyadSize = pfft::XY;
    } else {
      DyadSize = pfft::XYZ;
    }
  } else {
    ScalarSize = pfft::Scalar_dipoles;
    if (VerticalFilFlag == false) {
      DyadSize = pfft::XY_dipoles;
    } else {
      DyadSize = pfft::XYZ_dipoles;
    }
  } 

  setup_Pfft(DyadSize, ScalarSize, simuType_);
 
  if ((EMQS_flag==true) && (dipole_flag==false)) {
    Dipoles dipoles (layers);
    calcADyadic = CalcAfullwave3D(outerOperatorList, innerOperatorList, 0., dipole_flag,
				  dipoles, pfft::GXXc, simuType);
    DynamicGreenFunc kernel_XX(0., dipoles, pfft::GXXc);
    kernel_XX.set_shift(get_shift(ElementList));
    pfftDyadic.constructKernelDependent(kernel_XX, calcADyadic); 
    
    calcPScalar = CalcPfullwave3D(outerOperatorList, innerOperatorList, 0.,
				  dipole_flag, dipoles, pfft::GScalar, simuType);
    DynamicGreenFunc kernel_P(0., dipoles, pfft::GScalar);
    
    //kernel_P=DynamicGreenFunc(0., dipoles, pfft::GZZ_phi);
    kernel_P.set_shift(get_shift(ElementListP));
    pfftScalar.constructKernelDependent(kernel_P, calcPScalar);
  }


}


/******************************************************************* *
Assemble M matrix for nodal analysis
*******************************************************************/
void sysFormation::assemble_M_mat() 
{ 

  pfft::SpRowMat<double> Mf, Mfs, Mp, Mps, Mrs; 
  sparse sVs; 
  pfft::SpVec<double> Vs_temp;
  
  int err;
  
  const int Mf_nRow=(condStruc.mesh_mats).Mf_ptr.numRow(); 
  const int Mf_nCol=(condStruc.mesh_mats).Mf_ptr.numCol(); 
  const int Mfs_nRow=(condStruc.mesh_mats).Mfs_ptr.numRow(); 
  const int Mfs_nCol=(condStruc.mesh_mats).Mfs_ptr.numCol(); 
  const int Mp_nRow=(condStruc.mesh_mats).Mp_ptr.numRow(); 
  const int Mp_nCol=(condStruc.mesh_mats).Mp_ptr.numCol(); 
  const int Mps_nRow=(condStruc.mesh_mats).Mps_ptr.numRow(); 
  const int Mps_nCol=(condStruc.mesh_mats).Mps_ptr.numCol(); 
  const int Mrs_nRow=(condStruc.mesh_mats).Mrs_ptr.numRow(); 
  const int Mrs_nCol=(condStruc.mesh_mats).Mrs_ptr.numCol(); 
  const int Vs_nRow=(condStruc.mesh_mats).Vs_ptr.size();
  
 
 // const int M_nRow=Mf_nRow+Mfs_nRow+Mrs_nRow+Mp_nRow-Vs_nRow; const
  int M_nRow=Mf_nRow+Mfs_nRow+Mrs_nRow+Mp_nRow; 
  const int M_nCol=Mfs_nCol+Mrs_nCol+Mp_nCol;


  
  //copy into M_mat 
  /* 
     M=pfft::SpRowMat<double>(M_nRow, numFilaments+Mrs_nRow+numPanels); 
     insert_subMatrix(M, 0, 0, (condStruc.mesh_mats).Mf_ptr); 
     insert_subMatrix(M, Mf_nRow,0, (condStruc.mesh_mats).Mfs_ptr); 
     insert_subMatrix(M, Mf_nRow+Mfs_nRow, numFilaments, (condStruc.mesh_mats).Mrs_ptr);
     insert_subMatrix(M, Mf_nRow, numFilaments+Mrs_nRow, (condStruc.mesh_mats).Mps_ptr); 
     insert_subMatrix(M, Mf_nRow+Mfs_nRow+Mrs_nRow, numFilaments+Mrs_nRow, (condStruc.mesh_mats).Mp_ptr); */
  
  //copy into M_mat 
  M=pfft::SpRowMat<double>(M_nRow, numFilaments+numPanels+Mrs_nRow); 
  insert_subMatrix(M, 0, 0,(condStruc.mesh_mats).Mf_ptr); 
  insert_subMatrix(M, Mf_nRow,0,(condStruc.mesh_mats).Mfs_ptr); 
  insert_subMatrix(M, Mf_nRow+Mfs_nRow, numFilaments+numPanels, (condStruc.mesh_mats).Mrs_ptr); 
  insert_subMatrix(M, Mf_nRow, numFilaments, (condStruc.mesh_mats).Mps_ptr); 
  insert_subMatrix(M, Mf_nRow+Mfs_nRow+Mrs_nRow, numFilaments, (condStruc.mesh_mats).Mp_ptr);
  
  //debug_flag = 1; 
  if(debug_flag)
    output_to_file_real("M.mat", M); port_location_ = Mf_nRow+Mfs_nRow;
  
}


/******************************************************************* *
filling R, L and C matrices
*******************************************************************/
void sysFormation::fill_RMat() 
{ 
double R; filament<double> *f; int i;
 RVec = TNT::Vector<double> (numFilaments); 
 
 RMat=pfft::SpRowMat<double>(numFilaments, numFilaments); 
 for (f=condStruc.list_filaments.initialize(), i=0; f!=NULL;
      f=condStruc.list_filaments.iterate()) 
   {
     R=f->filLength()/(f->sigma()*f->filCrossA());
     RVec[i] = R; 
     RMat.insertElement(i,i,R); i++; 
   } 
 // debug_flag = 1; 
 if (debug_flag)
   output_to_file_real("R.mat", RMat);
 
}


void sysFormation::fill_LMat() 
{ 
  const filament<double> *f1, *f2;
  int i,j;
  
  LMat=pfft::SpRowMat<complex<double> >(numFilaments, numFilaments);
  for (f1=condStruc.list_filaments.initialize(), i=0; f1!=NULL;
       f1=condStruc.list_filaments.iterate()) 
    { 
      for (f2=condStruc.list_filaments.initialize(2),j=0; f2!=NULL;
	   f2=condStruc.list_filaments.iterate(2)) { 
	inductance.fill_L( LMat, i, j, f1, f2); 
	j++; 
      } 
      i++; 
    } 
  // debug_flag = 1; 
  if (debug_flag)
    output_to_file_complex("L.mat", LMat); 
}


void sysFormation::fill_PMat() 
{ 
  panel<double> *p1, *p2; int i,j;
  
  PMat=pfft::SpRowMat<complex<double> >(numPanels, numPanels); 
  for (p1=condStruc.list_panels.initialize(), i=0; p1!=NULL;
       p1=condStruc.list_panels.iterate()) { 
    for (p2=condStruc.list_panels.initialize(2), j=0; p2!=NULL;
	 p2=condStruc.list_panels.iterate(2)) { 
      chg.fill_p(PMat, i,j, p1, p2); 
      j++; 
    } 
    i++; 
  } 
  // debug_flag = 1; 
  if (debug_flag)
    output_to_file_complex("P.mat", PMat); 
}


/******************************************************************* *
Setting up system by providing freq info
*******************************************************************/
void sysFormation::setup_system(const Dipoles& dipoles, const double
freq) 
{ 
  setup_freq_info(freq); 
  setup_loads(zLoads); //if requiring full-wave or there is a substrate 
  //then compute L and P matrices for each freq point 
  if  (( dipole_flag==true)||(EMQS_flag==false)) {
    inductance.set_dipoles(dipoles, freq); 
    chg.set_dipoles(dipoles, freq);
    fill_PMat(); 
    cout<<"finished PMat"<<endl; 
    fill_LMat(); 
    cout<<"finished LMat"<<endl; } 
  setup_ZEM(zLoads); 
  cout<<"finshed setting up zem"<<endl; 
  M_mult_ZEM_mult_MT (mZEMmt, ZEM, M); 
  cout<<"finshed multiplying MZemMT"<<endl; 
  setup_RHS(zLoads); 
  cout<<"System matrix of size "<<mZEMmt.numRow()<< " by "<<mZEMmt.numCol()<<" has been setup successfully"<<endl;
  
}

void sysFormation::setup_freq_info (double freq) 
{ 
  omega=2*PI*freq;
  k0=layers.get_k_topmostLayer(); 
}

//void sysFormation::setup_loads(vector<complex<double> >& zLoads) {
//Ext_Node_Pair *ext; Resistance *resist;; complex<double> Z;
//
//  //condStruc.find_load_resistance(); //obtaining loads by examining
//  the resistance attached to //each port for
//  (resist=condStruc.list_resistances.initialize(); resist!=NULL;
//  resist=condStruc.list_resistances.iterate()) // { //rList=
//  ext->get_resistance_list(); //assuming all z values are R; no C or
//  L as loads Z=resist->get_r_value(); zLoads.push_back(Z); } //}

void sysFormation::setup_loads(vector<complex<double> >& zLoads) 
{
  Ext_Node_Pair *ext; vector<Resistance*> rList; complex<double> Z;

  condStruc.find_load_resistance(); 
//obtaining loads by examining the  resistance attached to 
//each port 
  for (ext=condStruc.list_ext_node_pairs.initialize(); ext!=NULL; ext=condStruc.list_ext_node_pairs.iterate()) { 
    rList=  ext->get_resistance_list(); 
    Z=0; 
    //assuming all z values are R; no C or L as loads 
    for (int i=0; i<rList.size(); i++)
      Z=Z+rList[i]->get_r_value(); zLoads.push_back(Z); 
  } 
}


void sysFormation::setup_ZEM(const vector<complex<double> >& zLoads) 
{
  pfft::SpRowMat<std::complex<double> > sL, Ps, loadMat, RsL;
  sL=pfft::SpRowMat<complex<double> >(numFilaments, numFilaments);
  Ps=pfft::SpRowMat<complex<double> >(numPanels, numPanels);
  ZEM=pfft::SpColMat<std::complex<double> >(M.numCol(),M.numCol()); 
  RsL=pfft::SpRowMat<std::complex<double> >(numFilaments, numFilaments);

  int numLoads=M.numCol()-(sL.numRow()+Ps.numRow());
  loadMat=pfft::SpRowMat<std::complex<double> >(numLoads, numLoads);
  

  matMultScalar(sL, LMat, IMAG*omega); matDivScalar(Ps, PMat,IMAG*omega); 
  matAdd(RsL, sL, RMat);

  
  
  for (int i=0; i<numLoads; i++) loadMat.insertElement(i,i,zLoads[i]);
  insert_subMatrix_spColMat(ZEM,0,0,RsL);
  insert_subMatrix_spColMat(ZEM, sL.numRow(), sL.numRow(), loadMat);
  insert_subMatrix_spColMat(ZEM, sL.numRow()+numLoads,
			    sL.numRow()+numLoads, Ps);
  
  // debug_flag = 1; 
  //if (debug_flag) 
  // output_to_file_complex("P.mat", PMat);
  
  FILE *fp, *fp2; fp = fopen ("ZEM_real.mat", "w"); 
  fp2 = fopen("ZEM_imag.mat", "w"); 
  for (size_t col=0; col<ZEM.numCol();col++) { 
    for (size_t row=0; row<ZEM[col].size(); row++) {
      fprintf(fp, "%i %i %e ", ZEM[col].index(row)+1, col+1,
	      real(ZEM[col].value(row) )); 
      fprintf(fp2, "%i %i %e ", ZEM[col].index(row)+1, col+1,  imag(ZEM[col].value(row)));
      fprintf(fp, "\n"); 
      fprintf(fp2, "\n"); } } 
  fclose(fp);
  fclose(fp2);
  
}

void sysFormation::setup_RHS(const vector<complex<double> >& zLoads) {
}


void sysFormation::M_mult_ZEM_mult_MT  (pfft::SpColMat<std::complex<double> >& mZEMmt, 
					const pfft::SpColMat<std::complex<double> >& ZEM, 
					const pfft::SpRowMat<double>& Mx) {
  
}


void sysFormation::solve_system_direct(vector<complex<double> >& Im,
                                       vector<complex<double> >& Ib,
                                       complex<double>& adm,
                                       complex<double>& Rsource) {
}


/******************************************************************* *
misc. auxillary functions
*******************************************************************/


/*****************************************************************/ 
// this can be moved to fastsub!! and a complete new file instead of systemFormation can be used

void sysFormation::Solver(const Dipoles& dipoles, const double freq,
                              const SimulationType simuType, const
                              FormulationType formulationType, const
                              SolverType solverTypeIn) 
{
  
  { //formulationType==FORMULATION_C 
    if (solverTypeIn==DIRECT_SOLVER)
      { 
	// call Xin's direct 
      } 
    else if (solverTypeIn==pFFT) 
      {
	pfftSolver(dipoles, freq, simuType); 
      }
    
  }				   
}


/********************************************************************************
**			IterativeSolver --
** 
*******************************************************************************/

void sysFormation::pfftSolver(const Dipoles& dipoles,
                                const double freq,
				const SimulationType simuType)
  
{
  bool CONVERGE;
  //setup_M_Elements();
  
  TNT::Vector<complex<double> > RHS(Mpfft.numRow());
  TNT::Vector<complex<double> > x0(Mpfft.numRow());
  complex<double> Zout;

  setup_freq_info(freq);

  if ((dipole_flag == true) || (simuType == FULLWAVE)) {
    calcADyadic = CalcAfullwave3D(outerOperatorList, innerOperatorList, k0, dipole_flag,
				  dipoles, pfft::GXXc, simuType);
    DynamicGreenFunc kernel_XX(k0, dipoles, pfft::GXXc);
    kernel_XX.set_shift(get_shift(ElementList));
    pfftDyadic.constructKernelDependent(kernel_XX, calcADyadic); 
    
    calcPScalar = CalcPfullwave3D(outerOperatorList, innerOperatorList, k0,
				  dipole_flag, dipoles, pfft::GScalar, simuType);
    DynamicGreenFunc kernel_P(0., dipoles, pfft::GScalar);
    
    //kernel_P=DynamicGreenFunc(0., dipoles, pfft::GZZ_phi);
    kernel_P.set_shift(get_shift(ElementListP));
    pfftScalar.constructKernelDependent(kernel_P, calcPScalar);
  }
    
  PreConditionerSetup();
  
  RHS = setup_RHS_InitialGuess();
  for (int i=0; i<Mpfft.numRow(); i++)
    x0[i]=0.;
  //	x0[i]=RHS[i];
  cout << "finished Preconditioner Setup" << endl;
  CONVERGE=gmres(this, x0, RHS,MpreCondMatMT, Mpfft.numRow(), 100, 1e-4);
  
  //cout << x0[port_location_] << endl;
  cout << x0[Mpfft.numRow()-2] << endl;
  cout << x0[Mpfft.numRow()-1] << endl;
  
  // cout << x0 << endl;
  // debug_flag = false;
  Zout = 4./(x0[Mpfft.numRow()-1]-x0[Mpfft.numRow()-2]);
    {
      FILE* output_file;
      output_file=fopen("output.mat", "a");
      fprintf(output_file, "%E %E %E\n ", freq, Zout.real(), Zout.imag());
      fclose(output_file);
    }
  
}

TNT::Vector<complex<double> > sysFormation::MatrixVectorProduct_C(TNT::Vector<complex<double> > RHS)
{
  
  
  TNT::Vector<complex<double> > AP_Vmpfft(RHS.size());
  //TNT::Vector<complex<double> > IbX(ElementListX.size());
  //TNT::Vector<complex<double> > IbY(ElementListY.size());
  //TNT::Vector<complex<double> > IbZ(ElementListZ.size());
  TNT::Vector<complex<double> > IbF(ElementList.size());
  TNT::Vector<complex<double> > IbX(ElementList.size());
  TNT::Vector<complex<double> > IbY(ElementList.size());
  TNT::Vector<complex<double> > IbZ(ElementList.size());
  TNT::Vector<complex<double> > IbP(ElementListP.size());
  
  TNT::Vector<complex<double> > Ib(MTpfft.numRow());
  
  matMultVec(Ib,MTpfft,RHS);
  

  for (int ip=0; ip<ElementListP.size(); ip++)
    IbP[ip]=Ib[ip+ElementList.size()];
  
  for (int i = 0; i<DirectionMatrix.size(); i++) {
    IbF[i] = Ib[i];
    IbX[i] = Ib[i] * abs(DirectionMatrix[i].x());
    IbY[i] = Ib[i] * abs(DirectionMatrix[i].y());
    IbZ[i] = Ib[i] * abs(DirectionMatrix[i].z());
  }

  
  TNT::Vector<complex<double> > Vb(ElementList.size() + ElementListP.size());
  TNT::Vector<complex<double> > y_F (ElementList.size());
  TNT::Vector<complex<double> > y_XX(ElementList.size());
  TNT::Vector<complex<double> > y_YY(ElementList.size());
  TNT::Vector<complex<double> > y_ZZ(ElementList.size());
  TNT::Vector<complex<double> > y_P(ElementListP.size());
  
  
  if (ElementListP.size() > 0) {
    pfftScalar.IEoperator(y_P, outerOperatorList[0], innerOperatorList[0],
			IbP);
  }
  
  if (ElementList.size() > 0) {
    pfftDyadic.IEoperatorDyadic(y_F, y_XX, y_YY, y_ZZ, 
			     outerOperatorList[0], innerOperatorList[0],
			     IbF, IbX, IbY, IbZ);
  }
  
  //cout << "finished step 2" << endl;
  //cout << y << endl;
  for (int i=0; i<ElementList.size(); i++) {
    Vb[i]=y_F[i]*IMAG*omega*1e-7;
  }
  
  for (int i = 0; i<DirectionMatrix.size(); i++) {
    if (DirectionMatrix[i].x() != 0)
      Vb[i] += (y_XX[i] * abs(DirectionMatrix[i].x())*IMAG*omega*1e-7);
    if (DirectionMatrix[i].y() != 0)
      Vb[i] += (y_YY[i] * abs(DirectionMatrix[i].y())*IMAG*omega*1e-7);
    if (DirectionMatrix[i].z() != 0)
      Vb[i] += (y_ZZ[i] * abs(DirectionMatrix[i].z())*IMAG*omega*1e-7);
  }
    
  
  for (int i=0; i<ElementListP.size(); i++) {
    Vb[ElementList.size()+i]
      =y_P[i]/(IMAG*omega*4.*M_PI*layers.get_epir_topmostlayer()*EPI0);
  }
  
  // cout << Vb << endl;
  matMultVec(AP_Vmpfft,Mpfft,Vb);
  //cout << AP_Vmpfft << endl; 
  TNT::Vector<complex<double> > Vr(ElementList.size());
  TNT::Vector<complex<double> > Ibr(ElementList.size());
  TNT::Vector<complex<double> > VrPadded(MTpfft.numRow());
 
 for (int i=0; i<ElementList.size();i++)
    Ibr[i]=Ib[i];
  matMultVec(Vr,RMat,Ibr);
  for (int i=0; i<ElementList.size();i++)
    VrPadded[i]=Vr[i];
  for (int i=0; i<ElementListP.size();i++)
    VrPadded[i+ElementList.size()]=0.;
  for (int i = numPanels + numFilaments; i < MTpfft.numRow(); i++)
    VrPadded[i] = RList[i - numPanels - numFilaments] * Ib[i];
  
  
  matMultVecAdd(AP_Vmpfft,Mpfft,VrPadded);
  
  return AP_Vmpfft;
  
}
/********************************************************
setup_RHS_InitialGuess();
*********************************************************************************/
TNT::Vector<complex<double> > sysFormation::setup_RHS_InitialGuess()
{
  
  TNT::Vector<std::complex<double> > RHSpfft(Mpfft.numRow());
  
  for (int i=0; i<Mpfft.numRow()-2; i++)
    //for (int i=0; i<Mpfft.numRow(); i++)
    RHSpfft[i]=(complex<double>(0.,0.));
  
  //  RHSpfft[port_location_]=(complex<double>(1.,0.));
  
  RHSpfft[Mpfft.numRow()-2]=(complex<double>(-1.0, 0.));
  RHSpfft[Mpfft.numRow()-1]=(complex<double>( 1.0, 0.));
  return RHSpfft;	
}





/******************************************************************************
	setup_M this function is used to reorder the rows and vectors of
******************************************************************************/
void sysFormation::setup_M_Elements(void)
{
  filament<double> *f1;
  panel<double> *p1, *p2;
  //int i,j;
  double f_x, f_y, f_z;
  int numfil_X=0, numfil_Y=0, numfil_Z=0;
  int Tnumfil_X=0, Tnumfil_Y=0, Tnumfil_Z=0, Tnumpanel_H=0, Tnumpanel_V=0;
  //	Mpfft = pfft::SpRowMat<double>(M.numRow()-2, M.numCol()-1);
  //	MTpfft = pfft::SpRowMat<double>(M.numCol()-1, M.numRow()-2);
  //	pfft::SpColMat<double>  MCol(M.numCol()-1, M.numRow());
  // cout << M.numCol() << endl;
  Mpfft = pfft::SpRowMat<double>(M.numRow(), M.numCol());
  MTpfft = pfft::SpRowMat<double>(M.numCol(), M.numRow());
  pfft::SpColMat<double>  MCol(M.numCol(), M.numRow());
  
  for (int i=0; i<M.numRow(); i++) {
    for (int j=0; j < M[i].size(); j++) {
      // if (M[i].index(j)<12)
      MCol.insertElement(i,M[i].index(j),M[i].value(j));
      // else if (M[i].index(j)>12)
      // 		MCol.insertElement(i,M[i].index(j)-1,M[i].value(j));
    }
  }
  
  
  for (f1=condStruc.list_filaments.initialize(); f1!=NULL;
       f1=condStruc.list_filaments.iterate()) {
    f_x= (f1->l_dir()).x();
    f_y= (f1->l_dir()).y();
    f_z= (f1->l_dir()).z();
    if (f_z != 0) 
      VerticalFilFlag = true;

    if (f_x != 0) Tnumfil_X++;
    else if (f_y != 0) Tnumfil_Y++;
    else if (f_z != 0) Tnumfil_Z++;
    else cout << " This is impossible:: filament without direction!! " << endl;
  }
  
  //DirectionMatrix = vector<point3D<double> >(Tnumfil_X + Tnumfil_Y + Tnumfil_Z); 
  
  for (p1=condStruc.list_panels.initialize(); p1!=NULL;
       p1=condStruc.list_panels.iterate()) 
    {
      if (p1->direction().x()!=0 || p1->direction().y()!=0) 
	Tnumpanel_H++;
      else if (p1->direction().z()!=0) 
	Tnumpanel_V++;
      else cout << p1->direction().x() << " " << p1->direction().y() << " " << p1->direction().z() << endl;
    }
  
  int globalCounter = 0;
  int idx = 0;
  for (f1=condStruc.list_filaments.initialize(); f1!=NULL;
       f1=condStruc.list_filaments.iterate()) {
    f_x= (f1->l_dir()).x();
    f_y= (f1->l_dir()).y();
    f_z= (f1->l_dir()).z();
    pfft::element ele(f1,numfil_X+numfil_Y+numfil_Z,-1,0);
    idx = numfil_X+numfil_Y+numfil_Z;
    if (f_x!=0) {
      for (int R_idx=0; R_idx < MCol[idx].size(); R_idx++) {
      	Mpfft.insertElement(MCol[idx].index(R_idx), globalCounter,  MCol[idx].value(R_idx));
      	MTpfft.insertElement(globalCounter, MCol[idx].index(R_idx), MCol[idx].value(R_idx));
      }
      numfil_X++;
      globalCounter++;
      ElementListX.push_back(ele);
      ElementList.push_back(ele);
      DirectionMatrix.push_back(f1->l_dir());
    } else if (f_y!=0) {
      for (int R_idx=0; R_idx < MCol[idx].size(); R_idx++) {
      	Mpfft.insertElement(MCol[idx].index(R_idx), globalCounter,  MCol[idx].value(R_idx));
      	MTpfft.insertElement(globalCounter, MCol[idx].index(R_idx), MCol[idx].value(R_idx));
      }
      numfil_Y++;
      globalCounter++;
      ElementListY.push_back(ele);
      ElementList.push_back(ele);
      DirectionMatrix.push_back(f1->l_dir());
    } else if (f_z!=0) {
      for (int R_idx=0; R_idx < MCol[idx].size(); R_idx++) {
	Mpfft.insertElement(MCol[idx].index(R_idx), globalCounter,  MCol[idx].value(R_idx));
	MTpfft.insertElement(globalCounter, MCol[idx].index(R_idx), MCol[idx].value(R_idx));
      }
      numfil_Z++;
      globalCounter++;
      ElementListZ.push_back(ele);
      ElementList.push_back(ele);
      DirectionMatrix.push_back(f1->l_dir());
    }
  }
  
  //cout << "number of source filaments is: " << Tnumfil_X+Tnumfil_Y+Tnumfil_Z << endl;
  
  
  int numPanel_H=0, numPanel_V=0;
  int idxP=0;
  int Tnumfil=Tnumfil_X+Tnumfil_Y+Tnumfil_Z;
  int idxV=Tnumfil+Tnumpanel_H;
  for (p1=condStruc.list_panels.initialize(); p1!=NULL;
       p1=condStruc.list_panels.iterate()) {
    pfft::element ele(p1,0,-2,0);
    idxP = Tnumfil+numPanel_H+numPanel_V;
    if ((p1->direction().x()!=0) || (p1->direction().y()!=0)) {
      for (int R_idx=0; R_idx < MCol[idxP].size(); R_idx++) {
      	Mpfft.insertElement(MCol[idxP].index(R_idx), Tnumfil+numPanel_H, MCol[idxP].value(R_idx));
      	MTpfft.insertElement(Tnumfil+numPanel_H, MCol[idxP].index(R_idx), MCol[idxP].value(R_idx));
      }
      numPanel_H++;
      ElementListPH.push_back(ele);
    }
    else if (p1->direction().z()!=0) {
      for (int R_idx=0; R_idx < MCol[idxP].size(); R_idx++) {
      	Mpfft.insertElement(MCol[idxP].index(R_idx), idxV+numPanel_V,  MCol[idxP].value(R_idx));
      	MTpfft.insertElement(idxV+numPanel_V, MCol[idxP].index(R_idx), MCol[idxP].value(R_idx));
      } 
      numPanel_V++;
      ElementListPV.push_back(ele);
    }
  }
  for (int i=0; i<ElementListPH.size(); i++)
    ElementListP.push_back(ElementListPH[i]);
  for (int i=0; i<ElementListPV.size(); i++)
    ElementListP.push_back(ElementListPV[i]);
  
  cout << "number of source panels is: " << ElementListP.size() << endl;
  // cout << "number of horizontal source panels is: " << numPanel_H << endl;
  // cout << "number of vertical source panels is: " << numPanel_V << endl;

  // output_to_file_real("Mpfft.mat", Mpfft);
  // output_to_file_real("MTpfft.mat", MTpfft);
  
  if (numFilaments+numPanels != Tnumfil + numPanel_V + numPanel_H) {
    cout << "something wrong here" << endl;
    exit(1);
  }
  
  if ((condStruc.mesh_mats).Mrs_ptr.numRow()) {
    cout << "how am i here" << endl;
    insert_subMatrix(Mpfft, (condStruc.mesh_mats).Mf_ptr.numRow()+(condStruc.mesh_mats).Mfs_ptr.numRow(), 
		     numFilaments+numPanels, (condStruc.mesh_mats).Mrs_ptr);
    insert_subMatrix(MTpfft, numFilaments+numPanels, 
		     (condStruc.mesh_mats).Mf_ptr.numRow()+(condStruc.mesh_mats).Mfs_ptr.numRow(), (condStruc.mesh_mats).Mrs_ptr);
  }
  
  //cout << Mpfft.numRow() << endl;
  //cout << Mpfft.numCol() << endl;
  
}


/*************************************************************************
* setup pfft this function serves to setup the elements and call
* all the different pfft routines with the appropriate initialization
**************************************************************************/
void sysFormation::setup_Pfft(
			      //const Dipoles& dipoles,
			      //const double freq,
			      const pfft::KernelSize DyadSize,
			      const pfft::KernelSize ScalarSize,
			      const SimulationType simuType)
  
{
  
  
  //ofstream fout("error.dat");
  size_t max_size = 1;
  
  
  size_t interpStencilSize =  max_size;
  size_t projectStencilSize = max_size;
  size_t directStencilSize = 2*max_size;

  cout << endl << "\t Setting Pfft::" << endl;
  
  outerOperatorList.push_back(pfft::D0_D0);
  innerOperatorList.push_back(pfft::D0_D0);
  
  if (ElementList.size() > 0) {
    calcA = CalcAfullwave3D(outerOperatorList, innerOperatorList);
    pfftDyadic = PfftA(ElementList, ElementList, calcA,
		       directStencilSize, projectStencilSize, interpStencilSize, 
		       DyadSize);
  }
  
  
  if (ElementListP.size() > 0) {
    
    calcP = CalcPfullwave3D(outerOperatorList, innerOperatorList);
    pfftScalar = PfftPanel(ElementListP, ElementListP, calcP,
			   directStencilSize, projectStencilSize, interpStencilSize, 
			   ScalarSize);
  }
  
  
}



void sysFormation::PreConditionerSetup(void)
{

  TNT::Vector<complex<double> > ans(Mpfft.numCol());
  //pfft::SpVec<std::complex<double> > ans(Mpfft.numCol());
  MpreCondMatMT=pfft::SpColMat<std::complex<double> >(Mpfft.numRow(),Mpfft.numRow());
  
  Resistance *resist;
  // TNT::Vector<complex<double> > External_Resistance(Mpfft.numRow());
  // MRMT = TNT::Vector<complex<double> > (Mpfft.numRow());
  
  int numResist = 0;
  for (resist=condStruc.list_resistances.initialize(); resist!=NULL;
       resist=condStruc.list_resistances.iterate()) {
    numResist++;
  }
  
  RList =  TNT::Vector<complex<double> > (numResist);
  int i = 0;
  for (resist=condStruc.list_resistances.initialize(); resist!=NULL;
       resist=condStruc.list_resistances.iterate()) {
    //cout <<  RList[i] << endl;
    RList[i] = (resist->get_r_value());
    i++;
  }  
  
  
  for (int i = 0; i < Mpfft.numRow(); i++) {
    ans[i] = 0.;
  }
  cout << "Setting up Preconditioner" << endl;
  for(size_t R_idx=0; R_idx < Mpfft.numRow(); R_idx++)
    {
      if ((((int)R_idx) % 100)==0)
	cout << "finished row # = " << R_idx  << endl;
      
      for (size_t i = 0; i < Mpfft[R_idx].size(); i++)
	{
	  size_t row = Mpfft[R_idx].index(i);
	  for (size_t col = 0; col < MTpfft[row].size(); col++)
	    {
	      if (row < ElementList.size())
		ans[MTpfft[row].index(col)] += Mpfft[R_idx].value(i) * MTpfft[row].value(col) 
		  * ((pfftDyadic.getPreCondMat())[row] * IMAG*omega*1e-7 + RVec[row]);
	      /* else if (row < ElementListX.size()+ElementListY.size())
		 ans[MTpfft[row].index(col)] += Mpfft[R_idx].value(i) * MTpfft[row].value(col) 
		 * (pfft_YY.getPreCondMat())[row-ElementListX.size()] * IMAG*omega*1e-7;
		 else if (row < ElementListX.size()+ElementListY.size()+ElementListZ.size())
		 ans[MTpfft[row].index(col)] += Mpfft[R_idx].value(i) * MTpfft[row].value(col) 
		 * (pfft_ZZ.getPreCondMat())[row-ElementListX.size()-ElementListY.size()] * IMAG*omega*1e-7;*/
	      else if (row < ElementListX.size()+ElementListY.size()+ElementListZ.size()+ElementListP.size())
		ans[MTpfft[row].index(col)] += Mpfft[R_idx].value(i) * MTpfft[row].value(col) 
		  * (pfftScalar.getPreCondMat())[row-ElementListX.size()-ElementListY.size()-ElementListZ.size()] 
		  / (IMAG*omega*4.*M_PI*layers.get_epir_topmostlayer()*EPI0);
	      else {
		cout << "how am i here" << endl;
		int ER_idx = row - numPanels - numFilaments;
		ans[MTpfft[row].index(col)] += Mpfft[R_idx].value(i) * MTpfft[row].value(col) 
		  * RList[ER_idx];
		cout <<  RList[ER_idx] << endl;
	      }
	    }
	}
      
      for (int i = 0; i < Mpfft.numRow(); i++) {
	if (ans[i] !=0.) 
	  MpreCondMatMT.insertElement(R_idx, i, ans[i]);
	ans[i]=0.;
      }
      
      
    }
}


double sysFormation::get_shift(vector<pfft::element> ElementList)
{
  double zmin=ElementList[0].zmin();
  for (int i=1; i<ElementList.size(); i++)
    zmin=min(zmin,ElementList[i].zmin());
  return zmin;
}

