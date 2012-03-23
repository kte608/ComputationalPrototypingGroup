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
#include "inductance.h"


void Inductance::init (const bool& EMQS_flag, const bool& dipole_flag,const double epi)
{
  isEMQS=EMQS_flag;
  isDipole=dipole_flag;
  epi0=epi;

  //only for the computation of dAdx and dAdy for the static case
  if ((isEMQS) && (isDipole)) 
    calcA_X.init(0, false, true);
}

void Inductance::set_dipoles(const Dipoles& dipoles_, const double freq)
{
  dipoles=&(dipoles_);
  K_num=sqrt((2.*PI*freq)*(2.*PI*freq)*MU0*epi0);

  if (isEMQS==false)
  {
    if (isDipole==false) //fullwave + no sub
      calcA_fullWave.init(K_num,true, false);
	else
		{ //fullwave+sub
      calcA_fullWave.init(K_num,true, false);
      calcA_X.init(K_num,true, true);
     // calcA_X.init(K_num,false, true);
		}
  }
}

void Inductance::fill_L(pfft::SpRowMat<complex<double> >& LMat,
                        const int i, const int j,
                        const filament<double> *fi,
                        const filament<double> *fj)
{
  double A_double=0.;
  complex<double> A_cx=0., dAdx=0., dAdy=0.;
  complex<double> A, A_ij, A_ji;
  
  double fi_x, fi_y, fi_z, fj_x, fj_y, fj_z;
  fi_x= (fi->l_dir()).x();
  fi_y= (fi->l_dir()).y();
  fi_z= (fi->l_dir()).z();
  fj_x= (fj->l_dir()).x();
  fj_y= (fj->l_dir()).y();
  fj_z= (fj->l_dir()).z();
  
  if ( isDipole==false)
    {
      if (isEMQS==true) //EMQS + no substrate
	{
	  if(i==j)
	    {
	      A_double=calcA_static.selfterm(fi);
	      //A_double=selfterm(fi->get_henry_fil_ptr());
	      LMat.insertElement(i,i,A_double);
	    }
	  else if (i>j)
	    {
	      A_double=calcA_static.mutual(fi, fj);
	      //A_double=mutual(fi->get_henry_fil_ptr(), fj->get_henry_fil_ptr());
	      LMat.insertElement(i,j,A_double);
	      LMat.insertElement(j,i,A_double);
	    }
	}
      else //full-wave +  no substrate
	{
	  if (i==j)
	    {
	      //calcA_fullWave(*fi,*fi, A_cx, dAdx, dAdy);
	      calcA_fullWave.selfterm(*fi, A_cx);
	      LMat.insertElement(i,i,A_cx);
	    }
	  else if (i>j)
	    {
	      A_cx = 0.0;
	      if( ( fi_x!=0 && fj_x!=0) || ( fi_y!=0 && fj_y!=0) || (fi_z!=0 && fj_z!=0))
		calcA_fullWave(*fi,*fj, A_cx, dAdx, dAdy);
	      LMat.insertElement(i,j,A_cx);
	      LMat.insertElement(j,i,A_cx);
	    }
	}
    }
  else
    {
      if (i==j)
	{
	  if ( fi_x!=0 || fi_y!=0)
	    {
	      A=get_selfterm_dipoles(fi,dipoles->coes_G_A_XX,dipoles->exps_G_A_XX);
	      A_cx=A_cx+A*(fi_x*fi_x+fi_y*fi_y);
	    }
	  if (fi_z!=0)
	    {
	      A=get_selfterm_dipoles(fi,dipoles->coes_G_A_ZZ,dipoles->exps_G_A_ZZ );
	      A_cx=A_cx+A*(fi_z*fi_z);
	    }
	  //cout<<"self"<<i<<" "<<i<<endl;
	  LMat.insertElement(i,i,A_cx);
	}
      else if (i>j)
	{
	  A_ij=A_ji=0.;
	  if( ( fi_x!=0 && fj_x!=0) || ( fi_y!=0 && fj_y!=0))  //G_xx
	    {
	      A=get_mutual_dipoles(fi,fj,dipoles->coes_G_A_XX,dipoles->exps_G_A_XX);
	      A_ij += A*(fi_x*fj_x+fi_y*fj_y);
	      A_ji += A*(fi_x*fj_x+fi_y*fj_y);
	    }
	  if ((fi_z!=0) && (fj_z!=0)) //G_zz
	    {
	      A=get_mutual_dipoles(fi,fj,dipoles->coes_G_A_ZZ,dipoles->exps_G_A_ZZ);
	      A_ij +=A*(fi_z*fj_z);
	      A_ji +=A*(fi_z*fj_z);
	    }
	  if ((fi_x!=0 && fj_z!=0) || (fj_x!=0 && fi_z!=0) || (fi_y!=0 && fj_z!=0) || (fj_y!=0 && fi_z!=0))
	    {
	      get_mutual_dipoles_asym(fi, fj, dipoles->coes_G_A_XZ,dipoles->exps_G_A_XZ, dAdx, dAdy);
	      
	      if ((fi_x!=0 && fj_z!=0) || (fj_x!=0 && fi_z!=0))
		{
		  A_ij +=dAdx*(fi_x*fj_z);
		  A_ji +=dAdx*(fj_x*fi_z);
		}
	      if ((fi_y!=0 && fj_z!=0) || (fj_y!=0 && fi_z!=0))
		{
		  A_ij +=dAdy*(fi_y*fj_z);
		  A_ji +=dAdy*(fj_y*fi_z);
		}
	    }
	  LMat.insertElement(i,j, A_ij);
	  LMat.insertElement(j,i, A_ji);
	}
    } //if  EMQS and dipole
}
////////////////////////

void Inductance::Myfill_L_REDUCED(pfft::SpRowMat<complex<double> >& LMat,
                        const int i, const int j,
                        const filament<double> *fi,
                        const filament<double> *fj)
{
  double A_double=0.;
  complex<double> A_cx=0., dAdx=0., dAdy=0.;
  complex<double> A, A_ij, A_ji;

  double fi_x, fi_y, fi_z, fj_x, fj_y, fj_z;
  fi_x= (fi->l_dir()).x();
  fi_y= (fi->l_dir()).y();
  fi_z= (fi->l_dir()).z();
  fj_x= (fj->l_dir()).x();
  fj_y= (fj->l_dir()).y();
  fj_z= (fj->l_dir()).z();

	if (i==j)
  {
    if (fi_x!=0)
      A_cx=get_selfterm_dipoles(fi,dipoles->coes_G_A_XX,dipoles->exps_G_A_XX );
    if (fi_y!=0)
      A_cx=get_selfterm_dipoles(fi,dipoles->coes_G_A_ZZ_MOD,dipoles->exps_G_A_ZZ );
    if (fi_z!=0)
      A_cx=get_selfterm_dipoles(fi,dipoles->coes_G_A_ZZ,dipoles->exps_G_A_ZZ );
    //cout<<"self"<<i<<" "<<i<<endl;
    LMat.insertElement(i,i,A_cx);
  }
  else if (i>j)
  {
    A_ij=A_ji=0.;
  	if ((fi_x!=0) && (fj_x!=0)) //G_xx
    {
      A=get_mutual_dipoles(fi,fj,dipoles->coes_G_A_XX,dipoles->exps_G_A_XX);
      A_ij +=A;
      A_ji +=A;
    }
  	if ((fi_y!=0) && (fj_y!=0)) //G_yy
    {
      A=get_mutual_dipoles(fi,fj,dipoles->coes_G_A_ZZ_MOD,dipoles->exps_G_A_ZZ);
      A_ij +=A;
      A_ji +=A;
    }
  	if ((fi_z!=0) && (fj_z!=0)) //G_zz
    {
      A=get_mutual_dipoles(fi,fj,dipoles->coes_G_A_ZZ,dipoles->exps_G_A_ZZ);
      A_ij +=A;
      A_ji +=A;
    }
    LMat.insertElement(i,j, A_ij);
    LMat.insertElement(j,i, A_ji);
  }
}



///////////////////////

void Inductance::Myfill_L(pfft::SpRowMat<complex<double> >& LMat,
                        const int i, const int j,
                        const filament<double> *fi,
                        const filament<double> *fj)
{
  double A_double=0.;
  complex<double> A_cx=0., dAdx=0., dAdy=0.;
  complex<double> A, A_ij, A_ji;

  double fi_x, fi_y, fi_z, fj_x, fj_y, fj_z;
  fi_x= (fi->l_dir()).x();
  fi_y= (fi->l_dir()).y();
  fi_z= (fi->l_dir()).z();
  fj_x= (fj->l_dir()).x();
  fj_y= (fj->l_dir()).y();
  fj_z= (fj->l_dir()).z();

	if (i==j)
  {
    if ( fi_x!=0 || fi_y!=0)
      A_cx=Myget_selfterm_dipoles(fi,dipoles->coes_G_A_XX,dipoles->exps_G_A_XX,
														dipoles->coes_G_A_XX_MOD,dipoles->exps_G_A_XX_MOD);
    if (fi_z!=0)
      A_cx=get_selfterm_dipoles(fi,dipoles->coes_G_A_ZZ,dipoles->exps_G_A_ZZ );
    //cout<<"self"<<i<<" "<<i<<endl;
    LMat.insertElement(i,i,A_cx);
  }
  else if (i>j)
  {
    A_ij=A_ji=0.;
    if( ( fi_x!=0 && fj_x!=0) || ( fi_y!=0 && fj_y!=0))  //G_xx
    {
      A=Myget_mutual_dipoles(fi,fj,dipoles->coes_G_A_XX,dipoles->exps_G_A_XX,
														dipoles->coes_G_A_XX_MOD,dipoles->exps_G_A_XX_MOD);
      A_ij += A;
      A_ji += A;
    }
  	if ((fi_z!=0) && (fj_z!=0)) //G_zz
    {
      A=get_mutual_dipoles(fi,fj,dipoles->coes_G_A_ZZ,dipoles->exps_G_A_ZZ);
      A_ij +=A;
      A_ji +=A;
    }
    if ((fi_x!=0 && fj_y!=0) || (fi_y!=0 && fj_x!=0) )  //G_xy
    {
      A_ij=Myget_mutual_dipoles_asym(fi, fj, dipoles->coes_G_A_XX_MOD,dipoles->exps_G_A_XX_MOD);
      A_ji=A_ij;
    }
    LMat.insertElement(i,j, A_ij);
    LMat.insertElement(j,i, A_ji);
  }
}

complex<double> Inductance::Myget_selfterm_dipoles( const filament<double> *f,
    const vector<complex<double> >& coes,
    const vector<complex<double> >& exps,
    const vector<complex<double> >& coes2,
    const vector<complex<double> >& exps2)
{
  double val_real;
  complex<double> val_real_cx=0., val_cx=0., val_cx_XX;
  complex<double> val_cx_YY=0., val_cx_XY=0.;
  complex<double> val_cx_temp=0., junk1=0., junk2=0.;
  const int numDipole=coes.size();
  const int numDipole2=coes2.size();
//  const int numDipole=1;
  int i;
  filament<complex<double> > cx_fil, cx_fil2;
  copy_fil(f, 0., cx_fil2);

    //calcA_fullWave(*f, *f, val_real_cx, junk1, junk2);
	for (i=0; i<numDipole; i++)
    {
			//cout << exps[i] << coes[i] << endl;
      get_cx_fil(f, exps[i], cx_fil);
      val_cx_temp=0.;
      calcA_X(cx_fil, cx_fil2, val_cx_temp, junk1, junk2);
      val_cx=val_cx+coes[i]*val_cx_temp;
    }
	if ((f->l_dir()).x()!=0) {
	for (i=0; i<numDipole2; i++)
    {
			//cout << exps[i] << coes[i] << endl;
      get_cx_fil(f, exps2[i], cx_fil);
      val_cx_XX=0.; val_cx_YY=0.; val_cx_XY=0.;
      calcA_X.XX_deriv(cx_fil, cx_fil2, 1, val_cx_XX, val_cx_YY, val_cx_XY);
      val_cx=val_cx+coes2[i]*val_cx_XX;
    }
	} else if ((f->l_dir()).y()!=0) {
	for (i=0; i<numDipole2; i++)
    {
      get_cx_fil(f, exps2[i], cx_fil);
      val_cx_XX=0.; val_cx_YY=0.; val_cx_XY=0.;
      calcA_X.XX_deriv(cx_fil, cx_fil2, 2, val_cx_XX, val_cx_YY, val_cx_XY);
			// cout << coes2[i] << coes2[i]*val_cx_YY << endl;
      val_cx=val_cx+coes2[i]*val_cx_YY;
    }
	}
	calcA_fullWave.selfterm(*f, val_real_cx);
	return val_real_cx+val_cx;
}

complex<double> Inductance::Myget_mutual_dipoles( const filament<double> *f1,
    const filament<double>* f2,
    const vector<complex<double> >& coes,
    const vector<complex<double> >& exps,
    const vector<complex<double> >& coes2,
    const vector<complex<double> >& exps2)
{
  double val_real;
  complex<double> val_cx=0., val_real_cx=0., val_cx_temp;
  complex<double> val_cx_XX=0., val_cx_YY=0., val_cx_XY=0.;
  complex<double> junk1=0., junk2=0.;
  const int numDipole=coes.size();
  const int numDipole2=coes2.size();
 // const int numDipole=1;
  int i;
  filament<complex<double> > cx_fil, cx_fil2;
  copy_fil(f1, 0., cx_fil2);


    ///calcA_fullWave(*f1, *f2, val_real_cx, junk1, junk2);
    for (i=0; i<numDipole; i++)
    {
      get_cx_fil(f2, exps[i], cx_fil);
      val_cx_temp=0.;
     /// calcA_X(cx_fil, *f1, val_cx_temp, junk1, junk2);
      calcA_X(cx_fil, cx_fil2, val_cx_temp, junk1, junk2);
      val_cx=val_cx+coes[i]*val_cx_temp;
    }
	if ((f1->l_dir()).x()!=0) {
		for (i=0; i<numDipole2; i++)
    {
			//cout << exps[i] << coes[i] << endl;
      get_cx_fil(f2, exps2[i], cx_fil);
      val_cx_temp=0.;
      calcA_X.XX_deriv(cx_fil, cx_fil2, 1, val_cx_XX, val_cx_YY, val_cx_XY);
      val_cx=val_cx+coes2[i]*val_cx_XX;
    }
	} else if ((f1->l_dir()).y()!=0) {
		for (i=0; i<numDipole2; i++)
    {
			//cout << exps[i] << coes[i] << endl;
      get_cx_fil(f2, exps2[i], cx_fil);
      val_cx_temp=0.;
      calcA_X.XX_deriv(cx_fil, cx_fil2, 2, val_cx_XX, val_cx_YY, val_cx_XY);
      val_cx=val_cx+coes2[i]*val_cx_YY;
    }
	}
    calcA_fullWave(*f1, *f2, val_real_cx, junk1, junk2);
    return val_real_cx+val_cx;
}

complex<double> Inductance::Myget_mutual_dipoles_asym(const filament<double> *f1,
    const filament<double>* f2,
    const vector<complex<double> >& coes,
    const vector<complex<double> >& exps)
{
  complex<double> val_cx_XX=0., val_cx_YY=0., val_cx_XY=0., val_cx=0.;
  filament<complex<double> > cx_fil, cx_fil2;
  copy_fil(f1, 0., cx_fil2);

  const int numDipole=coes.size();
 // const int numDipole=1;

    for (int i=0; i<numDipole; i++)
    {
      get_cx_fil(f2, exps[i], cx_fil);
      val_cx_XY=0.;
     /// calcA_X(cx_fil, *f1, junk, dIdx, dIdy);
      calcA_X.XX_deriv(cx_fil, cx_fil2, 3, val_cx_XX, val_cx_YY, val_cx_XY);
      val_cx=val_cx+coes[i]*val_cx_XY;
    }
	return val_cx;
}

//////////////////////


complex<double> Inductance::get_selfterm_dipoles( const filament<double> *f,
    const vector<complex<double> >& coes,
    const vector<complex<double> >& exps)
{
  double val_real;
  complex<double> val_real_cx=0., val_cx=0., val_cx_temp;
  complex<double> junk1=0., junk2=0.;
  const int numDipole=coes.size();
//  const int numDipole=1;
  int i;
  filament<complex<double> > cx_fil, cx_fil2;
  copy_fil(f, 0., cx_fil2);

  if (isEMQS==true)
  {
    val_real=calcA_static.selfterm(f);
    for (i=0; i<numDipole; i++)
    {
      get_cx_fil(f, exps[i], cx_fil);
      val_cx=val_cx+coes[i]*calcA_static.mutual(f, &cx_fil);
    }
    return val_real+val_cx;
  }
  else
  {
    //calcA_fullWave(*f, *f, val_real_cx, junk1, junk2);
	calcA_fullWave.selfterm(*f, val_real_cx);
	for (i=0; i<numDipole; i++)
    {			
			//cout << exps[i] << coes[i] << endl;
      get_cx_fil(f, exps[i], cx_fil);
      val_cx_temp=0.;
      calcA_X(cx_fil, cx_fil2, val_cx_temp, junk1, junk2);
      val_cx=val_cx+coes[i]*val_cx_temp;
    }
    return val_real_cx+val_cx;
  }
}


complex<double> Inductance::get_mutual_dipoles( const filament<double> *f1,
    const filament<double>* f2,
    const vector<complex<double> >& coes,
    const vector<complex<double> >& exps)
{
  double val_real;
  complex<double> val_cx=0., val_real_cx=0., val_cx_temp;
  complex<double> junk1=0., junk2=0.;
  const int numDipole=coes.size();
 // const int numDipole=1;
  int i;
  filament<complex<double> > cx_fil, cx_fil2;
  copy_fil(f1, 0., cx_fil2);


  if (isEMQS==true)
  {
    val_real=calcA_static.mutual(f1, f2);
    for (i=0; i<numDipole; i++)
    {
      get_cx_fil(f2, exps[i], cx_fil);
      complex<double> valcx(0.,0.);
      valcx= calcA_static.mutual(f1, &cx_fil);
      val_cx=val_cx+coes[i]*valcx;
    }
    return val_real+val_cx;
  }
  else
  {
    ///calcA_fullWave(*f1, *f2, val_real_cx, junk1, junk2);
    for (i=0; i<numDipole; i++)
    {
      get_cx_fil(f2, exps[i], cx_fil);
      val_cx_temp=0.;
     /// calcA_X(cx_fil, *f1, val_cx_temp, junk1, junk2);
      calcA_X(cx_fil, cx_fil2, val_cx_temp, junk1, junk2);
      val_cx=val_cx+coes[i]*val_cx_temp;
    }
    calcA_fullWave(*f1, *f2, val_real_cx, junk1, junk2);
    return val_real_cx+val_cx;
  }
}

void Inductance::get_mutual_dipoles_asym(const filament<double> *f1,
    const filament<double>* f2,
    const vector<complex<double> >& coes,
    const vector<complex<double> >& exps,
    complex<double>& dAdx, complex<double>& dAdy)
{
  complex<double>junk=0., dIdx=0., dIdy=0.;
  filament<complex<double> > cx_fil, cx_fil2;
  copy_fil(f1, 0., cx_fil2);

  const int numDipole=coes.size();
 // const int numDipole=1;

  if (isEMQS)
  {
    for (int i=0; i<numDipole; i++)
    {
      get_cx_fil(f2, exps[i], cx_fil);
      dIdx=dIdy=0.;
      calcA_X(cx_fil, cx_fil2, junk, dIdx, dIdy);
      dAdx=dAdx+coes[i]*dIdx;
      dAdy=dAdy+coes[i]*dIdy;
    }
  }
  else
  {
    for (int i=0; i<numDipole; i++)
    {
      get_cx_fil(f2, exps[i], cx_fil);
      dIdx=dIdy=0.;
     /// calcA_X(cx_fil, *f1, junk, dIdx, dIdy);
      calcA_X(cx_fil, cx_fil2, junk, dIdx, dIdy);
      dAdx=dAdx+coes[i]*dIdx;
      dAdy=dAdy+coes[i]*dIdy;
    }
  }
}


void Inductance::get_cx_fil(const filament<double> *f,
                            const complex<double> cx_location,
                            filament<complex<double> >& cx_fil)
{
  double subHeight=dipoles->substrateHeight;

  complex<double> x_val(f->loc0().x(),0);
  complex<double> y_val(f->loc0().y(),0);

  point3D<complex<double> > endPt0 (x_val, y_val,
                                    (2*subHeight-f->loc0().z())+IMAG*cx_location );


  complex<double> x_val_pt2(f->loc1().x(),0);
  complex<double> y_val_pt2(f->loc1().y(),0);
  point3D<complex<double> >  endPt1(x_val_pt2,y_val_pt2,
                                    (2*subHeight-f->loc1().z())+IMAG*cx_location) ;
  double width= f->filWidth();
  double height=f->filHeight();
  double sigma=f->sigma();

  cx_fil.initialize(endPt0, endPt1, width, height, sigma);
}


void Inductance::copy_fil(const filament<double> *f,
                            const complex<double> cx_location,
                            filament<complex<double> >& cx_fil)
{
  double subHeight=dipoles->substrateHeight;

  complex<double> x_val(f->loc0().x(),0);
  complex<double> y_val(f->loc0().y(),0);
  complex<double> z_val(f->loc0().z(),0);

  point3D<complex<double> > endPt0 (x_val, y_val, z_val);


  complex<double> x_val_pt2(f->loc1().x(),0);
  complex<double> y_val_pt2(f->loc1().y(),0);
  complex<double> z_val_pt2(f->loc1().z(),0);
  point3D<complex<double> >  endPt1(x_val_pt2,y_val_pt2,z_val_pt2) ;

  double width= f->filWidth();
  double height=f->filHeight();
  double sigma=f->sigma();

  cx_fil.initialize(endPt0, endPt1, width, height, sigma);
}
