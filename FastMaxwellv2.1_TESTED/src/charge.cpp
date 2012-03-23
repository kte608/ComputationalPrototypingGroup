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
#include "charge.h"
#include <math.h>
#include <stdio.h>

void Charge::init(const bool& EMQS_flag, const bool& dipole_flag, Structure& structure,
                  const double epi)
{

  //only initialize freq-independ variables
  isEMQS=EMQS_flag;
  isDipole=dipole_flag;
  condStruc=&(structure);
  epi0=epi;

  if (isEMQS)
  {
    int numPanels=condStruc->num_panels();
    calcp_static.init(0, false, true);
    
	calcp_static_cx.init(0, false, true);
  }
}


void Charge::set_dipoles(const Dipoles& dipoles_, const double freq)
{
  // set freq-dependent variables
  dipoles=&(dipoles_);
  K_num=sqrt((2.*PI*freq)*(2.*PI*freq)*MU0*epi0);

  if (isEMQS==false) //fullwave 
  {
    if (isDipole==false)
      calcp_fullWave.init(K_num,false, true);
    else {
      calcp_fullWave.init(K_num,false, true);
      calcp_fullWave_cx.init(K_num,false, true); }
  }
}

void Charge::fill_p(pfft::SpRowMat<complex<double> >& PMat,
                    const int i,const int j,
                    const panel<double>* p1,  const panel<double>* p2)
{

  complex<double> slp=0., val=0., val_temp=0., phi=0.,zn=0., phi_collocation=0., test=0.;

  if ( isDipole==false)
  {
    if (isEMQS==true)
    {
     	  calcp_static(*p2, p1->centroid(), slp, val, zn);
     	  val=val/(p2->area());
    //  val=calcp(chg_array[j], chg_array[i]->x,
      //          chg_array[i]->y, chg_array[i]->z, NULL);  //this already area normalized
			//cout << val << endl;

    }
    else
    {
      calcp_fullWave(*p2, p1->centroid(),slp,val, zn);
      val=val/p2->area();
    }

    PMat.insertElement(i,j,(val)/(4*M_PI*epi0));
  }
  else
  {
    double nx=p2->normal().x();
    double ny=p2->normal().y();
    double nz=p2->normal().z();
 	//cout<<"normals "<<nx<<" "<<ny<<" "<<nz<<endl;
	

		if ((nx !=0) || (ny!=0)) {
		phi = phi + get_phi_hor_dipoles(p2, p1->centroid(), dipoles->coes_G_PHI_H, dipoles->exps_G_PHI_H,i,j);
   // val=val + phi * (nx+ny);
    val=val + phi;
		}
    if (nz !=0) {
		phi = phi + get_phi_hor_dipoles(p2, p1->centroid(), dipoles->coes_G_PHI_H, dipoles->exps_G_PHI_H,i,j);
		//phi = phi + get_phi_vert_dipoles(p2, p1->centroid(), dipoles->coes_G_A_ZZ_MOD, dipoles->exps_G_A_ZZ,i,j);
   // val=val + phi * nz;
    val=val + phi;
		}
    PMat.insertElement(i,j,(val)/(4*M_PI*epi0));
  }
}

void Charge::Myfill_p(pfft::SpRowMat<complex<double> >& PMat,
                    const int i,const int j,
                    const panel<double>* p1,  const panel<double>* p2)
{

  complex<double> slp=0., val=0., val_temp=0., phi=0.,zn=0., phi_collocation=0., test=0.;

    double nx=p2->normal().x();
    double ny=p2->normal().y();
    double nz=p2->normal().z();
 	//cout<<"normals "<<nx<<" "<<ny<<" "<<nz<<endl;


		//phi = phi + get_phi_vert_dipoles(p2, p1->centroid(), dipoles->coes_G_A_ZZ_MOD, dipoles->exps_G_A_ZZ,i,j);
  if (i==j) {
		phi = phi + Myget_phi_vert_dipoles(p2, p1, dipoles->coes_G_A_ZZ_MOD, dipoles->exps_G_A_ZZ,i,j);
    val=val + phi;
    PMat.insertElement(i,j,(val)/(4*M_PI*epi0));
	}
  else if (i>j) {
		phi = phi + Myget_phi_vert_dipoles(p2, p1, dipoles->coes_G_A_ZZ_MOD, dipoles->exps_G_A_ZZ,i,j);
    val=val + phi;
    PMat.insertElement(i,j,(val)/(4*M_PI*epi0));
    PMat.insertElement(j,i,(val)/(4*M_PI*epi0));
	}
}


complex<double> Charge::get_phi_hor_dipoles(const panel<double>* p, const point3D<double>& eval,
    const vector<complex<double> >& coes,
    const vector<complex<double> >& exps, const int panel_i, const int panel_j)
{
  const int numDipole=coes.size();
 // const int numDipole=1;
  complex<double> val_real=0., val_cx=0., val_cx_temp=0., slp=0.,zn=0.;
  panel<complex<double> > cx_panel;

  if (isEMQS==true)
  {
    //      calcp_static(*p, eval, slp, val_real,zn);
    // 	 val_real=val_real/p->area();
    // cout << panel_j << panel_i << endl;
          calcp_static(*p, eval, slp, val_real,zn);
     	 val_real=val_real/p->area();
    //val_real=calcp(chg_array[panel_j], chg_array[panel_i]->x,
      //             chg_array[panel_i]->y, chg_array[panel_i]->z, NULL);
    // cout << val_real << endl;
    
	for (int i=0; i<numDipole; i++)
    {
	 get_cx_panel(p, exps[i], cx_panel);
      val_cx_temp=0.;
      calcp_static_cx(cx_panel, eval, slp, val_cx_temp,zn);
      val_cx_temp=val_cx_temp/p->area();
      val_cx=val_cx+coes[i]*val_cx_temp;
    }
    return val_real+val_cx;
  }
  else
  {
    calcp_fullWave(*p, eval, slp, val_real, zn);
    val_real=val_real/p->area();

    for (int i=0; i<numDipole; i++)
    {
			//cout << exps[i] << endl;
      get_cx_panel(p, exps[i], cx_panel);
      val_cx_temp=0.;
      calcp_fullWave_cx(cx_panel, eval, slp, val_cx_temp,zn);
      val_cx_temp=val_cx_temp/p->area();
			//cout << coes[i] << endl;
      val_cx=val_cx+coes[i]*val_cx_temp;
    }
    return val_real+val_cx;
  }
}

complex<double> Charge::get_phi_vert_dipoles(const panel<double>* p, const point3D<double>& eval,
    const vector<complex<double> >& coes,
    const vector<complex<double> >& exps, const int panel_i, const int panel_j)
{
  const int numDipole=coes.size();
 // const int numDipole=1;
  complex<double> val_real=0., val_cx=0., val_cx_temp=0., slp=0.,zn=0.;
  panel<complex<double> > cx_panel;

  if (isEMQS==true)
  {
     	  calcp_static(*p, eval, slp, val_real, zn);
     	  val_real=val_real/(p->area());
	
   // val_real=calcp(chg_array[panel_i], chg_array[panel_j]->x,
     //              chg_array[panel_j]->y, chg_array[panel_j]->z, NULL);

    for (int i=0; i<numDipole; i++)
    {

      get_cx_panel(p, exps[i], cx_panel);
      val_cx_temp=0.;
      calcp_static_cx(cx_panel, eval, slp, val_cx_temp, zn);
      val_cx_temp=val_cx_temp/p->area();
      val_cx=val_cx+coes[i]*val_cx_temp;
    }
    return val_real+val_cx;
  }
  else
  {
    calcp_fullWave(*p, eval, slp, val_real, zn);
    val_real=val_real/p->area();

    for (int i=0; i<numDipole; i++)
    {
      get_cx_panel(p, exps[i], cx_panel);
      val_cx_temp=0.;
      calcp_fullWave_cx(cx_panel, eval, slp, val_cx_temp, zn);
      val_cx_temp=val_cx_temp/p->area();
      val_cx=val_cx+coes[i]*val_cx_temp;
    }
    return val_real+val_cx;
  }
}



complex<double> Charge::Myget_phi_vert_dipoles(const panel<double>* p, const panel<double>* peval,
    const vector<complex<double> >& coes,
    const vector<complex<double> >& exps, const int panel_i, const int panel_j)
{
  const int numDipole=coes.size();
 // const int numDipole=1;
  complex<double> val_real=0., val_cx=0., val_cx_temp=0., slp=0.,zn=0., dlp=0.;
  panel<complex<double> > cx_panel;
  panel<double>  r_panel;
  size_t numPtsW=2;
  size_t numPtsH=2;
  vector<point3D<double > > OUTER_quad_pts;
  vector<double > OUTER_quad_weights;


  	OUTER_quad_pts.resize(numPtsW*numPtsH);
  	OUTER_quad_weights.resize(numPtsW*numPtsH);
    get_r_panel(peval, 0., r_panel);
		r_panel.panel_quadrature(numPtsW * numPtsH, OUTER_quad_pts, OUTER_quad_weights);
		for (int j=0; j<OUTER_quad_pts.size(); j++) {
			dlp=0;
	    calcp_fullWave(*p, OUTER_quad_pts[j], slp, dlp, zn);
			val_real=val_real+dlp*OUTER_quad_weights[j];
    }
		val_real=val_real/p->area()/peval->area();
		

    for (int i=0; i<numDipole; i++)
    {
      get_cx_panel(p, exps[i], cx_panel);
      val_cx_temp=0.;
			for (int j=0; j<OUTER_quad_pts.size(); j++) {
				dlp=0;
	      calcp_fullWave_cx(cx_panel, OUTER_quad_pts[j], slp, dlp, zn);
				val_cx_temp=val_cx_temp+dlp*OUTER_quad_weights[j];
    	}
			val_cx_temp=val_cx_temp/p->area()/peval->area();
      val_cx=val_cx+coes[i]*val_cx_temp;
    }
    return val_real+val_cx;
}


void Charge::get_cx_panel(const panel<double>* p,
                          const complex<double> cx_location,
                          panel<complex<double> >& cx_panel)
{
  double subHeight=dipoles->substrateHeight;
  vector<point3D<complex<double> > > cx_verts;

  
  cx_verts.resize(p->shape());

  for (int i=0; i<p->shape(); i++)
  {
    cx_verts[i]=point3D<complex<double> >(p->vertex(i).x(), p->vertex(i).y(),
                                          (2*subHeight-p->vertex(i).z())+IMAG*cx_location );
  }

  cx_panel.initialize(cx_verts);

}


void Charge::get_r_panel(const panel<double>* p,
                          const complex<double> cx_location,
                          panel<double> & r_panel)
{
  double subHeight=dipoles->substrateHeight;
  vector<point3D<double>  > r_verts;


  r_verts.resize(p->shape());

  for (int i=0; i<p->shape(); i++)
  {
    r_verts[i]=point3D<double> (p->vertex(i).x(), p->vertex(i).y(), p->vertex(i).z());
  }

  r_panel.initialize(r_verts);

}

/*charge *Charge::build_charges_list()
{
  charge *first_chg=NULL, *chg1=NULL, *chg2=NULL;
  panel<double>* p;
  int i;

  for (p=condStruc->list_panels.initialize(); p != NULL;
       p=condStruc->list_panels.iterate())
  {
    chg1 = new charge;
    if (chg2 != NULL) chg2->next = chg1;
    else first_chg=chg1;
    chg2 = chg1;


    for (i=0; i<4; i++)
    {
      chg1->corner[i][0] = (p->get_apex(i+1))->x();
      chg1->corner[i][1] = (p->get_apex(i+1))->y();
      chg1->corner[i][2] = (p->get_apex(i+1))->z();
    }

    chg1->shape = 4;
    chg1->next = NULL;
  }
  return first_chg;
}
*/
