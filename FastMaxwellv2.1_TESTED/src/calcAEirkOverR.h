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
//calcAEikrOverR handles computation of slp= \int_v \int V e^{-jkr}/r dr'dr
//dIdx= \int_v \int V d/dx(e^{-jkR}/R) dr'dr
// dIdy= \int_v \int V d/dy(e^{-jkR}/R) dr'dr
// r can be real or complex
// this program is a part of magnetic potential A computation

#ifndef CALCAEIKROVERR_H_
#define CALCAEIKROVERR_H_

#include "calcpForEikrOverR.h"
#include "calcaoneoverr.h"
#include "calcA_quadrature.h"
#include<iostream>


template <class T1, class T2>
class calcAEikrOverR
{
public:
  calcAEikrOverR(void){};
  calcAEikrOverR(const complex<double> K_num, const bool needSLP, const bool needDLP);
  void init(const complex<double> K_num, const bool needSLP, const bool needDLP);

  void selfterm (filament<T1> src,  complex<double>& slp);

  void operator() (filament<T1> src, filament<T2> eval,
                   complex<double>& slp, complex<double>& dIdx, complex<double>& dIdy);
  void XX_deriv (filament<T1> src, filament<T2> eval, int TYPE,
                   complex<double>& slp, complex<double>& dIdx, complex<double>& dIdy);

private:
  void compute_accelerated(filament<T1> src, filament<T2> eval, complex<double>& slp,
                           complex<double>& dIdx, complex<double>& dIdy);

  complex<double> K_;
  bool needSingleLayerPotential_;
  bool needDoubleLayerPotential_;
  bool needNothing_;
  calcpForEikrOverR<T1, T2> calcp_fw, calcp_XX;
  calcAQuadrature<T1, T2> calcA_quad;

};

template <class T1, class T2>
calcAEikrOverR<T1, T2>::calcAEikrOverR (
  const complex<double> K_num,
  const bool needSLP,
  const bool needDLP): K_(K_num), needSingleLayerPotential_(needSLP),
    needDoubleLayerPotential_(needDLP)
{
  needNothing_ = (!needSLP) && (!needDLP);
  if (! needNothing_)
  {
    calcp_fw.init(K_, needSLP, needDLP);
    calcp_XX.init(K_, needSLP, needDLP);
    calcA_quad.init(K_, needSLP, needDLP);
  }
}



template <class T1, class T2>
void calcAEikrOverR<T1, T2>::init (
  const complex<double> K_num,
  const bool needSLP,
  const bool needDLP)
{
  K_=(K_num);
  needSingleLayerPotential_=(needSLP);
  needDoubleLayerPotential_=(needDLP);

  needNothing_ = (!needSLP) && (!needDLP);
  if (! needNothing_)
  {
    calcp_fw.init(K_, needSLP, needDLP);
    calcp_XX.init(K_, needSLP, needDLP);
    calcA_quad.init(K_, needSLP, needDLP);
  }
}


template<class T1, class T2>
void calcAEikrOverR<T1, T2>::
operator()(filament<T1> src, filament<T2> eval,
           complex<double>& slp, complex<double>& dIdx,  complex<double>& dIdy)
{
  if (needNothing_) return;

  slp=0., dIdx=0., dIdy=0.;
//  calcA_quad(src, eval, slp, dIdx, dIdy);

  //   //determine if two filaments are the same
  //   if ( (src.loc0().x()==eval.loc0().x()) && (src.loc0().y()==eval.loc0().y())
  //        && (src.loc0().z()==eval.loc0().z()) && (src.loc1().x()==eval.loc1().x())
  //        && (src.loc1().y()==eval.loc1().y()) && (src.loc1().z()==eval.loc1().z())
  //        && (src.filHeight()==eval.filHeight()) && (src.filWidth()==eval.filWidth()))
  //   {
  //     //selfterm(src, slp);
  //     dIdx=0.;
  //     dIdy=0.;
  //   }
/*
	if (length(src.filCentroid()-eval.filCentroid())==0.) {
 	   selfterm(src, slp);
	}
else*/ {
	compute_accelerated(src, eval, slp, dIdx, dIdy);

	if ( needSingleLayerPotential_ )
	  slp=MU0*slp/(4*PI*src.filCrossA()*eval.filCrossA());
	if ( needDoubleLayerPotential_ )
	  {
	    dIdx=MU0*dIdx/(4*PI*src.filCrossA()*eval.filCrossA());
	    dIdy=MU0*dIdy/(4*PI*src.filCrossA()*eval.filCrossA());
	  }
  }
}



template<class T1, class T2>
void calcAEikrOverR<T1, T2>::compute_accelerated(filament<T1> src, filament<T2> eval,
    complex<double>& slp, complex<double>& dIdx, complex<double>& dIdy)
{
  size_t numPtsW=1;
  size_t numPtsH=1;
  size_t numPtsL=8;
  vector<point3D<T2> > quad_pts;
  vector<T2> quad_weights;
  vector<panel<T1> > facePanels;

  complex<double> slpLocal=0., dlpLocal=0., zn=0.;
  dIdx=slp=dIdy=0.;

  //obtain volume quadrature pts for the eval. filament--checked
  quad_pts.resize(numPtsW*numPtsH*numPtsL);
  quad_weights.resize(numPtsW*numPtsH*numPtsL);
  eval.volume_quadrature(numPtsW, numPtsH, numPtsL, quad_pts, quad_weights);

  
  //     cout<<"printing quad info"<<endl;
  //     for (int i=0; i< quad_pts.size(); i++)
  //       quad_pts[i].print_info();
  //     for (int i=0; i< quad_weights.size(); i++)
  //   	  cout<<quad_weights[i]<<endl;
  
  //obtain face panels for src. filament--checked
  facePanels=src.face_panels();
  //cout << facePanels.size() << endl;
  //   for (int i=0; i<facePanels.size(); i++)
  //     facePanels[i].print_info();
  
  
  for (int i=0; i<quad_pts.size(); i++)
    //  for (int i=0; i<1; i++)
    {
      for (int j=0; j<facePanels.size(); j++)
	// 	  for (int j=0; j<1; j++)
	{
	  // 		facePanels[j].print_info();
	  // 		quad_pts[i].print_info();
	  calcp_fw(facePanels[j], quad_pts[i], slpLocal, dlpLocal, zn);
	  //cout << slpLocal << dlpLocal << quad_weights[i] << facePanels[j].normal().x() << endl;
	  //cout<<"i "<<i+1<<" j "<<j+1<<" "<<slpLocal<<" "<<dlpLocal<<" "<<zn<<endl;
	  
	  if ( needSingleLayerPotential_ )
	    slp=slp+slpLocal*zn*quad_weights[i];
	  if ( needDoubleLayerPotential_ )
	    {
	      //cout<<dlpLocal<<" "<<facePanels[j].normal().x()<<" "<<facePanels[j].normal().y()<<endl;
	      dIdx=dIdx-dlpLocal*facePanels[j].normal().x()*quad_weights[i];
	      dIdy=dIdy-dlpLocal*facePanels[j].normal().y()*quad_weights[i];
	    }
	  //cout << "slpLocal" << slp <<"dlpLocal=" << dIdx << dIdy << endl;
	}
    }
}



template <class T1, class T2>
void calcAEikrOverR<T1, T2>::selfterm(filament<T1> src, complex<double>& slp)
{
  //   calcAOneOverR calcA_static;
  //   double wavelength=2*PI/abs(K_);
  //   int numDiv=(int) ceil(abs(src.filLength())/(wavelength/8))+10;
  //   cout<<numDiv<<endl;
  //   numDiv=MIN(numDiv,10);
  //   T1 dis;
  //
  //
  //   vector<filament<T1> > fils;
  //   fils.resize(numDiv);
  //   src.divFil(numDiv, fils);
  //   slp=0.;
  //
  //
  //   for (int i=0; i<numDiv; i++)
  //   {
  //     for (int j=0; j<numDiv; j++)
  //     {
  //       if (i==j)
  //       {
  //         slp+=calcA_static.selfterm(&fils[i]);
  //       }
  //       else
  //       {
  //         dis=length((fils[i].loc0()-fils[i].loc1())/2.+fils[i].loc0(),
  //                    (fils[j].loc0()-fils[j].loc1())/2.+fils[j].loc0());
  //         slp+=calcA_static.mutual(&fils[i], &fils[j])*exp(-IMAG*K_*dis);
  //       }
  //     }
  //   }

  slp=0.;
  calcA_quad.selfterm(src,slp);

}


template<class T1, class T2>
void calcAEikrOverR<T1, T2>::
XX_deriv (filament<T1> src, filament<T2> eval, int TYPE,
           complex<double>& dIdxx, complex<double>& dIdyy,  complex<double>& dIdxy)
{
  if (needNothing_) return;
  size_t numPtsW=2;
  size_t numPtsH=2;
  vector<panel<T1> > facePanels;
  panel<T1> temp;
  vector<point3D<T1> > OUTER_quad_pts, INNER_quad_pts;
  vector<T1> OUTER_quad_weights, INNER_quad_weights;
  vector<panel<T2> > EfacePanels;
	panel<T2> Etemp;

  complex<double> slpLocal=0., dlpLocal=0., zn=0., val=0., R=0.;
  complex<double> dIdyx=0., dIdXY=0., dIdYX;

  dIdxx=0., dIdyy=0., dIdxy=0.;
  
  OUTER_quad_pts.resize(numPtsW*numPtsH);
  OUTER_quad_weights.resize(numPtsW*numPtsH);
  
  EfacePanels=eval.face_panels();
  facePanels=src.face_panels();
  if (TYPE==1) {
    for (int k=0; k<EfacePanels.size() ; k++) {
      if (EfacePanels[k].normal().x()!=0.) {
	Etemp=EfacePanels[k];
	Etemp.panel_quadrature(numPtsW * numPtsH, OUTER_quad_pts, OUTER_quad_weights);
	for (int j=0; j<facePanels.size()  ; j++) {
	  if (facePanels[j].normal().x()!=0.) {
	    val=0.;
	    for (int i=0; i<OUTER_quad_pts.size(); i++) {
	      calcp_XX(facePanels[j], OUTER_quad_pts[i], slpLocal, dlpLocal, zn);
	      val=val+dlpLocal*OUTER_quad_weights[i];
	    }
	    dIdxx=dIdxx+val*facePanels[j].normal().x()*EfacePanels[k].normal().x();
	  } }
      } }
  } else if (TYPE==2) {
    for (int k=0; k<EfacePanels.size() ; k++) {
      if (EfacePanels[k].normal().y()!=0.) {
	Etemp=EfacePanels[k];
	Etemp.panel_quadrature(numPtsW * numPtsH, OUTER_quad_pts, OUTER_quad_weights);
	for (int j=0; j<facePanels.size()  ; j++) {
	  if (facePanels[j].normal().y()!=0.) {
	    val=0.;
	    for (int i=0; i<OUTER_quad_pts.size(); i++) {
	      calcp_XX(facePanels[j], OUTER_quad_pts[i], slpLocal, dlpLocal, zn);
	      val=val+dlpLocal*OUTER_quad_weights[i];
		}
	    dIdyy=dIdyy+val*facePanels[j].normal().y()*EfacePanels[k].normal().y();
	  } }
      } }
  } else if (TYPE==3) {
    for (int k=0; k<EfacePanels.size() ; k++) {
      if (EfacePanels[k].normal().y()!=0.) {
	Etemp=EfacePanels[k];
	Etemp.panel_quadrature(numPtsW * numPtsH, OUTER_quad_pts, OUTER_quad_weights);
	for (int j=0; j<facePanels.size()  ; j++) {
	  if (facePanels[j].normal().x()!=0.) {
	    val=0.;
	    for (int i=0; i<OUTER_quad_pts.size(); i++) {
	      calcp_XX(facePanels[j], OUTER_quad_pts[i], slpLocal, dlpLocal, zn);
	      val=val+dlpLocal*OUTER_quad_weights[i];
	    }
	    dIdxy=dIdxy+val*facePanels[j].normal().x()*EfacePanels[k].normal().y();
	  } }
      } }
  }
  
  dIdxx=MU0*dIdxx/(4*PI*src.filCrossA()*eval.filCrossA());
  dIdyy=MU0*dIdyy/(4*PI*src.filCrossA()*eval.filCrossA());
  dIdxy=MU0*dIdxy/(4*PI*src.filCrossA()*eval.filCrossA());
  
  
  //  for (int k=0; k<EfacePanels.size() ; k++) {
  //	if (EfacePanels[k].normal().x()!=0.) {
////	for (int kk=0; kk<EfacePanels[k].shape(); kk++) {
////  	cout << EfacePanels[k].vertex(kk).x() << EfacePanels[k].vertex(kk).y() << EfacePanels[k].vertex(kk).z() << endl;
////	}
//	Etemp=EfacePanels[k];
//	Etemp.panel_quadrature(numPtsW * numPtsH, OUTER_quad_pts, OUTER_quad_weights);
//	for (int j=0; j<facePanels.size()  ; j++) {
//	if (facePanels[j].normal().y()!=0.) {
//		val=0.;
//		for (int i=0; i<OUTER_quad_pts.size(); i++) {
//		//cout << OUTER_quad_pts[i].x()<< OUTER_quad_pts[i].y() << OUTER_quad_pts[i].z() << endl;
//		calcp_XX(facePanels[j], OUTER_quad_pts[i], slpLocal, dlpLocal, zn);
//		val=val+dlpLocal*OUTER_quad_weights[i];
//		}
//		cout << facePanels[j].normal().y()*EfacePanels[k].normal().x() << "YX" << dIdxy << endl;
//		dIdYX=dIdYX+val*facePanels[j].normal().y()*EfacePanels[k].normal().x();
//	} }
//	} }
//  for (int k=0; k<EfacePanels.size() ; k++) {
//	if (EfacePanels[k].normal().x()!=0.) {
//	for (int j=0; j<facePanels.size(); j++) {
//	if (facePanels[j].normal().y()!=0.) {
//		Etemp=EfacePanels[k];
//		Etemp.panel_quadrature(numPtsW * numPtsH, OUTER_quad_pts, OUTER_quad_weights);
//    temp=facePanels[j];
//		temp.panel_quadrature(numPtsW * numPtsH, INNER_quad_pts, INNER_quad_weights);
//		val=0.;
//		for (int i=0; i<OUTER_quad_pts.size(); i++) {
//			//cout << OUTER_quad_pts[i].x()<< OUTER_quad_pts[i].y() << OUTER_quad_pts[i].z() << endl;
//			for (int ii=0; ii<INNER_quad_pts.size(); ii++) {
//    		R=length(OUTER_quad_pts[i]-INNER_quad_pts[ii]);
//    		val=val
//+(exp(-1.*IMAG* K_*R)/R)*OUTER_quad_weights[i]*INNER_quad_weights[ii];
//			}
//		}
//	cout << facePanels[j].normal().y()*EfacePanels[k].normal().x() << "XY" << val << endl;
//	dIdxy=dIdxy+val*facePanels[j].normal().y()*EfacePanels[k].normal().x();
//		} }
//	} }
//
//  for (int k=0; k<EfacePanels.size() ; k++) {
//	if (EfacePanels[k].normal().y()!=0.) {
//	for (int j=0; j<facePanels.size(); j++) {
//	if (facePanels[j].normal().x()!=0.) {
//		Etemp=EfacePanels[k];
//		Etemp.panel_quadrature(numPtsW * numPtsH, OUTER_quad_pts, OUTER_quad_weights);
//    temp=facePanels[j];
//		temp.panel_quadrature(numPtsW * numPtsH, INNER_quad_pts, INNER_quad_weights);
//		val=0.;
//		for (int i=0; i<OUTER_quad_pts.size(); i++) {
//			//cout << OUTER_quad_pts[i].x()<< OUTER_quad_pts[i].y() << OUTER_quad_pts[i].z() << endl;
//			for (int ii=0; ii<INNER_quad_pts.size(); ii++) {
//    		R=length(OUTER_quad_pts[i]-INNER_quad_pts[ii]);
//    		val=val
//+(exp(-1.*IMAG* K_*R)/R)*OUTER_quad_weights[i]*INNER_quad_weights[ii];
//			}
//		}
//	cout << facePanels[j].normal().x()*EfacePanels[k].normal().y() << "YX" << val << endl;
//	dIdyx=dIdyx+val*facePanels[j].normal().x()*EfacePanels[k].normal().y();
//		} }
//	} }
////  for (int k=0; k<EfacePanels.size() ; k++) {
////	if (EfacePanels[k].normal().x()!=0.) {
////	for (int j=0; j<facePanels.size(); j++) {
////	if (facePanels[j].normal().y()!=0.) {
////		EfacePanels[k].panel_quadrature(numPtsW * numPtsH, OUTER_quad_pts, OUTER_quad_weights);
////  	facePanels[j].panel_quadrature(numPtsW * numPtsH, INNER_quad_pts, INNER_quad_weights);
////		val=0.;
////		for (int i=0; i<OUTER_quad_pts.size(); i++) {
////			for (int ii=0; ii<INNER_quad_pts.size(); ii++) {
////    		R=length(OUTER_quad_pts[i]-INNER_quad_pts[ii]);
////    		val=val
////+(exp(-1.*IMAG* K_*R)/R)*OUTER_quad_weights[i]*INNER_quad_weights[ii];
////			}
////		}
////	cout << facePanels[j].normal().y()*EfacePanels[k].normal().x() << "XY" << val << endl;
////	dIdxy=dIdxy+val*facePanels[j].normal().y()*EfacePanels[k].normal().x();
////		} }
////	} }
//	cout << "XIN" << dIdYX << "XIN" << dIdXY << "XY " << dIdxy << " YX" << dIdyx << endl;

}


//template<class T1, class T2>
//complex<double> calcAEikrOverR<T1, T2>::quad_panel(panel<T1>& src_p, panel<T1>& eval_p)
//
//{
//  vector<point3D<T1> > quad_pts;
//  vector<T1> quad_weights;
//  T1 R;
//  complex<double> val=0.;
//
//
//
//
//}

#endif
