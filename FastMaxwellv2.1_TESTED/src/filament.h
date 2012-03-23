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
#ifndef FILAMENT_H_
#define FILAMENT_H_

#include "point3D.h"
#include <cfloat>
#include <vector>
#include <iostream>
#include "math.h"
#include "const_def.h"
//#include "henry.h"
#include "panel.h"
#include "../pfft/vector3D.h"

using namespace std;
using namespace point_vec_operations;


template <class T>
class filament
{

public:
  typedef point3D<T> vector3D;

  filament(void) {};
  //filament(const point3D<T> &,
    //       const point3D<T> &,
      //     T, T, double, SEGMENT*, FILAMENT*);
  filament(const point3D<T> &,
           const point3D<T> &,
           T, T, double);

  //setting up functions that compute filament properties
  void setup_filament(void);
  void comp_axes_directions(void);
  void comp_vol_crossA(void);
  void initialize(const point3D<T> &,
                  const point3D<T> &,
                  T, T, double);

  //retrieval functions
  const T filLength (void) const {return filLength_;};
  const T filWidth (void) const {return filWidth_;};
  const T filHeight(void) const {return filHeight_;};
  const T filCrossA (void) const {return filCrossA_;};
  const T filVolume(void) const {return filVol_;};
  const double sigma(void) const {return filSigma_;};
  //SEGMENT *get_henry_seg_ptr (void) const { return henry_seg_ptr_; };
  //FILAMENT *get_henry_fil_ptr (void) const { return henry_fil_ptr_; };

  const point3D<T>* get_point(int num)const;
  const point3D<T> loc0(void) const {return endPt0_;};
  const point3D<T> loc1(void) const {return endPt1_;};
  const point3D<T> filCentroid(void) const {return centroid_;};

  const vector3D l_dir (void) const {return dirL_;};
  const vector3D w_dir (void) const {return dirW_;};
  const vector3D h_dir (void) const {return dirH_;};
  const vector<panel<T> >  face_panels(void) const {return facePanels;};



  void assign_number(int num)   {filIndex_=num;};
  int get_number () const { return filIndex_; };
  void set_endpt0(point3D<T> pt) {endPt0_=pt;};
  void set_endpt1(point3D<T> pt) {endPt1_=pt;};
  void set_width(T width) {filWidth_=width;};
  void set_height(T height) {filHeight_=height;};
  void set_sigma(double sigma) {filSigma_=sigma;};
  void set_length(T length_) {filLength_=length_;};

  void obtain_face_panels(void);
  void volume_quadrature(const size_t numPtsW,const size_t numPtsH, const size_t numPtsL,
                         vector<point3D<T> >& pts, vector<T>& weights);
  void divFil(const int numDiv, vector<filament<T> >& fils);
  void get_globel_vertices(vector< pfft::vector3D<T> >& lovalVerts);
 /* void get_cx_fil(const filament<double> f,
                  const complex<double> cx_location,
                  filament<complex<double> >& cx_fil, bool isReal);
   */
 private:
  //int dummy;
  void get_local_vertices(vector<point3D<T> >& localVerts);
  void get_local_quadraturePts(vector<point3D<T> >& fil_local_pts, 
			       vector<T>& fil_local_weights, 
			       const size_t numPtsW, const size_t numPtsH, 
			       const size_t numPtsL);
  
  void gaussQuadrature ( const size_t n, std::vector<double>& w,
                         std::vector<double>& x);
  
  point3D<T> endPt0_, endPt1_, centroid_;
  //double filWidth_, filHeight_,  filLength_;
  double filSigma_;
  T filWidth_, filHeight_, filLength_;
  vector3D dirL_, dirW_, dirH_;
  T filCrossA_, filVol_;
  int filIndex_;
  //SEGMENT *henry_seg_ptr_;
  //FILAMENT *henry_fil_ptr_;
  vector<panel<T> > facePanels;

};


/*
template <class T>
filament<T>::filament (const point3D<T> & endPt0,
                       const point3D<T> & endPt1,
                       T width, T height, double sigma,
                       SEGMENT *henry_seg_ptr, FILAMENT *henry_fil_ptr )
{
  endPt0_=endPt0;
  endPt1_=endPt1;
  filWidth_=width;
  filHeight_=height;
  filSigma_=sigma;
  setup_filament();
  henry_seg_ptr_=henry_seg_ptr;
  henry_fil_ptr_=henry_fil_ptr;
}
*/
template <class T>
filament<T>::filament (const point3D<T> & endPt0,
                       const point3D<T> & endPt1,
                       T width, T height, double sigma)
{
  endPt0_=endPt0;
  endPt1_=endPt1;
  filWidth_=width;
  filHeight_=height;
  filSigma_=sigma;
 // henry_seg_ptr_=NULL;
 // henry_fil_ptr_=NULL;
  setup_filament();
}

template<class T>
void filament<T>::initialize (const point3D<T> & endPt0,
                              const point3D<T> & endPt1,
                              T width, T height, double sigma)
{
  endPt0_=endPt0;
  endPt1_=endPt1;
  filWidth_=width;
  filHeight_=height;
  filSigma_=sigma;
//  henry_seg_ptr_=NULL;
//  henry_fil_ptr_=NULL;
  setup_filament();

}




template<class T>
const point3D<T>* filament<T>::get_point(int num)const
{
  if (num==1)
    return &endPt0_;
  else if (num==2)
    return &endPt1_;
}


template <class T>
void filament<T>::setup_filament(void)
{
  filCrossA_=filWidth_*filHeight_;
  filLength_=length(endPt0_, endPt1_);
  centroid_= (endPt1_-endPt0_)/2+endPt0_;
  filVol_=filWidth_*filLength_*filHeight_;
  dirL_= (endPt1_-endPt0_)/filLength_;

  point3D<double> GLOBAL_Z_AXIS (0,0,1.);
  dirW_= cross(GLOBAL_Z_AXIS, dirL_);
  if ((abs(dirW_.x()/filLength_) < EPS )
      && (abs(dirW_.y()/filLength_) < EPS))
  {
    dirW_.set_x(1.);
    dirW_.set_y(0.);
  }

  dirW_.normalize();

  dirH_=cross(dirW_, dirL_);
  dirH_*=-1;
  dirH_.normalize();

  obtain_face_panels();
}

template <class T>
void filament<T>::obtain_face_panels(void)
{
  vector<point3D<T> > vertices;

  get_local_vertices(vertices);
  for (int i=0; i<vertices.size(); i++){
	  vertices[i].transferLocalToGlobalCoord(endPt0_, dirW_, dirH_, dirL_);
  }
//   for (int i=0; i<vertices.size(); i++)
// 	  cout<<vertices[i].x()<<" "<<vertices[i].y()<<" "<<vertices[i].z()<<endl;
	  
  //cout<<vertices.size()<<endl;
  facePanels.resize(6);
  panel<T> facePanel1= panel<T> (vertices[3], vertices[1], vertices[0], vertices[2]);
  facePanels[0]=(facePanel1);
  panel<T> facePanel2= panel<T> (vertices[4], vertices[5], vertices[7], vertices[6]);
  facePanels[1]=(facePanel2);
  panel<T> facePanel3= panel<T> (vertices[5], vertices[4], vertices[0], vertices[1]);
  facePanels[2]=(facePanel3);
  panel<T> facePanel4= panel<T> (vertices[6], vertices[7], vertices[3], vertices[2]);
  facePanels[3]=(facePanel4);
  panel<T> facePanel5= panel<T> (vertices[0], vertices[4], vertices[6], vertices[2]);
  facePanels[4]=(facePanel5);
  panel<T> facePanel6= panel<T> (vertices[3], vertices[7], vertices[5], vertices[1]);
  facePanels[5]=(facePanel6);
}

template <class T>
		void filament<T>:: volume_quadrature(const size_t numPtsW, 
											 const size_t numPtsH, const size_t numPtsL,
                                     vector<point3D<T> >& pts, vector<T>& weights)
{

  get_local_quadraturePts(pts, weights, numPtsW, numPtsH, numPtsL);
  for (int i=0; i<pts.size(); i++)
	  pts[i].transferLocalToGlobalCoord(endPt0_, dirW_, dirH_, dirL_);
//     for (int i=0; i<pts.size(); i++)
//   	  cout<<i<<" "<<pts[i].x()<<" "<<pts[i].y()<<" "<<pts[i].z()<<endl;
  

}

template <class T>
void filament<T>:: get_local_quadraturePts(vector<point3D<T> >& fil_local_pts,
    vector<T >& fil_local_weights,
	const size_t numPtsW, const size_t numPtsH, const size_t numPtsL)
{
  vector<double> weightsW, posW, weightsL, posL, weightsH, posH;
  
  weightsW.resize(numPtsW);
  posW.resize(numPtsW);
  gaussQuadrature(numPtsW, weightsW, posW);
  
  weightsH.resize(numPtsH);
  posH.resize(numPtsH);
  if(numPtsW!=numPtsH){
	  gaussQuadrature(numPtsH, weightsH, posH);
  }	  else{
	  weightsH=weightsW;
	  posH=posW;
  }
  
  weightsL.resize(numPtsL);
  posL.resize(numPtsL);
  gaussQuadrature(numPtsL,weightsL,posL);

  T wStart=-filWidth_/2.;
  T hStart=-filHeight_/2.;
  T localL, localW, localH, weightL, weightW, weightH;

  int count=0;
  for (int i=0; i<numPtsL; i++)
  {
    localL=posL[i]*filLength_;
    weightL=weightsL[i]*filLength_;
    for (int j=0; j<numPtsW; j++)
    {
      localW=posW[j]*filWidth_+wStart;
      weightW=weightsW[j]*filWidth_;

      for (int k=0; k<numPtsH; k++)
      {
        localH=posH[k]*filHeight_+hStart;
        weightH=weightsH[k]*filHeight_;

        fil_local_pts[count]=(point3D<T>(localW, localH, localL));
        fil_local_weights[count]=(weightL*weightW*weightH);
		count++;
      }
    }
  }
}


template <class T>
void filament<T>:: get_local_vertices(vector<point3D<T> >& localVerts)
{
  T wStart=-filWidth_/2.;
  T hStart=-filHeight_/2.;
  T localX, localY, localZ;

  
  
  for (int i=0; i<2; i++)
  {
	  localZ=static_cast<double>(i)*filLength_;
    for (int j=0; j<2; j++)
    {
		localX=static_cast<double>(j)*filWidth_+wStart;
      for (int k=0; k<2; k++)
      {
		  localY=static_cast<double>(k)*filHeight_+hStart;
		
        localVerts.push_back(point3D<T>(localX, localY, localZ));
      }
    }
  }

}


template<class T>
void filament<T>::get_globel_vertices(vector<pfft::vector3D<T> >& localVerts)
{
	vector<point3D<T> > tmp;

	get_local_vertices(tmp);
  for (int i=0; i<tmp.size(); i++){
	  tmp[i].transferLocalToGlobalCoord(endPt0_, dirW_, dirH_, dirL_);
	//	cout << tmp[i].x() << " " << tmp[i].y() << " " << tmp[i].z() << endl;
    localVerts.push_back(pfft::vector3D<T>(tmp[i].x(), tmp[i].y(), tmp[i].z()));
  }
}





template<class T>
void filament<T>::divFil(const int numDiv, vector<filament<T> >& fils)
{

	
	T l_start=0., l_end=0.;
	T division= filLength_ /((double)numDiv);
	for (int i=0; i<numDiv; i++)
	{
 		point3D <T> endpt0 (0.,0.,l_start);
 		endpt0.transferGlobalToLocalCoord(point3D<T>(0.,0.,0.),dirW_, dirH_ , dirL_);
		//endpt0.print_info();
		
		l_end=l_start+division;
		//cout<<l_end<<endl;	
		point3D<T> endpt1(0.,0.,l_end);
		//endpt1.print_info();
		endpt1.transferGlobalToLocalCoord(point3D<T>(0.,0.,0.),dirW_, dirH_ , dirL_);
		//endpt1.print_info();
		
		l_start=l_end;
 		fils[i]=filament<T> (endpt0, endpt1, filWidth_, filHeight_, filSigma_);	
	}
}

template<class T>
void filament<T>::gaussQuadrature (
  const size_t n, // number of points
  std::vector<double>& w,  // weights
  std::vector<double>& x)   // points
{
  double x1 = 0;
  double x2 = 1;

  int m = (n+1)/2;
  double xm = .5 * (x2 + x1);
  double xl = .5 * (x2 - x1);
  double pp;

  for (int i = 1; i <= m; i++)
  {
    // Initial value of the ith root
    double z = cos(PI * (i-.25) / (n+.5));

    // Get the ith root by Newton method
    bool converge = false;
    do
    {
      // Evaluate Leg polynomial. p1 is the value at z
      double p1 = 1.0;
      double p2 = 0.0;
      for (int j=1; j <= n; j++)
      {
        double p3 = p2;
        p2 = p1;
        p1 = ((2.0*j-1.0) * z * p2 - (j-1.0) * p3) / j;
      }

      // pp is the derivative
      pp = n * (z*p1 - p2) / (z*z - 1.0);
      double z1 = z;
      z = z1 - p1/pp;

      if (abs(z-z1) < 1e1*DBL_EPSILON )
      {
        converge = true;
      }

    }
    while( !converge);

    x[i-1] = xm - xl*z;
    x[n+1-i-1] = xm + xl*z;
    w[i-1] = 2. * xl / ((1 - z*z) * pp * pp);
    w[n+1-i-1] = w[i-1];
  }
}


#endif
