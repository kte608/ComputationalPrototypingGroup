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
#ifndef PANEL_H_
#define PANEL_H_

#include "const_def.h"
#include "point3D.h"
#include <vector>
#include <complex>
#include <cfloat>
#include "../pfft/vector3D.h"

using namespace std;
using namespace point_vec_operations;


template <class T>
class panel
{
public:

  typedef  point3D<T> vector3D;

  panel(void){};

  //initialze vertices for a quad. panel
  panel(const point3D<T> &,
        const point3D<T> &,
        const point3D<T> &,
        const point3D<T> &
       );
  //initialze vertices for a triangular panel
  panel(const point3D<T> &,
        const point3D<T> &,
        const point3D<T> &
       );



  //setting up functions that compute panel properties
  void initialize(const vector<point3D<T> >&);
  void setup_panel(void);
  void comp_centroid_area(void);
  void comp_normal_tangent(void);
  void setdirection(vector3D direction) {direction_=direction;;}



  //retrieval functions
  int panel_index(void) const { return globalPanelIndex_;};
  const size_t shape(void) const { return vertex_.size();};
  const point3D<T>* get_apex(int num) const {return &( vertex_[num-1]); };
  const point3D<T> vertex(int num) const {return vertex_[num];};
  const vector3D tangent1(void) const { return tangent1_;};
  const vector3D tangent2(void) const {return tangent2_;};
  const vector3D normal(void) const {return normal_;};
  const vector3D direction(void) const {return direction_;};
  const point3D<T> centroid(void) const {return centroid_;};
  const T area(void) const {return area_; };
  void assign_number(int num)  {globalPanelIndex_=num; };
  int get_number(void)const {return globalPanelIndex_;};

  void shiftToLocalCoord ( const point3D<T>& localOrigin);
  void print_info();
  void gaussQuadrature(const size_t n,
                       std::vector<double>& w,
                       std::vector<double>& x);
  void panel_quadrature (const size_t numPts,
                         vector<point3D<T> >& pts, vector<T>& weights);
  void get_globel_vertices(vector< pfft::vector3D<T> >& lovalVerts);

private:
  T comp_triangle_area( const point3D<T>& ,
                        const point3D<T>& ,
                        const point3D<T>& );
  bool less_than (complex<double>, complex<double>);
  bool less_than (double, double);

  vector<point3D<T> > vertex_;
  int globalPanelIndex_;
  vector3D tangent1_;
  vector3D tangent2_;
  point3D<T> centroid_;
  T area_;
  vector3D normal_, direction_;

};




/*******************************************************
 * panel constructors -- for triangular and quad. panels
 *******************************************************/
template <class T>
panel<T>::panel (const point3D<T>& v1,
                 const point3D<T>& v2,
                 const point3D<T>& v3 )
{
  vertex_.resize(3);
  vertex_[0]=v1;
  vertex_[1]=v2;
  vertex_[2]=v3;


  setup_panel();
}

template <class T>
panel<T>::panel (const point3D<T>& v1,
                 const point3D<T>& v2,
                 const point3D<T>& v3,
                 const point3D<T>& v4)
{
  vertex_.resize(4);
  vertex_[0]=v1;
  vertex_[1]=v2;
  vertex_[2]=v3;
  vertex_[3]=v4;

  setup_panel();
}

template <class T>
void panel<T>::initialize(const vector<point3D<T> >& verts)
{
  vertex_.resize(verts.size());
  for (int i=0; i<verts.size(); i++)
    vertex_[i]=verts[i];
  setup_panel();
}

template <class T>
void panel<T>::setup_panel(void)
{
  comp_centroid_area();
  comp_normal_tangent();
}


template <class T>
void panel<T>::comp_centroid_area(void)
{
  T triangleArea;
  point3D<T> triangleCentroid;

  area_=0.;
  centroid_=point3D<T>(0.,0.,0.);
  for (size_t i=1; i<shape()-1; i++)
  {
    triangleArea=comp_triangle_area(vertex_[0], vertex_[i], vertex_[i+1]);
    triangleCentroid=((vertex_[0]+vertex_[i])+vertex_[i+1])/3;
    area_+=triangleArea;
    centroid_+=triangleCentroid*triangleArea;
  }

  if (area_==0.)
  {
    //pfft::errorMessage("element.cc:compCentroidAndArea",
    //"Panel area is zero! There must be a bug here!!");
  }
  else
  {
    centroid_ /= area_;
  }
}


template <class T>
T panel<T>::comp_triangle_area(
  const point3D<T>& v1,
  const point3D<T>& v2,
  const point3D<T>& v3)
{
  vector3D edge1=v2-v1;
  vector3D edge2=v2-v3;
  vector3D vec=crossProd(edge1, edge2);
  return 0.5*length(vec);
}


template <class T>
void panel<T>::comp_normal_tangent(void)
{
  if (shape()==4)
  {
    point3D<T> midPoint1 = (vertex_[0] + vertex_[1]) / 2.;
    point3D<T> midPoint2 = (vertex_[2] + vertex_[3]) / 2.;
    tangent1_ = midPoint2 - midPoint1;
    T edgeMidPointDist1_ = length(tangent1_);


    midPoint1 = (vertex_[0] + vertex_[3]) / 2.;
    midPoint2 = (vertex_[1] + vertex_[2]) / 2.;
    tangent2_ = midPoint2 - midPoint1;
    T edgeMidPointDist2_ = length(tangent2_);

    /* ((abs(edgeMidPointDist1_) < EPS) || (abs(edgeMidPointDist2_)< EPS) ) {
         std::cout << endl
    << "\t element.cc:compNormalAndTangent" << endl
    << "\t edgeMidPointDist is zero! There must be a bug here!!" 
    << endl;
      throw domain_error("edgeMidPointDist is zero! There must be a bug here!!");   
      
      }*/

  }
  else
  {
    tangent1_ = vertex_[1] - vertex_[0];
    tangent2_ = vertex_[1] - vertex_[2];
  }

  tangent1_.normalize();
  normal_ = crossProd(tangent1_, tangent2_);
  normal_.normalize();
  tangent2_ = crossProd(normal_, tangent1_);
  tangent2_.normalize();
}

template <class T>
void panel<T>::shiftToLocalCoord ( const point3D<T>& localOrigin)
{
  for(size_t ii=0; ii<vertex_.size(); ii++)
  {
    vertex_[ii].transferGlobalToLocalCoord(localOrigin, tangent1_, tangent2_, normal_);
  }
}

template<class T>
void panel<T>::panel_quadrature (const size_t numPts,
                                 vector<point3D<T> >& pts, vector<T>& weights)
{

	shiftToLocalCoord(centroid_);

  //note, panel is already shifted to local coord. system
  //perform quadrature on the already transformed panel
  T L1=length(vertex(3)-vertex(0));
  T L2=length(vertex(1)-vertex(0));

  T startX = vertex(0).x();
  if ( less_than(vertex(3).x() , vertex(0).x()))
    startX=  vertex(3).x();

  T startY = vertex(0).y();
  if ( less_than(vertex(1).y() ,  vertex(0).y()))
    startY=  vertex(1).y();

  double ratio =abs(L1/L2);
  int numQuadPtsY=(int)ceil(sqrt(1.*numPts/ratio));
  if ((numQuadPtsY % 2)!=0)
    numQuadPtsY++;
  int numQuadPtsX=(int)ceil(1.*numPts/numQuadPtsY);
  if ((numQuadPtsX % 2)!=0)
    numQuadPtsX++;


  vector<double> wY, posY;
  wY.resize(numQuadPtsY);
  posY.resize(numQuadPtsY);
  gaussQuadrature (numQuadPtsY, wY, posY);
  vector<double> wX, posX;
  wX.resize(numQuadPtsX);
  posX.resize(numQuadPtsX);
  gaussQuadrature (numQuadPtsX, wX, posX);

  T pX, weightX, pY, weightY;
  pts.resize(numQuadPtsY*numQuadPtsX);
  weights.resize(numQuadPtsY*numQuadPtsX);
  int count=0;
  for (int i=0; i< numQuadPtsX; i++)
  {
    pX=L1*posX[i]+startX;
    weightX=L1*wX[i];
    for (int j=0; j< numQuadPtsY; j++)
    {
      pY=L2*posY[j]+startY;
      weightY=L2*wY[j];
      pts[count]=point3D<T>(pX, pY, 0.);
      weights[count]=weightX*weightY;
      count++;
    }
  }

//  trasfer pts to global coords
  for (int i=0; i<pts.size(); i++)
    pts[i].transferLocalToGlobalCoord(centroid_, tangent1_, tangent2_, normal_);

}


template<class T>
void panel<T>::gaussQuadrature (
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

template<class T>
void panel<T>::get_globel_vertices(vector<pfft::vector3D<T> >& localVerts)
{

	for (int j=0; j<shape(); j++){
    localVerts.push_back(pfft::vector3D<T>(vertex(j).x(), vertex(j).y(), vertex(j).z()));
	}
}


template <class T>
void panel<T>::print_info()
{
  for(size_t ii=0; ii<vertex_.size(); ii++)
  {
    cout<<vertex_[ii].x()<<" "<<vertex_[ii].y()<<" "<<vertex_[ii].z()<<endl;
  }
}

template <class T>
bool panel<T>::less_than (complex<double> A, complex<double> B)
{
  return (real(A) < real(B));
}

template <class T>
bool panel<T>::less_than (double A, double B)
{
  return (A<B);
}


#endif
