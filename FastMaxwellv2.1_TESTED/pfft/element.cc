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

  ==========================================================================
*/

const static char cvsid[] = "$Id: element.cc,v 1.7 2003/02/11 03:10:26 zhzhu Exp $";

#include <fstream>
#include <string>
#include <stdexcept> //for exception handling
#include "element.h"
#include <cfloat> // for DBL_EPSILON
#include "utils.h" //for errorMessage()

using namespace std;
using namespace pfft;


element::element (
		  filament<double>* f,
		  const int name,
		  const int boundaryIndex,
		  const int boundaryType)
  : name_(name), boundaryIndex_(boundaryIndex), boundaryType_(boundaryType)
{
  fil_= f;
  fil_->get_globel_vertices(vertex_);
  boundaryIndex_=-1;
  setupElement();
}

element::element (
		  panel<double>* p,
		  const int name,
		  const int boundaryIndex,
		  const int boundaryType)
  : name_(name), boundaryIndex_(boundaryIndex), boundaryType_(boundaryType)
{
  panel_= p;
  panel_->get_globel_vertices(vertex_);
  boundaryIndex_=-2;
  setupElement();
}


/**********************************************************************
 * setupElement --
 **********************************************************************/
void
element::setupElement (
		       void)
{
  compCentroidAndArea();
  compBoundingSphereRadius();
}

/**********************************************************************
 * compCentroidAndArea --
 * The standard algorithm to calculate the area of a ploygon is to break it
 * up into triangles and sum up the area of each triangle.
 *
 * The standard algorithm to calculate the centroid, the center of gravity,
 * is to compute the centroid of each triangle and take the weighted average
 * of them. The weighting factor is the area of each triangle.
 *
 * For the sake of efficiency, these two functions are combined here.
 *
 * It is assumed here that the logically adjacent vertices are 
 * physically adjacent.
 * For example, vertex[0] is physically connected to vertex[1], 
 * and vertex[1] to vertex[2], etc
 **********************************************************************/
void
element::compCentroidAndArea (
			      void)
{
  double triangleArea;
  point3D triangleCentroid;
  point_vec_operations::point3D<double> centroid__;
  area_ = 0;
  centroid_ = vector3D<double>(0., 0., 0.);
  if (boundaryIndex_==-1) {
    centroid__=fil_->filCentroid();
    area_=fil_->filCrossA();
    centroid_.x()=centroid__.x();
    centroid_.y()=centroid__.y();
    centroid_.z()=centroid__.z();
  }
  else if (boundaryIndex_==-2) {
    centroid__=panel_->centroid();
    area_=panel_->area();
    centroid_.x()=centroid__.x();
    centroid_.y()=centroid__.y();
    centroid_.z()=centroid__.z();
  }
  else {
    pfft::errorMessage("element.cc:compCentroidAndArea",
		       "ele unknown shape! There must be a bug here!!");
  }    
}


/**********************************************************************
 * compBoundingSphereRadius --
 **********************************************************************/
void
element::compBoundingSphereRadius ( void)
{
  boundingSphereRadius_ = 0.;  zmin_=vertex_[0].z();
  for (size_t i=0; i<shape(); i++) {
    vector3D<double> vec = boundingSphereCenter() - vertex_[i];
    boundingSphereRadius_ = max(boundingSphereRadius_, length(vec));
    zmin_=min(zmin_,vertex_[i].z());
  }
}


/**********************************************************************
 * comparison --
 **********************************************************************/
bool
pfft::operator == (
		   const element& e1, 
		   const element& e2) 
{
  return (e1.shape() == e2.shape()) && (e1.vertex_ == e2.vertex_);
}

bool
pfft::operator != (
		   const element& e1, 
		   const element& e2) 
{
  return (e1.shape() != e2.shape()) || (e1.vertex_ != e2.vertex_);
}

/**********************************************************************
 * output --
 **********************************************************************/
ostream& 
pfft::operator << (
		   ostream& os, 
		   const element& e) 
{
  os << endl 
     << "name = " << e.name() << endl
     << "boundaryIndex = " << e.boundaryIndex() << endl
     << "shape = " << e.shape() << endl;

  for (size_t i=0; i<e.shape(); i++) 
    os << "vertex[" << i << "] = " << e.vertex(i); 

  os << "area = " << e.area() << endl
     << "centroid = " << e.centroid()
     << "boundingSphereRadius = " << e.boundingSphereRadius()
     << endl;
  return os;
}

/**********************************************************************
 * input --
 **********************************************************************/
istream&
pfft::operator >> (
		   istream& is,
		   element& e)
{
  string shapeKey;
  int shape;

  if (is) {
    is >> shapeKey >> e.boundaryIndex_;
    if ( (shapeKey=="3") || (shapeKey=="t") || (shapeKey=="T") ) 
      shape = 3;
    else if ( (shapeKey=="4") || (shapeKey=="q") || (shapeKey=="Q") ) 
      shape = 4;
    else
      throw domain_error("corrupted element shape data in element line");

    e.vertex_.clear();
    e.vertex_.reserve(shape);
    for (int i=0; i<shape; i++) {
      vector3D<double> vt;
      if ( !(is >> vt) ) 
	throw domain_error("corrupted vertex data in element line");
      e.vertex_.push_back(vt);
    }
    e.setupElement();
  }
  return is;
}


/**********************************************************************
 * findMinEdgeLength --
 **********************************************************************/
double
element::findMinEdgeLength (
			    void) const
{
  double minEdgeLength = 1e20;
  for (size_t i=0; i<vertex_.size(); i++) {
    size_t j = (i+1) % vertex_.size();
    double edgeLength = length(vertex_[i] - vertex_[j]);
    minEdgeLength = minEdgeLength < edgeLength ? minEdgeLength : edgeLength;
  }
  return minEdgeLength;
}

