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
#ifndef POINT3D_H
#define POINT3D_H
#include <math.h>
#include <iostream>
//#include "node.h"
//class Node;
using namespace std;

namespace point_vec_operations
{

  template <class T>
  class point3D
  {
  public:
    point3D(void) {};
    point3D(const T x, const T y, const T z): x_(x), y_(y), z_(z) {};
  //  point3D(const Node &n): x_(n.x()),y_(n.y()), z_(n.z()) {};

	typedef point_vec_operations::point3D<double> vector3D;


    const T x(void)const {return x_;};
    const T y(void)const {return y_;};
    const T z(void)const  {return z_;};
    const T abs(void) const;

    void set_x(const T x) { x_=x;};
    void set_y(const T y) { y_=y;};

    point3D<T>&  operator = (const point3D<T> & );
    point3D<T>&  operator += (const point3D<T> & );
    point3D<T>&  operator -= (const point3D<T> & );
    point3D<T>&  operator *= (const T &);
    point3D<T>&  operator /= (const T &);
    //int operator == (const point3D<T> & );
    //int paral ( const point3D<T>  );


    void transferLocalToGlobalCoord (
      const point3D<T>& origin,
      const point3D<T>& X,
      const point3D<T>& Y,
      const point3D<T>& Z);

    void transferGlobalToLocalCoord (
      const point3D<T>& origin,
      const point3D<T>& X,
      const point3D<T>& Y,
      const point3D<T>& Z);


    void set(T const& , T const& , T const& );
    void normalize(void);
    void print_info(void);

  private:
    T x_,y_,z_;
  };



  template <class T>
  void point3D<T>::set (T const& x, T const& y , T const& z)
  {
    x_=x;
    y_=y;
    z_=z;
  }

  template <class T>
  point3D<T>&  point3D<T>::operator = (point3D const& p2)
  {
    if (this==&p2) return *this;
    if (this != &p2)
    {
      x_=p2.x_;
      y_=p2.y_;
      z_=p2.z_;
      return *this;
    }
  }

  template <class T>
  point3D<T>&  point3D<T>::operator += (point3D const& p2 )
  {
    x_ +=p2.x_;
    y_ +=p2.y_;
    z_ += p2.z_;
    return *this;
  }

  template <class T>
  point3D<T>&  point3D<T>::operator -= (point3D const& p2 )
  {
    x_ -=p2.x_;
    y_ -=p2.y_;
    z_ -= p2.z_;
    return *this;
  }

  template <class T>
  point3D<T>&  point3D<T>::operator *= (T const& d)
  {
    x_  *= d;
    y_  *=d;
    z_  *=d;
    return *this;
  }

  template <class T>
  point3D<T>&  point3D<T>::operator /= (T const& d)
  {
    x_ /= d ;
    y_ /= d;
    z_ /= d;
    return *this;
  }


  template <class T>
  void point3D<T>::normalize(void)
  {
    T len=length(*this);
    if (len !=0.)
    {
      x_/=len;
      y_/=len;
      z_/=len;
    }
  }

  template <class T>
  const	T point3D<T>::abs() const
  {
    return sqrt(x_*x_ + y_*y_ + z_*z_);
  }

  template <class T1>
  point3D<T1> transferLocalToGlobalCoord (
    const point3D<T1>& vec,
    const point3D<T1>& origin,
    const point3D<T1>& X,
    const point3D<T1>& Y,
    const point3D<T1>& Z)
  {
    point3D<T1> ans(vec);
    ans.transferLocalToGlobalCoord(origin, X, Y, Z);
    return ans;
  }

  template <class T1>
  point3D<T1> transferGlobalToLocalCoord (
    const point3D<T1>& vec,
    const point3D<T1>& origin,
    const point3D<T1>& X,
    const point3D<T1>& Y,
    const point3D<T1>& Z)
  {
    point3D<T1> ans(vec);
    ans.transferGlobalToLocalCoord(origin, X, Y, Z);
    return ans;
  }

  template <class T>
  void point3D<T>::print_info(void)
  {
    std::cout<<x_<<" "<<y_<<" "<<z_<<std::endl;
  }

  //*********************************************************************
  // finish defining template class point3D
  // now we are defining functions in the namespace
  //*********************************************************************

  /**********************************************************************
  * transferGlobalToLocalCoord --
  * Vectors *this, origin, X, Y and Z are in global coordinate system.
  * This function find the local coordinates of vector *this.
  * This local coordinate system consists of origin, X, Y and Z. 
  **********************************************************************/
  template <class T1>
  void point3D<T1>::transferGlobalToLocalCoord (
    const point3D<T1>& origin,
    const point3D<T1>& X,
    const point3D<T1>& Y,
    const point3D<T1>& Z)
  {
    T1 x = (*this - origin) * X;
    T1 y = (*this  - origin) * Y;
    T1 z = (*this  - origin) * Z;

    x_ = x;
    y_ = y;
    z_ = z;
  }

  template <class T1>
  void point3D<T1>::transferLocalToGlobalCoord (
    const point3D<T1>& origin,
    const point3D<T1>& X,
    const point3D<T1>& Y,
    const point3D<T1>& Z)
  {
    point3D<T1> rotation(X.x(), Y.x(), Z.x());
    T1 x = (*this) * rotation + origin.x();

    rotation = point3D<T1>(X.y(), Y.y(), Z.y());
    T1 y = (*this) * rotation + origin.y();

    rotation = point3D<T1>(X.z(), Y.z(), Z.z());
    T1 z = (*this) * rotation + origin.z();

    x_ = x;
    y_ = y;
    z_ = z;
  }

  template <class T>
  point3D<T> operator - (const point3D<T>& v1, const point3D<T>& v2)
  {
    point3D<T> vt=v1;
    vt-=v2;
    return vt;
  }

  template <class T>
  point3D<T> operator + (const point3D<T>& v1, const point3D<T>& v2)
  {
    point3D<T> vt=v1;
    vt+=v2;
    return vt;
  }

  template<class T1, class T2>
  point3D<T1> operator / (const point3D<T1>& vec, const T2 & s)
  {
    point3D<T1> vt=vec;
    vt/=s;
    return vt;
  }

  template<class T1, class T2>
  point3D<T1> operator * (const point3D<T1>& vec, const T2 & s)
  {
    point3D<T1> vt=vec;
    vt*=s;
    return vt;
  }
  template<class T1, class T2>
  T1 operator * (const point3D<T1>& vec, const point3D<T2>& vec2)
  {
    return vec.x()*vec2.x()+vec.y()*vec2.y()+vec.z()*vec2.z();
  }

  template<class T>
  int operator == (const point3D<T>& p1, const point3D<T>& p2 )
  {
    return ((p1.x() == p2.x()) && (p1.y()==p2.y()) && (p1.z()== p2.z()));
  }




  template <class T>
  point3D<T> crossProd( const point3D<T>& v1,
                        const point3D<T>& v2)
  {
    T x = v1.y()*v2.z() - v1.z()*v2.y();
    T y = v1.z()*v2.x() - v1.x()*v2.z();
    T z = v1.x()*v2.y() - v1.y()*v2.x();
    return point3D<T>(x, y, z);

  }


  template <class T>
  T length(const point3D<T>& vec)
  {
    T len=vec.x()*vec.x()+vec.y()*vec.y()+vec.z()*vec.z();
    return sqrt(len);
  }


  template <class T>
  T length(const point3D<T>& vec1, const point3D<T>& vec2 )
  {
    point3D<T> vec=vec1-vec2;
    T len=vec.x()*vec.x()+vec.y()*vec.y()+vec.z()*vec.z();
    return sqrt(len);
  }

  template<class T1, class T2>
  point3D<T2>  cross ( const point3D<T1>& v1, const point3D<T2>& v2)
  {
    T2 x = v1.y()*v2.z() - v1.z()*v2.y();
    T2 y = v1.z()*v2.x() - v1.x()*v2.z();
    T2 z = v1.x()*v2.y() - v1.y()*v2.x();
    return point3D<T2>(x, y, z);
  }

  template<class T>
  T  angle ( const point3D<T>& p1,  const point3D<T>& p2)
  {
    point3D<T> p_null (0, 0, 0);

    if ((p1==p_null) || (p2==p_null))
      fprintf (stderr, "Warning: Point::get_angle_with() - operating on "
               "null vector\n");
    return acos(p1.x()*p2.x() + p1.y()*p2.y() + p1.z()*p2.z());
  }



  template <class T>
  int paral (  const point3D<T>&  p1,   const point3D<T>&  p2)
  {
    point3D<T> p_null (0, 0, 0);

    if ((p1==p_null) || (p2==p_null))
      fprintf (stderr, "Warning: Point::paral() - operating on null vector\n");
    point3D<T> temp=cross(p1, p2);
    if (temp.abs() == 0) return 1;
    else return 0;
  }


}
 //match namespace point_vec_operations

#endif

