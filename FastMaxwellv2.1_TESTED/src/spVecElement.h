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
#ifndef _SP_VEC_ELEMENT_H_
#define _SP_VEC_ELEMENT_H_




template <class T>
class SpVecElement
{
public:
  SpVecElement(void) : index_(0), value_(0) {};
  SpVecElement(size_t index, T value) : index_(index), value_(value) {};

  friend bool operator< <> (const SpVecElement<T>&, const SpVecElement<T>&);
  SpVecElement<T>& operator+= (const T s) { value_ += s; return *this; }
  SpVecElement<T>& operator-= (const T s) { value_ -= s; return *this; }
  SpVecElement<T>& operator*= (const T s) { value_ *= s; return *this; }
  SpVecElement<T>& operator/= (const T s) { value_ /= s; return *this; }

  const T value(void) const { return value_; }
  T& value(void) { return value_; }
  size_t index(void) const { return index_; }
  void set_value(T val){value_=val;};

private:
  size_t index_;
  T value_;

};

/**********************************************************************
  * output --
**********************************************************************/
template <class T>
std::ostream&
operator << (
  std::ostream& os,
  const SpVecElement<T>& ele)
{
  os << "(" << ele.index() << "," << ele.value() << ")";
  return os;
}

/**********************************************************************
		  * operator< --
**********************************************************************/
template<class T>
bool operator< (
  const SpVecElement<T>& ele1,
  const SpVecElement<T>& ele2)
{
  return (ele1.index_ < ele2.index_);
}




#endif

