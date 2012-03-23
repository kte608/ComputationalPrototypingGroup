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
#include "frequencysample.h"
#include "fastsub.h"
#include <cmath>
#include <stdexcept>
#include <iostream>

using namespace fastsub;
using namespace std;

/**********************************************************************
 * FrequencySample --
 **********************************************************************/
FrequencySample::FrequencySample (
  const double firstPoint,
  const double lastPoint,
  const int numSam,
  const SampleMethod samMethod)
    : firstPoint_(firstPoint), lastPoint_(lastPoint),
    numSam_(numSam), samMethod_(samMethod)
{
  firstPoint_ = max(100.,firstPoint);
  lastPoint_ = max(firstPoint_,lastPoint_);

  if (samMethod == LOG)
  {
    setupLogSamPoint();
  }
  else if (samMethod == LINEAR)
  {
    setupLinearSamPoint();
  }
  else
  {
    //errorMessage("frequencySample.cc : frequencySample",
     //            "Illegal frequency sampling method");
  }

  try
  {
    checkRange();
  }
  catch (const domain_error error)
  {
    cerr << error.what() << endl;
    //errorMessage("frequencySample.cc : frequencySample",
      //           "Illegal frequency sampling data");
  }

  if (numSam_ >= 3)
  {
    size_t midPointIndex = numSam_/2;
    middlePoint_ = samPointList_[midPointIndex];
  }
  else
  {
    middlePoint_ = firstPoint_;
  }
}

/**********************************************************************
 * setupLogSamPoint --
 **********************************************************************/
void
FrequencySample::setupLogSamPoint (
  void)
{
  double step;
  if (numSam_ >= 2)
  {
    step = (log10(lastPoint_) - log10(firstPoint_))/(numSam_-1);
    step = pow(10, step);
  }

  double f = firstPoint_;
  samPointList_.reserve(numSam_);
  samPointList_.push_back(f);
  for (int i = 1; i < numSam_; i++)
  {
    f *= step;
    samPointList_.push_back(f);
  }
}

/**********************************************************************
 * setupLinearSamPoint --
 **********************************************************************/
void
FrequencySample::setupLinearSamPoint (
  void)
{
  samPointList_.resize(numSam_);
  if (numSam_ >= 2)
  {
    double step = (lastPoint_ - firstPoint_) / (numSam_-1);
    for (size_t i = 0; i < numSam_; i++)
    {
      samPointList_[i] = firstPoint_ + i * step;
    }
  }
  else
  {
    samPointList_[0] = firstPoint_;
  }
}

/**********************************************************************
 * checkRange --
 **********************************************************************/
void
FrequencySample::checkRange (
  void)
{
  if (firstPoint_ < 0)
  {
    throw domain_error("Begin of the frequency range should be positive!");
  }

  if (lastPoint_ < 0)
  {
    throw domain_error("End of the frequency range should be positive!");
  }

  if (lastPoint_ < firstPoint_)
  {
    throw domain_error("Begin of the frequency range should be less than end of the frequency range!");
  }

  if ( ((firstPoint_ != lastPoint_) && (numSam_ == 1) )   ||
       ((firstPoint_ == lastPoint_) && (numSam_ != 1) ) )
  {
    throw domain_error("Inconsistant frequency range and sampling number!");
  }

  if (numSam_ <= 0)
  {
    throw domain_error("Number of frequency sampling points should be positive!");
  }

  if ( (numSam_ > 1) && (lastPoint_ < firstPoint_) )
  {
    throw domain_error("End of the frequency range should be larger than the begin the range!");
  }
}




