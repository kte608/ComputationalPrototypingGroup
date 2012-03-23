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
#ifndef FREQUENCYSAMPLE_H
#define FREQUENCYSAMPLE_H

#include <vector>

namespace fastsub
{

  enum SampleMethod { LOG, LINEAR, NUM_SAMPLE_METHOD };

  class FrequencySample
  {

  public:
    FrequencySample(const double firstPoint = 1e6,
                    const double lastPoint = 1e6,
                    const int numSam = 1,
                    const SampleMethod samMethod = LOG);
    double point(int i) const { return samPointList_[i]; }
    int numSam(void) const { return numSam_; }
    double firstPoint(void) const { return firstPoint_; }
    double lastPoint(void) const { return lastPoint_; }
    SampleMethod samMethod(void) const { return samMethod_; }
    double middlePoint(void) const { return middlePoint_; }

  private:
    double firstPoint_;
    double lastPoint_;
    double middlePoint_;
    int numSam_;
    SampleMethod samMethod_;
    std::vector<double> samPointList_;

    void setupLogSamPoint (void);
    void setupLinearSamPoint (void);
    void checkRange(void);

  };

} // namespace surf

#endif
