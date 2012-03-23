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
#ifndef LAYERENV_H
#define LAYERENV_H

#include <complex>
#include <vector>
#include "layer.h"


using namespace std;
class LayerEnv
{
  friend class Dipoles;

public:
  typedef complex<double> Complex;

  LayerEnv();
  void add_defined_layer(const Layer* newLayer);
  Layer* get_layer_from_name(char* name);
  Layer* get_layer_from_startHeight(double height);
  Layer* get_layer_from_index(int layerIndex);
  void compute_freq_related_info(double freq);
  int num_layers(void) const {return list_layers.size();};
  complex<double> get_k_topmostLayer(void) const;
  double get_epir_topmostlayer(void)const { return list_layers[0].epir();};
  double get_subHeight(void) const; 

private:
  vector<Layer> list_layers;
  vector<Layer>::iterator iter;
  bool layer_freq_flag;

};
#endif
