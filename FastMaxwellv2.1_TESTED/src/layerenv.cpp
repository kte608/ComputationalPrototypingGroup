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
#include "layerenv.h"
#include "const_def.h"
#include <math.h>


LayerEnv::LayerEnv(void)
{
  layer_freq_flag=false;
}

void LayerEnv::add_defined_layer(const Layer* newLayer)
{
  bool flag=false;

  if (strcmp(newLayer->name(), "top")==0)
  {
    iter=list_layers.begin();
    list_layers.insert(iter, *newLayer);
    flag=true;
  }

  if (flag==false)
  {
    for (iter=list_layers.begin(); iter!=list_layers.end(); iter++)
    {
      if ((*iter).startHeight() < newLayer->startHeight())
      {
        flag=true;
        break;
      }
    }
  }
  if (flag==false)
    list_layers.push_back(*newLayer);

  /*cout<<list_layers.size()<<" ";
  for (iter=list_layers.begin(); iter!=list_layers.end(); iter++)
  {
  	cout<<(*iter).name()<<" ";
  }*/
}


Layer* LayerEnv::get_layer_from_name(char *name)
{
  for (iter=list_layers.begin(); iter!=list_layers.end(); iter++)
    if (strcmp((*iter).name(), name)==0)
      return &(*iter);

  return NULL;

}

Layer* LayerEnv::get_layer_from_startHeight(double height)
{
  for (iter=list_layers.begin(); iter!=list_layers.end(); iter++)
  {
    if (fabs((*iter).startHeight()-height)<EPS)
      return  &(*iter);
  }
  return NULL;
}

Layer* LayerEnv::get_layer_from_index(int index)
{
	if (index>list_layers.size()-1)
		return NULL;
	else
		return &(list_layers[index]);
}

void LayerEnv::compute_freq_related_info(double freq)
{
  Complex k;
  double omega=2*PI*freq;

  for (iter=list_layers.begin(); iter!=list_layers.end(); iter++)
  {
    Complex epi (((*iter).epir())*EPI0,
                 -((*iter).sigma())/(omega));
    k=sqrt(omega*omega*MU0*epi);

    (*iter).set_epi(epi);
    (*iter).set_k(k);
  }
  layer_freq_flag=true;

}

complex<double> LayerEnv::get_k_topmostLayer(void) const {
	if (layer_freq_flag==true)
  		return list_layers[0].k();
	else{
		cout<<" Wave number for each layer hasn't been set yet"<<endl;
		return ZERO;
	}
}

double LayerEnv::get_subHeight(void) const {
	if (list_layers.size()>1)
		return list_layers[1].startHeight();
	else{
		cout<<"no substrate exists"<<endl;
		exit(1);
	}

}

