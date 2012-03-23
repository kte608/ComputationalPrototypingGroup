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
#include "internal_node.h"


Internal_Node::Internal_Node (int _type)
{
  type = _type;
  number = 0;
}

void Internal_Node::add_list_filaments_with_direction (
  List<filament<double> > *semi_list_filaments, int _direction)
{
  filament<double> *fil;

  for (fil = semi_list_filaments->initialize(); fil!=NULL;
       fil = semi_list_filaments->iterate() )
  {
    list_internal_node_filaments_entries.push_back(
      new Filament_Entry(fil, _direction) );
  }
}

void Internal_Node::add_resistance_with_direction(Resistance *resist,
    int _direction)
{
  list_internal_node_resistances_entries.push_back(
    new Resistance_Entry(resist, _direction) );
}

bool Internal_Node::is_resist_int_node(){
	if (type == IN_RESIST) return 1;
	return 0;
}

