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
#ifndef INTERNAL_NODE_H
#define INTERNAL_NODE_H

#include "const_def.h"
#include "panel.h"
#include "filament.h"
#include "list.h"

class Resistance;

struct Resistance_Entry
{
  Resistance_Entry (Resistance *_resistance_ptr, int _direction)
  {
    resistance_ptr = _resistance_ptr;
    direction = _direction;
  }
  Resistance *get_resist_ptr () const { return resistance_ptr; };
  int get_direction () const { return direction; };

  Resistance *resistance_ptr;
  int direction;
};


struct Filament_Entry
{
  Filament_Entry (filament<double> *_filament_ptr, int _direction){
	  filament_ptr = _filament_ptr;
	  direction = _direction;
  }
  filament<double> *get_fil_ptr () const { return filament_ptr; };
  int get_direction () const { return direction; };

  filament<double> *filament_ptr;
  int direction;
};



class Internal_Node
{
  friend class Structure;
  friend class conductor;

public:
  Internal_Node (int _type=IN_NORMAL);

  panel<double> * get_first_panel_ptr ()
  {return list_internal_node_panels.initialize(); };
  void set_number (const long int _number) { number = _number; };
  long int get_number () const { return number; };
  void add_list_filaments_with_direction (
    List<filament<double> > *semi_list_filaments, int _direction);
  void add_resistance_with_direction (Resistance *resist, int _direction);
  bool is_resist_int_node(void);


private:
  long int number;
  List <panel<double> > list_internal_node_panels;
  List <Resistance_Entry> list_internal_node_resistances_entries;
  List <Filament_Entry> list_internal_node_filaments_entries;
  int type;
};
#endif
