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
#ifndef EXT_NODE_PAIR_H
#define EXT_NODE_PAIR_H

#include "node.h"
#include "resistance.h"
#include <vector>

/**
@author Xin Hu
*/
class Ext_Node_Pair
{
public:
  Ext_Node_Pair (Node *first_node, Node *second_node, char *name=0);
  char *node_name (const int node_number) const;
  char *name () const { return port_name; };
  Node *get_node (const int node_number) const;
  std::vector<Resistance*> get_resistance_list(void)
  {return resistance_list;};
  void connect_to_resistance(Resistance* resist);
private:
  Node *node1, *node2;
  char *port_name;
  std::vector<Resistance*> resistance_list;
  
};

#endif


