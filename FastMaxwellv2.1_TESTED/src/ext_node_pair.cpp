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
#include "ext_node_pair.h"
#include <string.h>

Ext_Node_Pair::Ext_Node_Pair (Node *first_node, Node *second_node, char *name)
{
  node1 = first_node;
  node2 = second_node;
  if (name)
  {
    port_name = new char[strlen(name)+1];
    strcpy (port_name, name);
  }
  else
    port_name = NULL;
}

char *Ext_Node_Pair::node_name (const int node_number) const
{
  switch (node_number)
  {
  case 1:
    if (node1 != NULL)
      return node1->name_node();
    else
      return '\0';
    break;
  case 2:
    if (node1 != NULL)
      return node2->name_node();
    else
      return '\0';
    break;
  }
}

void Ext_Node_Pair::connect_to_resistance(Resistance* resist)
{
	resistance_list.push_back(resist);
}

Node *Ext_Node_Pair::get_node (const int node_number) const
{
  switch (node_number)
  {
  case 1:
    return node1;
    break;
  case 2:
    return node2;
    break;
  }
}

