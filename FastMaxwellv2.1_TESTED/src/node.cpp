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
#include "node.h"
#include <string>


Node::Node (char *node_name, double x_value, double y_value,
            double z_value, int type)
{
  x_=x_value;
  y_=y_value;
  z_=z_value;
  real_node_ptr = this;
  name_=new char [strlen(node_name)+1];
  strcpy(name_, node_name);
  type_=type;
}

Node* Node::real_node()
{
  if (real_node_ptr == this) return this;
  else return real_node_ptr->real_node();
}

void Node::add_int_node_ptr (Internal_Node *int_node_ptr) {

	if (is_real_node() == 0) {
		printf ("Errare humanum est: Node::add_int_node_ptr() - this is "
				"not a real_node\n");
		exit (-1);
	}
	list_node_internal_nodes.push_back (int_node_ptr); 
}
