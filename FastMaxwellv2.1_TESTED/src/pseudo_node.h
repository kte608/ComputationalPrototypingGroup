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
#ifndef PSEUDO_NODE_H
#define PSEUDO_NODE_H

#include "node.h"
#include <string.h>

/**
@author Xin Hu
*/
class Pseudo_Node
{
public:
  Pseudo_Node(char *node_name)
  {
    name = new char[strlen(node_name)+1];
    strcpy (name, node_name);
    defined_node_ptr_ = NULL;
  };
  void assign_defined_node (Node *_defined_node_ptr)
  { defined_node_ptr_ = _defined_node_ptr; };
  Node *get_defined_node () { return defined_node_ptr_; };
  char *name_node () { return name; };
private:
  char *name;
  Node *defined_node_ptr_;
};

#endif
