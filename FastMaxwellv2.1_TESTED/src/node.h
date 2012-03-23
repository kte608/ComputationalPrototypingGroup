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
#ifndef NODE_H
#define NODE_H

#include "list.h"
#include"conductor.h"
#include "internal_node.h"


class Node
{
public:
  Node(void) {};
  Node (char *, double, double, double, int);

  double x (void) const { return x_; };
  double y (void) const { return y_; };
  double z (void) const { return z_; };
  char* name_node(void) const {return name_;};

  void add_conductor(conductor* newCond )
  {
    list_node_conductors.push_back(newCond);
  }
  Node *real_node ();
  void assign_real_node (Node *node) { real_node_ptr = node; };
  List <conductor> *get_list_node_conductors_ptr ()
  { return &list_node_conductors; };
  List<Internal_Node> *get_list_internal_nodes ()
  { return &list_node_internal_nodes; };
  void set(double x, double y, double z) {x_=x; y_=y; z_=z;};
  void add_int_node_ptr (Internal_Node *int_node_ptr) ;
  int is_real_node () { return (real_node_ptr == this); };



private:
  double x_, y_, z_;
  char *name_;
  int type_;
  Node *real_node_ptr;
  List <conductor> list_node_conductors;
  List <Internal_Node> list_node_internal_nodes;

};

#endif
