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
#ifndef STRUCTURE_H
#define STRUCTURE_H

#include "list.h"
#include "node.h"
#include "pseudo_node.h"
#include "ext_node_pair.h"
#include "conductor.h"
//#include "henry.h"
//#include "incidence.h"
#include "fastsub.h"
//#include <vector>

#include <iostream>


using namespace fastsub;

struct mesh_matrices
{
 // Incidence_Matrix<filament<double> > *Mf_ptr, *Mfs_ptr;
 // Incidence_Matrix<panel<double> > *Mp_ptr, *Mps_ptr;
 // Incidence_Matrix<Resistance> *Mrs_ptr;
 // Source_Vector *Vs_ptr;
  
pfft::SpRowMat<double> Mf_ptr, Mfs_ptr, Mp_ptr, Mps_ptr, Mrs_ptr; 
pfft::SpVec<double> Vs_ptr;
};

class Structure
{
public:

  friend class sysFormation;
  friend class Charge;
  
  Structure(void);
  void add_defined_node (Node* newNode)
  { list_defined_nodes.push_back(newNode); }
  void add_conductor(conductor* newCond)
  { list_conductors.push_back(newCond); }
  void add_pseudo_node (Pseudo_Node *new_pseudo_node)
  { list_pseudo_nodes.push_back (new_pseudo_node); };
  void add_ext_node_pair (Ext_Node_Pair *new_ext)
  { list_ext_node_pairs.push_back (new_ext); };
  void add_resistance (Resistance *new_resist)
  { list_resistances.push_back(new_resist); };
  Node *get_node_from_name (char *name);
  Ext_Node_Pair *get_external_from_portname (char *name);
  Ext_Node_Pair *get_ext_node_pair_from_node_names(char *node1_name,
      char *node2_name);
  Resistance *get_resistance_from_resistance_name(char *name);
  Resistance *get_resistance_from_node_names(char *node1_name,
      char *node2_name);
  int num_panels(void)const {return list_panels.size();};
  int num_filaments(void) const {return list_filaments.size();};
  void discretize(SolverMethod);
  void find_load_resistance();

private:
  mesh_matrices mesh_mats;
  //vector<panel<double> > list_panels_vec;
  List<panel<double> > list_panels;
  List<filament<double> > list_filaments;
  //vector<filament<double> > list_filaments_vec;
  List<conductor> list_conductors;
  List <Internal_Node> global_list_internal_nodes;
  List<Node> list_defined_nodes;
  List<Pseudo_Node> list_pseudo_nodes;
  List<Ext_Node_Pair> list_ext_node_pairs;
  List <Resistance> list_resistances;
  //List<Resistance> list_load_resistances;
//  List <SEGMENT> list_henry_segments;

  //private funcs
  void prepare_near_conductors(void);
  void connect_conductors_meshes (void);
  void connect_resistances_meshes(void);
  void connect_sources_meshes (void);

};


#endif
