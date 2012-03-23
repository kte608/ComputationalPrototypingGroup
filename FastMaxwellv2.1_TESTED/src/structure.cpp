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
#include "structure.h"
#include "const_def.h"
#include "panel.h"
#include "filament.h"
#include "resistance.h"

#include <math.h>
#include <string.h>
#include <stdlib.h>



Structure::Structure ()
{
  // mesh_mats.Mf_ptr = NULL;
  // mesh_mats.Mp_ptr = NULL;
  // mesh_mats.Mfs_ptr = NULL;
  //  mesh_mats.Mps_ptr = NULL;
  // mesh_mats.Mrs_ptr = NULL;
  // mesh_mats.Vs_ptr = NULL;
    mesh_mats.Mf_ptr = pfft::SpRowMat<double>();
    mesh_mats.Mp_ptr = pfft::SpRowMat<double>();
    mesh_mats.Mfs_ptr = pfft::SpRowMat<double>();
    mesh_mats.Mps_ptr = pfft::SpRowMat<double>();
    mesh_mats.Mrs_ptr = pfft::SpRowMat<double>();
    mesh_mats.Vs_ptr = pfft::SpVec<double>();
}



void Structure::discretize (SolverMethod solve_method)
{
  conductor *conductor_ptr;
  List<panel<double> > *list_panels_conductor_ptr;
  List<filament<double> > *list_filaments_conductor_ptr;
  //List <SEGMENT> *list_henry_segments_conductor_ptr;
  List <Internal_Node> *list_internal_nodes_conductor_ptr;
  List <Internal_Node> *list_internal_nodes_resistance_ptr;
  // Incidence_Matrix<filament<double> > *Mf_conductor_ptr, *Mfs_conductor_ptr;
  // Incidence_Matrix<panel<double> > *Mp_conductor_ptr, *Mps_conductor_ptr;
  pfft::SpRowMat<double> Mf_conductor_ptr, Mfs_conductor_ptr, Mp_conductor_ptr, Mps_conductor_ptr; 
  
  panel<double> *p;
  filament<double> *f;
  Resistance *resist, *r;
  long int i;
  int e;
  
  std::cout<<"--> Stage: Discretization"<<endl;
  
  // create mesh or nodal matrices
  if (solve_method == CIRCUIT_SOLVE_MESH)
    {
      mesh_mats.Mf_ptr = pfft::SpRowMat<double>();
      mesh_mats.Mp_ptr = pfft::SpRowMat<double>();
      mesh_mats.Mfs_ptr = pfft::SpRowMat<double>();
      mesh_mats.Mps_ptr = pfft::SpRowMat<double>();
      mesh_mats.Mrs_ptr = pfft::SpRowMat<double>();
      mesh_mats.Vs_ptr = pfft::SpVec<double>();
    }
  
  //prepare conductors connected to a common node
  prepare_near_conductors ();
  
  //prepare each conductor and discretization
  conductor_ptr = list_conductors.initialize ();
  
  int global_mesh_line_p = 1;
  int global_mesh_line_f = 1;
  int global_mesh_line_s = 1;
  while (conductor_ptr != NULL)
    {
      conductor_ptr->discretize(solve_method,
				&list_panels_conductor_ptr,
				&list_filaments_conductor_ptr,
				//&list_henry_segments_conductor_ptr,
				&list_internal_nodes_conductor_ptr,
				&(mesh_mats.Mf_ptr), &(mesh_mats.Mp_ptr),
				&(mesh_mats.Mfs_ptr), &(mesh_mats.Mps_ptr),
				list_panels.size(),list_filaments.size(),
				&global_mesh_line_p, &global_mesh_line_f, 
				&global_mesh_line_s);
      
      list_panels.add_sub_list(list_panels_conductor_ptr);
      list_filaments.add_sub_list (list_filaments_conductor_ptr);
      //  list_henry_segments.add_sub_list (list_henry_segments_conductor_ptr);
      global_list_internal_nodes.add_sub_list (list_internal_nodes_conductor_ptr);
      
      
      
      delete list_panels_conductor_ptr;
      delete list_filaments_conductor_ptr;
      //delete list_henry_segments_conductor_ptr;
      
      
      if (solve_method == CIRCUIT_SOLVE_MESH)
	{
	  // process matrices to become lists, if they are not already
	  //   Mf_conductor_ptr->set_mode(MODE_LIST);
	  //   Mp_conductor_ptr->set_mode(MODE_LIST);
	  //   Mfs_conductor_ptr->set_mode(MODE_LIST);
   //   Mps_conductor_ptr->set_mode(MODE_LIST);

      // add mesh matrices to global mesh matrices
   //   mesh_mats.Mf_ptr->add_sub_matrix (Mf_conductor_ptr);
   //   mesh_mats.Mp_ptr->add_sub_matrix (Mp_conductor_ptr);
   //   mesh_mats.Mfs_ptr->add_sub_matrix (Mfs_conductor_ptr);
   //   mesh_mats.Mps_ptr->add_sub_matrix (Mps_conductor_ptr);
/*
      mesh_mats.Mf_ptr = (Mf_conductor_ptr);
      mesh_mats.Mp_ptr = (Mp_conductor_ptr);
      mesh_mats.Mfs_ptr = (Mfs_conductor_ptr);
      mesh_mats.Mps_ptr = (Mps_conductor_ptr);
      
      // delete conductor matrices
      delete Mf_conductor_ptr;
      delete Mp_conductor_ptr;
      delete Mfs_conductor_ptr;
      delete Mps_conductor_ptr;
      */
    }

    // each conductor gives a internal node pointer to real nodes
    // correspondent to each of his 'lateral' nodes. this is used later to
    // connect conductors, and to nodal analysis
    conductor_ptr->send_connect_info();

    // prepare next conductor
    conductor_ptr = list_conductors.iterate ();
  }

  // test to the number of meshes (before connecting conductors or sources)
  if (solve_method == CIRCUIT_SOLVE_MESH)
  {
    long int n_meshes_test, n_branches_test, n_junctions_test;

 //   n_meshes_test = mesh_mats.Mf_ptr->get_n_lines() + mesh_mats.Mp_ptr->get_n_lines() +
   //                 mesh_mats.Mfs_ptr->get_n_lines();
    n_meshes_test = mesh_mats.Mf_ptr.numRow() + mesh_mats.Mp_ptr.numRow() +
                    mesh_mats.Mfs_ptr.numRow();
    n_branches_test = list_filaments.size() + list_panels.size();
    n_junctions_test = global_list_internal_nodes.size() +1;
    if (n_meshes_test != (n_branches_test - n_junctions_test + 1) )
    {
      printf ("Errare humanum est: Structure::discretize() - "
              "number of meshes\n");
      exit (-1);
    }
  }

  // each resistance creates information for connections
  for (resist = list_resistances.initialize (); resist != NULL;
       resist = list_resistances.iterate () )
  {
    resist->create_connect_info(&list_internal_nodes_resistance_ptr);
    global_list_internal_nodes.add_sub_list (list_internal_nodes_conductor_ptr);
    // like in conductors, do not delete each resistance list of internal_nodes
  }

  // do remaining meshes, if input options specify MESH analysis
  if (solve_method == CIRCUIT_SOLVE_MESH)
  {
    // each resistance creates information for connections
    for (resist = list_resistances.initialize (); resist != NULL;
         resist = list_resistances.iterate () )
    {
      resist->create_connect_info(&list_internal_nodes_resistance_ptr);
      global_list_internal_nodes.add_sub_list (list_internal_nodes_conductor_ptr);
      // like in conductors, do not delete each resistance list of internal_nodes
    }

    // do meshes to connect conductors
    connect_conductors_meshes ();

    // do meshes to connect resistances
    connect_resistances_meshes ();

    // do meshes to connect sources
    connect_sources_meshes ();

    // process global matrices to become vectors
//    mesh_mats.Mf_ptr->set_mode (MODE_VECTOR);
//    mesh_mats.Mp_ptr->set_mode (MODE_VECTOR);
//    mesh_mats.Mfs_ptr->set_mode (MODE_VECTOR);
//    mesh_mats.Mps_ptr->set_mode (MODE_VECTOR);
//    mesh_mats.Mrs_ptr->set_mode (MODE_VECTOR);

//    printf ("Total number of Meshes: %ld\n", mesh_mats.Mf_ptr->get_n_lines() +
  //          mesh_mats.Mp_ptr->get_n_lines() + mesh_mats.Mfs_ptr->get_n_lines() );
    printf ("Total number of Meshes: %ld\n", mesh_mats.Mf_ptr.numRow() +
            mesh_mats.Mp_ptr.numRow() + mesh_mats.Mfs_ptr.numRow() );
  }


  for (p=list_panels.initialize(), i=1; p!=NULL;
       p=list_panels.iterate())
  {
    p->assign_number(i);
    i++;
  }

  for (f=list_filaments.initialize(), i=1; f!=NULL;
       f=list_filaments.iterate())
  {
    f->assign_number(i);
    i++;
  }

  for (r=list_resistances.initialize(), i=1; r!=NULL;
       r=list_resistances.iterate())
  {
    r->assign_number(i);
    i++;
  }
  printf ("Number of filaments: %ld\n", list_filaments.size());
  printf ("Number of panels: %ld\n", list_panels.size());

  // tell each incidence matrix how many column has - used later for
  // calculations
/*  if (solve_method == CIRCUIT_SOLVE_MESH)
  {
    mesh_mats.Mf_ptr->set_n_cols (list_filaments.size());
    mesh_mats.Mp_ptr->set_n_cols (list_panels.size());
    mesh_mats.Mfs_ptr->set_n_cols (list_filaments.size());
    mesh_mats.Mps_ptr->set_n_cols (list_panels.size());
    mesh_mats.Mrs_ptr->set_n_cols (list_resistances.size());
  } */


  // if (solve_method == CIRCUIT_SOLVE_MESH)
  //{
  //if   (mesh_mats.Mf_ptr.numCol() != (list_filaments.size())) 
  //std::cout << "error" << std::endl;
  //if  (mesh_mats.Mp_ptr.numCol() != (list_panels.size())) 
  //std::cout << "error" << std::endl;
  //if  (mesh_mats.Mfs_ptr.numCol() != (list_filaments.size()))
  //std::cout << "error" << std::endl;
  //if  (mesh_mats.Mps_ptr.numCol() != (list_panels.size()))
  //std::cout << "error" << std::endl;
  //if  (mesh_mats.Mrs_ptr.numCol() != (list_resistances.size()))
  //std::cout << "error" << std::endl;
  // }

  // output information about panel location (if specified)
  bool output_panels=true;
  if (output_panels == true)
  {
    FILE *fp = fopen ("panels.coords", "w");
    if (fp == NULL)
    {
      printf ("Error: Couldn't open file \"panels.coords\" to store panels "
              "coordinates info.\n");
      exit (-1);
    }
    for (p=list_panels.initialize(); p!=NULL;
         p=list_panels.iterate())
    {
      for (e=1; e<=4; e++)
        fprintf (fp, "[%e,%e,%e] ",
                 (p->get_apex(e))->x(),(p->get_apex(e))->y(),
                 (p->get_apex(e))->z());

      fprintf(fp, "\n");
    }
    fclose (fp);
    printf ("Panels information written to \"panels.coord\".\n");
  }

  // output information about filament location (if specified)
  // each column represent one filament: rows present x1, y1, z1,
  // x2, y2, z2, width, height, sigma,
  // widthdir[XX], widthdir[YY], widthdir[ZZ]
  bool output_filaments = true;
  if (output_filaments == true)
  {
    FILE *fp = fopen ("filaments.coords", "w");
    if (fp == NULL)
    {
      printf ("Error: Couldn't open file \"filaments.coords\" to store "
              "filaments location info.\n");
      exit (-1);
    }
    for (f=list_filaments.initialize(); f!=NULL;
         f=list_filaments.iterate())
    {
      for (e=1; e<=2; e++)
        fprintf(fp, "[%e,%e,%e] ", (f->get_point(e))->x(), (f->get_point(e))->y(),
                (f->get_point(e))->z());
      fprintf (fp, "%e ", f->filWidth());
      fprintf (fp, "%e ", f->filHeight());
      fprintf (fp, "%e ", f->sigma());
      fprintf (fp, "\n");
    }
    fclose (fp);
    printf ("Filaments information written to \"filaments.coord\".\n");
  }

  // output information about load resistances location (if specified)
  // each column represent one resistance: rows present x1, y1, z1,
  // x2, y2, z2, value
  bool output_load_resist = false;
  if (output_load_resist == true)
  {
    FILE *fp = fopen ("load_resist.coords", "w");
    if (fp == NULL)
    {
      printf ("Error: Couldn't open file \"load_resist.coords\" to store "
              "load resistances location info.\n");
      exit (-1);
    }
    for (r=list_resistances.initialize(); r!=NULL;r=list_resistances.iterate())
    {
      for (e=1; e<=2; e++)
        fprintf(fp, "[%e ,%e,%e] ",(r->get_node(e))->x(),(r->get_node(e))->y(),
                (r->get_node(e))->z());
      fprintf(fp, "%e  ", r->get_r_value());
      fprintf (fp, "\n");
    }

    fclose (fp);
    printf ("Load resistances information written to \"load_resist.coord\".\n");
  }

  cout<<"Discretization done"<<endl;
}


void Structure::connect_conductors_meshes ()
{
  Node *node;
  long int mesh_number;
  List<Internal_Node> *list_int_nodes;
  Internal_Node *int_node, *previous_int_node;
  int i;
  panel<double> *panel1_ptr, *panel2_ptr;

  // precaution
/*  if (mesh_mats.Mp_ptr == NULL)
  {
    printf ("Errare humanum est: Structure::connect_conductors_meshes()\n");
    exit (-1);
  }*/

  // set next mesh_number
//  mesh_number = (mesh_mats.Mp_ptr->get_n_lines()) + 1;
  mesh_number = (mesh_mats.Mp_ptr.numRow()) + 1;

  // for each defined node
  for (node=list_defined_nodes.initialize(); node != NULL;
       node=list_defined_nodes.iterate() )
  {
    // only real_nodes
    if ((node->is_real_node()) == 0) continue;

    // do meshes between pairs of panels
    // get list_internal_nodes of this real_node
    list_int_nodes = node->get_list_internal_nodes();

    // do meshes between pairs of panels
    i=0;
    // do not consider internal nodes of load resistances
    for (int_node = list_int_nodes->initialize();
         ((int_node != NULL) && (int_node->is_resist_int_node()) );
         int_node = list_int_nodes->iterate() );

    while (int_node != NULL)
    {
      if ((panel2_ptr = int_node->get_first_panel_ptr()) == NULL)
      {
        printf ("Error: Structure::connect_conductors_meshes() - there "
                "is any panel in this internal node.\n");
        exit(-1);
      }
      if (i != 0) {
        //mesh_mats.Mp_ptr->add_entry (0, mesh_number++, panel1_ptr, -1,
          //                           panel2_ptr, 1, NULL);
	mesh_mats.Mp_ptr.resize_insertElement(mesh_number, panel1_ptr->get_number(), -1);
	mesh_mats.Mp_ptr.insertElement(mesh_number-1, panel2_ptr->get_number(), 1);
	mesh_number++; }
      previous_int_node = int_node;
      panel1_ptr = panel2_ptr;

      i++;
      // do not consider internal nodes of load resistances
      for (int_node = list_int_nodes->iterate();
           ((int_node != NULL) && (int_node->is_resist_int_node()) );
           int_node = list_int_nodes->iterate() );
    }
  }
}

void Structure::connect_resistances_meshes()
{
  Resistance *resist;
  List<Internal_Node> *list_int_nodes;
  Internal_Node *int_node_ptr;
  panel<double> *panel1_ptr, *panel2_ptr;
  long int mesh_number_Mps, mesh_number_Mrs;

  // precaution
/*  if (mesh_mats.Mp_ptr == NULL)
  {
    printf ("Errare humanum est: Structure::connect_resistances_meshes()\n");
    exit (-1);
  }
*/
  // get number of already existing meshes (in Mp); number of meshes in Mrs is 0;
 // mesh_number_Mps = (mesh_mats.Mps_ptr->get_n_lines()) + 1;
  mesh_number_Mps = (mesh_mats.Mps_ptr.numRow()) + 1;
  mesh_number_Mrs = 0+1;

  for (resist=list_resistances.initialize(); resist != NULL;
       resist=list_resistances.iterate() )
  {

    // get panel connected to 1st node of this resistance
    if ((int_node_ptr = ((((resist->get_node(1))->real_node())->
                          get_list_internal_nodes())->initialize())) == NULL)
    {
      printf ("Error: Structure::connect_resistances_meshes() - nothing "
              "is connects to node %s of Resistance %s.\n",
              resist->node_name(1), resist->name() );
      exit(-1);
    }
    if ((panel1_ptr = int_node_ptr->get_first_panel_ptr()) == NULL)
    {
      printf ("Error: Structure::connect_resistances_meshes() - no panels "
              "on this internal node ?\n");
      exit(-1);
    }

    // get panel connected to 2nd node of this resistance
    if ((int_node_ptr = ((((resist->get_node(2))->real_node())->
                          get_list_internal_nodes())->initialize())) == NULL)
    {
      printf ("Error: Structure::connect_resistances_meshes() - nothing "
              "is connects to node %s of Resistance %s.\n",
              resist->node_name(2), resist->name() );
      exit(-1);
    }
    if ((panel2_ptr = int_node_ptr->get_first_panel_ptr()) == NULL)
    {
      printf ("Error: Structure::connect_resistances_meshes() - no panels "
              "on this internal node ?\n");
      exit(-1);
    }

    // add mesh in Mps and Mrs
    //mesh_mats.Mps_ptr->add_entry (0,mesh_number_Mps++, panel1_ptr, -1, panel2_ptr, 1, NULL);
    //mesh_mats.Mrs_ptr->add_entry (0, mesh_number_Mrs++, resist, 1, NULL);
	mesh_mats.Mps_ptr.resize_insertElement(mesh_number_Mps, panel1_ptr->get_number(), -1);
	mesh_mats.Mps_ptr.insertElement(mesh_number_Mps-1, panel2_ptr->get_number(), 1);
		mesh_number_Mps++;
	mesh_mats.Mrs_ptr.resize_insertElement(mesh_number_Mrs++, resist->get_number(),1);
  
  }
}

void Structure::connect_sources_meshes ()
{
  Ext_Node_Pair *ext_node_pair;
  List<Internal_Node> *list_int_nodes;
  Internal_Node *int_node_ptr;
  panel<double> *panel1_ptr, *panel2_ptr;
  long int mesh_number;
  // in case that one of the nodes are infinity
  int dont_add_source_to_1st_node=0, dont_add_source_to_2nd_node=0;


  // precaution
/*  if ((mesh_mats.Mp_ptr == NULL) || (mesh_mats.Vs_ptr == NULL))
  {
    printf ("Errare humanum est: Structure::connect_sources_meshes()\n");
    exit (-1);
  } */
  // get number of already existing meshes (in Mp);
 // mesh_number = (mesh_mats.Mp_ptr->get_n_lines()) + 1;
  mesh_number = (mesh_mats.Mp_ptr.numRow()) + 1;

  for (ext_node_pair=list_ext_node_pairs.initialize(); ext_node_pair != NULL;
       ext_node_pair=list_ext_node_pairs.iterate() )
  {
    if (ext_node_pair->get_node(1) == NULL)
      dont_add_source_to_1st_node=1;     // is infinity (don't add source)
    else
    {
      // get panel connected to 1st node of this ext_node_pair
      if ((int_node_ptr = ((((ext_node_pair->get_node(1))->real_node())->
                            get_list_internal_nodes())->initialize())) == NULL)
      {
        printf ("Error: Structure::connect_conductors_meshes() - nothing "
                "is connects to node %s of External Node %s.\n",
                ext_node_pair->node_name(1), ext_node_pair->name() );
        exit(-1);
      }
      if ((panel1_ptr = int_node_ptr->get_first_panel_ptr()) == NULL)
      {
        printf ("Error: Structure::connect_conductors_meshes() - no panels "
                "on this internal node ?\n");
        exit(-1);
      }
    }

    if (ext_node_pair->get_node(2) == NULL)
      dont_add_source_to_2nd_node=1;     // is infinity (don't add source)
    else
    {
      // get panel connected to 2nd node of this ext_node_pair
      if ((int_node_ptr = ((((ext_node_pair->get_node(2))->real_node())->
                            get_list_internal_nodes())->initialize())) == NULL)
      {
        printf ("Error: Structure::connect_conductors_meshes() - nothing "
                "is connects to node %s of External Node %s.\n",
                ext_node_pair->node_name(2), ext_node_pair->name() );
        exit(-1);
      }
      if ((panel2_ptr = int_node_ptr->get_first_panel_ptr()) == NULL)
      {
        printf ("Error: Structure::connect_conductors_meshes() - no panels "
                "on this internal node ?\n");
        exit(-1);
      }
    }

    if (dont_add_source_to_1st_node && dont_add_source_to_2nd_node)
    {
      printf ("Errare humanum est: Structure::connect_sources_meshes()\n");
      exit (-1);
    }

    /*if ((dont_add_source_to_1st_node == 0) && (dont_add_source_to_2nd_node == 0)){
      mesh_mats.Mp_ptr.resize_insertElement(mesh_number, panel1_ptr->get_number(), -1);
      mesh_mats.Mp_ptr.insertElement(mesh_number-1, panel2_ptr->get_number(), 1);
      mesh_number++;      
      } else */ 
    {
      // add mesh for 1st node
      if (dont_add_source_to_1st_node == 0)
	{
	  //  mesh_mats.Mp_ptr->add_entry (0, mesh_number, panel1_ptr, 1, NULL);
	  //  mesh_mats.Vs_ptr->add_entry (mesh_number++);  // mesh current has to have same
	  mesh_mats.Mp_ptr.resize_insertElement(mesh_number++, panel1_ptr->get_number(), 1);
	  //mesh_mats.Vs_ptr.push_back(panel1_ptr->get_number()/1.);	
	  
	}  // direction as source current
      
      // add mesh for 2nd node
      if (dont_add_source_to_2nd_node == 0)
	{
	  // mesh_mats.Mp_ptr->add_entry (0, mesh_number, panel2_ptr, 1, NULL);
	  // mesh_mats.Vs_ptr->add_entry (mesh_number++);  // mesh current has to have same
	  mesh_mats.Mp_ptr.resize_insertElement(mesh_number++, panel2_ptr->get_number(), 1);
	  //mesh_mats.Vs_ptr.push_back(panel2_ptr->get_number()/1.);	
	  
	}  // direction as source current
    }  
  }
}

void Structure::find_load_resistance()
{
  Ext_Node_Pair *ext_node_pair;
  Resistance *resist;
  char * nodeName1, *nodeName2;

  for (resist = list_resistances.initialize(); resist!=NULL;
       resist=list_resistances.iterate())
  {
    if ((resist->get_node(1) !=NULL) && (resist->get_node(2) !=NULL))
    {
      nodeName1=resist->node_name(1);
      nodeName2=resist->node_name(2);
      ext_node_pair=get_ext_node_pair_from_node_names(nodeName1, nodeName2);
      ext_node_pair->connect_to_resistance(resist);
    }
  }

}


void Structure::prepare_near_conductors()
{
  Node *n;
  conductor *s1, *s2, *chosen;
  int side, nep_tmp;
  double stretch;
  double h, additional_stretch;
  int sign1, sign2;
  
  for (n=list_defined_nodes.initialize(); n!=NULL;
       n=list_defined_nodes.iterate() )
    {
      List <conductor> *list_conductors_node_ptr;
      list_conductors_node_ptr = n->get_list_node_conductors_ptr();
      
      // operations necessary to this node
      if (list_conductors_node_ptr->size() > 1) {
	//test angle
	s1 = list_conductors_node_ptr->initialize();
	for (s2=list_conductors_node_ptr->iterate(); s2!=NULL;
	     s2=list_conductors_node_ptr->iterate())
	  {
	    //if ((paral((s1->vuw()),(s2->vuw())) +
	    //	 paral((s1->vuw()),(s2->vul())) +
	    // paral((s1->vuw()),(s2->vuh())) ) < 2)
	    //{
	    //fprintf (stderr, "Warning: Segments %s and %s share a common "
	    //	 "node not having a right\n\t\tangle between them.\n",
	    //	 s1->get_name(), s2->get_name());
	    //}
	  }
	
	
	// only 2 segments in the same direction
	s1 = list_conductors_node_ptr->initialize();
	s2 = list_conductors_node_ptr->iterate();
	if ((list_conductors_node_ptr->size() == 2) &&
	    paral(s1->vul(), s2->vul()) )
	  {
	    if (n==s1->node(1))
	      s1->set_nep1 (1);
	    else
	      s1->set_nep2 (1);
	    if (n==s2->node(1))
	      s2->set_nep1 (1);
	    else
	      s2->set_nep2 (1);
	  }
	
	// other cases
	else
	  {
	    // choose chosen
	    for (chosen=s1=list_conductors_node_ptr->initialize();
		 s1!=NULL;
		 s1=list_conductors_node_ptr->iterate() )
	      if (((s1->width()) * (s1->height())) >
		  ((chosen->width()) * (chosen->height()))) 
		chosen=s1;
	    
	    // list will contain only the non-chosen conductors
	    list_conductors_node_ptr->delete_element(chosen);
	    
	    // calculate stretch for chosen
	    stretch = 0;
	    nep_tmp = 0;
	    for (s1=list_conductors_node_ptr->initialize(); s1!=NULL;
		 s1=list_conductors_node_ptr->iterate() )
	      {
		if (paral(s1->vuw(), chosen->vul()))
		  if (((s1->width())/2) > stretch)
		    stretch = s1->width()/2;
		if (paral(s1->vuh(),chosen->vul()))
		  if (((s1->height())/2) > stretch)
		    stretch = s1->height()/2;
		if (paral(s1->vul(),chosen->vul()))
		  nep_tmp = 1;
	      }
	    if (n==chosen->node(1)) side = 1;
	    else side = 2;
	    
	    chosen->set_stretch (stretch, side);
	    if (n==chosen->node(1))
		chosen->set_nep1 (nep_tmp);
	    else
	      chosen->set_nep2 (nep_tmp);
	    
	    // operations on non-chosen conductors (will have negative stretch)
	    for (s1=list_conductors_node_ptr->initialize(); s1!=NULL;
		 s1=list_conductors_node_ptr->iterate() )
	      {
		if (n==s1->node(1))
		  side = 1;
		else
		  side = 2;
		
		if (paral(chosen->vuw(),s1->vul()))
		  s1->set_stretch (-(chosen->width())/2, side);
		else if (paral(chosen->vuh(), s1->vul()))
		  s1->set_stretch (-(chosen->height())/2, side);
		else if (paral (chosen->vul(), s1->vul()))
		  s1->set_stretch (-stretch, side);
		else
		  {
		    // Non-usual angle (different than 90 or 180 degrees).
		    // In this case, chosen is not stretch; there should be only one
		    // conductor to the same node; that other conductor will have
		    // a negative stretch of w_chosen/2*cos(angle between conds - PI/2)
		    if (list_conductors_node_ptr->size() > 2)
		      {
			printf ("Error: There should be only two conductors connected"
				" to node %s.\n", n->name_node());
			exit(-1);
		      }
		    // chosen is not stretch in this case
		      chosen->set_stretch (0, side);
		      if (n==chosen->node(1))
			chosen->set_nep1 (1);
		      else
			chosen->set_nep2 (1);
		      // fprintf (stderr, "Angle: %g\n", 180.0/PI*
		      //((chosen->get_vul()).get_angle_with(s1->get_vul())));
		      
		      
		      if (paral (chosen->vuh(), s1->vuh()))
			{
			  if (angle(s1->vul(), s2->vul()) > PI/2)
			    h =
			      cos((angle(chosen->vul(), s1->vul()))-PI/2)*
			      chosen->width() / 2;
			  else
			    {
			      h =
				sin(angle(chosen->vul(), s1->vul()))*
				chosen->width() / 2;
			      h = h + s1->width()/2 /
				tan (angle(chosen->vul(), s1->vul()));
			    }
			  if (fabs(h) > s1->length())
			    {
			      printf ("Errare humanum est: Structure::prepare_near_cond"
				      "uctors() - segment %s shortened %g of its length "
				      "(angle = %g).\n",
				      s1->get_name(), s1->length()/h, 180.0/PI*
				      (angle(chosen->vul(), s1->vul())));
			      exit(-1);
			    }
			  s1->set_stretch ((-h), side);
			}
		      
		      else if (paral(chosen->vuw(), s1->vuw()))
			{
			  if (angle(s1->vul(), s2->vul()) > PI/2)
			    h =
			      cos((angle(chosen->vul(), s1->vul()))-PI/2)*
			      chosen->height() / 2;
			  else
			    {
			      h =
				sin(angle (chosen->vul(), s1->vul()))*
				chosen->height() / 2;
			      h = h + s1->height()/2 /
				tan (angle(chosen->vul(), s1->vul()));
			    }
			  if (fabs(h) > s1->length())
			    {
			      printf ("Errare humanum est: Structure::prepare_near_cond"
				      "uctors() - segment %s shortened %g (angle = %g).\n",
				      s1->get_name(), s1->length()/h, 180.0/PI*
				      (angle(chosen->vul(), s1->vul())));
			      exit(-1);
			    }
			  s1->set_stretch ((-h), side);
			}
		      
		      else
			{
			  printf ("\t\tNon-usual angle between two conductors - length "
				  "directions have two be on the same plane.\n");
			  exit(-1);
			}
		      if (h!=0)
			fprintf (stderr, "\t\tSegment %s shortened 1/%4g of its "
                       "length (in terms of panels).\n",
				 s1->get_name(), s1->length()/h);
		    }
		  
		  if (n==s1->node(1))
		    s1->set_nep1 (1);
		  else
		    s1->set_nep2 (1);
		}
	    }
	}
  }
  
}


Node * Structure::get_node_from_name (char *name)
{
  Node *node;
  Pseudo_Node *pnode;

  node = list_defined_nodes.initialize ();
  while((node != NULL) && (strcmp(name, node->name_node()) != 0))
    node = list_defined_nodes.iterate ();


  if (node != NULL)
    return node;
  else
  {
    pnode = list_pseudo_nodes.initialize ();
    while ((pnode != NULL) && (strcmp(name, pnode->name_node()) != 0))
      pnode = list_pseudo_nodes.iterate ();
    if (pnode == NULL)
      return NULL;
    else
      return pnode->get_defined_node();
  }
}


Ext_Node_Pair *Structure::get_external_from_portname(char *portname)
{
  Ext_Node_Pair *ext;

  if (strcmp(portname," ") == 0)
    return NULL;

  ext = list_ext_node_pairs.initialize();
  for(; ext != NULL; ext = list_ext_node_pairs.iterate() )
    if (strcmp(portname,ext->name()) == 0)
      return ext;

  return NULL;
}

Ext_Node_Pair *Structure::get_ext_node_pair_from_node_names(char *node1_name,
    char *node2_name)
{
  Ext_Node_Pair *ext;
  if ((strcmp(node1_name," ") == 0) && (strcmp(node2_name," ") == 0))
    return NULL;
  for( ext = list_ext_node_pairs.initialize(); ext != NULL;
       ext = list_ext_node_pairs.iterate() )
    if ((strcmp(node1_name,ext->node_name(1)) == 0)&&
        (strcmp(node2_name,ext->node_name(2)) == 0))
      return ext;

  return NULL;
}

Resistance *Structure::get_resistance_from_resistance_name(char *resistance_name)
{
  Resistance *resist;

  if (strcmp(resistance_name,"") == 0)
    return NULL;

  resist = list_resistances.initialize();
  for(; resist != NULL; resist = list_resistances.iterate() )
    if (strcmp(resistance_name,resist->name()) == 0)
      return resist;

  return NULL;
}

Resistance *Structure::get_resistance_from_node_names(char *node1_name, char *node2_name)
{
  Resistance *resist;

  if ((strcmp(node1_name,"") == 0)||(strcmp(node2_name,"") == 0))
    return NULL;

  resist = list_resistances.initialize();
  for(; resist != NULL; resist = list_resistances.iterate() )
    if ((strcmp(node1_name,resist->node_name(1)) == 0) &&
        (strcmp(node2_name,resist->node_name(2)) == 0))
      return resist;

  return NULL;
}

