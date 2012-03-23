/*************************************************************************\

     Fastpep - Fast Parasitic Extraction Program for 3D Geometries

 File : resistance.cpp - class Resistance implementation
 Language : C++

 Date: August 12, 1997
\*************************************************************************/

#include "resistance.h"


/************************************\
 Method: Resistance::Resistance
 Description: constructor
 Arguments: 2 nodes and the (optional) name
 Return:
\************************************/
Resistance::Resistance (Node *first_node, Node *second_node, double r_value, char *name=0) {
  node1 = first_node;
  node2 = second_node;
  resistance_value = r_value;
  if (name) {
    resistance_name = new char[strlen(name)+1];
    strcpy (resistance_name, name);
  }
  else
    resistance_name = NULL;
  number = 0;
}


/************************************\
 Method: Resistance::node_name - public
 Description: return name of specified node
 Arguments: 
 Return: name of specified node
\************************************/
char *Resistance::node_name (const int node_number) const {
  switch (node_number) {
  case 1:
    return node1->name_node();
    break;
  case 2:
    return node2->name_node();
    break;
  }
}


/************************************\
 Method: Resistance::get_node - public
 Description: return pointer to specified node
 Arguments: 
 Return: pinterspecified node
\************************************/
Node *Resistance::get_node (const int node_number) const {
  switch (node_number) {
  case 1:
    return node1;
    break;
  case 2:
    return node2;
    break;
  }
}

/************************************\
 Method: Resistance::create_connect_info - public
 Description: create 2 internal nodes for each resistance.
              
              Each internal_node has a pointer to this resistance, with direction.
	      This is used later for nodal analysis.

              These internal_nodes will be added later to
	      global_list_of_internal_nodes.

	      This internal nodes are created here because load resistances
	      are not created durinf discreization.
	      This internal nodes are slitly different - don't have panels nor
	      filaments in its lists.
 Arguments: 
 Return: void
\************************************/
void Resistance::create_connect_info (List <Internal_Node>
			                 **list_internal_nodes_resistance_ptr) {
  Internal_Node *int_node;

  list_internal_nodes.push_back (new Internal_Node (IN_RESIST) );
  list_internal_nodes.push_back (new Internal_Node (IN_RESIST) );


  int_node = list_internal_nodes.initialize();
  int_node->add_resistance_with_direction (this, 1);

  int_node = list_internal_nodes.get_last_element();
  int_node->add_resistance_with_direction (this, -1);


  (node1->real_node())->
    add_int_node_ptr (list_internal_nodes.initialize());
  (node2->real_node())->
    add_int_node_ptr (list_internal_nodes.get_last_element());


  *list_internal_nodes_resistance_ptr = &list_internal_nodes;
}
