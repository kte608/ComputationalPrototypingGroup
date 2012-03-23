/*************************************************************************\
 
     Fastpep - Fast Parasitic Extraction Program for 3D Geometries
 
 File : resistance.h - class Resistance specification
 Language : C++
 
 Date: August 12, 1997
\*************************************************************************/
#ifndef RESISTANCE_H
#define RESISTANCE_H


#include "node.h"
#include "internal_node.h"

class Resistance
{
public:
  Resistance (Node *first_node, Node *second_node, double r_value, char *name);
  char *node_name (const int node_number) const;
  double get_r_value () const { return resistance_value; };
  char *name () const { return resistance_name; };
  Node *get_node (const int node_number) const;
  void assign_number (const long int _number) { number = _number; };
  long int get_number () const { return number; };
  void create_connect_info (List <Internal_Node>
                            **list_internal_nodes_resistance_ptr);
private:
  Node *node1, *node2;
  double resistance_value;
  char *resistance_name;
  List <Internal_Node> list_internal_nodes;
  long int number;
};

#endif
