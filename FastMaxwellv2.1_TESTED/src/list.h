/*************************************************************************\
 
     Fastpep - Fast Parasitic Extraction Program for 3D Geometries
 
 File : list.h - specification and implementation of template List
 Language : C++
 
 Date: June 17, 1997
\*************************************************************************/
#ifndef LIST_H
#define LIST_H


#include <stdio.h>
#include <stdlib.h>


// class template LNode
template <class Type>
class LNode
{
public:
  LNode <Type> *next;
  Type *element;
  LNode (LNode <Type> *next_ptr = 0, Type *elem = 0); // 0 = NULL
};


//
// class template List
//
template <class Type>
class List
{
private:
  long int n_elements;
  LNode <Type> *first_node;
  LNode <Type> *last_node;
  LNode <Type> *iterate_ptr;
  LNode <Type> *iterate_ptr_2;
  void destructor ();
public:
  List ();
  ~List () { destructor (); };
  void delete_all_list_elements ();
  long int size() const { return n_elements; };
  void push_back (Type *new_element);
  void insert_element_beginning (Type *new_element);
  void add_sub_list (List <Type> *sub_list);
  int delete_element (Type *del_element);
  Type * initialize  (int iter_number=1) ;
  Type * iterate (int iter_number=1)  ;
  Type * get_current_element (int iter_number=1);
  Type * get_last_element ();
};


template <class Type>
LNode <Type> ::LNode (LNode <Type> *next_ptr, Type *elem )
{
  next = next_ptr;
  element = elem;
}

template <class Type>
List <Type> ::List ()
{
  n_elements = 0;
  iterate_ptr = NULL;
  iterate_ptr_2 = NULL;
  first_node = NULL;
  last_node = NULL;
}


template <class Type>
void List <Type> :: destructor ()
{
  LNode <Type> *aux_node_ptr2;

  for (LNode <Type> *aux_node_ptr = first_node; aux_node_ptr; )
  {
    aux_node_ptr2 = aux_node_ptr;
    aux_node_ptr = aux_node_ptr->next;
    delete aux_node_ptr2;
  }
}


template <class Type>
void List <Type> :: delete_all_list_elements ()
{

  for (LNode <Type> *aux_node_ptr = first_node; aux_node_ptr; )
  {
    delete aux_node_ptr->element;
    aux_node_ptr = aux_node_ptr->next;
  }
}


template <class Type>
void List <Type> ::push_back (Type *new_element)
{
  LNode <Type> *temp = last_node;
  last_node = new LNode <Type> (NULL, new_element);
  if (temp)
    temp->next = last_node;
  else                              // e' o first element
    first_node = last_node;

  n_elements++;
}


template <class Type>
void List <Type> ::insert_element_beginning (Type *new_element)
{
  LNode <Type> *temp = first_node;
  first_node = new LNode <Type> (NULL, new_element);
  if (temp)
    first_node->next = temp;
  else                              // e' o first element
    last_node = first_node;

  n_elements++;
}


template <class Type>
void List <Type> ::add_sub_list (List <Type> *sub_list)
{
  Type *temp;

  temp = sub_list->initialize();
  while (temp)
  {
    push_back (temp);
    temp = sub_list->iterate();
  }
}



template <class Type>
int List <Type> ::delete_element (Type *del_element)
{
  LNode <Type> *node_ptr;
  LNode <Type> *node_ptr_aux;

  if (del_element == 0) return -1;
  node_ptr = first_node;
  if (first_node->element == del_element)
  {       // is in first node

    if (first_node->next)                       // there is more than one
      first_node = node_ptr->next;
    else
    {                                      // it's the only one
      last_node = NULL;
      first_node = NULL;
    }

    delete node_ptr;
    n_elements --;
  }

  while (node_ptr && node_ptr->next &&
         node_ptr->next->element != del_element)
    node_ptr = node_ptr->next;


  if (node_ptr && node_ptr->next
      && (node_ptr->next->element == del_element))
  {
    node_ptr_aux = node_ptr->next;
    node_ptr->next = node_ptr_aux->next;
    if (last_node == node_ptr_aux)
      last_node = node_ptr;
    /*        delete node_ptr_aux->element;    */
    delete node_ptr_aux;
    n_elements--;
  }
  else
    return -2;

  return 0;
}


template <class Type>
Type *List <Type> ::initialize (int iter_number) 
{
  LNode <Type> **iterate_ptr_ptr;

  if (iter_number==1) iterate_ptr_ptr=&iterate_ptr;
  else if (iter_number==2) iterate_ptr_ptr=&iterate_ptr_2;
  else
  {
    printf ("Errare humanum est: List<Type>::initialize()\n");
    exit(-1);
  }

  (*iterate_ptr_ptr) = first_node;
  if ((*iterate_ptr_ptr) == NULL)
    return NULL;
  else
    return (*iterate_ptr_ptr)->element;
}



template <class Type>
Type *List <Type> ::iterate (int iter_number) 
{
  LNode <Type> **iterate_ptr_ptr;

  if (iter_number==1) iterate_ptr_ptr=&iterate_ptr;
  else if (iter_number==2) iterate_ptr_ptr=&iterate_ptr_2;
  else
  {
    printf ("Errare humanum est: List<Type>::initialize()\n");
    exit(-1);
  }

  if (*iterate_ptr_ptr)
    *iterate_ptr_ptr = (*iterate_ptr_ptr)->next;
  if (*iterate_ptr_ptr) return (*iterate_ptr_ptr)->element;
  return NULL;
}


template <class Type>
Type *List <Type> ::get_current_element (int iter_number)
{
  LNode <Type> **iterate_ptr_ptr;

  if (iter_number==1) iterate_ptr_ptr=&iterate_ptr;
  else if (iter_number==2) iterate_ptr_ptr=&iterate_ptr_2;
  else
  {
    printf ("Errare humanum est: List<Type>::initialize()\n");
    exit(-1);
  }

  if (*iterate_ptr_ptr)
    return (*iterate_ptr)->element;
  return NULL;
}


template <class Type>
Type *List <Type> ::get_last_element ()
{
  if (last_node)
    return last_node->element;
  return NULL;
}

#endif
