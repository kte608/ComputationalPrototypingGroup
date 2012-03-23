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
#include "conductor.h"
#include "node.h"
#include "panel.h"
//#include "henry.h"
#include <iostream>
#include <math.h>
#include <stdlib.h>

conductor::conductor (char *_name, Node *_node1, Node *_node2,
                      double _width, double _height,
                      point_vec_operations::point3D<double> *_widthdir, double _sigma,
                      int _nl, int _fnw, int _fnh, int _pnw, int _pnh,
                      double _frw, double _frh, double _pfracl, double _pfracw, double _pfrach,
                      int _type, int _nep)
{
  name_=new char[strlen(_name)+1];
  strcpy(name_, _name);
  node1_=_node1;
  node2_=_node2;
  width_=_width;
  height_=_height;
  widthdir_=_widthdir;
  sigma_=_sigma;
  nl_=_nl;
  fnw_=_fnw;
  fnh_=_fnh;
  pnw_=_pnw;
  pnh_=_pnh;
  frw_=_frw;
  frh_=_frh;
  pfracl_=_pfracl;
  pfracw_=_pfracw;
  pfrach_=_pfrach;
  type_=_type;
  nep1_=_nep;
  nep2_=_nep;

  num_filaments_ = num_panels_ = 0;
  stretch1_ = stretch2_ = 0;

  // unit vectors along each of the 3 directions of the filament
  point3D<double> pn1 ((*node1_).x(),(*node1_).y(),(*node1_).z()), pn2 ((*node2_).x(),(*node2_).y(),(*node2_).z());

  vul_ = (pn2 - pn1);
  vul_ = vul_ / (vul_.abs());

  vuw_ = *widthdir_;
  vuw_ = vuw_ / (vuw_.abs());

  vuh_ = cross(vul_, vuw_);
  vuh_ = vuh_ / (vuh_.abs());

  // length
  length_ = (pn1 - pn2).abs();
  have_edge1_ = have_edge2_ = 1;
  have_2nd_edge1_ = have_2nd_edge2_ = 0;


}

conductor::~conductor()
{
  delete name_;
  delete widthdir_;
}


Node *conductor::node (int number)
{
  switch (number)
  {
  case 1:
    return node1_;
    break;
  case 2:
    return node2_;
    break;
  default:
    return NULL;
    break;
  }
}


void conductor::set_nep1 (int _nep)
{
  if ((_nep != 0) && (_nep !=1))
  {
    fprintf (stderr, "Error: set_nep1()\n");
    exit (-1);
  }
  nep1_ = _nep;
}



void conductor::set_nep2 (int _nep)
{
  if ((_nep != 0) && (_nep !=1))
  {
    fprintf (stderr, "Error: set_nep2()\n");
    exit (-1);
  }
  nep2_ = _nep;
}

void conductor::set_stretch (double _stretch, int side)
{
  if (side==1)
    stretch1_ = _stretch;
  else if (side==2)
    stretch2_ = _stretch;
  else
  {
    fprintf (stderr, "Error: set_stretch\n");
    exit (-1);
  }
}


void conductor::discretize(fastsub::SolverMethod solve_method,
                           List<panel<double> >** list_panels_conductor_ptr,
                           List<filament<double> >** list_filaments_conductor_ptr,
                           //List <SEGMENT> **list_henry_segments_conductor_ptr,
                           List <Internal_Node> **list_internal_nodes_conductor_ptr,
                           pfft::SpRowMat<double> *Mf_conductor_ptr,
                           pfft::SpRowMat<double> *Mp_conductor_ptr,
                           pfft::SpRowMat<double> *Mfs_conductor_ptr,
                           pfft::SpRowMat<double> *Mps_conductor_ptr,
			   int panels_start_index, int filaments_start_index,
			   int *globalmeshline_p, int *globalmeshline_f, int *globalmeshline_s)
{
  int ret;
  global_panels_index_ = panels_start_index;
  global_filaments_index_ = filaments_start_index;
  global_mesh_line_p = globalmeshline_p;
  global_mesh_line_f = globalmeshline_f;
  global_mesh_line_s = globalmeshline_s;
  //panel discretization
  ret=discretize_panels(solve_method, list_panels_conductor_ptr,Mp_conductor_ptr);
  if (ret == -1)
  {
    printf ("Errare humanum est: Segment::discretize () - panels\n");
    exit (-1);
  }


  //filament discretization
  ret = discretize_filaments (solve_method,
                              list_filaments_conductor_ptr,
                              //list_henry_segments_conductor_ptr,
                              Mf_conductor_ptr, Mfs_conductor_ptr,
                              Mps_conductor_ptr);
  if (ret == -1)
  {
    printf ("Errare humanum est: Segment::discretize () - filaments\n");
    exit (-1);
  }


  // set conductor variables
  num_filaments_ = (*list_filaments_conductor_ptr)->size();
  num_panels_ = (*list_panels_conductor_ptr)->size();

  // return list of internal_nodes - should not be deleted
  *list_internal_nodes_conductor_ptr = &list_internal_nodes;

}


int conductor::discretize_panels(fastsub::SolverMethod solve_method,
                                 List<panel<double> >** list_panels_conductor_ptr,
                                 pfft::SpRowMat<double> *Mp_conductor_ptr)
{
  double innerl, innerw, innerh;
  double pfracl_edge, pfracl_tip;
  double pfracl_edge1, pfracl_edge2, pfracl_tip1, pfracl_tip2;
  double pfracl_2nd_edge1, pfracl_2nd_edge2, pfracl_2nd_tip1, pfracl_2nd_tip2;
  point3D<double> p1, p2, p3;
  Internal_Node *int_node;
  panel<double> *pa, *panel1, *panel2;
  int section;
  int side;
  double l, w, h;
  int i, j, k;
  long int mesh_line;

  *list_panels_conductor_ptr = new List <panel<double> >;

  // correct frac, if necessary
  if (pnw_ <= 2) pfracw_ = 1;
  if (pnh_ <= 2) pfrach_ = 1;

  // size of inner panels (w and h)
  innerw = width_ / (pnw_-2 + 2*pfracw_);
  innerh = height_ / (pnh_-2 + 2*pfrach_);

  // length
  innerl = length_/nl_;
  pfracl_edge = pfracl_edge1 = pfracl_edge2 = pfracl_;
  pfracl_tip = pfracl_tip1 = pfracl_tip2 = 0.5-pfracl_;

  if ((stretch1_ < 0) && (stretch2_ < 0)) {
    innerl = (length_-fabs(stretch1_)-fabs(stretch2_))/nl_;
    have_edge1_ = 1;
    have_edge2_ = 1;
  }
  else if  (stretch1_ < 0) {
    innerl = (length_-fabs(stretch1_))/nl_;
    have_edge1_ = 1;
    have_edge2_ = 1;
  } 
  else if  (stretch2_ < 0) {
    innerl = (length_-fabs(stretch2_))/nl_;
    have_edge1_ = 1;
    have_edge2_ = 1;
  } 
  
  /*
  // in case that more segments are connect to node 1
  if (stretch1_ < 0)
  {
    if (fabs(stretch1_) < (innerl*pfracl_edge1))
    {
      have_edge1_ = 1;
      pfracl_edge1 = pfracl_edge1 - fabs(stretch1_)/innerl;
      //pfracl_edge1 = pfracl_edge1 + fabs(stretch1_)/innerl;
    }
    else
    {
      have_edge1_ = 0;
      pfracl_tip1 = 0.5 - fabs(stretch1_)/innerl;
      //pfracl_tip1 = pfracl_tip1 + fabs(stretch1_)/innerl;
      if (pfracl_tip1 <= 0)
      {
        fprintf (stderr, "Error: Cross section of Segments connected to "
                 "node %s is to big to be possible 2 or more segments "
                 "to share a common node. Reduce "
                 "number of sections on conductors connectd to this node. "
                 "You can also use ground-planes to model this connection.\n",
                 node1_->name_node ());
	exit (-1);
      }
    }
  }
  // in case that more segments are connect to node 2
  if (stretch2_ < 0)
  {
    if (fabs(stretch2_) < (innerl*pfracl_edge2))
    {
      have_edge2_ = 1;
      pfracl_edge2 = pfracl_edge2 - fabs(stretch2_)/innerl;
      //pfracl_edge2 = pfracl_edge2 + fabs(stretch2_)/innerl;
    }
    else
    {
      have_edge2_ = 0;
      pfracl_tip2 = 0.5 - fabs(stretch2_)/innerl;
      //pfracl_tip2 = pfracl_tip2 + fabs(stretch2_)/innerl;
      if (pfracl_tip2 <= 0)
      {
        fprintf (stderr, "Error: Cross section of Segments connected to "
                 "node %s is to big to be possible 2 or more segments "
                 "to share a common node. Reduce "
                 "number of sections on conductors connectd to this node. "
                 "You can also use ground-planes to model this connection.\n",
                 node2_->name_node ());
	exit (-1);
      }
    }
  }
  */

  // center of section
  p1.set (node(1)->x(), node(1)->y(), node(1)->z());
  if (stretch1_ < 0)
    p1 = p1 + (vul() * (-stretch1_));
    //p1 = p1 + (vul() * (stretch1_));

  // for each section
  for (section=1; section <= (nl_+3); section++)
  {
    // add internal node to list_internal_nodes and set int_node
    if ((section != 2) && (section != (nl_+3)))
    {
      list_internal_nodes.push_back (new(Internal_Node));
      if (section == 1)
        int_node=list_internal_nodes.initialize();
      else if ((section != 2) && (section != (nl_+3)))
        int_node=list_internal_nodes.iterate();
    }

    if ((section == 1) && (have_edge1_ == 0))
      continue;
    if ((section == (nl_+3)) && (have_edge2_ == 0))
      continue;
    if (section==1) l = innerl * pfracl_edge1;
    else if (section==(nl_+3)) l = innerl * pfracl_edge2;
    else if (section==2) l = innerl * pfracl_tip1;
    else if (section==(nl_+2)) l = innerl * pfracl_tip2;
    else l = innerl;

    //top and bottom sides
    p2 = (p1 + (vuh_*height_/2) - (vuw_*width_/2));
    for (i=0; i<2; i++)
   // for (i=0; i<1; i++)
    {
      for (j=1; j<=pnw_; j++)
      {     // first top then bottom
        if ((j==1) || (j==pnw_)) w = innerw*pfracw_;
        else w = innerw;
        pa = new panel<double> (point3D<double> (p2),
                                point3D<double> (p2+(vuw_*w)),
                                point3D<double> (p2+(vul_*l)+(vuw_*w)),
                                point3D<double> (p2+(vul_*l)));
	pa->setdirection(vul_);
	pa->assign_number(global_panels_index_+num_panels_);
	num_panels_++;
        if (pa->area()==0) {
          cout << p2.x() << " " << p2.y() << " " << p2.z() << endl;
        };
        (*list_panels_conductor_ptr)->push_back (pa);
	//cout<<pa->get_number()<<endl;
        (int_node->list_internal_node_panels).push_back(pa);
        p2 = (p2 + (vuw_*w));
      }
      p2 = (p1 - (vuh_*height_/2) - (vuw_*width_/2));
    }

	/*panel<double>* temp;
	for (temp=(int_node->list_internal_node_panels).initialize(); temp!=NULL;
			temp=(int_node->list_internal_node_panels).iterate())
	cout<<temp->get_number()<<endl;*/
	
    // left and right sides
    p2 = (p1 - (vuh_*height_/2) - (vuw_*width_/2));
    for (i=0; i<2; i++)
    {
      for (j=1; j<=pnh_; j++)
	{     // first left then right
	  if ((j==1) || (j==pnh_)) h = innerh*pfrach_;
	  else h = innerh;
	  
	  pa =  new panel<double> (point3D<double> (p2),
				   point3D<double> (p2+(vuh_*h)),
				   point3D<double> (p2+(vul_*l)+(vuh_*h)),
				   point3D<double>  (p2+(vul_*l)));
	  pa->setdirection(vul_);
	  pa->assign_number(global_panels_index_+num_panels_);
	  num_panels_++;
	  if (pa->area()==0) {
	    cout << p2.x() << " " << p2.y() << " " << p2.z() << endl;
	  };
	  (*list_panels_conductor_ptr)->push_back (pa);
	  (int_node->list_internal_node_panels).push_back (pa);
	  p2 = (p2 + (vuh_*h));
	}
      p2 = (p1 - (vuh_*height_/2) + (vuw_*width_/2));
    }
    
    // prepare p1 for next section
    p1 = (p1 + (vul_*l));
  }
  
  // conductor is stretched - add panels to beginning or end of segment
  if ((stretch1_ > 0) || (stretch2_ > 0))
    {
      // set have_2nd_edge1, pfracl_2nd_edge1 and pfracl_2nd_tip1
      if (stretch1_ > 0)
	{
	  if (stretch1_ <= (innerl * pfracl_edge1))
	    {
	      have_2nd_edge1_ = 0;
	      pfracl_2nd_tip1 = stretch1_/innerl;
	    }
	  else
	    {
	      have_2nd_edge1_ = 1;
	      pfracl_2nd_edge1 = pfracl_edge1;
	      pfracl_2nd_tip1 = stretch1_/innerl - pfracl_2nd_edge1;
	    }
	}
      
      // set have_2nd_edge2, pfracl_2nd_edge2 and pfracl_2nd_tip2
      if (stretch2_ > 0)
	{
	  if (stretch2_ <= (innerl * pfracl_edge2))
	    {
	      have_2nd_edge2_ = 0;
	      pfracl_2nd_tip2 = stretch2_/innerl;
	    }
	  else
	    {
	      have_2nd_edge2_ = 1;
	      pfracl_2nd_edge2 = pfracl_edge2;
	      pfracl_2nd_tip2 = stretch2_/innerl - pfracl_2nd_edge2;
	    }
	}
      
      // add panels for stretched conductors
      for (side=1; side<=2; side++)
	{
	  if (((side==1) && (stretch1_ <= 0)) || ((side==2) && (stretch2_ <= 0)))
	    continue;
	  if (side==1)
	    int_node = list_internal_nodes.initialize();
	  else if (side==2)
	    int_node = list_internal_nodes.get_last_element();
	  
	  // center of section
	  p1.set (node(side)->x(), node(side)->y(), node(side)->z());
	  if ((side==1) && (stretch1_ > 0))
	    p1 = p1 - (vul() * stretch1_);
	  
	  // for each section (edge and tip)
	  for (section=1; section <= 2; section++)
	    {
	      if (((section == 1) && (side==1) && (have_2nd_edge1_ == 0)) ||
		  ((section == 2) && (side==2) && (have_2nd_edge2_ == 0)))
		continue;
	      if (side==1)
		{
		  if (section==1) l = innerl * pfracl_2nd_edge1;
		  else l = innerl * pfracl_2nd_tip1;
		}
	      if (side==2)
		{
		  if (section==1) l = innerl * pfracl_2nd_tip2;
		  else l = innerl * pfracl_2nd_edge2;
		}
	      
	      // top and bottom sides
	      p2 = (p1 + (vuh_*height_/2) - (vuw_*width_/2));
	      for (i=0; i<2; i++)
		{
		  for (j=1; j<=pnw_; j++)
		    {     // first top then bottom
		      if ((j==1) || (j==pnw_)) w = innerw*pfracw_;
		      else w = innerw;
		      pa = new panel<double> (point3D<double> (p2),
					      point3D<double> (p2+(vuw_*w)),
					      point3D<double> (p2+(vul_*l)+(vuw_*w)),
					      point3D<double> (p2+(vul_*l)));
		      pa->setdirection(vul_);
		      pa->assign_number(global_panels_index_+num_panels_);
		      num_panels_++;
		      if (pa->area()==0) {
			cout << p2.x() << " " << p2.y() << " " << p2.z() << endl;
		      };
		      (*list_panels_conductor_ptr)->push_back (pa);
		      (int_node->list_internal_node_panels).push_back (pa);
		      p2 = (p2 + (vuw_*w));
		    }
		  p2 = (p1 - (vuh_*height_/2) - (vuw_*width_/2));
		}
	      // left and right sides
	      p2 = (p1 - (vuh_*height_/2) - (vuw_*width_/2));
	      for (i=0; i<2; i++)
		{
		  for (j=1; j<=pnh_; j++)
		    {     // first left then right
		      if ((j==1) || (j==pnh_)) h = innerh*pfrach_;
		      else h = innerh;
		      
		      pa = new panel<double> (point3D<double> (p2),
					      point3D<double>(p2+(vuh_*h)),
					      point3D<double>(p2+(vul_*l)+(vuh_*h)),
					      point3D<double>(p2+(vul_*l)));
		      pa->setdirection(vul_);
		      pa->assign_number(global_panels_index_+num_panels_);
		      num_panels_++;
		      if (pa->area()==0) {
			cout << p2.x() << " " << p2.y() << " " << p2.z() << endl;
		      };
		      (*list_panels_conductor_ptr)->push_back (pa);
		      (int_node->list_internal_node_panels).push_back (pa);
		      p2 = (p2 + (vuh_*h));
		    }
		  p2 = (p1 - (vuh_*height_/2) + (vuw_*width_/2));
		}
	      
	      // prepare p1 for next section
	      p1 = (p1 + (vul_*l));
	    }
	}
    }
  // cross section sides
  for (k=1; k<=2; k++)
    {
      if ((k==1) && nep1_) continue;
      if ((k==2) && nep2_) continue;
      if (((k==1) && (stretch1_<0)) || ((k==2) && (stretch2_<0)))
	{
	  fprintf (stderr, "Errare humanum est: Segment::discretize_panels()\n");
	  exit (-1);
	}
      if (k==1)
	int_node = list_internal_nodes.initialize();
      else if (k==2)
	int_node = list_internal_nodes.get_last_element();
      p1.set (node(k)->x(), node(k)->y(), node(k)->z());
      if ((k==1) && (stretch1_>0))
	p1 = p1 - (vul_*stretch1_);
      if ((k==2) && (stretch2_>0))
	p1 = p1 + (vul_*stretch2_);
      p2 = (p1 - (vuh_*height_/2) - (vuw_*width_/2));
      for (i=1; i<=pnw_; i++)
	{
	  if ((i==1) || (i==pnw_)) w = innerw*pfracw_;
	  else w = innerw;
	  p3 = p2;
	  for (j=1; j<=pnh_; j++)
	    {
	      if ((j==1) || (j==pnh_)) h = innerh*pfrach_;
	      else h = innerh;
	      
	      pa = new panel<double> (point3D<double> (p3),
				      point3D<double> (p3+(vuw_*w)),
				      point3D<double> (p3+(vuw_*w)+(vuh_*h)),
				      point3D<double> (p3+(vuh_*h)));
	      pa->setdirection(vul_);
	      pa->assign_number(global_panels_index_+num_panels_);
	      num_panels_++;
	      if (pa->area()==0) {
		cout << p2.x() << " " << p2.y() << " " << p2.z() << endl;
	      };
	      (*list_panels_conductor_ptr)->push_back (pa);
	      (int_node->list_internal_node_panels).push_back (pa);
	      p3 = (p3 + (vuh_*h));
	    }
	  p2 = (p2 + (vuw_*w));
	}
    }
  
  // built meshes through panels, if input_options specify MESH_analysis
  if (solve_method == fastsub::CIRCUIT_SOLVE_MESH)
    {
      //  (*Mp_conductor_ptr) = new Incidence_Matrix<panel<double> >;
      //(*Mp_conductor_ptr) = new Incidence_Matrix<panel<double> >;
      /*(*Mp_conductor_ptr) = pfft::SpRowMat<double>();*/
      // for each internal_node
      for (int_node = list_internal_nodes.initialize (), mesh_line=1;
	   int_node != NULL; int_node = list_internal_nodes.iterate () )
	{
	  panel1 = (int_node->list_internal_node_panels).initialize();
	  if (panel1 == NULL)
	    {
	      fprintf (stderr, "Errare humanum est: Segment::discretize_panels() - "
		       "any panel in this node?\n");
	      exit (-1);
	    }
	  // for each pair of panels of an internal node
	  for (panel2 = (int_node->list_internal_node_panels).iterate();
	       panel2 != NULL;
	       panel1 = panel2,
		 panel2 = (int_node->list_internal_node_panels).iterate() ) {
	    //(*Mp_conductor_ptr)->add_entry (0, mesh_line++, panel1, -1,
	    //                              panel2, 1, NULL);		   
	    Mp_conductor_ptr->resize_insertElement(*global_mesh_line_p, panel1->get_number(), -1);
	    Mp_conductor_ptr->insertElement((*global_mesh_line_p)-1, panel2->get_number(),1);
	    (*global_mesh_line_p)++;
	    mesh_line++;
	  }
	  
	}
    }
  else
    (*Mp_conductor_ptr) = pfft::SpRowMat<double>();
  //NULL;
  
  return 0;
}


int conductor::discretize_filaments (fastsub::SolverMethod solve_method,
                                     List <filament<double> > **list_filaments_conductor_ptr,
                                     //List <SEGMENT> **list_henry_segments_conductor_ptr,
                                     pfft::SpRowMat<double> *Mf_conductor_ptr,
                                     pfft::SpRowMat<double> *Mfs_conductor_ptr,
                                     pfft::SpRowMat<double> *Mps_conductor_ptr)
{
  Node **the_nodes;
  conductor **the_segs;
  List <filament<double> > *semi_list_filaments_ptr;
  List <filament<double> > *previous_semi_list_filaments_ptr=NULL;
  Internal_Node *int_node, *previous_int_node;
  filament<double> *filament1, *filament2;
  char name[30];
  double x,y,z,dx,dy,dz;
  int i,j,k,m;
  long int mesh_line, mesh_line_s;

  *list_filaments_conductor_ptr = new List <filament<double> >;
  //*list_henry_segments_conductor_ptr = new List <SEGMENT>;
  the_nodes = new Node* [nl_+1];
  the_segs = new conductor* [nl_];
  
  x = node1_->x();
  y = node1_->y();
  z = node1_->z();
  dx = (node2_->x() - node1_->x())/nl_;
  dy = (node2_->y() - node1_->y())/nl_;
  dz = (node2_->z() - node1_->z())/nl_;


  // divide this segment into (nl) sub-segments with (nl+1) nodes
  for(i = 0; i <= nl_; i++)
  {
    sprintf(name, "N%d", i);
    the_nodes[i] = new Node (name, x, y, z, NORMAL);

    if (i > 0)
    {
      sprintf(name, "E%d", i - 1);

      the_segs[i - 1] = new conductor (name, the_nodes[i - 1], the_nodes[i],
                                       width_, height_, widthdir_,
                                       sigma_, 1, fnw_,
                                       fnh_, pnw_, pnh_,frw_, frh_, 1, 1, 1, NORMAL, 0);
    }
    x += dx;
    y += dy;
    z += dz;
  }

  // there are now (nl+1) nodes and (nl) sub-segments ready for discretization;
  // discretize each sub-segment;
  // fill internal nodes, which have already pointers to panels;
  // built meshes - Mf, Mfs and Mps are created (if input_options
  // specify MESH_analysis).

  for(i=0, mesh_line=1, mesh_line_s=1, int_node=list_internal_nodes.initialize();
      i <= nl_ ;
      i++, int_node=list_internal_nodes.iterate() )
    {
      
      // discretize and add filaments of each conductor to this conductor list
      if (i != nl_)
	{
	  semi_list_filaments_ptr = call_assign_fil(the_segs[i]);
	  (*list_filaments_conductor_ptr)->add_sub_list (semi_list_filaments_ptr);
	}
      //fill internal_nodes list of filament pointers
      // add these filaments (left side of this node) to list in internal_node
      if (i!=0)
	int_node->add_list_filaments_with_direction( previous_semi_list_filaments_ptr, -1);
      
      // add filaments (right side of this node) to list in internal_node
      if (i!=nl_)
	int_node->add_list_filaments_with_direction (semi_list_filaments_ptr,1);
      
      if (solve_method == fastsub::CIRCUIT_SOLVE_MESH)
	{
	  if (i!=0)
	    {
	      filament1=previous_semi_list_filaments_ptr->initialize();
	      if (filament1 == NULL)
		{
		  fprintf (stderr, "Errare humanum est: Segment::discretize_"
			   "filaments() - any filament in this sub-segment?\n");
		  exit (-1);
		}
	      
	      //meshes through 1st filament and 1 panel in each of its nodes
	      //(*Mfs_conductor_ptr)->add_entry (0, mesh_line_s, filament1, 1, NULL);
	      //(*Mps_conductor_ptr)->add_entry (0, mesh_line_s++,
	      //			       (previous_int_node->list_internal_node_panels).initialize(),
	      //			       -1,
	      //			       (int_node->list_internal_node_panels).initialize(),
	      //			       1, NULL);
	      filament1->assign_number(global_filaments_index_+num_filaments_);
	      num_filaments_++;
	      Mfs_conductor_ptr->resize_insertElement(*global_mesh_line_s, filament1->get_number(), 1);
	      Mps_conductor_ptr->resize_insertElement(*global_mesh_line_s, 
						      (previous_int_node->list_internal_node_panels).initialize()->get_number(), -1);
	      Mps_conductor_ptr->insertElement((*global_mesh_line_s)-1, 
					       (int_node->list_internal_node_panels).initialize()->get_number(), 1);
	      (*global_mesh_line_s)++;
	      mesh_line_s++;
	      
	      
	      // meshes only through filaments
	      for (filament2=previous_semi_list_filaments_ptr->iterate();
		   filament2 != NULL;
		   filament1 = filament2,
		     filament2=previous_semi_list_filaments_ptr->iterate()) 
		{
		  //(*Mf_conductor_ptr)->add_entry (0, mesh_line++, filament1, 1,
		  //                              filament2, -1, NULL);
		  filament2->assign_number(global_filaments_index_+num_filaments_); 
		  num_filaments_++;
		  Mf_conductor_ptr->resize_insertElement(*global_mesh_line_f, filament1->get_number(), 1);
		  Mf_conductor_ptr->insertElement((*global_mesh_line_f)-1, filament2->get_number(),-1);
		  (*global_mesh_line_f)++;
		  mesh_line++; 
		}
	    }
	}
      // prepare next sub-segment
      delete previous_semi_list_filaments_ptr;
      previous_semi_list_filaments_ptr = semi_list_filaments_ptr;
      previous_int_node = int_node;
    }
  
  if (((*list_filaments_conductor_ptr)->size()) !=
      (nl_*fnw_*fnh_))
    {
      fprintf(stderr, "Incorrect number of computed filaments.\n");
      return -1;
    }
  
  /*
    for(i = 0; i <= nl_; i++) delete the_nodes[i];
    delete [] the_nodes;       // works ??
    for(i = 0; i < nl_; i++)  delete the_segs[i];
    delete [] the_segs;        // works ??
  */
  
  return 0;
  
}




List <filament<double> > * conductor::call_assign_fil (conductor *seg)
//, List <SEGMENT> *list_henry_segments_conductor_ptr)
{

 // charge chgdummy;
 // int intdummy=0, j;
  //SEGMENT *henry_seg;
  //FILAMENT *henry_fil;
  //NODES *henry_node1, *henry_node2;
  filament<double> *fil;
  List <filament<double> > *local_list_filaments;

  local_list_filaments = new List<filament<double> >;

/*  // create henry segment
  henry_seg = new SEGMENT;

  // add it to list of henry_segments of this conductor
  list_henry_segments_conductor_ptr->push_back(henry_seg);

  // fill SEGMENT
  henry_seg->widthdir = new double [3];
  (henry_seg->widthdir) [XX] = (seg->widthdir())->x();
  (henry_seg->widthdir) [YY] = (seg->widthdir())->y();
  (henry_seg->widthdir) [ZZ] = (seg->widthdir())->z();
  henry_seg->length = sqrt(
                        (seg->node(1)->x() - seg->node(2)->x())*
                        (seg->node(1)->x() - seg->node(2)->x())
                        + (seg->node(1)->y() - seg->node(2)->y())*
                        (seg->node(1)->y() - seg->node(2)->y())
                        + (seg->node(1)->z() - seg->node(2)->z())*
                        (seg->node(1)->z() - seg->node(2)->z())
                      );
  henry_node1 = new NODES;
  henry_node2 = new NODES;
  henry_node1->x = (seg->node(1))->x();
  henry_node1->y = (seg->node(1))->y();
  henry_node1->z = (seg->node(1))->z();
  henry_node2->x = (seg->node(2))->x();
  henry_node2->y = (seg->node(2))->y();
  henry_node2->z = (seg->node(2))->z();
  henry_seg->node[0] = henry_node1;
  henry_seg->node[1] = henry_node2;
  henry_seg->width = seg->width();
  henry_seg->height = seg->height();
  henry_seg->hinc = seg->fnh();
  henry_seg->winc = seg->fnw();
  henry_seg->r_height = seg->frh();
  henry_seg->r_width = seg->frw();
  henry_seg->sigma = seg->sigma();
  henry_seg->area = henry_seg->width * henry_seg->height;
  //henry_seg->type = seg->type();
  henry_seg->number = -1;
  henry_seg->name = seg->get_name();
  henry_seg->loops = NULL;
  henry_seg->is_deleted = 0;
  henry_seg->next = NULL;

  // call assignFil (don't connect (and don't use) fasthenry charges)

  assignFil (henry_seg, &intdummy, &chgdummy);*/
  local_list_filaments = assignFil (seg);

  // extract information: now there is a matrix of filaments in henry_seg

/*  for (j=0; j < henry_seg->num_fils; j++)
  {
    henry_fil = &(henry_seg->filaments[j]);

    point3D<double> p1 (henry_fil->x[0], henry_fil->y[0], henry_fil->z[0]);
    point3D<double> p2 (henry_fil->x[1], henry_fil->y[1], henry_fil->z[1]);
    fil = new filament<double> (p1, p2, henry_fil->width, henry_fil->height,
                                henry_seg->sigma, henry_seg, henry_fil);
    local_list_filaments->push_back (fil);

    // Free memory allocated in assignFil (NOT YET - DO IT after build L)

  }*/
  return local_list_filaments;
}


void conductor::send_connect_info ()
{
  (node1_->real_node())->
  add_int_node_ptr (list_internal_nodes.initialize());
  (node2_->real_node())->
  add_int_node_ptr (list_internal_nodes.get_last_element());
}



List<filament<double> >* conductor::assignFil(conductor* seg)
{
  int i,j,k;
  double x,y,z, delw, delh;
  double x0, x1, y0, y1, z0, z1;
  double hx, hy, hz;
  int temp, counter;
  int Hinc, Winc;
  double Hdiv, Wdiv;
  double fil_height, fil_width, seg_length;
  double h_from_edge, w_from_edge, min_height, min_width;
  filament<double> *filptr, *seg_filaments;
  int countfils;
  double wx, wy, wz, mag;  // direction of width
  Node *node0, *node1;
  int indices[4], row, col;
  double r_width, r_height;
  int num_fils;
  List <filament<double> > *local_list_filaments;
  local_list_filaments = new List<filament<double> >;
 
 
 
  Hinc = seg->fnh();
  Winc = seg->fnw();
  r_height = seg->frh();
  r_width = seg-> frw();
 
 
  double seg_widthdir[3];
  seg_widthdir [XX]= (seg->widthdir())->x();
  seg_widthdir [YY] = (seg->widthdir())->y();
  seg_widthdir [ZZ] = (seg->widthdir())->z();
 
  seg_length=sqrt( (seg->node(1)->x() - seg->node(2)->x())*
                   (seg->node(1)->x() - seg->node(2)->x())
                   + (seg->node(1)->y() - seg->node(2)->y())*
                   (seg->node(1)->y() - seg->node(2)->y())
                   + (seg->node(1)->z() - seg->node(2)->z())*
                   (seg->node(1)->z() - seg->node(2)->z()) );
 
 
 
  num_fils = Winc*Hinc;
  filament<double>* filament_list;
  filament_list = new filament<double> [num_fils];
  seg_filaments=(filament<double> *)filament_list;
 
 
  if (fabs(1.0 - r_height) < EPS)
    Hdiv = Hinc;
  else
  {
    temp = Hinc/2;
    Hdiv = 2*(1.0 - pow(r_height, (double)temp))/(1.0 - r_height);
    if (Hinc%2 == 1) Hdiv += pow(r_height, (double)(temp));
  }
 
  if (fabs(1.0 - r_width) < EPS)
    Wdiv = Winc;
  else
  {
    temp = Winc/2;
    Wdiv = 2*(1.0 - pow(r_width, (double)temp))/(1.0 - r_width);
    if (Winc%2 == 1) Wdiv += pow(r_width, (double)(temp));
  }
 
  node0 = seg->node(1);
  node1 = seg->node(2);
 
 
  // determine direction of width
  if (seg->widthdir_ != NULL)
  {
    wx = seg_widthdir[XX];
    wy = seg_widthdir[YY];
    wz = seg_widthdir[ZZ];
  }
  else
  {
    // default for width direction is in x-y plane perpendic to length
    // so do cross product with unit z
    wx = -(node1->y() - node0->y())*1.0;
    wy = (node1->x() - node0->x())*1.0;
    wz = 0;
    if ( fabs(wx/seg_length) < EPS && fabs(wy/seg_length) < EPS)
    {
      // if all of x-y is perpendic to length, then choose x direction
      wx = 1.0;
      wy = 0;
    }
    mag = sqrt(wx*wx + wy*wy + wz*wz);
    wx = wx/mag;
    wy = wy/mag;
    wz = wz/mag;
  }
 
  // height direction perpendicular to both
  hx = -wy*(node1->z() - node0->z()) + (node1->y() - node0->y())*wz;
  hy = -wz*(node1->x() - node0->x()) + (node1->z() - node0->z())*wx;
  hz = -wx*(node1->y() - node0->y()) + (node1->x() - node0->x())*wy;
  mag = sqrt(hx*hx + hy*hy + hz*hz);
  hx = hx/mag;
  hy = hy/mag;
  hz = hz/mag;
 
  filptr = seg_filaments;
  counter = 0;
 
  // this will fill the 'filament' array. It uses symmetry wrt z and y
  // it generates the four corner fils first, then the next one in...
  // 6/25/93 - added stuff to place fils in filament array so that adjacent
  //           fils are adjacent in the array.  This will make the meshes
  //           in M consist of fils that are near each other.
  h_from_edge = 0.0;
  min_height = 1.0/Hdiv;  // fil of smallest height
  for(i = 0; i < Hinc/2 || (Hinc%2 == 1 && i == Hinc/2); i++)
  {
 
    // height of the filaments for this row
    fil_height = (seg->height()/Hdiv)*pow(r_height, (double) i);
 
    if (i == 0)
      h_from_edge += min_height/2;
    else
      h_from_edge += min_height/2*pow(r_height,(double)(i-1))*(1+r_height);
 
    delh = (seg->height())*(0.5 - h_from_edge);
 
    if (delh < 0.0 && fabs(delh/seg->height()) > EPS)
      printf("uh oh, delh < 0. delh/height = %lg\n", delh/seg->height());
 
    w_from_edge = 0;
    min_width = 1.0/Wdiv;
 
    for(j = 0; j < Winc/2 || (Winc%2 == 1 && j == Winc/2); j++)
    {
      fil_width = (seg->width()/Wdiv)*pow(r_width, (double)j );
 
      if (j == 0)
        w_from_edge += min_width/2;
      else
        w_from_edge += min_width/2*pow(r_width,(double)(j-1))*(1+r_width);
 
      delw = (seg->width())*(0.5 - w_from_edge);
 
      if (delw < 0.0 && fabs(delw/seg->width()) > EPS)
        printf("uh oh, delw < 0. delw/width = %lg\n", delw/seg->width());
 
 
      countfils = 0;
      row = i;
      col = j;
      if (row%2 == 1)  col = (Winc - 1) - col;
      indices[countfils] = col + Winc*row;
      filptr = &(filament_list[indices[countfils]]);
 
      x0=node0->x() + hx*delh + wx*delw;
      x1=node1->x() + hx*delh + wx*delw;
      y0=node0->y() + hy*delh + wy*delw;
      y1=node1->y() + hy*delh + wy*delw;
      z0=node0->z() + hz*delh + wz*delw;
      z1=node1->z() + hz*delh + wz*delw;
      point3D<double> p1 (x0, y0, z0 );
      point3D<double> p2 (x1, y1, z1 );
      filptr= new filament<double> (p1, p2, fil_width, fil_height, seg->sigma());
      local_list_filaments->push_back(filptr);
      countfils++;
 
 
      // do symmetric element wrt y
      if(j != Winc/2)
      {
        row = i;
        col = (Winc - 1) - j;
        if (row%2 == 1)
          col = (Winc - 1) - col;
        indices[countfils] = col + Winc*row;
        filptr = &(filament_list[indices[countfils]]);
 
        x0 = node0->x() + hx*delh - wx*delw;
        x1 = node1->x() + hx*delh - wx*delw;
        y0 = node0->y() + hy*delh - wy*delw;
        y1 = node1->y() + hy*delh - wy*delw;
        z0 = node0->z() + hz*delh - wz*delw;
        z1 = node1->z() + hz*delh - wz*delw;
        point3D<double> p1 (x0, y0, z0 );
        point3D<double> p2 (x1, y1, z1 );
        filptr= new filament<double> (p1, p2, fil_width, fil_height, seg->sigma());
        local_list_filaments->push_back(filptr);
        countfils++;
      }
 
      // wrt z
      if(i != Hinc/2)
      {
        row = (Hinc - 1) - i;
        col = j;
        if (row%2 == 1)
          col = (Winc - 1) - col;
        indices[countfils] = col + Winc*row;
        filptr = &(filament_list[indices[countfils]]);
        x0 = node0->x() - hx*delh + wx*delw;
        x1 = node1->x() - hx*delh + wx*delw;
        y0 = node0->y() - hy*delh + wy*delw;
        y1 = node1->y() - hy*delh + wy*delw;
        z0 = node0->z() - hz*delh + wz*delw;
        z1 = node1->z() - hz*delh + wz*delw;
        point3D<double> p1 (x0, y0, z0 );
        point3D<double> p2 (x1, y1, z1 );
        filptr= new filament<double> (p1, p2, fil_width, fil_height, seg->sigma());
        local_list_filaments->push_back(filptr);
        filptr++;
        countfils++;
      }
 
      // wrt z and y
      if( i != Hinc/2 && j != Winc/2)
      {
        row = (Hinc - 1) - i;
        col = (Winc - 1) - j;
        if (row%2 == 1)
          col = (Winc - 1) - col;
        indices[countfils] = col + Winc*row;
        //cout<<"fourth "<<indices[countfils]<<endl;
        filptr = &(seg_filaments[indices[countfils]]);
 
        x0 = node0->x() - hx*delh - wx*delw;
        x1 = node1->x() - hx*delh - wx*delw;
        y0 = node0->y() - hy*delh - wy*delw;
        y1 = node1->y() - hy*delh - wy*delw;
        z0 = node0->z() - hz*delh - wz*delw;
        z1 = node1->z() - hz*delh - wz*delw;
        point3D<double> p1 (x0, y0, z0 );
        point3D<double> p2 (x1, y1, z1 );
        filptr= new filament<double> (p1, p2, fil_width, fil_height, seg->sigma());
        local_list_filaments->push_back(filptr);
        countfils++;
      }
    }
  }
 
  i = 0;
  while(i < Hinc*Winc)
    i++;
 
  if (i != Hinc*Winc)
  {
    fprintf(stderr, "Hey, not all filaments created in assignfil()! \n");
    exit(1);
  }
 
  return local_list_filaments;
}



