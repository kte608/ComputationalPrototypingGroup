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
#ifndef CONDUCTOR_H
#define CONDUCTOR_H

/**
@author Xin Hu
*/
#include "list.h"
#include "point3D.h"
#include "filament.h"
#include "panel.h"
#include "internal_node.h"
//#include "matrix.h"
//#include "incidence.h"
#include "fastsub.h"
#include "../pfft/spRowMat.h"

class Node;
class conductor
{
public:
  conductor (char *_name, Node *_node1, Node *_node2,
             double _width, double _height,
             point3D<double> *_widthdir, double _sigma,
             int _nl, int _fnw, int _fnh, int _pnw, int _pnh,
             double _frw, double _frh, double _pfracl, 
	     double _pfracw, double _pfrach, int _type,
             int _nep);
  ~conductor(void);

/*  void discretize(fastsub::SolverMethod, 
				  List<panel<double> >**, 
				  List<filament<double> >**,
                  List<SEGMENT >**, 
				  List <Internal_Node> **,
                  Incidence_Matrix<filament<double> > **,
                  Incidence_Matrix<panel<double> > **,
                  Incidence_Matrix<filament<double> > **,
                  Incidence_Matrix<panel<double> > **);*/

  void discretize(fastsub::SolverMethod, 
		List<panel<double> >**, 
		List<filament<double> >**,
                  //List<SEGMENT >**, 
				  List <Internal_Node> **,
                  pfft::SpRowMat<double> *,
                  pfft::SpRowMat<double> *,
                  pfft::SpRowMat<double> *,
                  pfft::SpRowMat<double> *,
		  int, int, int *, int *, int *);

  void send_connect_info();

  point3D<double> vul () { return vul_; };
  point3D<double> vuw () { return vuw_; };
  point3D<double> vuh () { return vuh_; };
  char* get_name() {return name_;};
  Node *node (int number);
  void set_nep1(int);
  void set_nep2(int);
  void set_stretch(double, int);
  const point3D<double>* widthdir(void) {return widthdir_;}
  const double width(void) const {return width_;};
  const double height(void) const {return height_;};
  const double length(void) const {return length_;};
  const double sigma(void) const {return sigma_;};
  const int nl(void) const {return nl_;};
  const int fnw(void) const {return fnw_;};
  const int fnh(void) const {return fnh_;};
  const int pnw(void) const {return pnw_;};
  const int pnh(void) const {return pnh_;};
  const double frw(void) const {return frw_;};
  const double frh(void) const {return frh_;};
  const double pfracl(void) const {return pfracl_;};
  const double pfracw(void) const {return pfracw_;};
  const double pfrach(void) const {return pfrach_;};
  const int nep1(void) const {return nep1_;};
  const int nep2(void) const {return nep2_;};


private:
  double width_, height_, length_, sigma_;
  double frw_, frh_, pfracl_, pfracw_, pfrach_;
  int nl_, fnw_, fnh_, pnw_, pnh_;
  int nep1_, nep2_;
  double stretch1_, stretch2_;
  int have_edge1_, have_edge2_;
  int have_2nd_edge1_, have_2nd_edge2_;
  char* name_;
  Node *node1_, *node2_, *node_;
  point3D<double> *widthdir_;
  int type_;
  point3D<double> vul_, vuw_, vuh_;
  int num_filaments_, num_panels_;
  int global_filaments_index_, global_panels_index_;
  int *global_mesh_line_p;
  int *global_mesh_line_f;
  int *global_mesh_line_s;
  List <Internal_Node> list_internal_nodes;


/*  int discretize_panels(fastsub::SolverMethod, 
						List<panel<double> >**, 
						Incidence_Matrix<panel<double> > **);
  int discretize_filaments(fastsub::SolverMethod,  
						   List<filament<double> >**, 
						   List <SEGMENT> **,
                           Incidence_Matrix<filament<double> > **,
                           Incidence_Matrix<filament<double> > **,
                           Incidence_Matrix<panel<double> > ** );
			   */
  int discretize_panels(fastsub::SolverMethod, 
					List<panel<double> >**, 
					pfft::SpRowMat<double> *);
  int discretize_filaments(fastsub::SolverMethod,  
					List<filament<double> >**, 
						   //List <SEGMENT> **,
                           pfft::SpRowMat<double> *,
                           pfft::SpRowMat<double> *,
                           pfft::SpRowMat<double> *);
  
  
 // List <filament<double> > * call_assign_fil (conductor *, List <SEGMENT> *);
  List <filament<double> > * call_assign_fil (conductor *);
  List <filament<double> > * assignFil(conductor *);

};

// extern "C"
// {
//   // functions of henry.c
//   charge *assignFil(SEGMENT *, int *, charge *);
// }

#endif
