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
#ifndef INITIALIZE_DISCRETIZATION_H
#define INITIALIZE_DISCRETIZATION_H

#include "defaults.h"
#include "structure.h"
#include "layerenv.h"

using namespace std;

class Initialize_discretization
{

public:
  Initialize_discretization (char*, Structure&, Defaults&, LayerEnv&);

private:
 	
  char* title_;
  int keep_;	
  double units_;
  
  int read_input_file(FILE*, Structure&, Defaults&, LayerEnv&);
  char* getaline(FILE*);
  char* getoneline(FILE*);
  char* plusline(FILE*);
  int notblankline(char *);
  void savealine();
  int dodot (char *,  Structure&, Defaults&);
  int changeunits (char *);
  int dodefaults (char *, Defaults&);
  int equivnodes (char *, Structure&);
  int addexternal (char *, Structure&);
  int addnode(char *, int,  Structure&, Defaults&);
  int addseg (char *, int,  Structure&, Defaults& );
  int addresistance (char *, Structure&);
  int addlayer(char*, LayerEnv&, Defaults&);
  int nothing (char *);
  

};

#endif
