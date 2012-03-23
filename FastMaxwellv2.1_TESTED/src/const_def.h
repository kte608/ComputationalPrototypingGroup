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
#ifndef CONST_DEF_H
#define CONST_DEF_H

#include <math.h>
#include <stdio.h>
#include <ctype.h>

/* For node and segment types */
#define NORMAL 0
#define GPTYPE 1
#define PSEUDO 2
#define EXTERNTYPE 3

//#define CIRCUIT_SOLVE_MESH 1
//#define CIRCUIT_SOLVE_NODAL 2

#define MAXLINE 1000
#define XX 0
#define YY 1
#define ZZ 2
#define EPS 1e-13

#define IN_NORMAL 0
#define IN_RESIST 1
//#ifndef PI
   #define PI 3.14159265358979323846
//#endif

#define EPI0  1e-9/(36*PI)
//#define EPI0 8.854187818E-12
#define MU0   4*PI*1e-7

#define MAXDEG 10

#ifndef MAX
#define MAX(A,B)  ( (A) > (B) ? (A) : (B) )
#endif

#define IMAG complex<double>(0.,1.)
#define  ZERO complex<double>(0.,0.)
#define  ONE complex<double>(1.,0.)
#define  TWO complex<double>(2.,0.)

/*#ifndef MIN
#define MIN(A,B)  ( (A) < (B) ? (A) : (B) )
#endif*/

void tolowercase (char *);

typedef char *sparse;

#define IN_NORMAL 0
#define IN_RESIST 1


#endif






