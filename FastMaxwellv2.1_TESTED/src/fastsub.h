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

#ifndef FASTSUB_H
#define FASTSUB_H

#include "frequencysample.h"
#include "layerenv.h"
//#include "sysformation.h"

namespace fastsub {

enum cmdLineOptions{
	INPUT_FILE_NAME			=	'i',
	SIMULATION_TYPE 		=	's',
	FREQUENCY_BEGIN 		=   'b',
	FREQUENCY_END  			= 	'e',
	FREQUENCY_NUM_STEPS		=	'n',
	FREQUENCY_SAM_METHOD	=	'f',
	BASIS_FUNC              =   'a',
	OUTPUT_FILE_NAME        =   'o',
	DEBUG_MODE              =   'd',
	HELP					=   'h',
	DISCRETIZE              =   'c'
};

enum SimulationType{EMQS, FULLWAVE};
enum SolverMethod{CIRCUIT_SOLVE_MESH};
enum FormulationType{FORMULATION_A, FORMULATION_C};
enum SolverType {DIRECT_SOLVER, pFFT};



void cmd_line_parsor(int argc, char* argv[], char **, char **,
				   SimulationType&, SolverMethod&, FrequencySample&, bool&, bool&, bool&);
					
void printUsage(void);
				   
}
#endif
