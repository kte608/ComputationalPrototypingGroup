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

#include "fastsub.h"
#include "structure.h"
#include "defaults.h"
#include "initialize_discretization.h"
#include "layerenv.h"
#include "generate_dipoles.h"
#include "sysformation.h"
#include "extractz.h"

#include <iostream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>


using namespace fastsub;
using namespace std;
using namespace pfft;

int main(int argc, char* argv[])
{
  char* inputFileName;
  char* outputFileName;
  FILE *fp, *fp2_r, *fp2_i, *fp3, *fp_layer;
  SimulationType simuType=EMQS;
  SolverMethod solMethod=CIRCUIT_SOLVE_MESH;
  FrequencySample freqSam;
  Defaults default_vals;
  Structure structure;
  LayerEnv layers;

  int numLayers=0;
  vector<std::complex<double> > Im, Ib;
  complex<double> adm, Rsource;
  bool output_flag, debug_flag, discretize_flag;
  
  
  //parse command line
  cmd_line_parsor(argc, argv, &inputFileName, &outputFileName, simuType,
                  solMethod, freqSam, output_flag, debug_flag, discretize_flag);
  //simuType=FULLWAVE;
  simuType=EMQS;
  //prep for discretization
  
  
  Initialize_discretization init (inputFileName, structure, default_vals, layers);
  //discretize
  structure.discretize(solMethod);
  
  
  numLayers=layers.num_layers();
  fp_layer = fopen ("layer_info.mat", "w");
  fprintf(fp_layer, "%e  %e %e \n", layers.get_layer_from_index(0)->epir(),
	  layers.get_layer_from_index(0)->sigma(), layers.get_layer_from_index(0)->startHeight());
  if(numLayers==2)
    fprintf(fp_layer, "%e  %e %e \n", layers.get_layer_from_index(1)->epir(),
	    layers.get_layer_from_index(1)->sigma(), layers.get_layer_from_index(1)->startHeight());
  else if (numLayers>2)
    {
      cout<<"only two layers can be handled right now."<<endl;
      exit(1);
    }
  fclose(fp_layer);
  
  //initiate sysFormation
  sysFormation sys(structure, layers, simuType, debug_flag);
  //exit(1);
  if(discretize_flag==true)
    exit(1);
  
      FILE* output_file;
      output_file=fopen("output.mat", "w");
      fclose(output_file);

  //perform a frequency sweep
  for (int i=0; i< freqSam.numSam(); i++)
    {
      cout<<"Solving at frequency "<< freqSam.point(i)<<endl;
      layers.compute_freq_related_info(freqSam.point(i));
      //generate dipole information
      Dipoles dipoles (layers);
      //sys.Iterate();
      sys.Solver(dipoles, freqSam.point(i),simuType,FORMULATION_C,pFFT);
      
      //extract solution
      //   extractZ extract(Ib,Im,adm, Rsource,freqSam.point(i));
      //   cout<<"adm "<<adm<<endl;
      
      
      //output solution
      //    if (output_flag==true)
      //    {
      //      //fprintf(fp, "%e ", freqSam.point(i));
      //      for (int j=0; j<Im.size(); j++)
      //      {
      //        if (imag(Im[j])<=0)
      //          fprintf(fp, "%e+%ej ", real(Im[j]), imag(Im[j]));
      //        else
      //          fprintf(fp, "%e+%ej ", real(Im[j]), imag(Im[j]));
      //      }
      //
      //      //fprintf(fp2, "%e ", freqSam.point(i));
      //      for (int j=0; j<Ib.size(); j++)
      //      {
      //        fprintf(fp2_r, "%e ", real(Ib[j]));
//        fprintf(fp2_i, "%e ", imag(Ib[j]));
//      }
//
//      fprintf(fp3, "F: %e \n", freqSam.point(i));
//      fprintf(fp3, "Z: %e+%ej \n", real(extract.Z()), imag(extract.Z()));
//      fprintf(fp3, "Y: %e+%ej \n", real(extract.Y()), imag(extract.Y()));
//      fprintf(fp3, "R: %e \n", extract.R());
//      fprintf(fp3, "X: %e \n", extract.X());
//      fprintf(fp3, "C: %e \n", extract.C());
//      fprintf(fp3, "L: %e \n", extract.L());
//	  fprintf(fp3, "Q: %e \n", extract.Q());
//      fprintf(fp3, "\n");
//
//      fclose (fp);
//      fclose (fp2_r);
//      fclose (fp2_i);
//      fclose (fp3);
//
//      printf ("\n");
//      printf ("Im vector will be written to Im.mat in local directory \n");
//      printf ("Ib vector will be written to Ib_real.mat and Ib_imag.mat in local directory \n");
//      printf ("Impedance info will be written to  %s\n", outputFileName);
//    }
    }
}





//parse the command line and the switches
void fastsub::cmd_line_parsor(int argc, char* argv[], char** fileName, char** outputFileName, SimulationType& simuType, SolverMethod& solMethod, FrequencySample& freqSam, bool& output_flag, bool& debug_flag, bool& discretize_flag)
{
  int option;
  double freqBegin=freqSam.firstPoint();
  double freqEnd=freqSam.lastPoint();
  bool set_freqBegin=false, set_freqEnd=false;
  int numSam=freqSam.numSam();
  SampleMethod samMethod=freqSam.samMethod();
  output_flag=false;
  debug_flag=false;
  discretize_flag=false;

  while ((option=getopt(argc, argv, "i:s:b:e:n:f:a:o:ci:hi:di:"))!=EOF)
  {
    switch (option)
    {
    case INPUT_FILE_NAME:
      *fileName=optarg;
      break;
    case SIMULATION_TYPE:
      simuType=static_cast<SimulationType>(atoi(optarg));
      break;
    case FREQUENCY_BEGIN:
      sscanf(optarg, "%lf", &freqBegin);
      set_freqBegin=true;
      break;
    case FREQUENCY_END:
      sscanf(optarg, "%lf", &freqEnd);
      set_freqEnd=true;
      break;
    case FREQUENCY_NUM_STEPS:
      numSam=atoi(optarg);
      break;
    case OUTPUT_FILE_NAME:
      output_flag=true;
      *outputFileName=optarg;
      break;
    case DEBUG_MODE:
      debug_flag=true;
      break;
    case DISCRETIZE:
      discretize_flag=true;
      break;
    case HELP:
      printUsage();
      exit(0);
      break;
    default:
      printUsage();
      exit(0);
    }
  }

  if ((set_freqBegin) && (!set_freqEnd))
    freqSam=FrequencySample(freqBegin, freqBegin, numSam, samMethod);
  else if   (!(set_freqBegin) && (set_freqEnd))
    freqSam=FrequencySample(freqEnd, freqEnd, numSam, samMethod);
  else
    freqSam=FrequencySample(freqBegin, freqEnd, numSam, samMethod);
}

void fastsub::printUsage()
{
  cout << endl
  << "\t Usage:  fastsub <options>" << endl
  << "\t -i <input structure file name>" << endl
  << "\t -b <begin of the frequncy range>"  << endl
  //<< "\t    default is " << freqSam.firstPoint() << "Hz" << endl
  //<< "\t -e <end of the frequncy range>" << endl
  //<< "\t    default is " << freqSam.lastPoint() << "Hz"  << endl
  //<< "\t -n <number of sampling points in the frequncy range>" << endl
  //<< "\t    default is " << freqSam.numSam() << endl
  //<< "\t -m <frequency sampling method>, 0:LOGRITHMIC, 1:LINEAR" << endl
  //<< "\t    default is " << freqSamMethodName[freqSam.samMethod()] << endl
  << "\t -o <output mesh file name>" << endl
  << "\t    Impedance info will be generated into the file in addition to Im.mat and Ib.mat" << endl
  << "\t    that contains the Ib and Im vectors, respectively." << endl
  << "\t -d debugging mode; outputs R, L and P matrices" <<endl
  << "\t -c discretization only "<<endl
  << "\t -s <simulation mode: 0=quasistatic; 1=fullwave>" <<endl
  << "\t -h Print this message and Exit" << endl <<endl
  << "\t  ex: ./fastsub -i tline.inp -b 1e6"<<endl;
}



