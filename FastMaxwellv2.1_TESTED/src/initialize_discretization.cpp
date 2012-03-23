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
#include "initialize_discretization.h"
#include "const_def.h"
#include "list.h"
#include "point3D.h"
#include "pseudo_node.h"
#include "ext_node_pair.h"
#include "resistance.h"
#include "layer.h"

#include <string>
#include <stdlib.h>
#include <math.h>
#include <fstream>


Initialize_discretization::Initialize_discretization(char* fileName,
    Structure& str, Defaults& defaults, LayerEnv& layers)
{
  keep_=0;
  units_=1.0;


  int ret;
  FILE* fp;
  cout<<"--> Stage:Reading input"<<endl;

  fp=fopen(fileName, "r");
  if (fp==NULL)
  {
    cout<<"couldn't open "<< fileName<<endl;
    exit(1);
  }

  cout<<"Reading from file "<< fileName<<endl;
  ret=read_input_file(fp, str, defaults, layers);
  if (ret)
  {
    printf ("\nError in input file. Program terminated without "
            "performing any analysis.\n");
    exit (-1);
  }

  fclose(fp);
  cout<<"Input processed successfully."<<endl;

}


int Initialize_discretization::read_input_file(FILE* fp, Structure& str, 
											   Defaults& defaults, LayerEnv& layers)
{
  char* line;
  int end=0, error=0;
  Layer* default_layer;

  line=getaline(fp);
  title_=new char [strlen(line)+1];
  strcpy(title_, line);

  while(!end)
  {
    line=getaline(fp);
    tolowercase(line);
    switch (line[0])
    {
    case '*':
      break;
    case '.':
      //cout<<"dot"<<endl;
      end = dodot(line + 1, str, defaults);
      break;
    case 'n':
      //cout<<"addnode"<<endl;
      end = addnode(line, NORMAL,str, defaults);
      break;
    case 'e':
      //cout<<"addseg"<<endl;
      end = addseg(line, NORMAL, str, defaults);
      break;
    case 'l':
      //cout<<"addlayer"<<endl;
      end = addlayer(line, layers, defaults);
      break;
    case 'r':
      //cout<<"addresistance"<<endl;
      end = addresistance(line, str);
      break;
    default:
      //cout<<"do nothing"<<endl;
      end = nothing(line);
      break;
    }
    if (end == 1)
    {
      error = 1;
      end = 0;
    }
  }

  if ((defaults.isepir()==1) && (defaults.isLayerSigma()==1))
    {
      default_layer=new Layer("top", defaults.epir(),
			      defaults.layerSigma(), HUGE_VAL);
      layers.add_defined_layer(default_layer);
    }
  else
    {
      cout<<"please enter values for epir and sigma"<<endl;
      end=1;
    } 
  
  if (end == 2) end = error;
  return end;
}


char* Initialize_discretization::getaline(FILE* fp)
{
  static char *all_lines = NULL;
  static int length = 0;
  char *line;
  int newlength;

  if (length == 0)
  {
    length = MAXLINE;
    all_lines = new char [length];
  }
  line = getoneline (fp);
  if (line == NULL)
  {
    fprintf(stderr, "Unexpected end of file\n");
    exit(1);
  }
  strcpy(all_lines, line);

  // concatenate any lines beginning with a '+' to all_lines
  while( (line = plusline(fp)) != NULL)
  {
    if ((newlength = strlen(all_lines) + strlen(line) + 1) > length)
    {
      if ( (all_lines = (char *) realloc(all_lines,
                                         MAX(newlength,length+MAXLINE)))
           == NULL )
      {
        fprintf(stderr,"couldn't get more space for a line. Needed %d chars\n",
                newlength);
        exit(1);
      }
      else
        length = MAX(newlength,length+MAXLINE);
    }
    strcat(all_lines, line);
  }
  return all_lines;
}


char *Initialize_discretization::getoneline (FILE *fp)
{
  static char line[MAXLINE] = { '\0' };
  char *retchar;

  if (keep_)
  {
    keep_ = 0;
    return line;

  }
  else
    do
    {
      retchar = fgets(line, MAXLINE, fp);
    }
    while(retchar != NULL && !notblankline(line));

  if (retchar != NULL && strlen(line) == MAXLINE - 1)
    fprintf(stderr,"Warning: line may be too long:\n%s\n",line);

  if (retchar == NULL)
    return NULL;
  else
    return line;
}


char *Initialize_discretization::plusline (FILE *fp)
{
  char *tmpline;

  tmpline = getoneline(fp);

  while(tmpline != NULL && tmpline[0] == '*')
    tmpline = getoneline(fp);

  if (tmpline == NULL)
    return NULL;
  else if (tmpline[0] == '+')
  {
    tmpline++;
    return tmpline;
  }
  else
  {
    savealine();
    return NULL;
  }
}


int Initialize_discretization::notblankline(char *string)
{
  while( *string!='\0' && isspace(*string))
    string++;

  if (*string == '\0') return 0;
  else return 1;
}

void Initialize_discretization::savealine ()
{
  if (keep_ != 0)
  {
    printf("already have one line stored\n");
    exit(1);
  }
  else
    keep_ = 1;
}

int Initialize_discretization::dodot (char *line, Structure& str, Defaults& defaults)
{
  int end;
  if (strncmp("uni",line, 3) == 0)
    end = changeunits(line);
  else if (strncmp("ext",line, 3) == 0)
    end = addexternal(line, str);
  else if (strncmp("fre",line, 3) == 0)
    cout<<"Need to take care of doing frequencies"<<endl;
  //end = choosefreqs(line, str);
  else if (strncmp("equ",line, 3) == 0)
    end = equivnodes(line, str);
  else if (strncmp("def",line, 3) == 0)
    end = dodefaults(line, defaults);
  else if (strncmp("end",line, 3) == 0)
    end = 2;
  else
  {
    printf("What is the following line?\n.%s\n",line);
    return 1;
  }
  return end;
}


int Initialize_discretization::changeunits (char *line)
{
  char unitname[20];
  double previous_units;

  previous_units = units_;

  if (sscanf(line, "%*s %s", unitname) != 1)
  {
    printf("couldn't read units in line:\n.%s\n", line);
    return 1;
  }
  else if (strncmp("mil",unitname, 3) == 0)
    units_ = 2.54e-5;                     // mils to meters
  else if (strncmp("in",unitname, 2) == 0)
    units_ = 0.0254;                      // inches to meters
  else if (strncmp("um",unitname, 2) == 0)
    units_ = 1e-6;                        // microns to meters
  else if (strncmp("mm",unitname, 2) == 0)
    units_ = 1e-3;                      // millimeters to meters
  else if (strncmp("cm",unitname, 2) == 0)
    units_ = 1e-2;                      // centimeters to m
  else if (strncmp("m",unitname, 1) == 0)
    units_ = 1.0;                      // meters to meters
  else if (strncmp("km",unitname, 1) == 0)
    units_ = 1e3;                      // kilometers to meters
  else
  {
    printf("unrecognizable unit name \n");
    return (1);
  }

  printf("all lengths multiplied by %lg to convert to meters\n", units_);
  return 0;
}



int Initialize_discretization::equivnodes (char *line, Structure& structure)
{

  List <Node> nlist;
  List <Pseudo_Node> pnlist;
  Node *nl, *node, *realnode;
  Pseudo_Node *pn;

  int skip;
  char name1[80];


  if (sscanf(line,"%*s%n",&skip) != 0)
  {
    printf("Hey, no fair\n");
    return 1;
  }
  line += skip;


  while(notblankline(line))
  {
    if (sscanf(line,"%s%n",name1,&skip) != 1)
    {
      printf("equivnodes: sscanf error on equiv line\n");
      return 1;
    }
    line += skip;

    node = structure.get_node_from_name(name1);
    if (node == NULL)
    {
      //node doesn't exist, add it to pseudo node list
      if (name1[0] != 'n')
        fprintf(stderr, "create_pn: Invalid node name %s. "
                "First letter should be 'n'\n", name1);
      pn = new Pseudo_Node (name1);
      pnlist.push_back (pn);
    }
    else
      nlist.push_back(node);
  }

  if (nlist.size() == 0)
  {
    fprintf(stderr, "equivnodes: No already defined nodes "
            "in .equiv statement.\n");
    exit(1);
  }

  // the node which all the others are to become 'equiv'ed to
  realnode = nlist.initialize();
  if(realnode==NULL) cout<<"NULL"<<endl;
  realnode = realnode->real_node();


  // assign the pseudo nodes to realnode
  for (pn = pnlist.initialize(); pn != NULL; pn = pnlist.iterate())
    pn->assign_defined_node(realnode);


  // make the others also equivalent
  for (nlist.initialize(), nl=nlist.iterate(); nl!= NULL; nl=nlist.iterate())
    nl->assign_real_node (realnode);


  for (pn = pnlist.initialize(); pn!=NULL; pn=pnlist.iterate())
    structure.add_pseudo_node (pn);


  return 0;
}


int Initialize_discretization::addexternal (char *line, Structure& structure)
{
  int skip, i;
  Node *node[2];
  static char name[80];
  static char names[2][80];
  static char portname[80];
  int err=0;
  Ext_Node_Pair *ext;

  if(sscanf(line, "%*s%n", &skip) != 0)
  {
    printf("Hey, no fair\n");
    return 1;
  }
  line += skip;

  // read the nodes which will be external
  for(i = 0; i < 2; i++)
  {
    if(sscanf(line, "%s%n",name,&skip) != 1)
    {
      printf("No node name to read for .external in line:\n%s", line);
      return 1;
    }
    line += skip;

    strcpy(names[i], name);

    // search for the named node in the list
    node[i] = structure.get_node_from_name(name);

    if ((node[i] == NULL) && (strcmp (name, "infinity")))
    {
      printf("addexternal: No node read in yet named %s\n",name);
      return 1;
    }

    if (strcmp (name, "infinity"))
      node[i] = node[i]->real_node();
  }

  if(sscanf(line, "%s%n",portname,&skip) != 1)
  {
    portname[0] = '\0';
  }
  else
  {
    line += skip;
  }
//
//  if ((ext = structure.get_external_from_portname(portname)) != NULL)
//  {
//    fprintf(stderr, "Error: Cannot name port between %s and %s as %s.\n"
//            "The name is already used for port between %s and %s\n",
//            names[0],names[1],portname,ext->node_name(1),ext->node_name(2));
//    err = 1;
//  }
//
  if (node[0] == node[1])
  {
    fprintf(stderr, "Error in .external: Nodes %s and %s are the "
            "same node. (.equiv'ed maybe?)\n", names[0], names[1]);
    err = 1;
  }
  structure.add_ext_node_pair(new Ext_Node_Pair (node[0], node [1], portname));
  return err;
}

int Initialize_discretization::dodefaults (char *line, Defaults& default_vals)
{
  int skip;
  double dumb;
  int dumbi;
  static char string[80];


  if (sscanf(line,"%*s%n",&skip) != 0)
  {
    printf("Hey, no fair\n");
    return 1;
  }
  line += skip;

  while(notblankline(line))
  {
    if (sscanf(line," x = %lf%n",&dumb,&skip) == 1)
    {
      default_vals.set_isx (1);
      default_vals.set_x (units_*dumb);
      line += skip;
    }
    else if (sscanf(line," y = %lf%n",&dumb,&skip) == 1)
    {
      default_vals.set_isy (1);
      default_vals.set_y( units_*dumb);
      line += skip;
    }
    else if (sscanf(line," z = %lf%n",&dumb,&skip) == 1)
    {
      default_vals.set_isz (1);
      default_vals.set_z (units_*dumb);
      line += skip;
    }
    else if (sscanf(line," h = %lf%n",&dumb,&skip) == 1)
    {
      default_vals.set_ish (1);
      default_vals.set_h ( units_*dumb);
      line += skip;
    }
    else if (sscanf(line," w = %lf%n",&dumb,&skip) == 1)
    {
      default_vals.set_isw (1);
      default_vals.set_w (units_*dumb);
      line += skip;
    }
    else if (sscanf(line," rho = %lf%n",&dumb,&skip) == 1)
    {
      default_vals.set_issigma (1);
      default_vals.set_sigma ( (1.0/dumb)/units_);
      line += skip;
    }
    else if (sscanf(line," sigma = %lf%n",&dumb,&skip) == 1)
    {
      default_vals.set_issigma (1);
      default_vals.set_sigma ( dumb/units_);
      line += skip;
    }
    else if (sscanf(line," nl = %d%n",&dumbi,&skip) == 1)
    {
      default_vals.set_isnl (1);
      default_vals.set_nl ( dumbi);
      line += skip;
    }
    else if (sscanf(line," fnh = %d%n",&dumbi,&skip) == 1)
    {
      default_vals.set_isfnh ( 1);
      default_vals.set_fnh (dumbi);
      line += skip;
    }
    else if (sscanf(line," fnw = %d%n",&dumbi,&skip) == 1)
    {
      default_vals.set_isfnw ( 1);
      default_vals.set_fnw ( dumbi);
      line += skip;
    }
    else if (sscanf(line," pnh = %d%n",&dumbi,&skip) == 1)
    {
      default_vals.set_ispnh ( 1);
      default_vals.set_pnh ( dumbi);
      line += skip;
    }
    else if (sscanf(line," pnw = %d%n",&dumbi,&skip) == 1)
    {
      default_vals.set_ispnw ( 1);
      default_vals.set_pnw ( dumbi );
      line += skip;
    }
    else if (sscanf(line," frw = %lf%n",&dumb,&skip) == 1)
    {
      default_vals.set_isfrw ( 1 );
      default_vals.set_frw (dumb );
      line += skip;
    }
    else if (sscanf(line," frh = %lf%n",&dumb,&skip) == 1)
    {
      default_vals.set_isfrh ( 1 );
      default_vals.set_frh ( dumb );
      line += skip;
    }
    else if (sscanf(line," pfracw = %lf%n",&dumb,&skip) == 1)
    {
      default_vals.set_ispfracw ( 1 );
      default_vals.set_pfracw ( dumb );
      line += skip;
    }
    else if (sscanf(line," pfrach = %lf%n",&dumb,&skip) == 1)
    {
      default_vals.set_ispfrach ( 1 );
      default_vals.set_pfrach ( dumb );
      line += skip;
    }
    else if (sscanf(line," pfracl = %lf%n",&dumb,&skip) == 1)
    {
      default_vals.set_ispfracl ( 1 );
      default_vals.set_pfracl ( dumb );
      line += skip;
    }
    else if ((sscanf(line," %s%n",string, &skip) == 1) &&
             (!(strcmp(string,"nep"))) )
    {
      default_vals.set_isnep ( 1 );
      default_vals.set_nep ( 1 );
      line += skip;
    }
    else if (sscanf(line, " layer_epir= %lf%n",&dumb,&skip) == 1)
    {
      default_vals.set_isepir(1);
      default_vals.set_epir(dumb);
      line+=skip;
    }
    else if (sscanf(line, " layer_sigma= %lf%n",&dumb,&skip) == 1)
    {
      default_vals.set_isLayerSigma(1);
      default_vals.set_layerSigma(dumb);
      line+=skip;
    }
    else
    {
      printf("don't know what piece this is in default:%s\n",line);
      return 1;
    }
  }  // end while(notblank...

  return 0;
}



int Initialize_discretization::addnode (char *line, int type, Structure& structure,
                                        Defaults& default_vals)
{
  double dumb;
  int skip;
  int xread = 0, yread = 0, zread = 0;
  static char name[80];
  Node* node;
  double nodex, nodey, nodez;

  // read name
  sscanf(line,"%s%n",name,&skip);
  line += skip;

  while(notblankline(line))
  {
    if (sscanf(line," x = %lf%n",&dumb,&skip) == 1)
    {
      nodex = units_*dumb;
      line += skip;
      xread = 1;
    }
    else if (sscanf(line," y = %lf%n",&dumb,&skip) == 1)
    {
      nodey = units_*dumb;
      line += skip;
      yread = 1;
    }
    else if (sscanf(line," z = %lf%n",&dumb,&skip) == 1)
    {
      nodez = units_*dumb;
      line += skip;
      zread = 1;
    }
    else
    {
      printf("don't know what piece this is in node %s: %s\n",name,line);
      return (1);
    }
  } // while not blank

  if (zread == 0)
  {
    if (default_vals.isz() == 1)
      nodez = default_vals.z();
    else
    {
      printf("Error: Z not given for node %s - discarded\n",name);
      return 1;
    }
  }
  if (yread == 0)
  {
    if (default_vals.isy() == 1)
      nodey = default_vals.y();
    else
    {
      printf("Error: Y not given for node %s - discarded\n",name);
      return 1;
    }
  }
  if (xread == 0)
  {
    if (default_vals.isx() == 1)
      nodex = default_vals.x();
    else
    {
      printf("Error: X not given for node %s - discarded\n",name);
      return 1;
    }
  }

  // make sure node isn't already defined
  if (structure.get_node_from_name (name) != NULL)
  {
    fprintf(stderr, "Node %s already defined. Cannot redefine it\n", name);
    exit(1);
  }

  node = new Node (name, nodex, nodey, nodez, type);

  structure.add_defined_node (node);

  return 0;
}


int Initialize_discretization::addseg (char *line, int type, Structure& structure,
                                       Defaults& default_vals)
{
  conductor* new_seg;
  double dumb;
  int skip, j, dumbi;
  int hread, wread, fnhread, fnwread, pnhread, pnwread,
  nlread, sigmaread, wxread,wyread, wzread, frwread,
  frhread, pfraclread, pfracwread, pfrachread, nepread;
  static char name[80], n_name[80], string[80];
  Node *anode;

  double *widthdir_tmp=NULL; // if width is not || to x-y plane and perpendicular
  // to the length, then this is 3 element vector in
  // in the direction of width
  point3D<double>* widthdir;
  double width, height;          //width and height to cross section
  int fnw, fnh;          // number of filament divisions in each dir (w and h)
  int pnw, pnh;          // number of panel divisions in each dir (w and h)
  int nl;                // number of sections
  Node *node[2];                 // nodes at the ends
  double sigma;                  // conductivity
  double frw, frh;               // for assignFil(), ratio of adjacent fils
  double pfracl, pfracw, pfrach; // fraction (of total size) of edge panels
  int nep;

  hread = wread = fnhread= fnwread = pnhread = pnwread=
          nlread = sigmaread = wxread= wyread=
          wzread = frwread=frhread = pfraclread = pfracwread=
                          pfrachread = nepread = 0 ;

  // read name
  sscanf(line,"e%s%n",name,&skip);
  line += skip;

  // read nodenames
  for (j = 0; j < 2; j++)
  {
    sscanf(line,"%s%n",n_name,&skip);
    line += skip;

    anode = structure.get_node_from_name(n_name);

    if (anode == NULL)
    {
      printf("Error: No node read in yet named %s\n",n_name);
      return 1;
    }
    else
    {
      node[j] = anode;
    }
  }

  while(notblankline(line))
  {
    if (sscanf(line," h = %lf%n",&dumb,&skip) == 1)
    {
      height = units_*dumb;
      line += skip;
      hread = 1;
    }
    else if (sscanf(line," w = %lf%n",&dumb,&skip) == 1)
    {
      width = units_*dumb;
      line += skip;
      wread = 1;
    }
    else if (sscanf(line," sigma = %lf%n",&dumb,&skip) == 1)
    {
      sigma = dumb/units_;
      line += skip;
      sigmaread = 1;
    }
    else if (sscanf(line," rho = %lf%n",&dumb,&skip) == 1)
    {
      sigma = (1.0/dumb)/units_;
      line += skip;
      sigmaread = 1;
    }
    else if (sscanf(line," fnh = %d%n",&dumbi,&skip) == 1)
    {
      fnh = dumbi;
      line += skip;
      fnhread = 1;
    }
    else if (sscanf(line," fnw = %d%n",&dumbi,&skip) == 1)
    {
      fnw = dumbi;
      line += skip;
      fnwread = 1;
    }
    else if (sscanf(line," pnh = %d%n",&dumbi,&skip) == 1)
    {
      pnh = dumbi;
      line += skip;
      pnhread = 1;
    }
    else if (sscanf(line," pnw = %d%n",&dumbi,&skip) == 1)
    {
      pnw = dumbi;
      line += skip;
      pnwread = 1;
    }
    else if (sscanf(line," nl = %d%n",&dumbi,&skip) == 1)
    {
      nl = dumbi;
      line += skip;
      nlread = 1;
    }
    else if (sscanf(line," wx = %lf%n",&dumb,&skip) == 1)
    {
      if (widthdir_tmp == NULL)
        widthdir_tmp = new double [3];
      widthdir_tmp[XX] = units_*dumb;
      line += skip;
      wxread = 1;
    }
    else if (sscanf(line," wy = %lf%n",&dumb,&skip) == 1)
    {
      if (widthdir_tmp == NULL)
        widthdir_tmp = new double [3];
      widthdir_tmp[YY] = units_*dumb;
      line += skip;
      wyread = 1;
    }
    else if (sscanf(line," wz = %lf%n",&dumb,&skip) == 1)
    {
      if (widthdir_tmp == NULL)
        widthdir_tmp = new double [3];
      widthdir_tmp[ZZ] = units_*dumb;
      line += skip;
      wzread = 1;
    }
    else if (sscanf(line," frw = %lf%n",&dumb,&skip) == 1)
    {
      frw = dumb;
      line += skip;
      frwread = 1;
    }
    else if (sscanf(line," frh = %lf%n",&dumb,&skip) == 1)
    {
      frh = dumb;
      line += skip;
      frhread = 1;
    }
    else if (sscanf(line," pfracw = %lf%n",&dumb,&skip) == 1)
    {
      pfracw = dumb;
      line += skip;
      pfracwread = 1;
    }
    else if (sscanf(line," pfrach = %lf%n",&dumb,&skip) == 1)
    {
      pfrach = dumb;
      line += skip;
      pfrachread = 1;
    }
    else if (sscanf(line," pfracl = %lf%n",&dumb,&skip) == 1)
    {
      pfracl = dumb;
      line += skip;
      pfraclread = 1;

    }
    else if ((sscanf(line," %s%n",string, &skip) == 1) &&
             (!(strcmp(string,"nep"))) )
    {
      nep = 1;
      line += skip;
      nepread = 1;
    }
    else
    {
      printf("Error:don't know what piece this is in seg %s: %s\n",name,line);
      return (1);
    }
  }  // while not blank


  if (hread == 0)
  {
    if (default_vals.ish() == 1)
      height = default_vals.h();
    else
    {
      printf("Error:H not given for seg %s\n",name);
      return 1;
    }
  }
  if (wread == 0)
  {
    if (default_vals.isw() == 1)
      width = default_vals.w();
    else
    {
      printf("Error:W not given for seg %s\n",name);
      return 1;
    }
  }

  if (sigmaread == 0)
  {
    if (default_vals.issigma() == 1)
      sigma = default_vals.sigma();
    else
    {
      printf("Error:Sigma or rho not given for seg %s\n",name);
      return 1;
    }
  }
  if (nepread == 0)
  {
    if (default_vals.isnep() == 1)
      nep = default_vals.nep();
    else
      nep = 0;
  }

  if (fnhread == 0)
  {
    if (default_vals.isfnh() == 1)
      fnh = default_vals.fnh();
    else
      fnh = 1;
  }
  if (fnwread == 0)
  {
    if (default_vals.isfnw() == 1)
      fnw = default_vals.fnw();
    else
      fnw = 1;
  }
  if (pnhread == 0)
  {          // nep!=0 ?  (I put here why ?)
    if (default_vals.ispnh() == 1)
      pnh = default_vals.pnh();
    else
      pnh = 1;
  }

  if (pnwread == 0)
  {          // nep!=0 ?  (I put here why ?)
    if (default_vals.ispnw() == 1)
      pnw = default_vals.pnw();
    else
      pnw = 1;
  }
  if (nlread == 0)
  {
    if (default_vals.isnl() == 1)
      nl = default_vals.nl();
    else
      nl = 1;
  }
  if (frwread == 0)
  {
    if (default_vals.isfrw() == 1)
      frw = default_vals.frw();
    else if (fnw > 2)
    {
      printf("Error:ratio of adjacent filaments along width "
             "not given for seg: %s\n",
             name);
      return 1;
    }
  }

  if (frhread == 0)
  {
    if (default_vals.isfrh() == 1)
      frh = default_vals.frh();
    else if (fnh > 2)
    {
      printf("Error:ratio of adjacent filaments along height "
             "not given for seg: %s\n", name);
      return 1;
    }
  }

  if (pfracwread == 0)
  {           // nep!=0 ?  (I put here why ?)
    if (default_vals.ispfracw() == 1)
      pfracw = default_vals.pfracw();
    else if (pnw > 2)
    {
      printf("Error:fraction of edge panels "
             "not given for seg: %s\n", name);
      return 1;
    }
  }

  if (pfrachread == 0)
  {          // nep!=0 ?  (I put here why ?)
    if (default_vals.ispfrach() == 1)
      pfrach = default_vals.pfrach();
    else if (pnh > 2)
    {
      printf("Error:fraction of edge panels "
             "not given for seg: %s\n", name);
      return 1;
    }
  }

  if (pfraclread == 0)
  {
    if (default_vals.ispfracl() == 1)
      pfracl = default_vals.pfracl();
    else if (nl > 2)
    {
      printf("Error:fraction of edge panels "
             "not given for seg: %s\n", name);
      return 1;
    }
  }

  if ((fnw < 1) || (fnh < 1))
  {
    printf("Error: number of filaments less than 1 for seg: %s\n",
           name);
    return 1;
  }

  if ((pnw < 1) || (pnh < 1))
  {
    printf("Error: number of panels less than 1 for seg: %s\n",
           name);
    return 1;
  }

  if (nl < 1)
  {
    printf("Error: number of sections less than 1 for seg: %s\n",
           name);
    return 1;
  }

  if ((frw < 1.0) || (frh < 1.0))
  {
    printf("Error: ratio of adjacent fils less than 1.0 for seg: %s\n",
           name);
    return 1;
  }

  if ((pfracw > 1.0) || (pfrach > 1.0))
  {
    printf("Error: fraction of edge panel greater than 1.0 for seg: %s\n",
           name);
    return 1;
  }

  if ((pfracw <= 0.0) || (pfrach  <= 0.0))
  {
    printf("Error: fraction of edge panel not positive for seg: %s\n",
           name);
    return 1;
  }

  if (pfracl >= 0.5)
  {
    printf("Error: fraction of edge panel (length) not less than 0.5 "
           "for seg: %s\n", name);
    return 1;
  }

  if (pfracl <= 0.0)
  {
    printf("Error: fraction of edge panel (length) not positive "
           "for seg: %s\n", name);
    return 1;
  }

  if ( widthdir_tmp == NULL)
  {

  //  point3D<double> pn1 (*(node[0])), pn2 (*(node[1]));
    point3D<double> pn1 ((*(node[0])).x(),(*(node[0])).y(),(*(node[0])).z()), pn2 ((*(node[1])).x(),(*(node[1])).y(),(*(node[1])).z());
    point3D<double> dir_z (0, 0, 1);
	if (cross(dir_z, (pn2-pn1)).abs() >EPS)	
		widthdir = new point3D<double> (cross(dir_z, (pn2-pn1)));
    else
    {
      point3D<double> dir_y (0, 1, 0);
	  widthdir=new point3D<double> (cross(dir_y, (pn2-pn1)));
	  
    }
    (*widthdir) =(* widthdir) / ((*widthdir).abs());

  }
  else
  {

    if ( (wxread + wyread + wzread) != 3)
    {
      printf("Error:not all of wx wy and wz specified for seg %s\n",name);
      return 1;
    }
    widthdir=new point3D<double> (widthdir_tmp[XX], widthdir_tmp[YY], widthdir_tmp[ZZ]) ;
    (*widthdir) = (*widthdir) / ((*widthdir).abs());

  }

  new_seg=new conductor (name, node[0], node[1], width, height, widthdir,
                         sigma, nl, fnw, fnh, pnw, pnh, frw, frh, pfracl, pfracw, pfrach, type, nep);

  structure.add_conductor (new_seg);


  // add this conductor to list_node_conductors of its two nodes
  node[0]->add_conductor(new_seg);
  node[1]->add_conductor(new_seg);

  return 0;
}

int Initialize_discretization::nothing (char *line)
{
  if (line[0] == '+')
    printf("Nothing to continue.\n%s\n",line);
  else
    printf("What is the following line? \n%s\n", line);
  return (1);
}


int Initialize_discretization::addresistance(char* line, Structure& structure)
{
  int skip, i;
  Node *node[2];
  static char name[80];
  static char names[2][80];
  static char resistance_name[80];
  double resistance_value;
  int err=0;
  Resistance *resist;

  // read name
  sscanf(line,"r%s%n",resistance_name,&skip);
  line += skip;

  // read the nodes which will define the resistance ends
  for(i = 0; i < 2; i++)
  {
    if(sscanf(line, "%s%n",name,&skip) != 1)
    {
      printf("No node name to read for resistance in line:\n%s", line);
      return 1;
    }
    line += skip;

    strcpy(names[i], name);

    // search for the named node in the list
    node[i] = structure.get_node_from_name(name);

    if(node[i] == NULL)
    {
      printf("Resistance: No node read in yet named %s\n",name);
      return 1;
    }

    node[i] = node[i]->real_node();
  }

  if (sscanf(line," %lf%n",&resistance_value,&skip) != 1)
  {
    printf("No value for resistance in line:\n%s", line);
    return 1;
  }
  if (resistance_value < 0.0)
  {
    printf("Resistance value cannot be negative, in line:\n%s", line);
    return 1;
  }
  line += skip;


//  if ((resist = structure.get_resistance_from_resistance_name(resistance_name))
//      != NULL)
//  {
//    fprintf(stderr, "Error: Cannot name resistance between %s and %s as %s.\n"
//            "The name is already used for resistance between %s and %s\n",
//            names[0], names[1], resistance_name,
//            resist->node_name(1), resist->node_name(2));
//    err = 1;
//  }

  if (node[0] == node[1])
  {
    fprintf(stderr, "Error in R: Nodes %s and %s are the "
            "same node. (.equiv'ed maybe?)\n", names[0], names[1]);
    err = 1;
  }

  structure.add_resistance (new Resistance (node[0], node [1],
                            resistance_value, resistance_name) );

  return err;
}

int Initialize_discretization::addlayer(char* line, LayerEnv& layers,
                                        Defaults& default_vals)
{
  static char name[80];
  int skip, epirRead=0, condRead=0, subHRead=0;
  double epir, sigma, subH, dumb;
  Layer* layeradd;

  //read name
  sscanf(line, "%s%n", name, &skip);
  line+=skip;

  while(notblankline(line))
  {
    if (sscanf(line, " layer_epir = %lf%n", &epir, &skip)==1)
    {
      line+=skip;
      epirRead=1;
    }
    else if (sscanf(line, " layer_sigma = %lf%n", &sigma, &skip)==1)
    {
      line+=skip;
      condRead=1;
    }
    else if (sscanf (line, " starth= %lf%n", &dumb, &skip)==1)
    {
      subH=units_*dumb;
      line+=skip;
      subHRead=1;
    }
    else
    {
      printf("don't know what piece this is in layer %s: %s\n", name, line);
      return 1;
    }
  } //matching while
  if (epirRead==0)
  {
    if (default_vals.isepir()==1)
      epir=default_vals.epir();
    else
    {
      printf("Error: epir is not given for layer %s\n", name);
      return 1;
    }
  }
  if(condRead==0)
  {
    if (default_vals.isLayerSigma()==1)
      sigma=default_vals.layerSigma();
    else
    {
      printf("Error: conductivity is not given for layer %s\n", name);
      return 1;
    }
  }
  if (subHRead==0)
  {
    printf("Error: subH is not given for layer %s\n", name);
    return 1;
  }

  //make sure the layer is not already defined by name or start height
  if ((layers.get_layer_from_name(name) !=NULL))
  {
    fprintf(stderr,"layer %s already defined.  Cannot redefine it\n", name);
    exit(1);
  }
  if ( (layers.get_layer_from_startHeight(subH) !=NULL))
  {
    fprintf(stderr,"layer with height %lf already defined.  Cannot redefine it\n", subH);
    exit(1);
  }


  layeradd=new Layer (name, epir, sigma, subH);
  layers.add_defined_layer(layeradd);

  return 0;
}
