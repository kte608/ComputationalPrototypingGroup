/* coilgen.c

 * Copyright (C) 2001 Claudio Girardi
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <ctype.h>		/* for isspace() */
#include <time.h>		/* for time() and ctime() */
#include <math.h>

#define VERSION "0.1alpha"

/* maximum line length for reading data files */
#define MAX_LINE_LENGTH 100

#ifndef PI
#ifndef M_PI
#define M_PI 3.141592654
#endif
#define PI M_PI
#endif


extern char *optarg;
extern int optind, opterr, optopt;

static char *fname = NULL;	/* alternate configuration file name */
static char *program_name;
static FILE *fin = NULL, *fout = NULL;
static int discr;
static float theta;

typedef struct _point_t {
  float x;
  float y;
  float z;
} point_t;

typedef struct _coil_data_t {
  float d;			/* coil diameter */
  float l;			/* coil length */
  float n;			/* number of turns */
  float wd;			/* wire diameter */
  point_t c;			/* coil center position */
} coil_data_t;

/* coils data are stored in a list */
typedef struct cln_t *cln_p;
typedef struct cln_t {		/* a coil list node */
  cln_p prev;			/* previous coil */
  coil_data_t coil;		/* coil data */
  cln_p next;			/* next coil */
} cln_t;


static void usage(void)
{
  fprintf(stderr,
	  "%s, by IN3OTD, %s, Version " VERSION "\n"
	  "Usage: %s [options] [outfile]\n"
/*	  "options: -c : load alternate coils description file\n" */
	  "options: -d : number of segments per turn\n"
	  "         -h : print this help\n"
	  "         -v : print version information\n",
	  program_name, __DATE__, program_name);
  fprintf(stderr, "\n");
  fprintf(stderr, "Mail bug reports and suggestions to <in3otd@qsl.net>\n");
  fprintf(stderr, "\n");
  exit(1);
}


static void version(void)
{
  fprintf(stderr,
	  /* State the canonical name for the program */
	  "coilgen " VERSION "\n"
	  "Copyright (C) 2001 Claudio Girardi <in3otd@qsl.net>\n"
	  "This program is distributed in the hope that it will be useful,\n"
	  "but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
	  "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
	  "GNU General Public License for more details.\n"
	  "\n");
}


static void err_die(char *err)
{
  fprintf(stderr, "%s: %s: %s\n", program_name, err, strerror(errno));
  exit(-1);
}


static void parse_command_line(int argc, char *argv[])
{
  int optc;			/* For getopt */

  /* Parse the options */
  while ((optc = getopt(argc, argv, "d:hv")) != -1) {
    switch (optc) {
    case 'a': /* for further releases... */
      theta = atof(optarg);
      break;
    case 'c': /* for further releases... */
      fname = strdup(optarg);
      break;
    case 'd':
      discr = atoi(optarg);
      break;
    case 'h':
      usage();
      exit(0);
      break;
    case 'v':
      version();
      exit(0);
      break;
    default:
      usage();
      exit(0);
    }
  }

  if (argc == optind + 1) {
    /* try to open the filename given */
    if (!(fout = fopen(argv[optind], "wb")))
      err_die(argv[optind]);
  } else {
    /* output to standard output */
    fout = stdout;
  }
}


/**
 * read_coils_file - read coils definition from file
 * @fin: input file pointer
 * @cln: pointer to a pointer to a coil list node
 *
 * read_coils_file() reads coils definition from the file pointed 
 * by @fin into a list structure pointed by @cln.
 */
static int read_coils_file(FILE * fin, cln_t ** cln)
{

  char buf[MAX_LINE_LENGTH], *line, *cp;
  int line_no = 0;
  cln_t *tmp_cln;
  int f_no;

  if (!cln)			/* invalid list node pointer */
    exit(-1);

  while (fgets(buf, MAX_LINE_LENGTH, fin)) {	/* read a whole line */
    line_no++;			/* number of line being read */
    /* trim leading whitespace */
    for (cp = buf; *cp && isspace((unsigned char) (*cp)); cp++);
    line = cp;
    /* remove comments */
    if ((cp = strchr(line, '#')) != 0)
      *cp = (char) 0;
    /* skip blank lines */
    if (!*line)			/* ...that is, if line == "\0" */
      continue;

    /* allocate a list node */
    tmp_cln = (cln_t *) malloc(sizeof(cln_t));
    tmp_cln->prev = tmp_cln->next = NULL;

    /* "parse" a line */
    f_no = sscanf(line, "%f%f%f%f%f%f%f", &tmp_cln->coil.d, &tmp_cln->coil.l, &tmp_cln->coil.n, &tmp_cln->coil.wd, &tmp_cln->coil.c.x, &tmp_cln->coil.c.y, &tmp_cln->coil.c.z);

    if (f_no == 7) {		/* all fields read */
#ifdef DEBUG
      printf("%f %f %f %f %f %f %f\n", tmp_cln->coil.d, tmp_cln->coil.l, tmp_cln->coil.n, tmp_cln->coil.wd, tmp_cln->coil.c.x, tmp_cln->coil.c.y, tmp_cln->coil.c.z);
#endif
      if (*cln) {
	(*cln)->next = tmp_cln;
	tmp_cln->prev = (*cln);
      }
      (*cln) = tmp_cln;
    } else {
      fprintf(stderr, "%s: error in input file at line %i!\n", program_name, line_no);
      exit(-1);
    }
  }
  return 0;
}


static int draw_coil(point_t p0, float radius, float pitch, float n_turns, int start_node, float wire_width, float wire_height)
{
  int i, N_nodes;
  float x0, y0, z0, x, y, z, r_corr, ind;

  N_nodes = discr * n_turns;

  ind = 0;			/* avoid warning about unused variable */
#ifdef DEBUG
  printf("Radius = %.3f Length = %.3f Pitch = %.3f #turns = %.3f S = %.3f\n", radius, n_turns * pitch, pitch, n_turns, sqrt((M_PI * n_turns * 2.0 * radius) * (M_PI * n_turns * 2.0 * radius) + M_PI * n_turns * pitch * n_turns * pitch));
  ind = 1e-03 * (4 * radius * radius * n_turns * n_turns) /
    (n_turns * pitch + 0.9 * radius);
  printf("Ind = %e uH\n", ind);
  printf("\n");
#endif

  /* correction factor, since we are approximating a circle with a polygon */
  r_corr = cos(M_PI / (2 * discr));
  r_corr *= r_corr;
  radius /= r_corr;

  fprintf(fout, "* Nodes\n");
  for (i = 0; i < N_nodes; i++) {
    x0 = p0.x + pitch * i / discr;
    y0 = p0.y + radius * sin(i * 2.0 * M_PI / discr);
    z0 = p0.z + radius * cos(i * 2.0 * M_PI / discr);

    x = x0 * cos(theta) - y0 * sin(theta);
    y = x0 * sin(theta) + y0 * cos(theta);
    z = z0;

    fprintf(fout, "N%i x=%.3f y=%.3f z=%.3f\n",
	    start_node + i, x, y, z);
  }

  fprintf(fout, "\n");
  fprintf(fout, "* Segments\n");
  for (i = 0; i < N_nodes - 1; i++) {
    fprintf(fout, "E%i N%i N%i w=%.3f h=%.3f\n",
	    start_node + i, start_node + i, start_node + i + 1,
	    wire_width, wire_height);
  }
  fprintf(fout, "\n");

  return (start_node + N_nodes - 1);
}


int main(int argc, char *argv[])
{

  float radius, pitch, n_turns, wire_width, wire_height;
  int start_node, end_node;
  point_t p0;
  cln_t *cln = NULL;
  time_t timep;


  /* Construct the name of the executable, without the directory part.  */
  program_name = strrchr(argv[0], '/');
  if (!program_name)
    program_name = argv[0];
  else
    ++program_name;

  /* Default values */
  theta = 0;
  discr = 16;

  parse_command_line(argc, argv);

  if (!(fin = fopen("coilrc", "rb")))
    err_die("coilrc");
  read_coils_file(fin, &cln);

  /* go to the first coil definition */
  while (cln->prev)
    cln = cln->prev;

  time(&timep);
  fprintf(fout, "* File generated by %s, %s", program_name, ctime(&timep));
  fprintf(fout, "\n");
  fprintf(fout, ".Units mm\n");
  fprintf(fout, ".default sigma=5.8e4\n");
  fprintf(fout, "\n");

  start_node = 0;

  while (cln) {			/* for every coil... */
    /* convert dimensions to millimeters */
    wire_height = wire_width = cln->coil.wd * 1000.0;
    p0.x = (cln->coil.c.x - cln->coil.l / 2.0) * 1000.0;
    p0.y = cln->coil.c.y * 1000.0;
    p0.z = cln->coil.c.z * 1000.0;
    radius = (cln->coil.d / 2.0) * 1000.0;
    n_turns = cln->coil.n;
    pitch = (cln->coil.l / n_turns) * 1000.0;

    end_node = draw_coil(p0, radius, pitch, cln->coil.n, start_node, wire_width, wire_height);

    fprintf(fout, "\n");
    fprintf(fout, "* define the 'port'\n");
    fprintf(fout, ".external N%i N%i\n", start_node, end_node);
    fprintf(fout, "\n");

    cln = cln->next;
    start_node = end_node + 1;

  }

  fprintf(fout, "* frequency range\n");
  fprintf(fout, "*   the frequency is 1/(2*pi), so the resulting reactance\n");
  fprintf(fout, "*   is numerically equal to the inductance\n");
  fprintf(fout, ".freq fmin=0.1592 fmax=0.1592 ndec=1\n");
  fprintf(fout, "\n");
  fprintf(fout, "* The end\n");
  fprintf(fout, ".end\n");

  if (fout != stdout) fclose(fout);

  exit(0);

}
