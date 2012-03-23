/*!\page LICENSE LICENSE
 
Copyright (C) 2003 by the Board of Trustees of Massachusetts Institute of Technology, hereafter designated as the Copyright Owners.
 
License to use, copy, modify, sell and/or distribute this software and
its documentation for any purpose is hereby granted without royalty,
subject to the following terms and conditions:
 
1.  The above copyright notice and this permission notice must
appear in all copies of the software and related documentation.
 
2.  The names of the Copyright Owners may not be used in advertising or
publicity pertaining to distribution of the software without the specific,
prior written permission of the Copyright Owners.
 
3.  THE SOFTWARE IS PROVIDED "AS-IS" AND THE COPYRIGHT OWNERS MAKE NO
REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED, BY WAY OF EXAMPLE, BUT NOT
LIMITATION.  THE COPYRIGHT OWNERS MAKE NO REPRESENTATIONS OR WARRANTIES OF
MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT THE USE OF THE
SOFTWARE WILL NOT INFRINGE ANY PATENTS, COPYRIGHTS TRADEMARKS OR OTHER
RIGHTS. THE COPYRIGHT OWNERS SHALL NOT BE LIABLE FOR ANY LIABILITY OR DAMAGES
WITH RESPECT TO ANY CLAIM BY LICENSEE OR ANY THIRD PARTY ON ACCOUNT OF, OR
ARISING FROM THE LICENSE, OR ANY SUBLICENSE OR USE OF THE SOFTWARE OR ANY
SERVICE OR SUPPORT.
 
LICENSEE shall indemnify, hold harmless and defend the Copyright Owners and
their trustees, officers, employees, students and agents against any and all
claims arising out of the exercise of any rights under this Agreement,
including, without limiting the generality of the foregoing, against any
damages, losses or liabilities whatsoever with respect to death or injury to
person or damage to property arising from or out of the possession, use, or
operation of Software or Licensed Program(s) by LICENSEE or its customers.
 
*/


#include "mulGlobal.h"
int calcp_calls; 
int mv_prods; 
extern int **dirflags; 

extern int num_calcp_panel, num_calcp_grid; 
main(argc, argv)
int argc;
char *argv[];
{
  int ttliter, i, j, num_cond;
  charge *chglist, *nq, *input_problem();
  ssystem *sys, *mulInit();
  double **capmat, dirtimesav, mulsetup, initalltime, ttlsetup, ttlsolve;
  double relperm;
  int autmom, autlev, numMom, numLev;
  int gsize; 

  extern int fulldirops, fullPqops;
  extern int num_dummy_panels, num_dielec_panels; 
  extern int num_both_panels, num_cond_panels, up_size, eval_size;
  extern char *title, *ps_file_base;
  extern long memcount;
  extern double prectime, conjtime, dirtime, multime, uptime, downtime;
  extern double evaltime, lutime, fullsoltime, prsetime;
  extern char *kill_name_list;

  Name *name_list;

#if DUMPPS == ON || DUMPPS == ALL
  char filename[BUFSIZ];
#endif

#if CAPVEW == ON
  extern char **argvals;
  extern int argcnt, m_, q_, dd_;
  extern double ***axes;
#endif

  charge *cp, *cp_t; 
  double pot, x, y, z, theta; 
  FILE *fp, *fopen(); 
  double calcp(), ecalcp(); 
  int time; 

#if (PRECOND !=NONE) 
  fprintf(stderr, "Preconditioner not implemented.\nRecompile code with PRECOND = NONE in mulGlobal.h \n"); 
  exit(0); 
#endif 

  /* initialize memory and time counters, etc. */
  fulldirops = fullPqops = 0;
  prectime = conjtime = dirtime = multime = uptime = downtime = 0.0;
  evaltime = lutime = fullsoltime = mulsetup = 0.0;
  memcount = 0;
  mv_prods = 0; 

 num_calcp_panel =  num_calcp_grid = 0; 
  calcp_calls = 0; 

  CALLOC(title, BUFSIZ, char, ON, AMSC);

  /* initialize defaults, etc */
  autmom = autlev = ON;
  relperm = 1.0;
  argvals = argv;
  argcnt = argc;
  CALLOC(axes, 10, double **, ON, AMSC);
  for(i = 0; i < 10; i++) {
    CALLOC(axes[i], 2, double *, ON, AMSC);
    for(j = 0; j < 2; j++) {
      CALLOC(axes[i][j], 3, double, ON, AMSC);
    }
  }

  gsize = 3; 
  Rcoef = 0.0; 
  /* get the list of all panels in the problem */
  /* - many command line parameters having to do with the postscript
       file dumping interface are passed back via globals (see mulGlobal.c) */
  chglist = input_problem(argv, argc, &autmom, &autlev, &relperm, 
			  &numMom, &numLev, &name_list, &num_cond, &gsize, &Rcoef);

#if CAPVEW == ON
  /* if no fastcap run is to be done, just dump the psfile */
  if(m_) {
#if MATPIC == ON 
    dump_struct(chglist, NULL); 
#else 
    if(!q_) get_ps_file_base(argv, argc);
    dump_ps_geometry(chglist, NULL, 0, dd_);
#endif  
    exit(0); 
  }

#endif 

  starttimer;
  sys = mulInit(autlev, numLev, numMom, chglist, gsize);  /* Set up cubes, charges. */
  stoptimer;
  initalltime = dtime;

  numLev = sys->depth;

  sys->num_cond = num_cond;
  sys->cond_names = name_list;

  fprintf(stdout, "\nINPUT SUMMARY\n");

#if CMDDAT == ON
  fprintf(stdout, "  Expansion order: %d\n", numMom);
  fprintf(stdout, "  Number of partitioning levels: %d\n", numLev);
  fprintf(stdout, "  Overall permittivity factor: %.3g\n", relperm);
  fprintf(stdout, "  Grid Level :%d\n", sys->gridsize+1); 
#endif

  /* Figure out number of panels and conductors. */
  eval_size = up_size = num_dummy_panels = num_dielec_panels = 0;
  num_both_panels = num_cond_panels = 0;
  for(nq = chglist; nq != NULL; nq = nq->next) {
    if(nq->dummy) num_dummy_panels++;
    else if(nq->surf->type == CONDTR) {
      num_cond_panels++;
    }
    else if(nq->surf->type == DIELEC) num_dielec_panels++;
    else if(nq->surf->type == BOTH) num_both_panels++;
  }
  if (num_dielec_panels > 0) { 
    fprintf(stderr, "\nDielectric panels not implemented.   Check input files.\n"); 
    exit(0); 
  }

  up_size = num_cond_panels + num_both_panels + num_dielec_panels;
  eval_size = up_size + num_dummy_panels;

#if DISSRF == OFF
  fprintf(stdout, "Title: `%s'\n", title);
#endif
  fprintf(stdout, "  Total number of panels: %d\n", up_size);
  fprintf(stdout, "    Number of conductor panels: %d\n", num_cond_panels);
  fprintf(stdout, "    Number of dielectric interface panels: %d\n", 
	  num_dielec_panels);
  fprintf(stdout, 
	  "    Number of thin conductor on dielectric interface panels: %d\n", 
	  num_both_panels);
  /*fprintf(stdout, "  Number of extra evaluation points: %d\n", 
	  num_dummy_panels);*/
  fprintf(stdout, "  Number of conductors: %d\n", num_cond);

#if NAMDAT == ON
  dumpCondNames(stdout, name_list);
#endif

  if(num_both_panels > 0) {
    fprintf(stderr, 
	    "Thin cond panels on dielectric interfaces not supported\n");
    exit(0);
  }

#if CKCLST == ON
  fprintf(stdout, "Checking panels...");
  if(has_duplicate_panels(stdout, chglist)) {
    fprintf(stdout, "charge list has duplicates\n");
    exit(-1);
  }
  fprintf(stdout, "no duplicates\n");
#endif

#if MULDAT == ON
  dumpMulSet(sys, numLev, numMom);
#endif
  fflush(stdout);

#if DUMPPS == ON || DUMPPS == ALL
  strcpy(filename, "psmat.ps");
  dump_ps_mat(filename, 0, 0, eval_size, eval_size, argv, argc, OPEN);
#endif
  
  mulMatDirect(sys);		/* Compute the direct part matrices. */

#ifndef NOGRID 
  mulMatGrid(sys);  		/* Build the grid. */
#endif 

#if DIRSOL == OFF		/* with DIRSOL just want to skip to solve */

#if PRECOND == BD
  starttimer;
  bdmulMatPrecond(sys);
  stoptimer;
  prsetime = dtime;		/* preconditioner set up time */
#endif

#if PRECOND != NONE
  printf("Preconditioner not implemented for grid code. \n"); 
  abort(); 
#endif 

#if PRECOND == OL
  starttimer;
  olmulMatPrecond(sys);
  stoptimer;
  prsetime = dtime;		/* preconditioner set up time */
#endif

#if DMPREC == ON
  dump_preconditioner(sys, chglist, 1);	/* dump prec. and P to matlab file */
#endif

#if DPSYSD == ON
  dissys(sys);
#endif

#if CKDLST == ON
  chkList(sys, DIRECT);
#endif

#endif				/* DIRSOL == OFF */
  dumpnums(ON, eval_size, up_size); /* save num/type of pot. coeff calcs */

  dumpnums2(); 
  printf("  calcp() evals, panel: %6d (%6d self) grid: %6d \n", num_calcp_panel, up_size, num_calcp_grid); 
  dirtimesav = dirtime;		/* save direct matrix setup time */
  dirtime = 0.0;		/* make way for direct solve time */

#if DIRSOL == OFF

#if DUMPPS == ON
  dump_ps_mat(filename, 0, 0, eval_size, eval_size, argv, argc, CLOSE);
#endif

#if DISSYN == ON
  dumpSynop(sys);
#endif

#if DMTCNT == ON
  dumpMatBldCnts(sys);
#endif

#endif				/* DIRSOL == ON */

#if OPCNT == ON
  clearCounters();   /* Clears the operation counters. */
#endif 

  fprintf(stdout, "\nITERATION DATA");
  ttliter = capsolve(&capmat, sys, chglist, eval_size, up_size, num_cond,
		     name_list);

#ifndef FAKESPHERE
#if MKSDAT == ON		/* dump symmetrized, 4 pi eps scaled matrix */
  mksCapDump(capmat, num_cond, relperm, &name_list);
#endif
#endif 

#if TIMDAT == ON 
  ttlsetup = initalltime + dirtimesav + mulsetup;
  multime = uptime + downtime + evaltime;
  ttlsolve = dirtime + multime + prectime + conjtime;

  fprintf(stdout, "\nTIME AND MEMORY USAGE SYNOPSIS\n");
#endif

#ifdef OTHER
  if(TIMDAT == ON) {
    fprintf(stdout, 
	    "Warning: compilation with OTHER flag gives incorrect times\n");
  }
#endif


#if TIMDAT == ON
  fprintf(stdout, "Total time: %g\n", ttlsetup + ttlsolve);
  fprintf(stdout, "  Total setup time: %g\n", ttlsetup);
  fprintf(stdout, "    Direct matrix setup time: %g\n", dirtimesav);
  fprintf(stdout, "    Multipole matrix setup time: %g\n", mulsetup);
  fprintf(stdout, "    Initial misc. allocation time: %g\n", initalltime);
  fprintf(stdout, "  Total iterative P*q = psi solve time: %g\n", ttlsolve);
  fprintf(stdout, "    P*q product time, direct part: %g\n", dirtime);
  fprintf(stdout, "    Total P*q time, multipole part: %g\n", multime);
  fprintf(stdout, "      Upward pass time: %g\n", uptime);
  fprintf(stdout, "      Downward pass time: %g\n", downtime);
  fprintf(stdout, "      Evaluation pass time: %g\n", evaltime);
  fprintf(stdout, "    Preconditioner solution time: %g\n", prectime);
  fprintf(stdout, "    Iterative loop overhead time: %g\n", conjtime);
  fprintf(stdout, "Total MV-products: %d\n", mv_prods); 
  fprintf(stdout, "MV-product time: %g\n", ttlsolve / (double) mv_prods); 

  if(DIRSOL == ON) {		/* if solution is done by Gaussian elim. */
    fprintf(stdout,"\nTotal direct, full matrix LU factor time: %g\n",lutime);
    fprintf(stdout,"Total direct, full matrix solve time: %g\n",fullsoltime);
    fprintf(stdout, "Total direct operations: %d\n", fulldirops);
  }
  else if(EXPGCR == ON) {	/* if solution done iteratively w/o multis */
    fprintf(stdout,"\nTotal A*q operations: %d (%d/iter)\n", 
	    fullPqops, fullPqops/ttliter);
  }

  fprintf(stdout, "Total memory allocated: %d kilobytes ", memcount/1024);
  uallocEfcy(memcount);

  fprintf(stdout, "  Q2M  matrix memory allocated: %7.d kilobytes\n",
	  memQ2M/1024);
  memcount = memQ2M;
  fprintf(stdout, "  Q2L  matrix memory allocated: %7.d kilobytes\n",
	  memQ2L/1024);
  memcount += memQ2L;
  fprintf(stdout, "  Q2P  matrix memory allocated: %7.d kilobytes\n",
	  memQ2P/1024);
  memcount += memQ2P;
  fprintf(stdout, "  L2L  matrix memory allocated: %7.d kilobytes\n",
	  memL2L/1024);
  memcount += memL2L;
  fprintf(stdout, "  M2M  matrix memory allocated: %7.d kilobytes\n",
	  memM2M/1024);
  memcount += memM2M;
  fprintf(stdout, "  M2L  matrix memory allocated: %7.d kilobytes\n",
	  memM2L/1024);
  memcount += memM2L;
  fprintf(stdout, "  M2P  matrix memory allocated: %7.d kilobytes\n",
	  memM2P/1024);
  memcount += memM2P;
  fprintf(stdout, "  L2P  matrix memory allocated: %7.d kilobytes\n",
	  memL2P/1024);
  memcount += memL2P;
  fprintf(stdout, "  Q2PD matrix memory allocated: %7.d kilobytes\n",
	  memQ2PD/1024);
  memcount += memQ2PD;


  fprintf(stdout, "  CUB structure memory allocated: %7.d kilobytes\n",
	  memCUB/1024);
  memcount += memCUB;
  fprintf(stdout, "  CHG structure memory allocated: %7.d kilobytes\n",
	  memCHG/1024);
  memcount += memCHG;

  fprintf(stdout, "  Miscellaneous mem. allocated: %7.d kilobytes\n",
	  memMSC/1024);
  memcount += memMSC;
  fprintf(stdout, "  Total memory (check w/above): %7.d kilobytes\n",
	  memcount/1024);
#endif
  printf("%d calcp() calls\n", calcp_calls); 

}

dump_struct(chglist,qv)
charge *chglist; 
double *qv; 
{

/* write a MATLAB readable file */ 
  charge *cp; 
  double *x, *y, *z, *q; 
  int index; 
  FILE *fp, *fopen(); 
  int type; 
  int size, i; 
  double *c; 

  size = 0; 
  for (cp = chglist; cp != NULL; cp = cp->next) {
    size++; 
  }
  fp = fopen("panels.mat", "w"); 

  x = (double *) calloc( 4*size, sizeof(double)); 
  y = (double *) calloc( 4*size, sizeof(double)); 
  z = (double *) calloc( 4*size, sizeof(double)); 
  q = (double *) calloc( 4*size, sizeof(double)); 
  c = (double *) calloc( 4*size, sizeof(double)); 

  /* do the triangles */ 
  index = 0; 
  for (cp = chglist; cp != NULL; cp = cp->next) {
    if (cp->shape == 3) {
      for (i=0; i<cp->shape; i++) {
        x[3*index+i] = cp->corner[i][0] * cp->X[0] + cp->corner[i][1]*cp->Y[0] + 
          cp->corner[i][2] * cp->Z[0] + cp->x; 

        y[3*index+i] = cp->corner[i][0] * cp->X[1] + cp->corner[i][1]*cp->Y[1] + 
          cp->corner[i][2] * cp->Z[1] + cp->y; 

        z[3*index+i] = cp->corner[i][0] * cp->X[2] + cp->corner[i][1]*cp->Y[2] + 
          cp->corner[i][2] * cp->Z[2] + cp->z; 
        
      }
      if (qv != NULL)
        q[index] = qv[cp->index]; 
      else 
        q[index] = 0.0; 
      c[index] = cp->cond; 
      index++; 
    }
  }
  
  
#ifdef sun 
  type = 1000; 
#else 
  type = 0000; 
#endif 
  if (index > 0) {
    savemat(fp, type, "xt", 3, index, 0, x, NULL);
    savemat(fp, type, "yt", 3, index, 0, y, NULL);
    savemat(fp, type, "zt", 3, index, 0, z, NULL);
    savemat(fp, type, "qt", 3, index, 0, q, NULL);
    savemat(fp, type, "ct", 1, index, 0, q, NULL);
  }
  
  /* now the quads */ 

  index = 0; 
  for (cp = chglist; cp != NULL; cp = cp->next) {
    if (cp->shape == 4) {
      for (i=0; i<cp->shape; i++) {
        x[4*index+i] = cp->corner[i][0] * cp->X[0] + cp->corner[i][1]*cp->Y[0] + 
          cp->corner[i][2] * cp->Z[0] + cp->x; 

        y[4*index+i] = cp->corner[i][0] * cp->X[1] + cp->corner[i][1]*cp->Y[1] + 
          cp->corner[i][2] * cp->Z[1] + cp->y; 

        z[4*index+i] = cp->corner[i][0] * cp->X[2] + cp->corner[i][1]*cp->Y[2] + 
          cp->corner[i][2] * cp->Z[2] + cp->z; 
        
      }
      if (qv != NULL) 
        q[index] = qv[cp->index]; 
      else 
        q[index] = 0.0; 
      c[index] = cp->cond; 
      index++; 
    }
  }
  if (index > 0) {
    savemat(fp, type, "xq", 4, index, 0, x, NULL);
    savemat(fp, type, "yq", 4, index, 0, y, NULL);
    savemat(fp, type, "zq", 4, index, 0, z, NULL);
    savemat(fp, type, "qq", 4, index, 0, q, NULL);
    savemat(fp, type, "cq", 1, index, 0, c, NULL);
  }

  fclose(fp); 

}
