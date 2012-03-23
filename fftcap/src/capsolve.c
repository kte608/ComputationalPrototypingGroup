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

/* This routine takes the cube data struct and computes capacitances. */
int capsolve(capmat, sys, chglist, size, real_size, numconds, name_list)
double ***capmat;		/* pointer to capacitance matrix */
ssystem *sys;
charge *chglist;
Name *name_list;
int size, numconds, real_size;	/* real_size = total #panels, incl dummies */
{
  int i, cond, iter, maxiter = MAXITER, ttliter = 0;
  charge *nq;
  double *q, *p, *r, *ap;
  double **bp, **bap;
  extern double fullsoltime;
  surface *surf;
  extern ITER *kill_num_list, *kinp_num_list;
  char *getConductorName();
  extern double iter_tol;

#if CAPVEW == ON
  extern ITER *qpic_num_list;
  extern int q_;
  extern int dd_;
#endif

#if DIRSOL == ON
  extern double *trimat, *sqrmat; /* globals in blkDirect.c */
#endif

  FILE *fp, *fopen(); 

  /* Allocate space for the capacitance matrix. */
  CALLOC((*capmat), numconds+1, double*, ON, AMSC);
  for(i=1; i <= numconds; i++)  {
    CALLOC((*capmat)[i], numconds+1, double, ON, AMSC);
  }

  /* Allocate space for cg vectors , r=residual and p=projection, ap = Ap. */
  CALLOC(q, size+1, double, ON, AMSC);
  CALLOC(r, size+1, double, ON, AMSC);

#if DIRSOL != ON		/* too much to allocate if not used */
  /* allocate for gcr accumulated basis vectors (moved out of loop 30Apr90) */
  fflush(stdout);		/* so header will be saved if crash occurs */

  CALLOC(bp, maxiter+1, double*, ON, AMSC);
  CALLOC(bap, maxiter+1, double*, ON, AMSC);

  /* moved inside of gcr to save memory 22OCT90
  for(i=0; i <= maxiter; i++) {
    CALLOC(bp[i], size+1, double, ON, AMSC);
    CALLOC(bap[i], size+1, double, ON, AMSC);
  } */
#endif

  /* ualloc_verify(); */

  /* P is the "psuedo-charge" for multipole. Ap is the "psuedo-potential". */
  p = sys->q;
  ap = sys->p;

  /* Loop through all the conductors. */
  for(cond=1; cond <= numconds; cond++) {
    
    /* skip conductors in the -rs and the -ri kill list */
    if(want_this_iter(kill_num_list, cond)
       || want_this_iter(kinp_num_list)) continue;

    fprintf(stdout, "\nStarting on column %d (%s)\n", cond, 
	    getConductorName(cond, &name_list));
    fflush(stdout);

    /* Set up the initial residue vector and charge guess. */
    for(i=1; i <= size; i++) r[i] = q[i] = 0.0;
    i = 0;
    for(nq = chglist; nq != NULL; nq = nq->next) {
      if(nq->cond == cond && !nq->dummy 
	 && (nq->surf->type == CONDTR || nq->surf->type == BOTH)) 
        r[nq->index] = 1.0;
    }

    if((iter = gmres(sys,q,p,r,ap,bp,bap,size,maxiter,iter_tol,chglist)) 
       > maxiter) {
      fprintf(stderr, "NONCONVERGENCE AFTER %d ITERATIONS\n", maxiter);
      exit(0);
    }
    ttliter += iter;

#if DMPCHG == LAST
    fprintf(stdout, "\nPanel charges, iteration %d\n", iter);
    dumpChgDen(stdout, q, chglist);
    fprintf(stdout, "End panel charges\n");
#endif

#if CAPVEW == ON
    /* dump shaded geometry file if only if this column picture wanted */
    /* (variable names are messed up - iter list now is list of columns) */
    if(want_this_iter(qpic_num_list, cond) || (q_ && qpic_num_list == NULL)) {
      dump_ps_geometry(chglist, q, cond, dd_);    
    }
#endif

    /* Calc cap matrix entries by summing up charges over each conductor. */
    /* use the permittivity ratio to get the real surface charge */
    /* NOT IMPLEMENTED: fancy stuff for infinitessimally thin conductors */
    /* (once again, permittivity data is poorly organized, lots of pointing) */
    for(i=1; i <= numconds; i++) (*capmat)[i][cond] = 0.0;
    for(nq = chglist; nq != NULL; nq = nq->next) {
      if(nq->dummy || (surf = nq->surf)->type != CONDTR) continue;
      (*capmat)[nq->cond][cond] += surf->outer_perm * q[nq->index];
    }

#if RAWDAT == ON
    if(ITRDAT == OFF) fprintf(stdout, "\n");
    fprintf(stdout, "cond=%d iters=%d\n", cond, iter);

    for(i=1; i <= numconds; i++) {
      fprintf(stdout, "c%d%d=%g  ", i, cond, (*capmat)[i][cond]);
      if(i % 4 == 0) fprintf(stdout, "\n");
    }
    fprintf(stdout, "\n\n");
#endif

#if ITRDAT == ON && RAWDAT == OFF
    fprintf(stdout, "%d iterations\n", iter);
#endif
}
  fflush(stdout);

  for (i=1; i<=size; i++) sys->q[i] = q[i]; 

  return(ttliter);
}


/* 
  Preconditioned(possibly) Generalized Minimum Residual. 
  */
int gmres(sys, q, p, r, ap, bv, bh, size, maxiter, tol, chglist)
ssystem *sys;
double *q, *p, *ap, *r, tol;
double **bv, **bh;
int size, maxiter;
charge *chglist;
{
  int iter, i, j;
  double rnorm, norm, maxnorm=10.0;
  double beta, hi, hip1, length;
  extern double conjtime, prectime;
#if EXPGCR == ON
  extern double *sqrmat;
#endif

  static double *c=NULL, *s=NULL, *g=NULL, *y=NULL;

  starttimer;

  /* Allocation first time through. */
  if(c == NULL) {
    CALLOC(c, size+1, double, ON, AMSC);
    CALLOC(s, size+1, double, ON, AMSC);
    CALLOC(g, size+1, double, ON, AMSC);
    CALLOC(y, size+1, double, ON, AMSC);
  }
  
  /* Set up v^1 and g^0. */
  INNER(rnorm, r, r, size);
  rnorm = sqrt(rnorm);
  for(i=1; i <= size; i++) {
    p[i] = r[i] / rnorm;
    g[i] = 0.0;
  }
  g[1] = rnorm;

  stoptimer;
  conjtime += dtime;

#if ITRDAT == ON
  fprintf(stdout, "||res|| = %g\n", rnorm); /* initial guess residual norm */
#endif
  
  for(iter = 1; (iter <= maxiter) && (rnorm > tol) && (iter <=size); iter++) {
    
    starttimer;
    /* allocate the back vectors if they haven't been already */
    if(bv[iter] == NULL) {
      CALLOC(bv[iter], size+1, double, ON, AMSC);
      CALLOC(bh[iter], iter+2, double, ON, AMSC);
    }
    
    /* Save p as the v{iter}. */
    for(i=1; i <= size; i++) bv[iter][i] = p[i];
    
    stoptimer;
    conjtime += dtime;

    /* Form Av{iter}. */
    computePsi(sys, p, ap, size, chglist);

    starttimer;
    
    /* Initialize v^{iter+1} to Av^{iter}. */
    for(i=1; i <= size; i++) p[i] = ap[i];
    
    /* Make v^{iter+1} orthogonal to v^{i}, i <= iter. */
    for(j=1; j <= iter; j++) {
      INNER(hi, ap, bv[j], size);
      for(i=1; i <= size; i++) p[i] -= hi * bv[j][i];
      bh[iter][j] = hi;
    }
    
    /* Normalize v^{iter+1}. */
    INNER(norm, p, p, size);    
    norm = sqrt(norm);
    if (norm > 0.0) 
      for(i=1; i <= size; i++) p[i] /= norm;
    bh[iter][iter+1] = norm;
    
    /* Apply rotations to new h column. */
    for(i=1; i < iter; i++) {
      hi = bh[iter][i];
      hip1 = bh[iter][i+1];
      bh[iter][i] = c[i] * hi - s[i] * hip1;
      bh[iter][i+1] = c[i] * hip1 + s[i] * hi;
    }
    
    /* Compute new rotations. */
    hi = bh[iter][iter];
    hip1 = bh[iter][iter+1];
    length = sqrt(hi * hi + hip1 * hip1);
    c[iter] = hi/length;
    s[iter] = -hip1/length;
    
    /* Apply new rotations. */
    bh[iter][iter] = c[iter] * hi - s[iter] * hip1;
    bh[iter][iter+1] = c[iter] * hip1 + s[iter] * hi;
    /* ASSERT(g[iter+1] == 0); WHY IS THIS HERE ??? */
    hi = g[iter];
    g[iter] = c[iter] * hi;
    g[iter+1] = s[iter] * hi;    
    
    rnorm = ABS(g[iter+1]);

    stoptimer;
    conjtime += dtime;
    
#if ITRDAT == ON
    fprintf(stdout, "||res|| = %g\n", rnorm);
#else
    fprintf(stdout, "%d ", iter);
    if((iter) % 15 == 0 && iter != 0) fprintf(stdout, "\n");
#endif
    fflush(stdout);
  }
  /* Decrement from the last increment. */
  iter--;

  starttimer;
  
  /* Compute solution, note, bh is bh[col][row]. */
  for(i=1; i <= iter; i++) y[i] = g[i];
  for(i = iter; i > 0; i--) {
    y[i] /=  bh[i][i];
    for(j = i-1; j > 0; j--) {
      y[j] -= bh[i][j]*y[i];
    }
  }
  for(i=1; i <= size; i++) {
    q[i] = 0.0;
    for(j=1; j <= iter; j++) {
      q[i] += y[j] * bv[j][i];
    }
  }

  stoptimer;
  conjtime += dtime;
  
#if PRECOND != NONE
  /* Undo the preconditioning to get the real q. */
  starttimer;
  for(i=1; i <= size; i++) {
    p[i] = q[i];
    ap[i] = 0.0;
  }
  mulPrecond(sys, PRECOND);
  for(i=1; i <= size; i++) {
    q[i] = p[i];
  }
  stoptimer;
  prectime += dtime;
#endif

  if(rnorm > tol) {
    fprintf(stdout, "\ngmres: WARNING exiting without converging\n");
  }
  return(iter);
}

extern int mv_prods; 
/* 
ComputePsi computes the potential from the charge vector, or may
include a preconditioner.  It is assumed that the vectors for the
charge and potential have already been set up and that the potential
vector has been zeroed.  ARBITRARY VECTORS CAN NOT BE USED.
*/
computePsi(sys, q, p, size, chglist)
ssystem *sys;
double *q, *p;
int size;
charge *chglist;
{
  extern double dirtime, uptime, downtime, evaltime, prectime;
  extern int real_size;
  int i,j;
  FILE *fp, *fopen(); 
  double *atmp; 
  static double *ptmp = NULL; 
  double dtmp, rtmp, rms; 

  ASSERT(p == sys->p);
  ASSERT(q == sys->q);
  for(i=1; i <= size; i++) p[i] = 0;

#if PRECOND != NONE
  starttimer;
  mulPrecond(sys, PRECOND);
  stoptimer;
  prectime += dtime;
#endif

  starttimer;
  mulDirect(sys); 

  stoptimer;

  dirtime += dtime;

#ifndef NOGRID 
  starttimer; 
  mulGrid(sys);   
  stoptimer; 
#endif 

  evaltime += dtime; 

#if DMPCHG == LAST
  fprintf(stdout, "\nPanel potentials divided by areas\n");
  dumpChgDen(stdout, p, chglist);
  fprintf(stdout, "End panel potentials\n");
#endif

#if OPCNT == ON
  printops();
#endif				/* OPCNT == ON */

  mv_prods++; 

}

void pinv(A, rows, cols)
double **A; 
int rows, cols; 
{
  
  char jobu[1], jobvt[1]; 
  int m[1], n[1], lda[1], lwork[1], info[1]; 
  int i,j,k; 
  double *dptr; 
  int ldu[1], ldvt[1]; 
  double *U, *V, *dwork; 
  double tol; 
  double *svals, cond; 
  double tmp; 
  int nsv, rank; 
  FILE *fp, *fopen(); 
  
  nsv = MIN(rows,cols); 
  dptr = (double *) calloc(rows*cols+1, sizeof(double)); 
  U = (double *) calloc(rows*rows+1, sizeof(double)); 
  V = (double *) calloc(cols*cols+1, sizeof(double)); 
  svals = (double *) calloc(rows + cols, sizeof(double)); 

  if (rows > cols) {
    printf("pinv code not set up for rows >  cols\n"); 
    exit(0); 
  }

/*  printf("pinv: %dx%d\n", rows, cols); */ 
  
  jobu[0] = 'A'; 
  jobvt[0] = 'A'; 
  m[0] = rows; 
  n[0] = cols; 
  
  lwork[0] = 10 * 8 * (rows + cols); 
  dwork = (double *) calloc(lwork[0], sizeof(double)); 

  ldu[0] = m[0]; 
  ldvt[0] = n[0]; 

  for (i=0; i < rows; i++) {
    for (j=0; j < cols; j++) {
      dptr[i + j * rows] = A[i][j]; 
    }
  }

  { 
    int job; 
    double * evec; 
    evec = (double *) calloc(cols + 1, sizeof(double)); 
    job = 11; 
    dsvdc_(dptr, &m[0], &m[0], &n[0], svals, evec, U, &ldu[0], V, &ldvt[0], dwork, &job, &info[0]); 
    free(evec); 
  }

  tol = 1.0e-15 * svals[0] * MAX(rows,cols); 

  rank = nsv; 
  for (i=0; i< nsv; i++) {
    if (svals[i] < tol) { 
      rank = i;  
      break; 
    }
  }

  if (rank < nsv) { 
    fprintf(stderr, "Warning: \n"); 
    fprintf(stderr, "pinv: rank is %d vs. %d\n", rank, nsv); 
  }

  tol = 1.0e-10 * svals[0]; 
  for (i=0; i< nsv; i++) {
    if (svals[i] < tol)
      svals[i] = 0.0; 
    else 
      svals[i] = 1.0 / svals[i]; 
  }

   
  /* U -> U^T */ 
  for (i=0; i< rows; i++) {
    for (j=0; j< i; j++) {
      tmp = U[i+j*rows]; 
      U[i+j*rows] = U[j+i*rows]; 
      U[j+i*rows] = tmp; 
    }
  }
  
  for (i=0; i< rows; i++) {
    /* loop over rows i */ 
    for (j=0; j < rows; j++) {
      /* loop over cols j */ 
      /* multiply U^T by inverse singular vals */ 
      U[i + j * rows] = U[i+j*rows] * svals[i]; 
    }
  }

  /* a[i][j] = sum_k b[i][k] c[k][j] */ 
  for (i=0; i< cols; i++) {
    for (j=0; j< rows; j++) {
      A[i][j] = 0.0; 
      for (k=0; k < rows; k++) {
        A[i][j] += V[i + k * cols] * U[k + j * rows]; 
      }
    }
  }


  printf("pinv:(%d,%d) \n", rows, cols); 
  if (info[0] > 0) {
    printf("%d superdiagonals did not converge\n", info[0]); 
  }
}

