/**********************************************************************
 * 
 * hmm_util.c
 * 
 * copyright (c) 2001-9, Karl W Broman
 * modified from hmm_main.c by Brian S Yandell and Laura M Shannon (c) 2011
 *
 *     This program is free software; you can redistribute it and/or
 *     modify it under the terms of the GNU General Public License,
 *     version 3, as published by the Free Software Foundation.
 * 
 *     This program is distributed in the hope that it will be useful,
 *     but without any warranty; without even the implied warranty of
 *     merchantability or fitness for a particular purpose.  See the GNU
 *     General Public License, version 3, for more details.
 * 
 *     A copy of the GNU General Public License, version 3, is available
 *     at http://www.r-project.org/Licenses/GPL-3
 * 
 * C functions for the R/qtl package
 *
 * Contains: init_stepf, stepfc, forward, backward, golden
 *
 * These are used in hmm_bcsft to simplify calculations.
 * They could be used in hmm_main.c
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h> 
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include "hmm_main.h"
#include "util.h"


void printArrayDouble(int numIndices, double *array) {
	int i;
	char string[100];
	for (i = 0; i < numIndices; i++) {
		sprintf(string, "%f\t", array[i]);
		Rprintf(string);
	}
	Rprintf("\n");
}

void printArrayInt(int numIndices, int *array) {
	int i;
	char string[100];
	for (i = 0; i < numIndices; i++) {
		sprintf(string, "%d\t", array[i]);
		Rprintf(string);
	}
	Rprintf("\n");
}

void init_stepf(double *rf, double *rf2, int n_gen, int n_mar, int *cross_scheme, 
		double stepf(int, int, double, double, int *),
		double **probmat)
{
  Rprintf("Starting init_stepf()\n");

  int j,obs1,obs2,tmp1;
  char verboseString[100];
    
  for(j=0; j<n_mar; j++) {
    for(obs2=1; obs2<=n_gen; obs2++) {
      tmp1 = ((obs2 * (obs2 - 1)) / 2) - 1;
      for(obs1=1; obs1<=obs2; obs1++)
	probmat[j][obs1 + tmp1] = stepf(obs1, obs2, rf[j], rf2[j], cross_scheme); /*rf2[j] is labeled as a "junk" input argument in prob_bcsftb()"*/
    }
  }
  Rprintf("probmat[marker][genotype] now contains:\n"); /*These are the probabilities of true genotypes at each marker? */
  for(j = 0; j < n_mar; j++) {
	printArrayDouble(10, probmat[j]);
  } 
}

void init_stepf_exHet(double *rf, double *rf2, int n_gen, int n_mar, int *cross_scheme, 
		double stepf(int, int, double, double, int *, double *),
		double **probmat, double *het)
{
  int j,obs1,obs2,tmp1;
  
  for(j=0; j<n_mar; j++) {
    for(obs2=1; obs2<=n_gen; obs2++) {
      tmp1 = ((obs2 * (obs2 - 1)) / 2) - 1;
      for(obs1=1; obs1<=obs2; obs1++)
	probmat[j][obs1 + tmp1] = stepf(obs1, obs2, rf[j], rf2[j], cross_scheme, het);
    }
  }
}

double stepfc(int obs1, int obs2, int mar, double **probmat)
{
  int tmp1;
  
  /* make obs1 <= obs2 */
  if(obs1 > obs2) {
    tmp1 = obs2;
    obs2 = obs1;
    obs1 = tmp1;
  }
  tmp1 = ((obs2 * (obs2 - 1)) / 2) - 1;
  return(probmat[mar][obs1 + tmp1]);
}

void forward_prob(int i, int n_mar, int n_gen, int curpos, int *cross_scheme, double error_prob,
	     int **Geno, double **probmat, double **alpha,
	     double initf(int, int *), 
	     double emitf(int, int, double, int *))
{
  char verboseString[100]; 
  Rprintf("Starting forward_prob()\n");
   	/* forward equations */

  /* Note: true genotypes coded as 1, 2, ...
     but in the alpha's and beta's, we use 0, 1, ... */

  int j,v,v2;
  double errortol,salpha;

  /* initialize alpha */
  /* curpos = -1: use error_prob always */
  /* curpos >= 0: use TOL except when j == curpos, then use error_prob */

  errortol = error_prob;
  if(curpos > 0) errortol = TOL;
  for(v=0; v<n_gen; v++) { //initf here is init_bcsftb; emitf is emit_bcsftb.
    alpha[v][0] = initf(v+1, cross_scheme) + emitf(Geno[0][i], v+1, errortol, cross_scheme);
    Rprintf("alpha[v][0]:");
    sprintf(verboseString, "%f", alpha[v][0]);
    Rprintf(verboseString);
    Rprintf("\t");
  }
  Rprintf("\n");
  if(curpos == 0) errortol = TOL;

  for(j=1; j<n_mar; j++) {
    if(curpos == j) errortol = error_prob;
    
    for(v=0; v<n_gen; v++) {
      salpha = alpha[0][j-1] + stepfc(1, v+1, j-1, probmat);
      
      for(v2=1; v2<n_gen; v2++)
	salpha = addlog(salpha, alpha[v2][j-1] + stepfc(v2+1, v+1, j-1, probmat));
      
      alpha[v][j] = salpha + emitf(Geno[j][i], v+1, errortol, cross_scheme);
    }
    if(curpos == j) errortol = TOL;
  }
}
void backward_prob(int i, int n_mar, int n_gen, int curpos, int *cross_scheme, double error_prob,
		   int **Geno, double **probmat, double **beta,
		   double initf(int, int *), 
		   double emitf(int, int, double, int *))
{
  Rprintf("Starting backward_prob()\n");
  /* backward equations */

  /* Note: true genotypes coded as 1, 2, ...
     but in the alpha's and beta's, we use 0, 1, ... */

  /* could divide this into forward and backward, then use forward second time
     in est_map */

  int j2,v,v2;
  double errortol,sbeta;

  /* initialize alpha and beta */
  for(v=0; v<n_gen; v++)
    beta[v][n_mar-1] = 0.0;

  /* curpos = -1: use error_prob always */
  /* curpos >= 0: use TOL except when j2+1 == curpos, then use error_prob */
  errortol = error_prob;
  if(curpos >= 0) errortol = TOL;

  for(j2=n_mar-2; j2>=0; j2--) {
    if(curpos == j2+1) errortol = error_prob;
    
    for(v=0; v<n_gen; v++) {
      sbeta = beta[0][j2+1] + stepfc(v+1, 1, j2, probmat) +
	emitf(Geno[j2+1][i], 1, errortol, cross_scheme);
      
      for(v2=1; v2<n_gen; v2++) {
	sbeta = addlog(sbeta, beta[v2][j2+1] + stepfc(v+1, v2+1, j2, probmat) +
		       emitf(Geno[j2+1][i], v2+1, errortol, cross_scheme));
      }
      beta[v][j2] = sbeta;
    }
    if(curpos == j2+1) errortol = TOL;
  }
}

void calc_probfb(int i, int n_mar, int n_gen, int curpos, double **alpha, double **beta,
		 double ***Genoprob)
{
  int j,v,j0,jmax;
  double s;

  j0 = 0;
  jmax = n_mar;
  if(curpos >= 0) {
    j0 = curpos;
    jmax = j0 + 1;
  }

  /* calculate genotype probabilities */
  for(j=j0; j<jmax; j++) {
    s = Genoprob[0][j][i] = alpha[0][j] + beta[0][j];
    for(v=1; v<n_gen; v++) {
      Genoprob[v][j][i] = alpha[v][j] + beta[v][j];
      s = addlog(s, Genoprob[v][j][i]);
    }
    for(v=0; v<n_gen; v++) 
      Genoprob[v][j][i] = exp(Genoprob[v][j][i] - s);
  }
}

double golden_search(double *countmat, int n_gen, int maxit, double tol, int *cross_scheme,
	      double comploglik(double, int, double *, int *))
{
  /* Golden section search. */
  /* en.wikipedia.org/wiki/Golden_section_search */

  static double resphi = 0.0;
  double x[4],y[4];
  int iter;

  if(resphi == 0.0)
    resphi = 1.5 - sqrt(5.0) / 2.0;
  
  x[0] = 0.0;
  x[2] = 1.0;
  y[0] = comploglik(0.0, n_gen, countmat, cross_scheme);
  y[2] = comploglik(0.5, n_gen, countmat, cross_scheme);

  if(y[2] < y[0]) {
    x[1] = x[0];
    x[0] = x[2];
    x[2] = x[1];
    y[1] = y[0];
    y[0] = y[2];
    y[2] = y[1];
  }
  
  x[1] = x[0] + resphi * (x[2] - x[0]);
  y[1] = comploglik(x[1], n_gen, countmat, cross_scheme);

  /* x0 and x2 are the current bounds; the minimum is between them.
   * x1 is the center point, which is closer to x0 than to x2. */

  for(iter=0; iter<maxit; iter++) {
    /* Create a new possible center in the area between x1 and x2, closer to x1. */
    x[3] = x[1] + resphi * (x[2] - x[1]);

    /* Evaluate termination criterion */
    if(fabs(x[2] - x[0]) < tol)
      break;
 
    y[3] = comploglik(x[3], n_gen, countmat, cross_scheme);

    if(y[3] >= y[1]) {
      x[0] = x[1];
      x[1] = x[3];
      y[0] = y[1];
      y[1] = y[3];
    }
    else {
      x[2] = x[0];
      x[0] = x[3];
      y[2] = y[0];
      y[0] = y[3];
    }
  }
  /* handle boundary situations cleanly */
  if((x[0] == 0.0 && y[0] >= y[1]) || (x[2] == 0.0 && y[2] >= y[1])) return(0.0);
  if((x[0] == 1.0 && y[0] >= y[1]) || (x[2] == 1.0 && y[2] >= y[1])) return(1.0);

  x[1] = (x[2] + x[0]) / 2.0;
  /* make negative if does not converge */
  if(iter >= maxit)
    x[1] = - x[1];
  return(x[1]);
}
 
double golden_search_exHet(double *countmat, int n_gen, int maxit, double tol, int *cross_scheme,
	      double comploglik(double, int, double *, int *, double *), double *het, FILE *outFile)
{
  /* Golden section search. */
  /* en.wikipedia.org/wiki/Golden_section_search */

  static double resphi = 0.0;
  double x[4],y[4];
  int iter;
  double rf;
  char verboseString[100];

  printArrayDouble(15, countmat);

  // 1:2:1, t=2 r = 0.1

  countmat[0] = 202.5;
  countmat[1] = 90;
  countmat[2] = 410;
  countmat[3] = 5;
  countmat[4] = 90;
  countmat[5] = 202.5;

  //Excess heterozygosity, h = 0.75, r = 0.1
  countmat[0] = 86.2601;
  countmat[1] = 75.3498;
  countmat[2] = 674.65044;
  countmat[3] = 2.12988;
  countmat[4] = 75.3498;
  countmat[5] = 86.2601;

  //1:2:1 t=8 r=0.1
  //
  //Excess heterozygosity t=8, h=0.75, r = 0.1
  //
  countmat[0] = 317.118;
  countmat[1] = 69.8702;
  countmat[2] = 63.61365;
  countmat[3] = 162.4108;
  countmat[4] = 69.8702;
  countmat[5] = 317.118;

  //Excess heterozygosity t=3, h=0.75, r = 0.1
  //
  countmat[0] = 153.173;
  countmat[1] = 107.3472;
  countmat[2] = 455.1529;
  countmat[3] = 23.8062;
  countmat[4] = 107.3472;
  countmat[5] = 153.173;

  //Excess heterozygosity t=4, h=0.75, r = 0.1
  //
  countmat[0] = 204.919;
  countmat[1] = 114.806;
  countmat[2] = 307.0685;
  countmat[3] = 53.4808;
  countmat[4] = 114.806;
  countmat[5] = 204.919;

  //Excess heterozygosity t=5, h=0.75, r = 0.1
  //
  countmat[0] = 244.824;
  countmat[1] = 109.2422;
  countmat[2] = 207.16382;
  countmat[3] = 84.7042;
  countmat[4] = 109.2422;
  countmat[5] = 244.824;

  //Excess heterozygosity t=6, h=0.75, r = 0.1
  //
  countmat[0] = 275.519;
  countmat[1] = 97.5414;
  countmat[2] = 139.76345;
  countmat[3] = 114.116;
  countmat[4] = 97.5414;
  countmat[5] = 275.519;

  //Excess heterozygosity t=7, h=0.75, r = 0.1
  //
  countmat[0] = 299.076;
  countmat[1] = 83.6872;
  countmat[2] = 94.2913;
  countmat[3] = 140.182;
  countmat[4] = 83.6872;
  countmat[5] = 299.076;
  
  //Excess heterozygosity t=7, h=0.75, r = 0.1
  //
  countmat[0] = 299.076;
  countmat[1] = 83.6872;
  countmat[2] = 94.2913;
  countmat[3] = 140.182;
  countmat[4] = 83.6872;
  countmat[5] = 299.076;
  
  //Excess heterozygosity t=7, h=0.75, r = 0.01
  //
  countmat[0] = 397.122;
  countmat[1] = 11.82784;
  countmat[2] = 166.150652;
  countmat[3] = 15.94934;
  countmat[4] = 11.82784;
  countmat[5] = 397.122;
  
  //Excess heterozygosity t=3, h=0.75, r = 0.01
  //
  countmat[0] = 211.277;
  countmat[1] = 12.7473;
  countmat[2] = 549.75316;
  countmat[3] = 2.19788;
  countmat[4] = 12.7473;
  countmat[5] = 211.277;
 
  //Excess heterozygosity t=3, h=0.75, r = 0.001
  //
  countmat[0] = 217.993;
  countmat[1] = 1.29659;
  countmat[2] = 561.2031247;
  countmat[3] = 0.216846;
  countmat[4] = 1.29659;
  countmat[5] = 217.993;
  
  //Excess heterozygosity t=3, h=0.75, r = 0.25
  //
  countmat[0] = 89.3508;
  countmat[1] = 200.538;
  countmat[2] = 361.9633;
  countmat[3] = 58.2612;
  countmat[4] = 200.538;
  countmat[5] = 89.3508;
  
  printArrayDouble(15, countmat);

  for (rf = 0; rf <= 0.5; rf = rf+0.001) {
	  fprintf(outFile, "%f\t%f\n", rf, comploglik(rf, n_gen, countmat, cross_scheme, het));
	  //sprintf(verboseString, "%f\t%f\n", rf, comploglik(rf, n_gen, countmat, cross_scheme, het));
	  //Rprintf(verboseString);
  }

  if(resphi == 0.0)
    resphi = 1.5 - sqrt(5.0) / 2.0;

  x[0] = 0.0;
  x[2] = 1.0;
  y[0] = comploglik(0.0, n_gen, countmat, cross_scheme, het);
  y[2] = comploglik(0.5, n_gen, countmat, cross_scheme, het);

  if(y[2] < y[0]) {
    x[1] = x[0];
    x[0] = x[2];
    x[2] = x[1];
    y[1] = y[0];
    y[0] = y[2];
    y[2] = y[1];
  }
  
  x[1] = x[0] + resphi * (x[2] - x[0]);
  y[1] = comploglik(x[1], n_gen, countmat, cross_scheme, het);

  /* x0 and x2 are the current bounds; the minimum is between them.
   * x1 is the center point, which is closer to x0 than to x2. */

  for(iter=0; iter<maxit; iter++) {
    /* Create a new possible center in the area between x1 and x2, closer to x1. */
    x[3] = x[1] + resphi * (x[2] - x[1]);

    /* Evaluate termination criterion */
    if(fabs(x[2] - x[0]) < tol)
      break;
 
    y[3] = comploglik(x[3], n_gen, countmat, cross_scheme, het);

    if(y[3] >= y[1]) {
      x[0] = x[1];
      x[1] = x[3];
      y[0] = y[1];
      y[1] = y[3];
    }
    else {
      x[2] = x[0];
      x[0] = x[3];
      y[2] = y[0];
      y[0] = y[3];
    }
  }
  /* handle boundary situations cleanly */
  if((x[0] == 0.0 && y[0] >= y[1]) || (x[2] == 0.0 && y[2] >= y[1])) return(0.0);
  if((x[0] == 1.0 && y[0] >= y[1]) || (x[2] == 1.0 && y[2] >= y[1])) return(1.0);

  x[1] = (x[2] + x[0]) / 2.0;
  /* make negative if does not converge */
  if(iter >= maxit)
    x[1] = - x[1];
  return(x[1]);
}

/* end of hmm_util.c */
