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
#include "hmm_bcsft.h"


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

void generateCountmat(double *countmat, double rf, int t, double *het) {
	double transpr[10];
	Rprintf("Within generate countmat\n");
	
  int k;
  for(k=0; k<10; k++)
    transpr[k] = 0.0;
    
  double r = rf;
  double B_11, B_12, B_14, B_22, B_23;
  double u, d;
  double h = *het;

  double hpowt, h2;
  double r2, r3, r4, r5, u2, u3, u4, d2;
  int i;
  char text[200];
    
  r = rf;
  if ((r > 0.4999999) && (r < 0.5000001)) {
	  r = 0.4999999;
  }

  if (r < 0.0000001) {
	  r = 0.0000001;
  }

  //u=(2.0*h)/(2.0-2.0*h);

  hpowt = pow(h, t);
  h2 = pow(h, 2.0);

  r2 = pow(r, 2.0);
  r3 = pow(r, 3.0);
  r4 = pow(r, 4.0);
  r5 = pow(r, 5.0);

  u =  -(2.0*h*r - r + sqrt((r2 - 2.0*h*r + h)*(2.0*h*r - 2.0*r - h + r2 + 1.0)) - 2.0*h*r2 + r2)/(h + 2.0*r - 2.0*h*r + 2.0*h*r2 - 2.0*r2 - 1.0);

  u2 = pow(u, 2.0);
  u3 = pow(u, 3.0);
  u4 = pow(u, 4.0);

  d = 2.0*(pow((1.0-r),2))+8.0*u*r*(1-r)+ 2.0*(r2)+ 2.0*(u2)*((pow((1.0-r),2.0))+(r2));

  d2 = pow(d, 2.0);
  
  B_11 = (1.0/2.0)*((2.0*(- 8.0*r3*u3 + 8.0*r3*u2 + 12.0*r2*u3 - 12.0*r2*u2 - 2.0*d*r2*u + d*r2 - 4.0*r*u3 + 8.0*r*u2 + 2.0*d*r*u - 2.0*d*r - 2.0*u2 + d))/((d + 4.0*r*u2 - 2.0*u2)*(- 4.0*r2*u2 + 4.0*r*u2 - 2.0*u2 + d)) - (d*pow((-(4.0*r*u2 - 2.0*u2)/d),t))/(2.0*(d*u2 + 4.0*r*u4 - 2.0*u4)) + (4.0*hpowt*(- u*r2 + u*r))/(h*(4.0*r2*u2 - 4.0*r*u2 + 2.0*u2 - d*h)) + (d*pow(((4.0*r2*u2 - 4.0*r*u2 + 2.0*u2)/d),t)*(16.0*r2*u2 - 16.0*r3*u2 + 8.0*r4*u2 - d*h - 8.0*r*u2 + 2.0*u2 - 4.0*d*r2*u + 2.0*d*h*r + 4.0*d*r*u - 2.0*d*h*r2 - 4.0*d*h*r*u + 4.0*d*h*r2*u))/(2.0*(2.0*r2*u2 - 2.0*r*u2 + u2)*(32.0*r2*u4 - 32.0*r3*u4 + 16.0*r4*u4 + d2*h - 2.0*d*u2 - 16.0*r*u4 + 4.0*u4 + 4.0*d*r*u2 - 4.0*d*r2*u2 - 2.0*d*h*u2 + 4.0*d*h*r*u2 - 4.0*d*h*r2*u2)));
  
  B_14 = (1.0/2.0)*((d*pow((-(4.0*r*u2 - 2.0*u2)/d),t))/(2.0*(d*u2 + 4.0*r*u4 - 2.0*u4)) + (2.0*(- 8.0*r3*u3 + 12.0*r2*u3 - 2.0*d*r2*u + d*r2 - 4.0*r*u3 + 2.0*d*r*u))/((d + 4.0*r*u2 - 2.0*u2)*(- 4.0*r2*u2 + 4.0*r*u2 - 2.0*u2 + d)) + (4.0*hpowt*(- u*r2 + u*r))/(h*(4.0*r2*u2 - 4.0*r*u2 + 2.0*u2 - d*h)) + (d*pow(((4.0*r2*u2 - 4.0*r*u2 + 2.0*u2)/d),t)*(16.0*r2*u2 - 16.0*r3*u2 + 8.0*r4*u2 - d*h - 8.0*r*u2 + 2.0*u2 - 4.0*d*r2*u + 2.0*d*h*r + 4.0*d*r*u - 2.0*d*h*r2 - 4.0*d*h*r*u + 4.0*d*h*r2*u))/(2.0*(2.0*r2*u2 - 2.0*r*u2 + u2)*(32.0*r2*u4 - 32.0*r3*u4 + 16.0*r4*u4 + d2*h - 2.0*d*u2 - 16.0*r*u4 + 4.0*u4 + 4.0*d*r*u2 - 4.0*d*r2*u2 - 2.0*d*h*u2 + 4.0*d*h*r*u2 - 4.0*d*h*r2*u2)));

  B_12 = (1.0/4.0)*((4.0*d*(- r2 + r)*pow(((4.0*r2*u2 - 4.0*r*u2 + 2.0*u2)/d),t))/((2.0*u*r2 - 2.0*u*r + u)*(4.0*r2*u2 - 4.0*r*u2 + 2.0*u2 - d*h)) - (8.0*hpowt*(- u*r2 + u*r))/(h*(4.0*r2*u2 - 4.0*r*u2 + 2.0*u2 - d*h)));
  
  B_22 = (d*pow(((4.0*r2*u2 - 4.0*r*u2 + 2.0*u2)/d),t))/(4.0*(2.0*r2*u2 - 2.0*r*u2 + u2)) - (d*pow((-(4.0*r*u2 - 2.0*u2)/d),t))/(4.0*(2.0*r*u2 - u2));
  
  B_23 = (d*pow(((4.0*r2*u2 - 4.0*r*u2 + 2.0*u2)/d),t))/(4.0*(2.0*r2*u2 - 2.0*r*u2 + u2)) + (d*pow((-(4.0*r*u2 - 2.0*u2)/d),t))/(4.0*(2.0*r*u2 - u2));

  transpr[0] = B_11;   /* AABB */
  transpr[1] = B_12;   /* AABb aaBb*/
  transpr[2] = B_14;   /* AAbb */
  transpr[3] = B_22;   /* AaBb */
  transpr[4] = B_23;   /* AabB */
  transpr[5] = B_11;   /* aabb */
  transpr[6] = B_12;   /* AaBB Aabb*/
  

  //sprintf(text, "%s\t%f\t%f\t%f\n", "Marginal probabilities 7, 8, 9: ", transpr[7], transpr[8], transpr[9]);
  //Rprintf(text);
  /* marginal probabilities for a single marker from the joint probability function*/
  transpr[7] = transpr[0] + transpr[1] + transpr[2];    /* AA */
  transpr[7] = log(transpr[7]);				 
  transpr[9] = transpr[7];				/* aa */

  transpr[8] = transpr[1] + transpr[3] + transpr[4] + transpr[1];
  transpr[8] = log(transpr[8]);
	
		
  countmat[0] = transpr[0] * 1000.0;
  countmat[1] = transpr[1] * 1000.0 * 2.0;
  countmat[2] = transpr[3] * 1000.0 + transpr[4] * 1000.0;
  countmat[3] = transpr[2] * 1000.0 * 2.0;
  countmat[4] = transpr[6] * 1000.0 * 2.0;
  countmat[5] = transpr[5] * 1000.0;
  return;
}

double golden_search_exHet(double *countmat, int n_gen, int maxit, double tol, int *cross_scheme,
	      double comploglik(double, int, double *, int *, double *), double *het, FILE *outFile)
{
  /* Golden section search. */
  /* en.wikipedia.org/wiki/Golden_section_search */
  warning("This is golden search is only meant for making a plot. This will not work for general purposes (it doesn't even use the input marker data).\n");
  static double resphi = 0.0;
  double x[4],y[4];
  int iter;
  double rf;
  double *pfixedHet;
  double fixedHet;
  fixedHet = 0.648;
  pfixedHet = &fixedHet;
  char verboseString[100];

  printArrayDouble(15, countmat);

  for (rf = 0; rf <= 0.5; rf = rf+0.001) {
  	generateCountmat(countmat, rf, cross_scheme[1], pfixedHet);
  	fprintf(outFile, "%f\t", rf);

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
  	if((x[0] == 0.0 && y[0] >= y[1]) || (x[2] == 0.0 && y[2] >= y[1])) {
  		fprintf(outFile, "%f\n", 0.0); //return(0.0);
		continue;
  	}
  	if((x[0] == 1.0 && y[0] >= y[1]) || (x[2] == 1.0 && y[2] >= y[1])) {
  		fprintf(outFile, "%f\n", 1.0); //return(1.0);
  		continue;
	}
  	x[1] = (x[2] + x[0]) / 2.0;
  	/* make negative if does not converge */
  	if(iter >= maxit)
    		x[1] = - x[1];
    		
  	fprintf(outFile, "%f\n", x[1]);
  }
  return(x[1]);
}

/* end of hmm_util.c */
