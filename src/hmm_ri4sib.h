/**********************************************************************
 * 
 * hmm_ri4sib.h
 * 
 * copyright (c) 2009, Karl W Broman
 *
 * last modified Apr, 2009
 * first written Apr, 2009
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
 * Contains: init_ri4sib, emit_ri4sib, step_ri4sib, step_special_ri4sib,
 *           calc_genoprob_ri4sib, calc_genoprob_special_ri4sib,
 *           argmax_geno_ri4sib, sim_geno_ri4sib,
 *           est_map_ri4sib, nrec2_ri4sib, logprec_ri4sib, est_rf_ri4sib,
 *           marker_loglik_ri4sib, calc_pairprob_ri4sib, 
 *           errorlod_ri4sib, calc_errorlod_ri4sib
 *
 * These are the init, emit, and step functions plus
 * all of the hmm wrappers for 4-way RIL by sib mating.
 *
 * Genotype codes:    1-4
 * "Phenotype" codes: 0=missing; otherwise binary 1-15, with bit i
 *                    indicating SNP compatible with parent i
 *
 **********************************************************************/

double init_ri4sib(int true_gen, int *cross_scheme);

double emit_ri4sib(int obs_gen, int true_gen, double error_prob, int *cross_scheme);
  
double step_ri4sib(int gen1, int gen2, double rf, double junk, int *cross_scheme);

double step_special_ri4sib(int gen1, int gen2, double rf, double junk, int *cross_scheme);

void calc_genoprob_ri4sib(int *n_ind, int *n_mar, int *geno, 
			  double *rf, double *error_prob, double *genoprob);
  
void calc_genoprob_special_ri4sib(int *n_ind, int *n_mar, int *geno, 
				  double *rf, double *error_prob, double *genoprob);

void argmax_geno_ri4sib(int *n_ind, int *n_pos, int *geno,
			double *rf, double *error_prob, int *argmax);

void sim_geno_ri4sib(int *n_ind, int *n_pos, int *n_draws, int *geno, 
		     double *rf, double *error_prob, int *draws);

void est_map_ri4sib(int *n_ind, int *n_mar, int *geno, double *rf, 
		    double *error_prob, double *loglik, int *maxit, 
		    double *tol, int *verbose);

/* expected no. recombinants */
double nrec2_ri4sib(int obs1, int obs2, double rf, int *cross_scheme);

/* log [joint probability * 8] */
double logprec_ri4sib(int obs1, int obs2, double rf, int *cross_scheme);

void est_rf_ri4sib(int *n_ind, int *n_mar, int *geno, double *rf, 
		   int *maxit, double *tol);

void marker_loglik_ri4sib(int *n_ind, int *geno,
			  double *error_prob, double *loglik);

void calc_pairprob_ri4sib(int *n_ind, int *n_mar, int *geno, 
			  double *rf, double *error_prob, 
			  double *genoprob, double *pairprob);

double errorlod_ri4sib(int obs, double *prob, double error_prob);

void calc_errorlod_ri4sib(int *n_ind, int *n_mar, int *geno, 
			  double *error_prob, double *genoprob, 
			  double *errlod);

/* end of hmm_ri4sib.h */
