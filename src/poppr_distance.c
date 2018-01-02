/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# This software was authored by Zhian N. Kamvar and Javier F. Tabima, graduate 
# students at Oregon State University; Jonah C. Brooks, undergraduate student at
# Oregon State University; and Dr. Nik Grünwald, an employee of USDA-ARS.
#
# Permission to use, copy, modify, and distribute this software and its
# documentation for educational, research and non-profit purposes, without fee, 
# and without a written agreement is hereby granted, provided that the statement
# above is incorporated into the material, giving appropriate attribution to the
# authors.
#
# Permission to incorporate this software into commercial products may be
# obtained by contacting USDA ARS and OREGON STATE UNIVERSITY Office for 
# Commercialization and Corporate Development.
#
# The software program and documentation are supplied "as is", without any
# accompanying services from the USDA or the University. USDA ARS or the 
# University do not warrant that the operation of the program will be 
# uninterrupted or error-free. The end-user understands that the program was 
# developed for research purposes and is advised not to rely exclusively on the 
# program for any reason.
#
# IN NO EVENT SHALL USDA ARS OR OREGON STATE UNIVERSITY BE LIABLE TO ANY PARTY 
# FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
# LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, 
# EVEN IF THE OREGON STATE UNIVERSITY HAS BEEN ADVISED OF THE POSSIBILITY OF 
# SUCH DAMAGE. USDA ARS OR OREGON STATE UNIVERSITY SPECIFICALLY DISCLAIMS ANY 
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE AND ANY STATUTORY 
# WARRANTY OF NON-INFRINGEMENT. THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS"
# BASIS, AND USDA ARS AND OREGON STATE UNIVERSITY HAVE NO OBLIGATIONS TO PROVIDE
# MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS. 
#
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <stdio.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R.h>
#include <R_ext/Utils.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
int perm_count;

SEXP pairwise_covar(SEXP pair_vec);
SEXP pairdiffs(SEXP freq_mat);
SEXP permuto(SEXP perm);
SEXP bruvo_distance(SEXP bruvo_mat, SEXP permutations, SEXP alleles, SEXP m_add, SEXP m_loss, SEXP old_model);
double bruvo_dist(int *in, int *nall, int *perm, int *woo, int *loss, int *add, int old_model);
void swap(int *x, int *y);  
void permute(int *a, int i, int n, int *c);
int fact(int x);
double mindist(int perms, int alleles, int *perm, double **dist);
void genome_add_calc(int* genos,
	int perms,
	int alleles,
	int* perm,
	double** dist,
	int zeroes,
	int* zero_ind,
	int curr_zero,
	int miss_ind,
	int* replacement,
	int* replaced_alleles,
	int* facts,
	int inds,
	int curr_ind,
	double* genome_add_sum,
	int* tracker,
	int old_model);
int multinomial_coeff(int* ARR, int n, int* facts);
void genome_loss_calc(int *genos, int nalleles, int *perm_array, int woo, 
		int *loss, int *add, int *zero_ind, int curr_zero, int zeroes, 
		int miss_ind, int curr_allele, double *genome_loss_sum, 
		int* replacements, int* facts, int *loss_tracker, int old_model);
/*
 * UNUSED FUNCTIONS
void fill_short_geno(int *genos, int nalleles, int *perm_array, int *woo, 
		int *loss, int *add, int zeroes, int *zero_ind, int curr_zero, 
		int miss_ind, int *replacement, int inds, int curr_ind, double *res, 
		int *tracker);
void print_distmat(double** dist, int* genos, int p);
*/		
		

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Calculates the root product of pairwise comparisons of each of the variances of
each locus.

Input: A vector of length n where n is the number of loci in the population.
Output: A vector of length n*(n-1)/2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
SEXP pairwise_covar(SEXP pair_vec)
{
	int I;
	int i;
	int j;
	int count;
	SEXP Rout;
	I = length(pair_vec);
	PROTECT(pair_vec = coerceVector(pair_vec, REALSXP));
	PROTECT(Rout = allocVector(REALSXP, (I*(I-1)/2) ));
	count = 0;
	for(i = 0; i < I-1; i++)
	{
		R_CheckUserInterrupt();
		for(j = i+1; j < I; j++)
		{
			REAL(Rout)[count++] = sqrt(REAL(pair_vec)[i] * REAL(pair_vec)[j]);
		}
	}
	UNPROTECT(2); // pair_vec; Rout
	return Rout;
}
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Calculates the absolute yes or no distance of the alleles at a locus. 

Input: An n x m matrix where n is the number of individuals, and m is the number
of alleles at a single locus. 

Output: A vector of length n*(n-1)/2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
SEXP pairdiffs(SEXP freq_mat)
{
	int rows;
	int cols;
	int i;
	int j;
	int k;
	int val;
	int count;
	SEXP Rout;
	SEXP Rdim;
	Rdim = getAttrib(freq_mat, R_DimSymbol);
	rows = INTEGER(Rdim)[0];
	cols = INTEGER(Rdim)[1];
	int* inmat = INTEGER(freq_mat);
	count = 0;
	PROTECT(Rout = allocVector(INTSXP, rows*(rows-1)/2));

	for(i = 0; i < rows-1; i++)
	{
		R_CheckUserInterrupt();
		for(j = i+1; j < rows; j++)
		{
			val = 0;
			for (k = 0; k < cols; k++)
			{
				if (inmat[i + k*rows] == NA_INTEGER || inmat[j + k*rows] == NA_INTEGER)
				{
					val = 0;
					break;
				}
				val += abs(inmat[i + k*rows] - inmat[j + k*rows]);
			}
			INTEGER(Rout)[count++] = val;
		}
	}
	UNPROTECT(1);
	return Rout;
}
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
permuto will return a vector of all permutations needed for bruvo's distance.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
SEXP permuto(SEXP perm)
{
	int permutations;
	int i;
	int per;
	SEXP Rval;
	/* 
		IMPORTANT: INITIALIZE THE COUNTER. THE POINTER IS NOT RELEASED FROM
		MEMORY OTHERWISE.
	*/
	perm_count = 0;
	perm = coerceVector(perm, INTSXP);
	per = INTEGER(perm)[0];
	int *allele_array;
	allele_array = R_Calloc(per, int);
	permutations = fact(per)*per;
	for(i = 0; i < per; i++)
	{
		allele_array[i] = i;
	}
	PROTECT(Rval = allocVector(INTSXP, permutations));
	permute(allele_array, 0, per-1, INTEGER(Rval));
	UNPROTECT(1);
	R_Free(allele_array);
	return Rval;
}
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Calculates Bruvo's distance over a matrix of individuals by loci. This is
calcluated regardless of ploidy. For more information, see Bruvo et al. 2006

Parameters:
bruvo_mat - a matrix of individuals by loci, one column per allele.
permutations - a vector of indeces for permuting the number of alleles. 
alleles - the ploidy of the population. 
m_loss - an indicator for the genome loss model
m_add - an indicator for the genome addition model

Returns:

A matrix of n*(n-1)/2 rows and the same number of columns as bruvo_mat

Notes:
Currently bruvo_mat should not contain NAs. The way to deal with them, as they
are based off of microsatellites, is to replace them with zero. The main
distance algorithm in turn will return a distance of 100 for any individuals
with missing data. In the wrapping R function, 100s will be converted to NAs
and then the average over all loci will be taken. 
 
Tests:
 
All tests for Bruvo's distance are based off of the published values. The 
following R code will test that Bruvo's distance is working properly:
 
# Published data from Bruvo et al 2004 (with extra zeroes)
tg <- structure(c(0L, 0L, 0L, 20L, 20L, 24L, 23L, 26L, 24L, 43L), 
                .Dim = c(2L, 5L), .Dimnames = list(c("1", "2"), NULL))
tg1 <- tg

# Comparison function
amat <- function(a1, a2) 1 - 2^(-abs(a1 - a2))

# Calculate distance between two samples
brvo <- function(allmat){
  n     <- ncol(allmat)
  rows  <- seq(n)
  pairs <- expand.grid(allmat[1, ], allmat[2, ]) 
  mat   <- matrix(apply(pairs, 1, function(x) amat(x[1], x[2])), n, n)
  # Infinite model accounted for here
  mat[allmat[1, ] == "0", ] <- 1
  mat[, allmat[2, ] == "0"] <- 1
  perms <- matrix(.Call("permuto", n, PACKAGE = "poppr") + 1, nrow = n) 
  walks <- apply(perms, 2, function(i) mean(mat[matrix(c(rows, i), ncol = 2)]))
  min(walks)
}

# Calculate distance between two samples with different models
Rbruvo <- function(dat, add = TRUE, loss = TRUE, old_model = FALSE){
  dat <- t(apply(dat, 1, sort))
  dat <- dat[, colSums(dat, na.rm = TRUE) > 0]
  if (!add & !loss){
    return(brvo(dat))
  } 
  nzeroes <- rowSums(dat == 0)
  zeroes  <- which.max(nzeroes)
  repls   <- seq(nzeroes[zeroes])
  zlocat  <- dat[zeroes, ] == 0
  mloss   <- 0
  madd    <- 0
  if (loss){
    donor    <- which.min(nzeroes)
    combs    <- lapply(repls, function(i) dat[donor, ]) 
    combs    <- do.call(expand.grid, combs)
    combs    <- t(apply(combs, 1, sort))
    combs    <- if (old_model) unique(combs) else combs
    combs    <- if (nrow(combs) == 1) t(combs) else combs
    losslist <- vector(mode = "list", length = nrow(combs))
    for (i in seq(nrow(combs))){
      losslist[[i]] <- dat
      losslist[[i]][zeroes, repls] <- combs[i, ]
    }
    mloss <- sum(vapply(losslist, brvo, numeric(1)))/length(losslist)
  }
  if (add){
    combs   <- lapply(repls, function(i) dat[zeroes, !zlocat]) 
    combs   <- do.call(expand.grid, combs)
    combs   <- t(apply(combs, 1, sort))
    combs   <- if (old_model) unique(combs) else combs
    combs   <- if (nrow(combs) == 1) t(combs) else combs
    addlist <- vector(mode = "list", length = nrow(combs))
    for (i in seq(nrow(combs))){
      addlist[[i]] <- dat
      addlist[[i]][zeroes, zlocat] <- combs[i, ]
    }
    madd <- sum(vapply(addlist, brvo, numeric(1)))/length(addlist)
  }
  
  (mloss + madd)/(sum(c(mloss, madd) > 0L))
}

# Function to run poppr's version of Bruvo's distance
requireNamespace("poppr")
poppr_bruvo <- function(dat, add = TRUE, loss = TRUE, old_model = FALSE){
  dat <- t(apply(dat, 1, sort))
  dat <- dat[, colSums(dat, na.rm = TRUE) > 0]
  ploid <- ncol(dat)
  .Call("bruvo_distance",
        dat,
        .Call("permuto", ploid, PACKAGE = "poppr"),
        ploid,
        add,
        loss,
        old_model,
        PACKAGE = "poppr")
}

# Function to run all the models and return a vector
run_models <- function(dat, FUN = "Rbruvo", ...){
  models <- expand.grid(add = 0:1 == 1L, loss = 0:1L == 1L)
  res    <- apply(models, 1, function(x) do.call(FUN, c(list(dat), x, ...)))
  names(res) <- c("infinite", "add", "loss", "combined")
  res
}

# works for published values
all.equal(run_models(tg, "poppr_bruvo"), run_models(tg, "Rbruvo"))

# works with recursion
tg1[1, ] <- c(0, 0, 0, 0, 24)
all.equal(run_models(tg1, "poppr_bruvo"), run_models(tg1, "Rbruvo"))
all.equal(run_models(tg1, "poppr_bruvo", old_model = TRUE), run_models(tg1, "Rbruvo", old_model = TRUE))

tg2 <- tg1
tg2[1, -(1:3)] <- c(102, 104)/2
tg2[2, -1] <- c(104,104,106,110)/2
all.equal(run_models(tg1, "poppr_bruvo"), run_models(tg1, "Rbruvo"))
all.equal(run_models(tg1, "poppr_bruvo", old_model = TRUE), run_models(tg1, "Rbruvo", old_model = TRUE))

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
SEXP bruvo_distance(SEXP bruvo_mat, SEXP permutations, SEXP alleles, SEXP m_add, SEXP m_loss, SEXP old_model)
{
	int rows;   // number of rows
	int cols;   // number of columns
	int ploidy; // maximum ploidy
	int P;      // The number of factorial combinations of alleles.
	int loss;   // indicator for genome loss model.
	int add;    // indicator for genome addition model.
	int* perm;  // pointer to permutation vector.
	int* pmat;  // pointer to pair of samples.
	
	// indices ------------------------------
	int allele; // allele index
	int i;      // sample one
	int j;      // sample two
	int locus;  // index for the first column of a locus
	int clm;    // index for the column shift in the incoming matrix
	int count;  // counter for the output
	
	// R objects ------------------------------
	SEXP Rdim;        // dimensions of the bruvo_mat
	SEXP Rval;        // output vector
	SEXP pair_matrix; // temporary array to store two samples at a single locus
	
	// Initialization ------------------------------
	count = 0;
	P = length(permutations);
	Rdim = getAttrib(bruvo_mat, R_DimSymbol);
	rows = INTEGER(Rdim)[0];
	cols = INTEGER(Rdim)[1];
	ploidy = INTEGER(coerceVector(alleles, INTSXP))[0];
	loss = asLogical(m_loss);
	add = asLogical(m_add);
	PROTECT(bruvo_mat = coerceVector(bruvo_mat, INTSXP));
	perm = INTEGER(coerceVector(permutations, INTSXP));
	PROTECT(Rval = allocMatrix(REALSXP, rows*(rows-1)/2, cols/ploidy));
	PROTECT(pair_matrix = allocVector(INTSXP, 2*ploidy));
	pmat = INTEGER(pair_matrix);
	
	for(locus = 0; locus < cols; locus += ploidy)
	{
		for(i = 0; i < rows - 1; i++)
		{
			R_CheckUserInterrupt(); // in case the user wants to quit
			for(allele = 0; allele < ploidy; allele++) 
			{
				clm = (allele + locus)*rows;
				pmat[allele] = INTEGER(bruvo_mat)[i + clm];
			}
			for(j = i + 1; j < rows; j++)
			{
				for(allele = 0; allele < ploidy ; allele++)
				{
					clm = (allele + locus)*rows;
					pmat[allele + ploidy] = INTEGER(bruvo_mat)[j + clm];
				}
				REAL(Rval)[count++] = bruvo_dist(pmat, &ploidy, perm, &P, &loss, &add, asInteger(old_model));
			}
		}
	}
	UNPROTECT(3); // bruvo_mat; Rval; pair_matrix
	return Rval;
}

/*==============================================================================
================================================================================
*	Internal C Functions
================================================================================
==============================================================================*/

/*	The algorithm for the permutation function is modified from:
*	http://www.geeksforgeeks.org/archives/767 */

/* Function to swap values at two pointers */
void swap (int *x, int *y)
{
	int temp;
	temp = *x;
	*x = *y;
	*y = temp;
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
*	Function to print permutations of string
*		This function takes four parameters:
*		1. String
*		2. Starting index of the string
*		3. Ending index of the string. 
*		4. pointer to array of size n*n! 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void permute(int *a, int i, int n, int *c) 
{
	int j;
	if (i == n)
	{
		/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		*	'a' will be the array containing the numeric sequence to be
		*	permuted. It will be reshuffled into a new pattern each
		*	time it reaches this control structure. To place the value
		*	into the array 'c', the pointer for a needs to be incremented
		*	over all its elements.
		~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
		perm_count += n+1;
		int ind = perm_count;
		for(j = n; j >= 0; j--)
		{
			c[--ind] = *(a + j);
		}
	}
	else
	{
		for (j = i; j <= n; j++)
		{
			swap((a + i), (a + j));
			permute(a, i + 1, n, c);
			swap((a + i), (a + j)); //backtrack
		}
	}
}

/* A factorial function for calculating permutations */
int fact(int x)
{
	int f = 1;
	int u;
	for (u = x; u > 1; u--)
	{
		f *= u;
	}
	return f;
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	Bruvo's Distance implementing Addition, Loss, and Infinite models for
	imputation of missing data.

	NOTE: The input needs to be divided by the repeat length beforehand for this
	to work. 

	in: a matrix of two individuals
	out: a double value that will be the output from bruvo's distance.
	n: number of individuals(2)
	nall: number of alleles
	perm: a vector from the permn function in R
	woo: p * p!
	loss: TRUE/FALSE: impute under genome loss model.
	add: TRUE/FALSE: impute under genome addition model. 
==============================================================================*/
double bruvo_dist(int *in, int *nall, int *perm, int *woo, int *loss, int *add, int old_model)
{
	// R_CheckUserInterrupt();
	int i; 
	int j; 
	int counter = 0;    // counter used for building arrays
	int n = 2;          // number of individuals
	int p = *nall;      // number of alleles
	int w = *woo;       // number of permutations = p * p!
	int loss_indicator = *loss; // 1 if the genome loss model should be used
	int add_indicator = *add;   // 1 if the genome addition model should be used

	int* genos;    // array to store the genotypes
	genos = R_Calloc(2*p, int);

	int zerocatch[2]; 	// 2 element array to store the number of missing 
					  	// alleles in each genotype. 
	 
	int** zero_ind;     // array to store the indices of the missing data for
						// each genotype
	zero_ind = R_Calloc(2, int*);
	zero_ind[0] = R_Calloc(p, int);
	zero_ind[1] = R_Calloc(p, int);

	int zerodiff;      	// used to check the amount of missing data different 
				       	// between the two genotypes

	double minn = 100; 	// The minimum distance 

	// reconstruct the genotype table.
	zerocatch[0] = 0;
	zerocatch[1] = 0;
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < p; j++)
		{
			// Catch missing data here.
			if (in[counter] == 0)
			{
				if (zerocatch[i] == p - 1)
				{
					goto finalsteps;
				}
				zerocatch[i] += 1;
				// Rprintf("[%d]", zerocatch[i] - 1);
				zero_ind[i][zerocatch[i] - 1] = j;
			}
			genos[i*p + j] = in[counter++];
			// Rprintf("%d\t", genos[i*p + j]);
		}
		// Rprintf("\n");
	}
	// Rprintf("\n");
	zerodiff = abs(zerocatch[0] - zerocatch[1]);
	/*==========================================================================
	* Removing superfluous zeroes from the data. This is in the case that both
	* of the genotypes contain one or more zeroes.
	*
	* smaller and larger refer to the size of the genotypes.
	==========================================================================*/
	if (zerocatch[0] > 0 && zerocatch[1] > 0)
	{
		int smaller = 0;
		int larger = 1;
		int reduction = 0; 
		int i;
		int j;
		int zero_counter;
		int *perm_array;
		int *new_geno;
		int *new_alleles;

		if (zerodiff == 0)
		{
			reduction = p - zerocatch[0];
		}
		else
		{
			if (zerocatch[0] < zerocatch[1])
			{
				smaller = 1;
				larger = 0;
			}
			reduction = p - (zerocatch[smaller] - zerodiff);
		}
	
		new_alleles = R_Calloc(reduction, int);
		for (i = 0; i < reduction; i++)
		{
			new_alleles[i] = i;
		}
		w = fact(reduction) * reduction;
		perm_array = R_Calloc(w, int);
		perm_count = 0;
		permute(new_alleles, 0, reduction - 1, perm_array);
		new_geno = R_Calloc(reduction*n, int);
		counter = 0;
		for (i=0; i < n; i++)
		{
			zero_counter = zerocatch[larger];
			for (j = 0; j < p; j++)
			{
				if (genos[i*p + j] == 0 && zero_counter > 0)
				{
					zero_counter--;
				}
				else
				{
					new_geno[counter++] = genos[i*p + j];
				}
			}
		}
	
		minn = bruvo_dist(new_geno, &reduction, perm_array, &w, 
								&loss_indicator, &add_indicator, old_model);
		R_Free(perm_array);
		R_Free(new_geno);
		R_Free(new_alleles);
		goto finalsteps;
	}

	double** dist; // array to store the distance
	int* facts;    // array to store factorials
	int idx_plus;  // i + 1 index for factorial calculation
	dist  = R_Calloc(p, double*);
	facts = R_Calloc(p, int);
	// Construct distance matrix of 1 - 2^{-|x|}.
	// This is constructed column by column. 
	// Genotype 1: COLUMNS
	// Genotype 2: ROWS
	for(i = 0; i < p; i++)
	{
		idx_plus = i + 1;
		facts[i] = (idx_plus == 1) ? idx_plus : idx_plus * facts[i - 1];
	
		dist[i]  = R_Calloc(p, double);
		for(j = 0; j < p; j++)
		{
			dist[i][j] = 1 - pow(2, -abs(genos[0*p + j] - genos[1*p + i]));
		}
	}

	// Rprintf("DISTANCE:");
	// print_distmat(dist, genos, p);

	if (zerocatch[0] > 0 || zerocatch[1] > 0)
	{
		int ind;    
		int miss_ind = 1;
		int z;
		int tracker = 0;
		int loss_tracker = 0;
		int comparison_factor = 1;
		double genome_add_sum = 0;
		double genome_loss_sum = 0;

		if (zerocatch[0] > 0) // The columns contain the zero value
		{
			miss_ind = 0;
		}
		ind = zero_ind[miss_ind][0];

		/*======================================================================
		*	INFINITE MODEL
		*	Infinite model will simply replace the distance of the comparisons
		*	containing the missing allele to 1.
		======================================================================*/
		if(loss_indicator != 1 && add_indicator != 1)
		{
			// Rprintf("TO INFINITY!\n");
			for (z = 0; z < zerocatch[miss_ind]; z++)
			{
				ind = zero_ind[miss_ind][z];
				if (zerocatch[0] == 0)
				{
					for (j = 0; j < p; j++)
					{
						dist[ind][j] = 1;
					}
				}
				else
				{
					for (j = 0; j < p; j++)
					{
						dist[j][ind] = 1;
					}						
				}
			}
			goto finalcalc;
		}
		/*======================================================================
		*	GENOME ADDITION MODEL
		*	Genome Addition model uses the observed values of the short 
		*	genotype for the replacement allele. This is achieved by simply
		*	shifting the columns or rows of the distance matrix and 
		*	recalculating the minimum distance. 
		======================================================================*/

		if (add_indicator == 1)
		{
			int Nobs; // Number of observed alleles in the shorter genotype
			Nobs = p - zerocatch[miss_ind];	
			int* replacements;
			int* replaced_alleles;
			replacements = R_Calloc(Nobs, int);
			replaced_alleles = R_Calloc(zerocatch[miss_ind], int);
			int short_counter = 0;
			// fill the replacements array with the indices of the replacement
			// distances.
			for (i = 0; i < p; i++)
			{
				if (genos[miss_ind*p + i] > 0)
				{
					// replaced_alleles[short_counter] = genos[miss_ind*p + i];
					replacements[short_counter++] = i;
				}
			}	
			// Rprintf("ADD!\n");
			// print_distmat(dist, genos, p);
			for (i = 0; i < Nobs; i++)
			{
				genome_add_calc(genos, w, p, perm, dist, zerocatch[miss_ind],
					zero_ind[miss_ind], 0, miss_ind, replacements, 
					replaced_alleles, facts, Nobs, i, &genome_add_sum, &tracker, old_model);
				// Rprintf("current add sum = %.6f\n", genome_add_sum);
			}
			R_Free(replacements);
			R_Free(replaced_alleles);
		}
		/*======================================================================
		*	GENOME LOSS MODEL
		*	Genome Loss model uses the alleles from the larger genotype to
		*	reconstruct the allelic state of the smaller. This means that
		*	they need to be replaced and passed through the function again.
		======================================================================*/
		if (loss_indicator == 1)
		{
			// Rprintf("LOSS!\n");
			int* loss_replacement;
			loss_replacement = R_Calloc(zerocatch[miss_ind], int);
			for (i = 0; i < p; i++)
			{
				genome_loss_calc(genos, p, perm, w, &loss_indicator, 
					&add_indicator, zero_ind[miss_ind], 0, zerocatch[miss_ind], 
					miss_ind, i, &genome_loss_sum, loss_replacement, facts, 
					&loss_tracker, old_model);
			}
			R_Free(loss_replacement);
		}
		if (tracker == 0)
		{
			tracker = 1;
		}
		if (loss_tracker == 0)
		{
			loss_tracker = 1;
		}
		genome_loss_sum = genome_loss_sum/loss_tracker;
		// Rprintf("LOSS SUM: %.6f\n", genome_loss_sum/p);
		genome_add_sum = genome_add_sum/tracker;
		// Rprintf("ADD SUM: %.6f\n", genome_add_sum/p);
		comparison_factor = loss_indicator + add_indicator;
		minn = (genome_add_sum + genome_loss_sum)/(p*comparison_factor);
	}
	else 
	{
		finalcalc: minn = mindist(w, p, perm, dist)/p;
	}
	for (i = 0; i < p; i++)
	{
		R_Free(dist[i]);
	}
	R_Free(dist);
	R_Free(facts);
finalsteps: 
	R_Free(genos);
	R_Free(zero_ind[0]);
	R_Free(zero_ind[1]);
	R_Free(zero_ind);
	return minn;
}

// Calculate the multinomial coefficient
//
// @param ARR an unordered integer array
// @param n the length of that array
// @param facts an array of factorials that is at least of length n. This will
//        be of the form 1, 2, 6, 24, 20, ... and is equivalent to the following
//        in R: exp(cumsum(log(1:n)))
//
// This will return the multinomial coefficient, which the number of ways an 
// array can be permuted to create unique, ordered arrays. It takes the form of
//
// 	     n!
// 	------------
// 	k1!k2!...km!
// 
// Where n is the number of elements in the array and km is the number of 
// indices of the mth category. An example of the multinomial coefficient in
// action is the number of ways a triploid heterozygous genotype can be permuted:
//
// 	GENOTYPE: A B C
//
//	  3!     6
//	------ = -
//	1!1!1!   1
//
// Or a partial heterozygote:
//
//	GENOTYPE: A B B
//
//	  3!     6
//	------ = -
//	 1!2!    2
//
// In terms of this program, this is used in the genotype addition/loss models
// where we want to sort through all possible combinations of the alleles at 
// the ambiguous genotypes. Because the computation of the distance does not 
// depend on the order, instead of recomputing this (which takes n^2) steps, we
// simply multiply it by this coefficient.
//
int multinomial_coeff(int* ARR, int n, int* facts) {
	// Deal with the simple cases (1 or 2 missing alleles) to avoid higher cost
	// than simply repeating the measure.
	if (n == 1)
	{
		return(1);
	}
	if (n == 2) // diploid expansion 1AA 2AB 1BB
	{
		if (ARR[0] == ARR[1])
		{
			return(1);
		}
		return(2);
	}
	int res;
	int i;
	int count = 1;
	int* TEMP;
	TEMP  = R_Calloc(n, int);
	for (i = 0; i < n; i++)
	{
		TEMP[i] = ARR[i];
	}
	R_qsort_int(TEMP, 1, n); // Note: this takes the arguments of what indices to sort, 1-indexed
	int previous_value = TEMP[0];
	res = count;
	for (i = 1; i < n; i++)
	{
		if (TEMP[i] != previous_value)
		{
			res *= facts[count - 1];
			// Rprintf("%d count: %d (%d)\n", TEMP[i - 1], count, res);
			count = 1;
		}
		else
		{
			count++;
		}
		previous_value = TEMP[i];

	}
	res *= facts[count - 1];
	// Rprintf("%d count: %d (%d)\n", TEMP[n - 1], count, res);
	res = facts[n - 1]/res; // n!/a!b!c!
	R_Free(TEMP);
	return(res);
}

/*==============================================================================
* 	GENOME ADDITION MODEL
*	This will replace the portions of the distance matrix generated via the 
*	missing values with columns (or rows) of observed values in a combinatorial 
*	way. This will involve recursion.
* 
*	Arguments:
*	---------
*	perms - number of permutations (passed to mindist)
*	alleles - number of maximum alleles (also passed to mindist)
*	*perm - permutation array (passed to mindist)
*	**dist - distance matrix to be manipulated (also passed to mindist)
*	
*	zeroes - number of zeroes present in the shorter genotype.
*	*zero_ind - array containing indices of the zero values of the short geno
*	curr_zero - the index of the current zero index for zero_ind
*	miss_ind - the index for the genotype with missing data. Necessary for
*				determining rows or columns
*	*replacement - array containing indices of replacement genotypes.
*   *replaced_alleles - array to contain the replacements
*	inds - number of replacement genotypes.
*	curr_ind - the index of the current replacement genotype.
*
*	*genome_add_sum - pointer to the total number of the genome addition model.
*	*tracker - pointer to counter for the number of calculations for addition. 
*	old_model - a binary indicator of whether or not to calculate the old model
*	            previously implemented. 
==============================================================================*/
void genome_add_calc(int* genos,
	int perms,
	int alleles,
	int* perm,
	double** dist,
	int zeroes,
	int* zero_ind,
	int curr_zero,
	int miss_ind,
	int* replacement,
	int* replaced_alleles,
	int* facts,
	int inds,
	int curr_ind,
	double* genome_add_sum,
	int* tracker,
	int old_model)
{
	// R_CheckUserInterrupt();
	int i;
	int j;
	int mult;
	//==========================================================================
	// Part 1: fill one row/column of the matrix.
	//==========================================================================
	if(miss_ind > 0)
	{
		for (j = 0; j < alleles; j++)
		{
			dist[zero_ind[curr_zero]][j] = dist[replacement[curr_ind]][j];
		}
	}
	else
	{
		for (j = 0; j < alleles; j++)
		{
			dist[j][zero_ind[curr_zero]] = dist[j][replacement[curr_ind]];	
		}
	}
	// Here we are gathering the array of replaced alleles for the binomial
	// expansion calculation. 
	replaced_alleles[curr_zero] = curr_ind;
	//==========================================================================
	// Part 2: Iterate through the rest of the possible combinations.
	//
	// The first for loop iterates through all possible individuals.
	// The first if loop will check if there are any more slots to be filled.
	// if there aren't, then the minimum distance will be calculated on the
	// matrix as it stands and then the sum will be returned. 
	//==========================================================================
	for (i = curr_ind; i < inds; i++)
	{
		if (curr_zero < zeroes - 1)
		{
			genome_add_calc(genos, perms, alleles, perm, dist, zeroes, zero_ind, 
				++curr_zero, miss_ind, replacement, replaced_alleles, facts, inds, i, genome_add_sum, 
				tracker, old_model);
			if (curr_zero == zeroes - 1)
			{
				return;
			}
		}
		else
		{
			// for (int z = 0; z < zeroes; z++)
			// {
			// 	Rprintf("%d ", replaced_alleles[z]);
			// }
			mult = (old_model) ? 1 : multinomial_coeff(replaced_alleles, zeroes, facts);
			// Rprintf("Multiplier: %d\n", mult);
			*genome_add_sum += mindist(perms, alleles, perm, dist) * mult;
			*tracker += 1 * mult;
			if (zeroes == 1 || i == inds - 1)
			{
				return;
			}
		}
		curr_zero--;
		
	}
	return;
}

/*==============================================================================
*	Genome Loss Model
*	Replace missing alleles in the shorter genotype with all possible
*	combinations of alleles in the larger genotype and recall bruvo_dist. There
*	are choose((n+k-1), k) possible combinations where n is the number of
*	alleles in the larger genotype and k is the number of missing alleles in the
*	shorter genotype. 
*	
*	Arguments:
*	---------
*	PASSED TO BRUVO_DIST:
*	*genos - genotype array
*	nalleles - number of maximum alleles.
*	*perm_array - permutation array.
*	*woo - nalleles * nalleles!
*	*loss - genome loss model indicator.
*	*add - genome addition indicator.
*
*	UNIQUE TO THIS FUNCTION:
*	*zero_ind - array containing indices of the zero values of the short geno
*	curr_zero - the index of the current zero index for zero_ind
*	zeroes - number of zeroes present in the shorter genotype.
*	miss_ind - the index for the genotype with missing data.
*	curr_allele - the current index for the replacement alleles of the full geno
*
*	*genome_loss_sum - pointer to the total number of the genome loss model.
*	*replacements - a vector containing the replacement alleles for multinomial
*	                coefficient calculation
*	*facts - a vector containing the factorial calculations for the multinomial
*	*loss_tracker - pointer to counter for the number of calculations for loss.
*	old_model - a binary indicator of whether or not to calculate the old model
*	            previously implemented. 
==============================================================================*/
void genome_loss_calc(int *genos, int nalleles, int *perm_array, int woo, 
		int *loss, int *add, int *zero_ind, int curr_zero, int zeroes, 
		int miss_ind, int curr_allele, double *genome_loss_sum, 
		int* replacements, int* facts, int *loss_tracker, int old_model)
{
	// R_CheckUserInterrupt();
	int i; 
	int full_ind;
	full_ind = 1 + (0 - miss_ind);
	int to_replace = miss_ind*nalleles;
	int donor      = full_ind*nalleles;
	int mult;

	genos[to_replace + zero_ind[curr_zero]] = genos[donor + curr_allele];
	replacements[curr_zero] = curr_allele; 

	for (i = curr_allele; i < nalleles; i++)
	{
		if (curr_zero < zeroes - 1)
		{
			genome_loss_calc(genos, nalleles, perm_array, woo, loss, add, 
				zero_ind, ++curr_zero, zeroes, miss_ind, i, genome_loss_sum, 
				replacements, facts, loss_tracker, old_model);
			if (curr_zero == zeroes - 1)
			{
				return;
			}
		}
		else
		{
			mult = (old_model) ? 1 : multinomial_coeff(replacements, zeroes, facts);

			// for (int z = 0; z < zeroes; z++)
			// {
			// 	Rprintf("%d ", replacements[z]);
			// }
			// Rprintf("\t Multiplier: %d\n", mult);

			*genome_loss_sum += bruvo_dist(genos, &nalleles, perm_array, 
				&woo, loss, add, old_model)*nalleles*mult;
			*loss_tracker += 1 * mult;
			if (zeroes == 1 || i == nalleles - 1)
			{
				return;
			}
		}
		curr_zero--;
	}
	return;
}

/*==============================================================================
 * 
 * This function is currently unused.
 * 
 * 
 * Notes for fill_short_geno: This will act much in the same way as
 * genome_loss_calc, except it will fill the shorter genotype with all possible
 * combinations of that genotype before sending it through bruvo_dist with
 * one full genotype. 
 *
 * Things that need to be set before running this:
 * - replacement is an array of the non-missing alleles from the shorter 
 *   genotype.
 * - inds is the number of non-missing alleles.
 * - *res will be minn
 * - *tracker will count the number of iterations this goes through in order
 *   to get an average. 
==============================================================================
void fill_short_geno(int *genos, int nalleles, int *perm_array, int *woo, 
		int *loss, int *add, int zeroes, int *zero_ind, int curr_zero, 
		int miss_ind, int *replacement, int inds, int curr_ind, double *res, 
		int *tracker)
{
	// R_CheckUserInterrupt();
	int i; //full_ind;
	genos[miss_ind*nalleles + zero_ind[curr_zero]] = 
		genos[miss_ind*nalleles + replacement[curr_ind]];
	for (i = curr_ind; i < inds; i++)
	{
		if (curr_zero < zeroes - 1)
		{
			fill_short_geno(genos, nalleles, perm_array, woo, loss, add, zeroes, 
				zero_ind, ++curr_zero, miss_ind, replacement, inds, i, res, 
				tracker);
			if (curr_zero == zeroes - 1)
			{
				return;
			}
		}
		else
		{
			*res += bruvo_dist(genos, &nalleles, perm_array, woo, loss, 
						add);
			*tracker += 1;
			if (zeroes == 1 || i == nalleles - 1)
			{
				return;
			}
		}
		curr_zero--;
	}
	return;
}
*/


/*
// Multiset coefficient: fact(n+k-1)/(fact(k)*fact(n-1))
*/

double mindist(int perms, int alleles, int *perm, double **dist)
{
	int i, j;
	int w = perms;
	int p = alleles;
	int counter = 0;
	double res = 0;
	double minn = 100;

	for(i = 0; i < w; i += p)
	{
		for(j = 0; j < p; j++)
		{
			if (j == 0)
			{
				res = dist[*(perm + counter)][j];
				counter++;
				if(res > minn)
				{
					j = p;
					counter = i + w/p;
					i = counter;
				}				
			}
			else
			{
				res += dist[*(perm + counter)][j];
				counter++;
				if(j < p-1 && res > minn)
				{				
					counter += (p-j-1);
					j = p;
				}
			}
		}
		/*	Checking if the new calculated distance is smaller than the smallest
		 distance seen. */
		if ( res < minn )
		{
			minn = res;
		}
	}
	return minn;
}

/* Helper function to print a distance matrix
 * mainly used for debugging

void print_distmat(double** dist, int* genos, int p)
{
	int i;
	int j;

	Rprintf("\n \t");
	for (i = 0; i < p; i++)
	{
		Rprintf("%d\t", genos[i]);
	}
	Rprintf("\n");
	for (i = 0; i < p; i++)
	{
		Rprintf("%d\t", genos[i+p]);
		for (j = 0; j < p; j++)
		{
			Rprintf("%.4f\t", dist[i][j]);
		}
		Rprintf("\n");
	}
	return;
}

*/

