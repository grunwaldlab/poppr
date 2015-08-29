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
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>

// Include openMP if the compiler supports it
#ifdef _OPENMP
#include <omp.h>
#endif

// Assumptions:
//  All genotypes have the same number of SNPs available.
//  All SNPs are diploid or haploid, depending on the function.
//  All SNPs are binary, stored as single bits in hexadecimal characters.

/*

Zygosity struct
===============

A struct used for storing information about sample genotype in 8 locus chunks.

*/
struct zygosity
{
  char c1;  // An 8 locus fragment of one chromosome 
  char c2;  // The corresponding chromosomal pairing of the 8 loci
            // 1 : dominant
            // 0 : recessive

  // The following indicate presence or absence of 
  char ch;  // Heterozygous sites
  char cd;  // Homozygous dominant sites
  char cr;  // Homozygous recessive sites
};

/*

Locus struct
============

A struct used for tallying the total number of dominant and recessive alleles,
heterozygous sites, and total number of genotypes sampled across a set of
samples.

*/

struct locus
{
  double d;  // Number of dominant alleles found at this locus across genotypes
  double r;  // Number of recessive alleles found at this locus across genotypes
  double h;  // Number of genotypes that are heterozygous at this locus
  double n;  // Number of genotypes that contributed data to this struct  
};


SEXP bitwise_distance_haploid(SEXP genlight, SEXP missing, SEXP requested_threads);
SEXP bitwise_distance_diploid(SEXP genlight, SEXP missing, SEXP differences_only, SEXP requested_threads);
SEXP association_index_haploid(SEXP genlight, SEXP missing, SEXP requested_threads);
SEXP association_index_diploid(SEXP genlight, SEXP missing, SEXP differences_only, SEXP requested_threads);
SEXP get_pgen_matrix_genind(SEXP genind, SEXP freqs, SEXP pops, SEXP npop);
SEXP get_pgen_matrix_genlight(SEXP genlight, SEXP window);
void fill_Pgen(double *pgen, struct locus *loci, int interval, SEXP genlight);
void fill_loci(struct locus *loc, SEXP genlight);
void fill_zygosity(struct zygosity *ind);
char get_similarity_set(struct zygosity *ind1, struct zygosity *ind2);
int get_zeros(char sim_set);
int get_difference(struct zygosity *z1, struct zygosity *z2);
int get_distance(struct zygosity *z1, struct zygosity *z2);
int get_distance_custom(char sim_set, struct zygosity *z1, struct zygosity *z2);

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Calculates the pairwise differences between samples in a genlight object. The
distances represent the number of sites between individuals which differ.

Input: A genlight object containing samples of haploids.
       A boolean representing whether missing data should match (TRUE) or not.
       An integer representing the number of threads that should be used.
Output: A distance matrix representing the number of differences between each sample.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
SEXP bitwise_distance_haploid(SEXP genlight, SEXP missing, SEXP requested_threads)
{
  // This function calculates the raw genetic distance between samples in 
  // a genlight object. The general flow of this function is as follows:
    // Define and initialize variables
    // Retrieve data from R objects passed in as arguments
    // Prepare for multithreading if compiled to do so
    // Initialize variables needed for the next loop
    // Loop through every genotype/sample in the genlight object, call each on i:
      // Retrieve genetic data and information for sample i from R objects
      // initialize multi threading for the next loop, if compiled to do so
      // Loop through every genotype up to and not including i to cover all pairings:
        // Retrieve data for sample j
        // Prepare to correct for missing data in both i and j
        // Loop through each chunk of 8 loci of both i and j, call each chunk k:
          // Find the locations in chunk k in which i and j are the same
          // Correct for missing data by forcing a match (or not match)
          // Update the output matrix with the distances found in this chunk.
    // Fill the final R return object and return it.


  SEXP R_out;               // output matrix (n x n)
  SEXP R_gen_symbol;        // gen slot
  SEXP R_chr_symbol;        // snp slot
  SEXP R_nap_symbol;        // NA.posi slot
  SEXP R_gen;               // 
  int num_gens;             // number of genotypes
  SEXP R_chr1_1;
  SEXP R_chr2_1;
  SEXP R_nap1;
  SEXP R_nap2;
  int nap1_length;
  int nap2_length;
  int chr_length;
  int next_missing_index_i;
  int next_missing_index_j;
  int next_missing_i;
  int next_missing_j;
  int missing_match;
  int num_threads;
  char mask;
  int i;
  int j;
  int k;

  int** distance_matrix;
  int cur_distance;
  char tmp_sim_set;

  // These variables and function calls are used to access elements of the
  // genlight object. ie, R_gen_symbol is being set up as an equivalent to the
  // @gen accessor for genlights.
   
  PROTECT(R_gen_symbol = install("gen")); // Used for accessing the named
                                          // elements of the genlight object
  PROTECT(R_chr_symbol = install("snp"));
  PROTECT(R_nap_symbol = install("NA.posi"));  

  // This will be a LIST of type LIST:RAW
  // Set R_gen to genlight@gen, a vector of genotypes in the genlight object 
  // stored as SNPbin objects.
  R_gen = getAttrib(genlight, R_gen_symbol); 
  
  // Set num_gens to contain the total number of samples/genotypes within the
  // genlight object
  num_gens = XLENGTH(R_gen);  

  // Set up and initialize the matrix for storing total distance between each
  // pair of genotypes
  PROTECT(R_out = allocVector(INTSXP, num_gens*num_gens));
  distance_matrix = R_Calloc(num_gens,int*);
  for(i = 0; i < num_gens; i++)
  {
    distance_matrix[i] = R_Calloc(num_gens,int);
  }

  #ifdef _OPENMP
  {
    // Set the number of threads to be used in each omp parallel region
    if(INTEGER(requested_threads)[0] == 0)
    {
      num_threads = omp_get_max_threads();
    }
    else
    {
      num_threads = INTEGER(requested_threads)[0];
    }
    omp_set_num_threads(num_threads);
  }
  #else
  {
    num_threads = 1;
  }
  #endif

  next_missing_index_i = -1;
  next_missing_index_j = -1;
  next_missing_i = -1;
  next_missing_j = -1;
  mask = 0;
  tmp_sim_set = 0;
  nap1_length = 0;
  nap2_length = 0;
  chr_length = 0;
  missing_match = asLogical(missing);

  // Loop through every genotype 
  for(i = 0; i < num_gens; i++)
  {
    // Set R_chr1_1 to be genlight@gen[[i]]@snp[[1]], aka a raw list
    // representing the entire first set of chromosomes in this genotype
    R_chr1_1 = VECTOR_ELT(getAttrib(VECTOR_ELT(R_gen,i),R_chr_symbol),0);

    // Set chr_length to represent the number of 8 locus chunks in this genotype
    chr_length = XLENGTH(R_chr1_1);

    // Set R_nap1 to be genlight@gen[[i]]@NA.posi aka a vector of integers each
    // representing the location of a locus with missing data
    R_nap1 = getAttrib(VECTOR_ELT(R_gen,i),R_nap_symbol);

    // Set nap1_length to be the number of loci containing missing data in this
    // sample
    
    nap1_length = XLENGTH(R_nap1);
    
    // Loop through every other genotype

    // The threading code splits this loop into pieces and divides it between
    // all created threads Care must be taken to ensure that no "shared" or
    // unlisted variables are being written to by more than one thread at a
    // time. Private variables can be accessed without worry, but have overhead
    // to create for each thread.

    #ifdef _OPENMP
    #pragma omp parallel for schedule(guided) \
      private(j,cur_distance,R_chr2_1,R_nap2,next_missing_index_j,next_missing_j,next_missing_index_i,next_missing_i,\
              tmp_sim_set, k, mask, nap2_length) \
      shared(R_nap1, nap1_length, i, distance_matrix)
    #endif

    for(j = 0; j < i; j++)
    {
      cur_distance = 0;
      // These will be arrays of type RAW
      // Set R_chr2_1 to be genlight@gen[[j]]@snp[[1]], aka a raw list 
      // representing the entire first set of chromosomes in this genotype
      R_chr2_1 = VECTOR_ELT(getAttrib(VECTOR_ELT(R_gen,j),R_chr_symbol),0);

      // Set R_nap2 to be genlight@gen[[j]]@NA.posi aka a vector of integers 
      // each representing the location of a locus with missing data
      R_nap2 = getAttrib(VECTOR_ELT(R_gen,j),R_nap_symbol);

      // Set nap2_length to be the number of loci containing missing data in 
      // this sample
      nap2_length = XLENGTH(R_nap2);

      // Set up the initial values for tracking missing data
      if(nap2_length > 0)
      { 
        // New sample j with missing, so we start tracking at the first missing 
        // allele
        next_missing_index_j = 0;

        // Get the index of the first missing allele, compensating for R's base 
        // 1 indexing.
        next_missing_j = (int)INTEGER(R_nap2)[next_missing_index_j] - 1;
      }
      else
      {
        next_missing_index_j = -1; 
        next_missing_j = -1;
      }
      if(nap1_length > 0)
      {
        // Since we'll be going back over the first genotype again, we need to 
        // go back through its missing alleles
        next_missing_index_i = 0;
        next_missing_i = (int)INTEGER(R_nap1)[next_missing_index_i] - 1; 
      }
      else
      {
        next_missing_index_i = -1;
        next_missing_i = -1;
      }
      // Finally, loop through all the chunks of SNPs for these genotypes
      for(k = 0; k < chr_length; k++) 
      {
        // Find the sites where the two individuals differ. 1's for match, 0's
        // for different 
        // ~ is the bitwise unary NOT operator, flipping all 1's
        // to 0's and vice versa 

        // ^ is the bitwise binary XOR operator, comparing two sets of bits and
        // making a new set with 1's everywhere the two sets were different, and
        // 0's where they were the same So this sets tmp_sim_set to have 1's at
        // any locus in which sample i and j were the same

        tmp_sim_set = ~((char)RAW(R_chr1_1)[k] ^ (char)RAW(R_chr2_1)[k]);

        // Check for missing values and force them to match
        while(next_missing_index_i < nap1_length && next_missing_i < (k+1)*8 && next_missing_i >= (k*8) )
        {
          // Handle missing bit in both samples with a mask
          
          // A mask is a set of bits designed specifically for use with bitwise
          // operators in order to get a desired result. In this case, we are
          // making an 8 bit mask with 1's wherever data is missing, and 0's
          // everywhere else.

          // << is a binary bitwise operator that shifts all the bits in the
          // left hand side to the left by the number of bits specified by the
          // right hand side. In this case we take 0000001 and move the 1 to the
          // left until it overlaps the next missing data, which, since we are
          // processing this in 8 loci chunks, will be next_missing_i modulo 8.

          mask = 1 << (next_missing_i%8); 

          if(missing_match)
          {
            // |= is a binary bitwise operator that performs bitwise OR and
            // stores it in the left hand side.

            // Using OR with a mask forces all bits that are 1's in the mask to
            // be 1's in the final product as well, regardless of what it was in
            // the original. In other words, this overwrites tmp_sim_set to show
            // a match anywhere the mask has a 1, which is anywhere we found
            // missing data.

            tmp_sim_set |= mask; // Force the missing bit to match
          }
          else
          {
            // &= is the binary bitwise operator that performs bitwise AND and
            // stores it in the left hand side.

            // Using AND with a mask forces all bits that are 0's in the mask to
            // be 0's in, the final product as well, regardless of what it was
            // in the original. In this case we are negating mask, so ~mask has
            // 0's everywhere it used to have 1's, ie where missing data was
            // found, and thus forces tmp_sim_set to show a difference wherever
            // missing data was found.

            tmp_sim_set &= ~mask; // Force the missing bit to not match
          }

          // Find the index of the next missing value in sample i.
          next_missing_index_i++; 
          next_missing_i = (int)INTEGER(R_nap1)[next_missing_index_i] - 1;
        }       

        // Repeat for j
        while(next_missing_index_j < nap2_length && next_missing_j < (k+1)*8 && next_missing_j >= (k*8))
        {
          // Handle missing bit in both samples with mask
          // See above for detailed comments on this process.
          mask = 1 << (next_missing_j%8);

          if(missing_match)
          {
            tmp_sim_set |= mask; // Force the missing bit to match
          }
          else
          {
            tmp_sim_set &= ~mask; // Force the missing bit to not match
          }

          next_missing_index_j++;
          next_missing_j = (int)INTEGER(R_nap2)[next_missing_index_j] - 1;
        }       

        // Add the distance from these 8 loci into the total between these two 
        // genotypes
        cur_distance += get_zeros(tmp_sim_set);
      }

      // Store the distance between these two genotypes in the distance matrix
      // Note that this could be a conflict between threads since
      // distance_matrix is shared However, since each iteration of this loop
      // will have a different value for j and the same value for i, no two
      // threads will ever have the same (i,j) combination, nor will any threads
      // (i,j) be another threads (j,i), since j < i for all threads.
      
      distance_matrix[i][j] = cur_distance;
      distance_matrix[j][i] = cur_distance;
    } // End parallel
  } 

  // Fill the output matrix
  for(i = 0; i < num_gens; i++)
  {
    for(j = 0; j < num_gens; j++)
    {
      INTEGER(R_out)[i + j*num_gens] = distance_matrix[i][j];
    }
  }

  for(i = 0; i < num_gens; i++)
  {
    R_Free(distance_matrix[i]);
  }
  R_Free(distance_matrix);
  UNPROTECT(4); 
  return R_out;
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Calculates the pairwise differences between samples in a genlight object. The
distances represent the number of sites between individuals which differ in
zygosity (if differences_only is TRUE) or the number of differing alleles (0 for
matching zygosity, 1 for the distance between a heterozygote and a homozygote,
and 2 for the distance between differing homozygotes).

Input: A genlight object containing samples of diploids.
       A boolean representing whether missing data should match (TRUE) or not.
       A boolean representing whether distance (FALSE) or differences (TRUE)
          should be returned.
       An integer representing the number of threads that should be used.
Output: A distance matrix representing the distance between each sample in the
          genlight object.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
SEXP bitwise_distance_diploid(SEXP genlight, SEXP missing, SEXP differences_only, SEXP requested_threads)
{
  // This function calculates the raw genetic distance between samples in 
  // a genlight object. The general flow of this function is as follows:
    // Define and initialize variables
    // Retrieve data from R objects passed in as arguments
    // Prepare for multithreading if compiled to do so
    // Initialize variables needed for the next loop
    // Loop through every genotype/sample in the genlight object, call each on i:
      // Retrieve genetic data and information on both chromosomes for sample i from R objects
      // initialize multi threading for the next loop, if compiled to do so
      // Loop through every genotype up to and not including i to cover all pairings:
        // Retrieve data for both chromosomes of sample j
        // Prepare to correct for missing data in both i and j
        // Loop through each chunk of 8 loci of both i and j, call each chunk k:
          // Use the zygosity struct to find the locations in chunk k in which i and j are the same
          // Correct for missing data by forcing a match (or not match)
          // Update the output matrix with the distances found in this chunk.
    // Fill the final R return object and return it.

  SEXP R_out;
  SEXP R_gen_symbol;
  SEXP R_chr_symbol;  
  SEXP R_nap_symbol;
  SEXP R_gen;
  int num_gens;
  SEXP R_chr1_1;
  SEXP R_chr1_2;
  SEXP R_chr2_1;
  SEXP R_chr2_2;
  SEXP R_nap1;
  SEXP R_nap2;
  int nap1_length;
  int nap2_length;
  int chr_length;
  int next_missing_index_i;
  int next_missing_index_j;
  int next_missing_i;
  int next_missing_j;
  int missing_match;
  int only_differences;
  int num_threads;
  char mask;
  struct zygosity set_1;
  struct zygosity set_2;
  int i;
  int j;
  int k;

  int** distance_matrix;
  int cur_distance;
  char tmp_sim_set;

  // These variables and function calls are used to access elements of the genlight object.
  // ie, R_gen_symbol is being set up as an equivalent to the @gen accessor for genlights.  
  PROTECT(R_gen_symbol = install("gen")); // Used for accessing the named elements of the genlight object
  PROTECT(R_chr_symbol = install("snp"));
  PROTECT(R_nap_symbol = install("NA.posi"));  

  // This will be a LIST of type LIST:RAW
  // Set R_gen to genlight@gen, a vector of genotypes in the genlight object stored as SNPbin objects.
  R_gen = getAttrib(genlight, R_gen_symbol); 
  // Set num_gens to contain the total number of samples/genotypes within the genlight object
  num_gens = XLENGTH(R_gen);

  // Set up and initialize the matrix for storing total distance between each pair of genotypes
  PROTECT(R_out = allocVector(INTSXP, num_gens*num_gens));
  distance_matrix = R_Calloc(num_gens,int*);
  for(i = 0; i < num_gens; i++)
  {
    distance_matrix[i] = R_Calloc(num_gens,int);
  }

  #ifdef _OPENMP
  {
    // Set the number of threads to be used in each omp parallel region
    if(INTEGER(requested_threads)[0] == 0)
    {
      num_threads = omp_get_max_threads();
    }
    else
    {
      num_threads = INTEGER(requested_threads)[0];
    }
    omp_set_num_threads(num_threads);
  }
  #else
  {
    num_threads = 1;
  }
  #endif

  next_missing_index_i = 0;
  next_missing_index_j = 0;
  next_missing_i = -1;
  next_missing_j = -1;
  mask = 0;
  tmp_sim_set = 0;
  nap1_length = 0;
  nap2_length = 0;
  chr_length = 0;
  missing_match = asLogical(missing);
  only_differences = asLogical(differences_only);

  // Loop through every genotype 
  for(i = 0; i < num_gens; i++)
  {
    // Set R_chr1_1 to be genlight@gen[[i]]@snp[[1]], 
    // aka a raw list representing the entire first set of chromosomes in this genotype
    R_chr1_1 = VECTOR_ELT(getAttrib(VECTOR_ELT(R_gen,i),R_chr_symbol),0); 
    // Set R_chr1_2 to be genlight@gen[[i]]@snp[[2]], 
    // aka a raw list representing the entire second set of chromosomes in this genotype
    R_chr1_2 = VECTOR_ELT(getAttrib(VECTOR_ELT(R_gen,i),R_chr_symbol),1); 
    // Set chr_length to represent the number of chunks (8 loci) in this genotype
    // Note that the length of R_chr*_* should all be identical.
    chr_length = XLENGTH(R_chr1_1);
    // Set R_nap1 to be genlight@gen[[i]]@NA.posi
    // aka a vector of integers each representing the location of a locus with missing data
    R_nap1 = getAttrib(VECTOR_ELT(R_gen,i),R_nap_symbol); // Vector of the indices of missing values 
    // Set nap1_length to be the number of loci containing missing data in this sample
    nap1_length = XLENGTH(R_nap1);
    // Loop through every other genotype
    // The threading code splits this loop into pieces and divides it between all created threads
    // Care must be taken to ensure that no "shared" or unlisted variables are being written to
    // by more than one thread at a time. Private variables can be accessed without worry, but have
    // overhead to create for each thread.
    #ifdef _OPENMP
    #pragma omp parallel for schedule(guided) \
      private(j,cur_distance,R_chr2_1,R_chr2_2,R_nap2,next_missing_index_j,next_missing_j,next_missing_index_i,next_missing_i,\
              set_1,set_2,tmp_sim_set, k, mask, nap2_length) \
      shared(R_nap1, nap1_length, i, distance_matrix)
    #endif
    for(j = 0; j < i; j++)
    {
      cur_distance = 0;
      // These will be arrays of type RAW
      // Set R_chr2_1 to be genlight@gen[[j]]@snp[[1]], 
      // aka a raw list representing the entire first set of chromosomes in this genotype
      R_chr2_1 = VECTOR_ELT(getAttrib(VECTOR_ELT(R_gen,j),R_chr_symbol),0); // Set 1
      // Set R_chr2_2 to be genlight@gen[[j]]@snp[[2]], 
      // aka a raw list representing the entire second set of chromosomes in this genotype
      R_chr2_2 = VECTOR_ELT(getAttrib(VECTOR_ELT(R_gen,j),R_chr_symbol),1); // Set 2
      // Set R_nap2 to be genlight@gen[[j]]@NA.posi
      // aka a vector of integers each representing the location of a locus with missing data
      R_nap2 = getAttrib(VECTOR_ELT(R_gen,j),R_nap_symbol); // Vector of the indices of missing values
      // Set nap2_length to be the number of loci containing missing data in this sample
      nap2_length = XLENGTH(R_nap2);

      // Set up the initial values for tracking missing data
      if(nap2_length > 0)
      {
        // New sample j with missing, so we start tracking at the first missing allele
        next_missing_index_j = 0; 
        // Get the index of the first missing allele, compensating for R's base 1 indexing.
        next_missing_j = (int)INTEGER(R_nap2)[next_missing_index_j] - 1; 
      }
      else
      {
        next_missing_index_j = -1;
        next_missing_j = -1;
      }
      if(nap1_length > 0)
      {
        // Since we'll be going back over the first genotype again, we need to go back through its missing alleles
        next_missing_index_i = 0;
        next_missing_i = (int)INTEGER(R_nap1)[next_missing_index_i] - 1; //Compensate for Rs 1 based indexing
      }
      else
      {
        next_missing_index_i = -1;
        next_missing_i = -1;
      }
      // Finally, loop through all the chunks of SNPs for these genotypes
      for(k = 0; k < chr_length; k++) 
      {
        // Store both sets of the next 8 bits (the k*8th through k*8+7th loci) of sample i into set_1
        set_1.c1 = (char)RAW(R_chr1_1)[k];
        set_1.c2 = (char)RAW(R_chr1_2)[k];
        // Then call the fill_zygosity function to calculate each locus's zygosity
        fill_zygosity(&set_1);

        // Store both sets of the next 8 bits (the k*8th through k*8+7th loci) of sample j into set_1
        set_2.c1 = (char)RAW(R_chr2_1)[k];
        set_2.c2 = (char)RAW(R_chr2_2)[k];
        // Then call the fill_zygosity function to calculate each locus's zygosity
        fill_zygosity(&set_2);

        // get_similarity_set returns 8 bits with 1's where the two sets shared zygosity and 0's elsewhere
        tmp_sim_set = get_similarity_set(&set_1,&set_2);

        // Check for missing values and force them to match
        while(next_missing_index_i < nap1_length && next_missing_i < (k+1)*8 && next_missing_i >= (k*8) )
        {
          // Handle missing bit in both sets and both samples with a mask
          // A mask is a set of bits designed specifically for use with bitwise operators
            // in order to get a desired result. In this case, we are making an 8 bit mask
            // with 1's wherever data is missing, and 0's everywhere else.
          // << is a binary bitwise operator that shifts all the bits in the left hand side
            // to the left by the number of bits specified by the right hand side. In this
            // case we take 0000001 and move the 1 to the left until it overlaps the next
            // missing data, which, since we are processing this in 8 loci chunks, will be
            // next_missing_i modulo 8.
          mask = 1 << (next_missing_i%8); 
          if(missing_match)
          {
            // |= is a binary bitwise operator that performs bitwise OR and stores it in
              // the left hand side.
            // Using OR with a mask forces all bits that are 1's in the mask to be 1's in
              // the final product as well, regardless of what it was in the original.
              // In other words, this overwrites tmp_sim_set to show a match anywhere
              // the mask has a 1, which is anywhere we found missing data.
            tmp_sim_set |= mask; // Force the missing bit to match
          }
          else
          {
            // &= is the binary bitwise operator that performs bitwise AND and stores it in
              // the left hand side.
            // Using AND with a mask forces all bits that are 0's in the mask to be 0's in,
              // the final product as well, regardless of what it was in the original.
              // In this case we are negating mask, so ~mask has 0's everywhere it used to
              // have 1's, ie where missing data was found, and thus forces tmp_sim_set to
              // show a difference wherever missing data was found.
            tmp_sim_set &= ~mask; // Force the missing bit to not match
          }
          // Find the index of the next missing value in sample i.
          next_missing_index_i++; 
          next_missing_i = (int)INTEGER(R_nap1)[next_missing_index_i] - 1;
        }       
        // Repeat for j
        while(next_missing_index_j < nap2_length && next_missing_j < (k+1)*8 && next_missing_j >= (k*8))
        {
          // Handle missing bit in both sets and both samples with a mask
          // See above for detailed comments on this process.
          mask = 1 << (next_missing_j%8);
          if(missing_match)
          {
            tmp_sim_set |= mask; // Force the missing bit to match
          }
          else
          {
            tmp_sim_set &= ~mask; // Force the missing bit to not match
          }
          next_missing_index_j++;
          next_missing_j = (int)INTEGER(R_nap2)[next_missing_index_j] - 1;
        }       

        // Add the distance from this word into the total between these two genotypes
        if(only_differences)
        {
          // If we only want to count sites with differing zygosity, we just need to
          // count the zeros in the similarity set and add that to the total distance.
          cur_distance += get_zeros(tmp_sim_set);
        }
        else
        {
          // If we want the actual distance, we need a more complicated analysis of 
          // the sets. See get_distance_custom for details.
          cur_distance += get_distance_custom(tmp_sim_set,&set_1,&set_2);
        }
      }
      // Store the distance between these two genotypes in the distance matrix
      // Note that this could be a conflict between threads since distance_matrix is shared
      // However, since each iteration of this loop will have a different value for j and the
      // same value for i, no two threads will ever have the same (i,j) combination, nor will
      // any threads (i,j) be another threads (j,i), since j < i for all threads.
      distance_matrix[i][j] = cur_distance;
      distance_matrix[j][i] = cur_distance;
    } // End parallel
  } 

  // Fill the output matrix
  for(i = 0; i < num_gens; i++)
  {
    for(j = 0; j < num_gens; j++)
    {
      INTEGER(R_out)[i + j*num_gens] = distance_matrix[i][j];
    }
  }

  for(i = 0; i < num_gens; i++)
  {
    R_Free(distance_matrix[i]);
  }
  R_Free(distance_matrix);
  UNPROTECT(4); 
  return R_out;
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Calculates the index of association of a genlight object of haploids.

Input: A genlight object containing samples of diploids.
       A boolean representing whether or not missing values should match. 
       An integer representing the number of threads to be used.
Output: The index of association for this genlight object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
SEXP association_index_haploid(SEXP genlight, SEXP missing, SEXP requested_threads)
{
  // This function calculates the index of association for samples in 
  // a genlight object. The general flow of this function is as follows:
    // Define and initialize variables
    // Retrieve data from R objects passed in as arguments
    // Fill chunk_matrix with all the genetic data in the genlight object for quick access
    // Prepare for multithreading if compiled to do so
    // Initialize variables needed for the next loop
    // initialize multi threading for the next loop, if compiled to do so
    // Loop through every chunk of 8 loci in the genlight object, call each on i:
      // Loop through every genotype in the genlight object, call each one j:
      // Prepare the genetic data for j stored in chunk_matrix for comparison
      // Prepare to correct for missing data in j
      // Loop through every genotype from j on to cover all pairings, call this k:
        // Prepare the genetic data for k stored in chunk_matrix for comparison
        // Find locations at which j and k differ
        // Correct for missing data in j and k
        // Loop through each locus in this chunk, called x:
          // Sum the distance between j and k at each locus in this chunk in its
          // corresponding element in M with the previous distances
          // Do the same for M2 but summing in the squared distance at that locus.
    // Retrieve the distance matrix for this genlight object from the bitwise_distance functions
    // Loop over this distance matrix:
      // Sum the distances between all pairs of samples into D
      // Sum the squared distances between all pairs of samples into D2
    // Calculate N choose 2, store as Nc2
    // Calculate the observed variance, Vo = (D2 - D*D/Nc2) / Nc2
    // Calculate the variance at each loci using vars[i] = (M2[i] - (M[i]*M[i])/Nc2) / Nc2
    // Calculate the expected variance by summing over vars: Ve = sumi(vars[i])
    // Calculate the denominator of the final formula by summing over all pairs of loci:
      // denom = 2*sumij(sqrt(vars[i]*vars[j])), such that all combinations of i and j are covered once
    // Store and return the final value for the index of association, (Ve - Vo)/denom

  SEXP R_out;
  SEXP R_gen_symbol;
  SEXP R_chr_symbol;  
  SEXP R_nap_symbol;
  SEXP R_nloc_symbol; // For accessing the number of SNPs in each genotype
  SEXP R_gen;
  int num_gens;
  int num_chunks;
  int chunk_length;
  int num_loci;
  SEXP R_chr1_1;
  SEXP R_chr2_1;
  SEXP R_nap1;
  SEXP R_nap2;
  SEXP R_dists;
  SEXP R_nloc;
  int nap1_length;
  int nap2_length;
  int chr_length;
  int next_missing_index_i;
  int next_missing_index_j;
  int next_missing_i;
  unsigned char missing_mask_i; // These store a map of missing values in the current locus
  unsigned char missing_mask_j; // where 1 is a missing value and 0 is non-missing
  int next_missing_j;
  int missing_match;
  int only_differences;
  int num_threads;
  char mask;
  struct zygosity set_1;
  struct zygosity set_2;
  int i;
  int j;
  int k;
  int x;

  double* vars; // Variance at each locus
  double* M;  // Sum of distances at each locus
  double* M2; // Sum of squared distances at each locus
  int D;   // Sum of distances between each sample
  int D2;  // Sum of squared distances between each sample
  double Vo; // Observed variance
  double Ve; // Expected variance
  double Nc2;  // num_gens choose 2
  double denom; // The denominator for the index of association function

  char** chunk_matrix;
  int cur_distance;
  unsigned char Sn;             // Used for temporary bitwise calculations
  unsigned char offset;
  unsigned char val; 


  // These variables and function calls are used to access elements of the genlight object.
  // ie, R_gen_symbol is being set up as an equivalent to the @gen accessor for genlights.  
  PROTECT(R_gen_symbol = install("gen")); 
  PROTECT(R_chr_symbol = install("snp"));
  PROTECT(R_nap_symbol = install("NA.posi"));  
  PROTECT(R_nloc_symbol = install("n.loc"));

  // This will be a LIST of type LIST:RAW
  // Set R_gen to genlight@gen, a vector of genotypes in the genlight object stored as SNPbin objects.
  R_gen = getAttrib(genlight, R_gen_symbol); 
  // Set num_gens to contain the total number of samples/genotypes within the genlight object
  num_gens = XLENGTH(R_gen);

  // Set R_chr1_1 to be genlight@gen[[0]]@snp[[1]], 
  // aka a raw list representing the entire first set of chromosomes in this genotype
  // Fetching this here expressly to calculate the number of chunks present in each sample
  R_chr1_1 = VECTOR_ELT(getAttrib(VECTOR_ELT(R_gen,0),R_chr_symbol),0); 
  num_chunks = XLENGTH(R_chr1_1);

  // Set R_nloc to be genlight@n.loc, 
  // aka an integer representing the number of loci in each sample
  R_nloc = getAttrib(genlight,R_nloc_symbol);
  num_loci = INTEGER(R_nloc)[0];

  // Prepare and allocate the output matrix
  PROTECT(R_out = allocVector(REALSXP, 1));
  // Prepare and allocate a matrix to store the SNPbin data from all samples
  // so that we don't need to fetch them over and over.
  chunk_matrix = R_Calloc(num_gens,char*);
  for(i = 0; i < num_gens; i++)
  {
    chunk_matrix[i] = R_Calloc(num_chunks,char);
    // Set R_chr1_1 to be the entire SNPbin object for this sample
    R_chr1_1 = VECTOR_ELT(getAttrib(VECTOR_ELT(R_gen,i),R_chr_symbol),0); 
    for(j = 0; j < num_chunks; j++)
    {
      // Set each element of chunk_matrix[i] to be one chunk (8 loci) of the sample
      chunk_matrix[i][j] = (char)RAW(R_chr1_1)[j];
    }
  }

  // This should be 8 for all chunks. If this assumption is wrong things will fail.
  chunk_length = 8;

  // Allocate memory for various arrays and matrices.
  vars = R_Calloc(num_chunks*chunk_length, double);
  M = R_Calloc(num_chunks*chunk_length, double);
  M2 = R_Calloc(num_chunks*chunk_length, double);

  #ifdef _OPENMP
  {
    // Set the number of threads to be used in each omp parallel region
    if(INTEGER(requested_threads)[0] == 0)
    {
      num_threads = omp_get_max_threads();
    }
    else
    {
      num_threads = INTEGER(requested_threads)[0];
    }
    omp_set_num_threads(num_threads);
  }
  #else
  {
    num_threads = 1;
  }
  #endif

  next_missing_index_i = 0;
  next_missing_index_j = 0;
  next_missing_i = -1;
  next_missing_j = -1;
  missing_mask_i = 0;
  missing_mask_j = 0;
  mask = 0;
  nap1_length = 0;
  nap2_length = 0;
  chr_length = 0;
  missing_match = asLogical(missing);

  // Loop through all SNP chunks
  #ifdef _OPENMP
  #pragma omp parallel for schedule(guided) \
    private(i,j,k,x,R_chr1_1,R_chr2_1,R_nap1,R_nap2,Sn,offset,val,\
            next_missing_index_j,next_missing_j,next_missing_index_i,next_missing_i,missing_mask_i,missing_mask_j,\
            set_1,set_2, mask, nap1_length, nap2_length) \
    shared(M, M2, missing_match, only_differences, R_gen, R_nap_symbol,\
           num_gens, num_loci, num_chunks, chunk_length, chunk_matrix)
  #endif
  for(i = 0; i < num_chunks; i++)
  {
    // Loop through all samples
    for(j = 0; j < num_gens; j++)
    {
      // Fill c1 with the next chunk
      set_1.c1 = chunk_matrix[j][i]; // TODO: There is likely a better way to do this for haploids
      // Prepare for correcting missing data
      // Set R_nap1 to be genlight@gen[[j]]@NA.posi
      // aka a vector of integers each representing the location of a locus with missing data
      R_nap1 = getAttrib(VECTOR_ELT(R_gen,j),R_nap_symbol); // Vector of the indices of missing values
      // Set nap1_length to be the number of loci containing missing data in this sample
      nap1_length = XLENGTH(R_nap1);
      next_missing_index_i = 0;
      next_missing_i = 0;
      missing_mask_i = 0;
      // Find the index of each missing value at or before the current locus chunk
      while(next_missing_index_i < nap1_length && (int)INTEGER(R_nap1)[next_missing_index_i]-1 < ((i+1)*8))
      {
        next_missing_i = (int)INTEGER(R_nap1)[next_missing_index_i]-1;
        // If this missing value is in the current locus chunk
        if(next_missing_i >= (i*8) && next_missing_i < ((i+1)*8))
        {
          // Store its inner-chunk location in the mask creating a mask with 1's 
          // in the location of any locus with missing data within this chunk of 8 loci.
          missing_mask_i |= (1<<(next_missing_i%8));
        }
        next_missing_index_i++;
      }
      // Loop through the rest of the genotypes
      for(k = j+1; k < num_gens; k++)
      {
        // Fill c1 with the next chunk
        set_2.c1 = chunk_matrix[k][i]; // TODO: Better way for haploids?

        // Get the locations at which these two chunks differ
        // ^ is the bitwise binary XOR operator, comparing two sets of bits and making a new
          // set with 1's everywhere the two sets were different, and 0's where they were the same
        // So this makes Sn a set with 1's anywhere set_1 and set_2 are different, or in other words
          // at each locus in this chunk at which sample j and sample k have differing alleles.
        Sn = (set_1.c1 ^ set_2.c1);

        // Prepare for correcting missing data
        // Set R_nap2 to be genlight@gen[[k]]@NA.posi
        // aka a vector of integers each representing the location of a locus with missing data
        R_nap2 = getAttrib(VECTOR_ELT(R_gen,k),R_nap_symbol); // Vector of the indices of missing values
        nap2_length = XLENGTH(R_nap2);
        next_missing_index_j = 0;
        next_missing_j = 0;
        missing_mask_j = 0;
        // Find the index of each missing value at or before the current locus chunk
        while(next_missing_index_j < nap2_length && (int)INTEGER(R_nap2)[next_missing_index_j]-1 < ((i+1)*8))
        {
          next_missing_j = (int)INTEGER(R_nap2)[next_missing_index_j]-1;
          // If this missing value is in the current locus chunk
          if(next_missing_j >= (i*8) && next_missing_j < ((i+1)*8))
          {
            // Store its inner-chunk location in the mask creating a mask with 1's 
            // in the location of any locus with missing data within this chunk of 8 loci.
            missing_mask_j |= (1<<(next_missing_j%8));
          }
          next_missing_index_j++;
        }
        // Use the missing data masks to correct for missing data
        if(missing_match)
        {
          // &= is the binary bitwise operator that performs bitwise AND and stores it in
            // the left hand side.
          // Using AND with a mask forces all bits that are 0's in the mask to be 0's in,
            // the final product as well, regardless of what it was in the original.
            // In this case we are negating mask, so ~mask has 0's everywhere it used to
            // have 1's, ie where missing data was found, and thus forces Sn to
            // show a match wherever missing data was found.
          Sn &= ~missing_mask_i; // Force the missing bit to 0 since Sn is matching-not
          Sn &= ~missing_mask_j; // Force the missing bit to 0 since Sn is matching-not
        }
        else
        {
          // |= is a binary bitwise operator that performs bitwise OR and stores it in
            // the left hand side.
          // Using OR with a mask forces all bits that are 1's in the mask to be 1's in
            // the final product as well, regardless of what it was in the original.
            // In other words, this overwrites Sn to show a difference anywhere
            // the mask has a 1, which is anywhere we found missing data.
          Sn |= missing_mask_i; // Force the missing bit to 1 since SN is matching-not
          Sn |= missing_mask_j; // Force the missing bit to 1 since SN is matching-not
        }
        // Loop through 0 through chunk_length, call this x
        for(x = 0; x < chunk_length; x++)
        {
          // Find the offset for the next bit/locus of interest 
          offset = x;
          // Create mask to retrieve just that bit/locus by shifting the 1 in 00000001
            // to the position of the current locus of interest.
          mask = 1<<offset;
          // mask&Sn will zero out every bit other than one we are interest in, then
            // left shifting it with >>offset pushes that bit all the way to the right.
            // In other words, if Sn was a 1 at bit x, val is set to 1. If Sn was a 0
            // at that bit, val is set to 0.
          val = (mask&Sn)>>offset;
          // Update M[current locus] with val 
          M[x + i*chunk_length] += (double)val;
          // Update M2[current locus] with val squared
          M2[x + i*chunk_length] += (double)val*(double)val;
        }
      }
    }
  }

  // Get the distance matrix from bitwise_distance
  R_dists = bitwise_distance_haploid(genlight, missing, requested_threads);

  // Calculate the sum and squared sum of distances between samples
  D = 0;
  D2 = 0;
  #ifdef _OPENMP
  #pragma omp parallel for schedule(guided) reduction(+ : D,D2) private(i,j) 
  #endif
  for(i = 0; i < num_gens; i++)
  {
    for(j = 0; j < i; j++)
    {
      D += INTEGER(R_dists)[i + j*num_gens];
      D2 += INTEGER(R_dists)[i + j*num_gens]*INTEGER(R_dists)[i + j*num_gens];
    }
  }

  // Calculate (num_gens choose 2), which will always be (n*n-n)/2 
  Nc2 = (num_gens*num_gens - num_gens)/2.0;
  // Calculate the observed variance using D and D2
  // Preceding a variable with a datatype, ie (double) forces the computer
    // to treat the variable as that data type for that instance. This is needed
    // here to prevent the occasional implicit typecasting error that was causing
    // this to perform integer division instead of floating point division.
  Vo = ((double)D2 - ((double)D*(double)D)/Nc2) / Nc2;

  // Calculate and fill a vector of variances
  for(i = 0; i < num_loci; i++)
  {
    vars[i] = (M2[i] - (M[i]*M[i])/Nc2) / Nc2;
  }
  // Calculate the expected variance
  Ve = 0;
  for(i = 0; i < num_loci; i++)
  {
    Ve += vars[i];
  }

  // Calculate the denominator for the index of association
  denom = 0;
  #ifdef _OPENMP
  #pragma omp parallel for schedule(guided) reduction(+ : denom) private(i, j)
  #endif
  for(i = 0; i < num_loci; i++)
  {
    for(j = i+1; j < num_loci; j++)  
    {
      if(i != j)
      {
        denom += sqrt(vars[i]*vars[j]);
      }
    }
  }
  denom = 2 * denom;

  // Calculate and store the index of association
  REAL(R_out)[0] = (Vo - Ve) / denom;

  for(i = 0; i < num_gens; i++)
  {
    R_Free(chunk_matrix[i]);
  }
  R_Free(chunk_matrix);
  R_Free(vars);
  R_Free(M);
  R_Free(M2);
  UNPROTECT(5); 
  return R_out;

}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Calculates the index of association of a genlight object of diploids.

Input: A genlight object containing samples of diploids.
       A boolean representing whether or not missing values should match. 
       A boolean representing whether distances or differences should be counted.
       An integer representing the number of threads to be used.
Output: The index of association for this genlight object over the specified loci
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
SEXP association_index_diploid(SEXP genlight, SEXP missing, SEXP differences_only, SEXP requested_threads)
{
  // This function calculates the index of association for samples in 
  // a genlight object. The general flow of this function is as follows:
    // Define and initialize variables
    // Retrieve data from R objects passed in as arguments
    // Fill chunk_matrix with all the genetic data in the genlight object for quick access
    // Prepare for multithreading if compiled to do so
    // Initialize variables needed for the next loop
    // initialize multi threading for the next loop, if compiled to do so
    // Loop through every chunk of 8 loci in the genlight object, call each on i:
      // Loop through every genotype in the genlight object, call each one j:
      // Prepare the genetic data for j stored in chunk_matrix for comparison
      // Prepare to correct for missing data in j
      // Loop through every genotype from j on to cover all pairings, call this k:
        // Prepare the genetic data for k stored in chunk_matrix for comparison
        // Find locations at which j and k differ
        // Correct for missing data in j and k
        // Loop through each locus in this chunk, called x:
          // Sum the distance between j and k at each locus in this chunk in its
          // corresponding element in M with the previous distances
          // Do the same for M2 but summing in the squared distance at that locus.
    // Retrieve the distance matrix for this genlight object from the bitwise_distance functions
    // Loop over this distance matrix:
      // Sum the distances between all pairs of samples into D
      // Sum the squared distances between all pairs of samples into D2
    // Calculate N choose 2, store as Nc2
    // Calculate the observed variance, Vo = (D2 - D*D/Nc2) / Nc2
    // Calculate the variance at each loci using vars[i] = (M2[i] - (M[i]*M[i])/Nc2) / Nc2
    // Calculate the expected variance by summing over vars: Ve = sumi(vars[i])
    // Calculate the denominator of the final formula by summing over all pairs of loci:
      // denom = 2*sumij(sqrt(vars[i]*vars[j])), such that all combinations of i and j are covered once
    // Store and return the final value for the index of association, (Ve - Vo)/denom

  SEXP R_out;
  SEXP R_gen_symbol;
  SEXP R_chr_symbol;  
  SEXP R_nap_symbol;
  SEXP R_nloc_symbol; // For accessing the number of SNPs in each genotype
  SEXP R_gen;
  int num_gens;
  int num_chunks;
  int chunk_length;
  int num_loci;
  SEXP R_chr1_1;
  SEXP R_chr1_2;
  SEXP R_chr2_1;
  SEXP R_chr2_2;
  SEXP R_nap1;
  SEXP R_nap2;
  SEXP R_dists;
  SEXP R_nloc;
  int nap1_length;
  int nap2_length;
  int chr_length;
  int next_missing_index_i;
  int next_missing_index_j;
  int next_missing_i;
  int next_missing_j;
  unsigned char missing_mask_i; // These store a map of missing values in the current locus
  unsigned char missing_mask_j; // where 1 is a missing value and 0 is non-missing
  int missing_match;
  int only_differences;
  int num_threads;
  char mask;
  struct zygosity set_1;
  struct zygosity set_2;
  int i;
  int j;
  int k;
  int x;

  double* vars; // Variance at each locus
  double* M;  // Sum of distances at each locus
  double* M2; // Sum of squared distances at each locus
  int D;   // Sum of distances between each sample
  int D2;  // Sum of squared distances between each sample
  double Vo; // Observed variance
  double Ve; // Expected variance
  double Nc2;  // num_gens choose 2
  double denom; // The denominator for the index of association function

  char** chunk_matrix;
  unsigned char Sn;             // Used for temporary bitwise calculations
  unsigned char Hnor;
  unsigned char Hs;
  unsigned char offset;
  unsigned char val; 


  // These variables and function calls are used to access elements of the genlight object.
  // ie, R_gen_symbol is being set up as an equivalent to the @gen accessor for genlights.  
  PROTECT(R_gen_symbol = install("gen")); // Used for accessing the named elements of the genlight object
  PROTECT(R_chr_symbol = install("snp"));
  PROTECT(R_nap_symbol = install("NA.posi"));  
  PROTECT(R_nloc_symbol = install("n.loc"));

  // This will be a LIST of type LIST:RAW
  // Set R_gen to genlight@gen, a vector of genotypes in the genlight object stored as SNPbin objects.
  R_gen = getAttrib(genlight, R_gen_symbol); 
  // Set num_gens to contain the total number of samples/genotypes within the genlight object
  num_gens = XLENGTH(R_gen);

  // Set R_chr1_1 to be genlight@gen[[0]]@snp[[1]], 
  // aka a raw list representing the entire first set of chromosomes in this genotype
  // Fetching this here expressly to calculate the number of chunks present in each sample
  R_chr1_1 = VECTOR_ELT(getAttrib(VECTOR_ELT(R_gen,0),R_chr_symbol),0); 
  num_chunks = XLENGTH(R_chr1_1);

  // Set R_nloc to be genlight@n.loc, 
  // aka an integer representing the number of loci in each sample
  R_nloc = getAttrib(genlight,R_nloc_symbol);
  num_loci = INTEGER(R_nloc)[0];

  // Prepare and allocate the output matrix
  PROTECT(R_out = allocVector(REALSXP, 1));
  // Prepare and allocate a matrix to store the SNPbin data from all samples
  // so that we don't need to fetch them over and over.
  chunk_matrix = R_Calloc(num_gens*2,char*);
  for(i = 0; i < num_gens; i++)
  {
    chunk_matrix[i*2] = R_Calloc(num_chunks,char);
    chunk_matrix[i*2+1] = R_Calloc(num_chunks,char);
    // Set R_chr1_1 and R_chr1_2 to be the entire SNPbin object for this sample
    R_chr1_1 = VECTOR_ELT(getAttrib(VECTOR_ELT(R_gen,i),R_chr_symbol),0);
    R_chr1_2 = VECTOR_ELT(getAttrib(VECTOR_ELT(R_gen,i),R_chr_symbol),1); 
    for(j = 0; j < num_chunks; j++)
    {
      // Set each element of chunk_matrix[i*2] to be one chunk (8 loci) of the sample i
        // and each element of chunk_matrix[i*2+1] to be the second set of that same chunk 
        // in order to store the full diploid sample for each sample
      chunk_matrix[i*2][j] = (char)RAW(R_chr1_1)[j];
      chunk_matrix[i*2+1][j] = (char)RAW(R_chr1_2)[j];
    }
  }

  // This should be 8 for all chunks. If this assumption is wrong things will fail.
  chunk_length = 8;

  // Allocate memory for various arrays and matrices.
  vars = R_Calloc(num_chunks*chunk_length, double);
  M = R_Calloc(num_chunks*chunk_length, double);
  M2 = R_Calloc(num_chunks*chunk_length, double);

  #ifdef _OPENMP
  {
    // Set the number of threads to be used in each omp parallel region
    if(INTEGER(requested_threads)[0] == 0)
    {
      num_threads = omp_get_max_threads();
    }
    else
    {
      num_threads = INTEGER(requested_threads)[0];
    }
    omp_set_num_threads(num_threads);
  }
  #else
  {
    num_threads = 1;
  }
  #endif

  next_missing_index_i = 0;
  next_missing_index_j = 0;
  next_missing_i = -1;
  next_missing_j = -1;
  missing_mask_i = 0;
  missing_mask_j = 0;
  mask = 0;
  nap1_length = 0;
  nap2_length = 0;
  chr_length = 0;
  missing_match = asLogical(missing);
  only_differences = asLogical(differences_only);

  // Loop through all SNP chunks
  #ifdef _OPENMP
  #pragma omp parallel for schedule(guided) \
    private(i,j,k,x,R_chr1_1,R_chr1_2,R_chr2_1,R_chr2_2,R_nap1,R_nap2,Sn,Hnor,Hs,offset,val,\
            next_missing_index_j,next_missing_j,next_missing_index_i,next_missing_i,missing_mask_i,missing_mask_j,\
            set_1,set_2, mask, nap1_length, nap2_length) \
    shared(M, M2, missing_match, only_differences, R_gen, R_nap_symbol,\
           num_gens, num_loci, num_chunks, chunk_length, chunk_matrix)
  #endif
  for(i = 0; i < num_chunks; i++)
  {
    // Loop through all samples
    for(j = 0; j < num_gens; j++)
    {
      // Fill c1 and c2 with the next two chunks, representing a pair
      set_1.c1 = chunk_matrix[j*2][i];
      set_1.c2 = chunk_matrix[j*2+1][i]; // This chunk is the chromosonal pair of the first
      // Call fill_zygosity to fill out the rest of the fields
      fill_zygosity(&set_1);
      // Prepare for correcting missing data
      // Set R_nap1 to be genlight@gen[[j]]@NA.posi
      // aka a vector of integers each representing the location of a locus with missing data
      R_nap1 = getAttrib(VECTOR_ELT(R_gen,j),R_nap_symbol); // Vector of the indices of missing values
      // Set nap1_length to be the number of loci containing missing data in this sample
      nap1_length = XLENGTH(R_nap1);
      next_missing_index_i = 0;
      next_missing_i = 0;
      missing_mask_i = 0;
      // Find the index of each missing value at or before the current locus chunk
      while(next_missing_index_i < nap1_length && (int)INTEGER(R_nap1)[next_missing_index_i]-1 < ((i+1)*8))
      {
        next_missing_i = (int)INTEGER(R_nap1)[next_missing_index_i]-1;
        // If this missing value is in the current locus chunk
        if(next_missing_i >= (i*8) && next_missing_i < ((i+1)*8))
        {
          // Store its inner-chunk location in the mask creating a mask with 1's 
          // in the location of any locus with missing data within this chunk of 8 loci.
          missing_mask_i |= (1<<(next_missing_i%8));
        }
        next_missing_index_i++;
      }
      // Loop through the rest of the genotypes
      for(k = j+1; k < num_gens; k++)
      {
        // Fill c1 and c2 with the next two chunks, representing a pair
        set_2.c1 = chunk_matrix[k*2][i];
        set_2.c2 = chunk_matrix[k*2+1][i]; // This chunk is the chromosonal pair of the first
        // Call fill_zygosity to fill out the rest of the fields
        fill_zygosity(&set_2);

        // Get the similarity set between the two, then negate it. Store as Sn
        Sn = ~get_similarity_set(&set_1,&set_2);
        // OR the two H fields together, negate it, and store that as Hnor
        Hnor = ~(set_1.ch | set_2.ch);
        // AND together Hnor and Sn to get a set of differing homozygote sites, store as Hs
        Hs = Sn & Hnor;

        // Prepare for correcting missing data
        R_nap2 = getAttrib(VECTOR_ELT(R_gen,k),R_nap_symbol); // Vector of the indices of missing values
        nap2_length = XLENGTH(R_nap2);
        next_missing_index_j = 0;
        next_missing_j = 0;
        missing_mask_j = 0;
        // Find the index of each missing value at or before the current locus chunk
        while(next_missing_index_j < nap2_length && (int)INTEGER(R_nap2)[next_missing_index_j]-1 < ((i+1)*8))
        {
          next_missing_j = (int)INTEGER(R_nap2)[next_missing_index_j]-1;
          // If this missing value is in the current locus chunk
          if(next_missing_j >= (i*8) && next_missing_j < ((i+1)*8))
          {
            // Store its inner-chunk location in the mask creating a mask with 1's 
            // in the location of any locus with missing data within this chunk of 8 loci.
            missing_mask_j |= (1<<(next_missing_j%8));
          }
          next_missing_index_j++;
        }
        // Use the missing data masks to correct for missing data
        if(missing_match)
        {
          // &= is the binary bitwise operator that performs bitwise AND and stores it in
            // the left hand side.
          // Using AND with a mask forces all bits that are 0's in the mask to be 0's in,
            // the final product as well, regardless of what it was in the original.
            // In this case we are negating mask, so ~mask has 0's everywhere it used to
            // have 1's, ie where missing data was found, and thus forces Sn and Hs to
            // show a match wherever missing data was found.
          Sn &= ~missing_mask_i; // Force the missing bit to 0 since Sn is matching-not
          Sn &= ~missing_mask_j; // Force the missing bit to 0 since Sn is matching-not
          Hs &= ~missing_mask_i; // Hs is differing homozygote sites, so missing should be 0
          Hs &= ~missing_mask_j; // Hs is differing homozygote sites, so missing should be 0
        }
        else
        {
          // |= is a binary bitwise operator that performs bitwise OR and stores it in
            // the left hand side.
          // Using OR with a mask forces all bits that are 1's in the mask to be 1's in
            // the final product as well, regardless of what it was in the original.
            // In other words, this overwrites Sn to show a difference anywhere
            // the mask has a 1, which is anywhere we found missing data.
          Sn |= missing_mask_i; // Force the missing bit to 1 since SN is matching-not
          Sn |= missing_mask_j; // Force the missing bit to 1 since SN is matching-not
          Hs |= (~set_2.ch & missing_mask_i);// Hs should be 1 if set2.ch is 0 at this bit
                                  // This makes missing data a distance of 1 from a heterozygote site
                                  // and 2 from either homozygote site.
          Hs |= (~set_1.ch & missing_mask_j);// Hs should be 1 if set2.ch is 0 at this bit
                                  // This makes missing data a distance of 1 from a heterozygote site
                                  // and 2 from either homozygote site.
          // In other words, missing_match == FALSE will return the largest possible distance
            // between two samples at a locus given what, if any, data is known.

        }
        // Loop through 0 through chunk_length, call this x
        for(x = 0; x < chunk_length; x++)
        {
          // Find the offset for the next bit/locus of interest 
          offset = x;
          // Create mask to retrieve just that bit/locus by shifting the 1 in 00000001
            // to the position of the current locus of interest.
          mask = 1<<offset;
          // mask&Sn will zero out every bit other than one we are interest in, then
            // left shifting it with >>offset pushes that bit all the way to the right.
            // In other words, if Sn was a 1 at bit x, val is set to 1. If Sn was a 0
            // at that bit, val is set to 0.
          val = (mask&Sn)>>offset;
          // If only_differences is false, ie we want distance not differences
          if(!only_differences)
          {
            // We need to add 1 to val if Hs was 1 at this locus, meaning we add one
              // if the two samples were differing homozygotes at that site.
            val += (mask&Hs)>>offset;
          }
          // Update M[current locus] with val 
          M[x + i*chunk_length] += (double)val;
          // Update M2[current locus] with val squared
          M2[x + i*chunk_length] += (double)val*(double)val;
        }
      }
    }
  }

  // Get the distance matrix from bitwise_distance
  R_dists = bitwise_distance_diploid(genlight, missing, differences_only, requested_threads);

  // Calculate the sum and squared sum of distances between samples
  D = 0;
  D2 = 0;
  x = 0;
  #ifdef _OPENMP
  #pragma omp parallel for schedule(guided) reduction(+ : D,D2) private(i,j) 
  #endif
  for(i = 0; i < num_gens; i++)
  {
    for(j = 0; j < i; j++)
    {
      D += INTEGER(R_dists)[i + j*num_gens];
      D2 += INTEGER(R_dists)[i + j*num_gens]*INTEGER(R_dists)[i + j*num_gens];
    }
  }

  // Calculate C(num_gens,2), which will always be (n*n-n)/2 
  Nc2 = (num_gens*num_gens - num_gens)/2.0;
  // Calculate the observed variance using D and D2
  Vo = ((double)D2 - ((double)D*(double)D)/Nc2) / Nc2;

  // Calculate and fill a vector of variances
  for(i = 0; i < num_loci; i++)
  {
    vars[i] = (M2[i] - (M[i]*M[i])/Nc2) / Nc2;
  }
  // Calculate the expected variance
  Ve = 0;
  for(i = 0; i < num_loci; i++)
  {
    Ve += vars[i];
  }

  // Calculate the denominator for the index of association
  denom = 0;
  #ifdef _OPENMP
  #pragma omp parallel for schedule(guided) reduction(+ : denom) private(i, j)
  #endif
  for(i = 0; i < num_loci; i++)
  {
    for(j = i+1; j < num_loci; j++)  
    {
      if(i != j)
      {
        denom += sqrt(vars[i]*vars[j]);
      }
    }
  }
  denom = 2 * denom;

  // Calculate and store the index of association
  REAL(R_out)[0] = (Vo - Ve) / denom;

  for(i = 0; i < num_gens*2; i++)
  {
    R_Free(chunk_matrix[i]);
  }
  R_Free(chunk_matrix);
  R_Free(vars);
  R_Free(M);
  R_Free(M2);
  UNPROTECT(5); 
  return R_out;

}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Calculates and returns a matrix of Pgen values for each genotype and loci in
the genind or genclone  object.

Input: A genind or genclone object containing samples of diploids.
       A frequency matrix constructed in R with makefreq(genind2genpop(genind)).
       A vector of population indices for all samples.
       An integer specifying the number of populations.
Output: A matrix containing the Pgen value of each genotype at each locus.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
SEXP get_pgen_matrix_genind(SEXP genind, SEXP freqs, SEXP pops, SEXP npop)
{
  SEXP R_out;
  SEXP R_tab_symbol;
  SEXP R_loc_symbol;
  SEXP R_pop_symbol;
  SEXP R_ploidy_symbol;
  SEXP R_tab;
  PROTECT(R_tab_symbol = install("tab")); // Used for accessing the named elements of the genind object
  PROTECT(R_loc_symbol = install("loc.n.all"));
  PROTECT(R_ploidy_symbol = install("ploidy"));
  double* pgens;
  int* indices;
  int* ploidy;
  int num_gens;
  int num_loci;
  int num_alleles;
  int num_pops;
  int index;
  int missing_last;
  int pop;
  int size;

  R_tab = getAttrib(genind, R_tab_symbol); 
  num_gens = INTEGER(getAttrib(R_tab, R_DimSymbol))[0];
  num_alleles = INTEGER(getAttrib(R_tab, R_DimSymbol))[1];
  num_loci = XLENGTH(getAttrib(genind, R_loc_symbol));
  num_pops = INTEGER(npop)[0];
  ploidy = INTEGER(getAttrib(genind, R_ploidy_symbol));
  missing_last = 0;
  size = num_gens*num_loci;
  indices = R_Calloc(size*2, int);
  PROTECT(R_out = allocMatrix(REALSXP, num_gens, num_loci));
  pgens = REAL(R_out);

  // Fill indices with the indices inside freqs of alleles found in each genotype
  // Note that R_tab is column-major ordered due to being passed from R
  for(int i = 0; i < num_gens; i++)
  {
    index = 0;
    for(int j = 0; j < num_alleles; j++)
    {
      if(INTEGER(R_tab)[i + j*num_gens] == NA_INTEGER)
      {
        // Add missing value sentinels into the index array once for each missing locus
        if(missing_last == 0)
        {
          indices[i*num_loci*2 + index] = -1;;
          indices[i*num_loci*2 + index+1] = -1;
          index += 2;
          missing_last = 1;
        }
      }
      else if(INTEGER(R_tab)[i + j*num_gens] == 2)
      {
        indices[i*num_loci*2 + index] = j;
        indices[i*num_loci*2 + index + 1] = j;
        index += 2;
        missing_last = 0;
      }
      else if(INTEGER(R_tab)[i+ j*num_gens] == 1)
      {
        indices[i*num_loci*2 + index] = j;
        index += 1;
        missing_last = 0;
      }
    }
  } 

  for(int i = 0; i < num_gens; i++)
  {
    pop = INTEGER(pops)[i] - 1;
    index = 0;
    for(int j = 0; j < num_loci; j++)
    {
      if(indices[i*num_loci*2 + index] == -1) // And therefore also [~ index+1] == -1
      {
        // Set the pgen value of this genotype at this locus to missing.
        pgens[i + j*num_gens] = NA_REAL;
      }
      else if (ploidy[i] == 1)
      {
        // If the ploidy is 1, the probability of a genotype at that locus 
        // is the allele frequency
        pgens[i + j*num_gens] = log(REAL(freqs)[pop + indices[i*num_loci*2 + index]*num_pops]);
      }
      else
      {
        pgens[i + j*num_gens] = log(REAL(freqs)[pop + indices[i*num_loci*2 + index]*num_pops]) + log(REAL(freqs)[pop + indices[i*num_loci*2 + index+1]*num_pops]);
        // Account for both permutations of heterozygous loci
        if (indices[i*num_loci*2 + index] != indices[i*num_loci*2 + index + 1])
        {
          pgens[i + j*num_gens] += log(2);
        }
      }
      index += 2;
    }
  }
  R_Free(indices);
  UNPROTECT(4);
  return R_out;
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Calculates and returns a matrix of Pgen values for the given genlight object.

Input: A genlight object containing samples of diploids.
Output: A matrix containing the Pgen value of each locus in each genotype in the 
        genlight object.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
SEXP get_pgen_matrix_genlight(SEXP genlight, SEXP window)
{

  // TODO: Accept an interval, then use either 8 or num_loci as default when calling from R, never 0

  SEXP R_out;
  SEXP R_gen_symbol;
  SEXP R_loc_symbol;
  SEXP R_pop_symbol;
  PROTECT(R_gen_symbol = install("gen")); // Used for accessing the named elements of the genlight object
  PROTECT(R_loc_symbol = install("n.loc"));
  PROTECT(R_pop_symbol = install("pop"));
  struct locus* loci; 
  double *pgens;
  int num_gens;
  int num_loci;
  int num_sets;
  int num_pops;
  int interval;
  int size;
  // Set the interval to calculate Pgen over every 8 loci
  interval = INTEGER(window)[0];

  num_gens = XLENGTH(getAttrib(genlight, R_gen_symbol));
  num_loci = INTEGER(getAttrib(genlight, R_loc_symbol))[0];
  num_sets = ceil((double)num_loci/(double)interval); // Number of sets of loci for which pgen values should be computed
  size = num_gens*num_sets;
  pgens = R_Calloc(size, double);
  PROTECT(R_out = allocVector(REALSXP,size));

  // Find the number of populations by taking the max over all genotypes
  num_pops = 0;
  for(int i = 0; i < num_gens; i++)
  {
    num_pops = (INTEGER(getAttrib(genlight, R_pop_symbol))[i] > num_pops) ? INTEGER(getAttrib(genlight, R_pop_symbol))[i] : num_pops; 
  }

  // Allocate memory for the array of locus struct
  loci = R_Calloc(num_loci*num_pops,struct locus);
  // Call fill_loci to get allelic frequency information
  fill_loci(loci, genlight);

  fill_Pgen(pgens,loci,interval,genlight);    

  for(int i = 0; i < num_gens; i++)
  {
    for(int j = 0; j < num_sets; j++)
    {
      // Transpose the matrix into column-major ordering before returning to R
      REAL(R_out)[i + j*num_gens] = (pgens[i*num_sets + j]);
    }
  }

  R_Free(loci);
  R_Free(pgens);
  UNPROTECT(4);
  return R_out;
}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Fills an array of doubles with the Pgen value associated with each individual
found in the genlight object. These values represent the probability of each
individual having been produced via random mating of the population, as estimated
by the samples present in the genlight object

Input: A pointer to an array of doubles to be filled.
        NOTE: This array MUST have a length equal to the number of genotypes found
        in the genlight object times the number of loci divided by specified interval
        ie, num_gens*ceil((double)num_loci/(double)interval)
       A pointer to an array of allelic frequencies filled with fill_loci
       THe number of loci which should be considered in each Pgen value 
       A genlight object from which the individual genotypes can be obtained.
Output: None. Fills in the array of doubles with the log of the Pgen value of each 
        individual genotype in the genlight object.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void fill_Pgen(double *pgen, struct locus *loci, int interval, SEXP genlight)
{

  double log_product;
  int num_loci;
  int num_gens;
  int chr_length;
  int group;
  int num_groups;
  int pop;
  struct locus* loc;
  struct zygosity zyg;
  SEXP R_gen;
  SEXP R_gen_symbol;
  SEXP R_chr_symbol;
  SEXP R_loc_symbol;
  SEXP R_nap_symbol;
  SEXP R_pop_symbol;
  SEXP R_chr1;
  SEXP R_chr2;
  SEXP R_nap;
  int nap_length;
  int next_missing_index;
  int next_missing;

  PROTECT(R_gen_symbol = install("gen")); // Used for accessing the named elements of the genlight object
  PROTECT(R_loc_symbol = install("n.loc"));
  PROTECT(R_chr_symbol = install("snp"));
  PROTECT(R_nap_symbol = install("NA.posi"));  
  PROTECT(R_pop_symbol = install("pop"));

  R_gen = getAttrib(genlight, R_gen_symbol);
  num_gens = XLENGTH(R_gen);
  num_loci = INTEGER(getAttrib(genlight, R_loc_symbol))[0];
  num_groups = ceil((double)num_loci/(double)interval);
  nap_length = 0;
  next_missing_index = 0;
  next_missing = -1;
  
  // for each genotype
  for(int i = 0; i < num_gens; i++)
  {
    R_chr1 = VECTOR_ELT(getAttrib(VECTOR_ELT(R_gen,i),R_chr_symbol),0);
    R_chr2 = VECTOR_ELT(getAttrib(VECTOR_ELT(R_gen,i),R_chr_symbol),1);
    chr_length = XLENGTH(R_chr1);
    R_nap = getAttrib(VECTOR_ELT(R_gen,i),R_nap_symbol); // Vector of the indices of missing values 
    nap_length = XLENGTH(R_nap);
    if(nap_length > 0)
    {
      next_missing_index = 0;
      next_missing = (int)INTEGER(R_nap)[next_missing_index] - 1; //Compensate for Rs 1 based indexing
    }
    else
    {
      next_missing_index = -1;
      next_missing = -1;
    }
    // Get the population of this genotype
    pop = INTEGER(getAttrib(genlight, R_pop_symbol))[i] - 1;
    // restart the log_product and group numbers
    log_product = 0;
    group = 0;
    // for each locus group
  
    for(int j = 0; j < chr_length; j++)
    {
      // Find zygosity for each site in this locus
      zyg.c1 = (char)RAW(R_chr1)[j];
      zyg.c2 = (char)RAW(R_chr2)[j];
      fill_zygosity(&zyg);

      // pseudo code for the next loop
      // for each bit in locus
        // if genotype is heterozygous at this site
          // add one to heterzygotes
          // product *= dominants * recessives / 4n^2 + (-1*dominants*recessives/4n^2 * (1 - h*4n / (2*dominants*recessives)))
          // product *= h/2n
        // else the genotype is homozygous dominant
          // product *= dominants*dominants / 4n^2 + ((dominants/2n) * (1 - dominants/2n) * (1 - h*4n / (2*dominants*recessives)))
          // product *= (rd + dh - 2nh)/(2rn)
        // else the genotype is homozygous recessive
          // product *= recessives*recessives / 4n^2 + ((recessives/2n) * (1 - recessives/2n) * (1 -  h*4n / (2*dominants*recessives)))
          // product *= (rd + rh - 2nh)/(2dn)
      // end for
  
      for(int k = 0; k < 8 && j*8+k < num_loci; k++)
      {
        // Skip any missing data
        if(next_missing_index < nap_length && next_missing != j*8+k)
        {
          loc = &loci[pop*num_loci + j*8 + k]; // loc[pop][bit]
          if(((zyg.ch >> k) & 1) == 1)
          {
            // handle each heterozygous site in this locus
            // log_product += log(loc->h) - log(loc->n); // - log(2) + log(2)  // Pgen(f)
            log_product += log(loc->d) + log(loc->r) - (log(2*(loc->n))+log(2*(loc->n))) + log(2);  // Pgen
          }
          else if(((zyg.cd >> k) & 1) == 1)
          {
            // handle each homozygous dominant site in this locus
            // log_product += log( (loc->r)*(loc->d) + (loc->d)*(loc->h) - (2*(loc->n))*(loc->h) ) - log( 2*(loc->r)*(loc->n) ); // Pgen(f)
            log_product += log(loc->d) + log(loc->d) - (log(2*(loc->n))+log(2*(loc->n))); // Pgen
          }
          else if(((zyg.cr >> k) & 1) == 1)
          {
            // handle each homozygous recessive site in this locus
            // log_product += log( (loc->r)*(loc->d) + (loc->r)*(loc->h) - (2*(loc->n))*(loc->h) ) - log( 2*(loc->d)*(loc->n) ); // Pgen(f)
            log_product += log(loc->r) + log(loc->r) - (log(2*(loc->n))+log(2*(loc->n))); // Pgen
          }
          else
          {
            // Error
          }
        }
        else if(next_missing == j*8+k)
        {
          // Update the missing data variables in order to skip the next portion of missing data
          next_missing_index++; 
          next_missing = (int)INTEGER(R_nap)[next_missing_index] - 1;
        }
        if((group+1)*interval == j*8+k+1 || j*8+k+1 == num_loci)
        {
          pgen[i*num_groups + group] = log_product;
          group++;
          log_product = 0;
        }
      }
    }
  }

  UNPROTECT(5);

}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Fills an array of struct locus objects based on allelic frequencies found in
the provided genlight object.

Input: A pointer to an array of locus objects to be filled.
        NOTE: This array MUST have a length equal to the number of loci found
        in each genotype in the genlight object.
       A genlight object from which the alleles and loci can be gathered.
Output: None. Fills in the allelic frequencies and other information found in
        each locus struct.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void fill_loci(struct locus *loc, SEXP genlight)
{ 
  // ~Pseudo code~
  // for each genotype in genlight
    // for each locus group in genotype
      // zyg.c1 = locus[0]
      // zyg.c2 = locus[1]
      // fill_zygosity(&zyg) // Find zygosity at each site in locus
      // for each bit in locus group 
        // if not missing
          // // Exactly one of the next three lines will add 1 (or 2) to its respective counter
            // Add 1 to h if this is a heterozygous site
          // loc[locus*8+bit].h += (zyg.ch & (1<<bit)) >> bit // Or is that 1 << 7-bit ?
            // Add 1 to dominant if this is heterozygous or 2 if this is homozygous dominant
          // loc[locus*8+bit].d += ((zyg.ch & (1<<bit)) >> bit) + 2 * (zyg.cd & (1<<bit)) 
            // Add 1 to recessive if this is heterozygous or 2 if this is homozygous recessive
          // loc[locus*8+bit].r += ((zyg.ch & (1<<bit)) >> bit) + 2 * (zyg.cr & (1<<bit)) 
            // Add 1 to the number of contributing genotypes regardless of what else was added
          // loc[locus*8+bit].n += 1

  struct zygosity zyg;
  SEXP R_gen_symbol;
  SEXP R_loc_symbol;
  SEXP R_chr_symbol;  
  SEXP R_nap_symbol;
  SEXP R_pop_symbol;
  SEXP R_gen;
  int num_gens;
  int num_loci;
  SEXP R_chr1_1;
  SEXP R_chr1_2;
  SEXP R_nap1;
  int nap1_length;
  int chr_length;
  int next_missing_index_i;
  int next_missing_i;
  int i;
  int bit;
  int byte;
  int pop;

  PROTECT(R_gen_symbol = install("gen")); // Used for accessing the named elements of the genlight object
  PROTECT(R_loc_symbol = install("n.loc"));
  PROTECT(R_chr_symbol = install("snp"));
  PROTECT(R_nap_symbol = install("NA.posi"));  
  PROTECT(R_pop_symbol = install("pop"));

  // This will be a LIST of type LIST:RAW
  R_gen = getAttrib(genlight, R_gen_symbol);
  num_gens = XLENGTH(R_gen);
  num_loci = INTEGER(getAttrib(genlight, R_loc_symbol))[0];
  
  next_missing_index_i = 0;
  next_missing_i = -1;
  nap1_length = 0;
  chr_length = 0;

  // Loop through every genotype 
  for(i = 0; i < num_gens; i++)
  {
    R_chr1_1 = VECTOR_ELT(getAttrib(VECTOR_ELT(R_gen,i),R_chr_symbol),0); // Chromosome 1
    R_chr1_2 = VECTOR_ELT(getAttrib(VECTOR_ELT(R_gen,i),R_chr_symbol),1); // Chromosome 2
    chr_length = XLENGTH(R_chr1_1);
    R_nap1 = getAttrib(VECTOR_ELT(R_gen,i),R_nap_symbol); // Vector of the indices of missing values 
    nap1_length = XLENGTH(R_nap1);
    if(nap1_length > 0)
    {
      next_missing_index_i = 0;
      next_missing_i = (int)INTEGER(R_nap1)[next_missing_index_i] - 1; //Compensate for Rs 1 based indexing
    }
    else
    {
      next_missing_index_i = -1;
      next_missing_i = -1;
    }
    // Get the population of this genotype
    pop = INTEGER(getAttrib(genlight, R_pop_symbol))[i] - 1;
    // Loop through all the chunks of SNPs in this genotype
    for(byte = 0; byte < chr_length; byte++) 
    {
      zyg.c1 = (char)RAW(R_chr1_1)[byte];
      zyg.c2 = (char)RAW(R_chr1_2)[byte];
      fill_zygosity(&zyg);
      
      // TODO: Find a way to parallelize either here or the locus loop.
      //       Unfortunately, next_missing_i complicates this.
      for(bit = 0; bit < 8 && byte*8+bit < num_loci; bit++)
      {
        // if not missing
        if(next_missing_i != byte*8+bit)
        { 
          // The following lines will add 1 (or 2) to whichever value(s) need to be increased, and 0 to the others
          loc[pop*num_loci + byte*8+bit].h += (zyg.ch >> bit) & 1;
          loc[pop*num_loci + byte*8+bit].d += ((zyg.ch >> bit) & 1) + 2 * ((zyg.cd >> bit) & 1); 
          loc[pop*num_loci + byte*8+bit].r += ((zyg.ch >> bit) & 1) + 2 * ((zyg.cr >> bit) & 1); 
          loc[pop*num_loci + byte*8+bit].n += 1;
        }
        else
        {
          if(next_missing_index_i < nap1_length-1)
          {
            next_missing_index_i++;
            next_missing_i = (int)INTEGER(R_nap1)[next_missing_index_i] - 1;
          }
        }
      }
    }
  } 

  UNPROTECT(5);
}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Calculates the zygosity at each location of a given section. The zygosity struct
must have c1 and c2 filled before calling this function.

Input: A zygosity struct pointer with c1 and c2 filled for the given section.
Output: None. Fills the ch, cd, and cr values in the provided struct to indicate
        heterozygous, homozygous dominant, and homozygous recessive sites respectively,
        where 1's in each string represent the presence of that zygosity and 0's
        represent a different zygosity.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void fill_zygosity(struct zygosity *ind)
{
  ind->ch = ind->c1 ^ ind->c2;  // Bitwise Exclusive OR provides 1's only at heterozygous sites
  ind->cd = ind->c1 & ind->c2;  // Bitwise AND provides 1's only at homozygous dominant sites
  ind->cr = ~(ind->c1 | ind->c2); // Bitwise NOR provides 1's only at homozygous recessive sites
}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Finds the locations at which two samples have differing zygosity.

Input: Two zygosity struct pointers with c1 and c2 filled, each representing the
       same section from two samples.
Output: A char representing a binary string of differences between the
        two samples in the given section. 0's represent a difference in zygosity
        at that location and 1's represent matching zygosity.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
char get_similarity_set(struct zygosity *ind1, struct zygosity *ind2)
{

  char s;   // The similarity values to be returned
  char sx;  // 1's wherever both are heterozygous
  char sa;  // 1's wherever both are homozygous dominant
  char sn;  // 1's wherever both are homozygous recessive

  sx = ind1->ch & ind2->ch;
  sa = ind1->cd & ind2->cd;
  sn = ind1->cr & ind2->cr;

  s = sx | sa | sn; // 1's wherever both individuals share the same zygosity  

  return s;
}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Calculates the number of zeros in a binary string. Used by get_difference to
find the number of differences between two samples in a given section.

Input: A char representing a binary string of differences between samples
       where 1's are matches and 0's are differences.
Output: The number of zeros in the argument value, representing the number of 
        differences between the two samples in the location used to generate
        the sim_set.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int get_zeros(char sim_set)
{
  int zeros = 0;
  int tmp_set = sim_set;
  int digits = 8; //sizeof(char);
  int i;

  for(i = 0; i < digits; i++)
  {
    if(tmp_set%2 == 0)
    {
      zeros++;
    }
    tmp_set = tmp_set >> 1; // Drop the rightmost digit and shift the others over 1
  }

  return zeros;
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Counts the number of differences between two partially filled zygosity structs.
c1 and c2 must be filled in both structs prior to calling this function.

Input: Two zygosity struct pointers representing the two sections to be compared.
       c1 and c2 must be filled in both structs.
Output: The number of locations in the given section that have differing zygosity
        between the two samples.
        cx, ca, and cn will be filled in both structs as a byproduct of this function.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int get_difference(struct zygosity *z1, struct zygosity *z2)
{
  int dif = 0;
  fill_zygosity(z1);
  fill_zygosity(z2);
  dif = get_zeros(get_similarity_set(z1,z2));

  return dif;
}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Counts the distance between two partially filled zygosity structs.
c1 and c2 must be filled in both structs prior to calling this function.

Input: Two zygosity struct pointers representing the two sections to be compared.
       c1 and c2 must be filled in both structs.
Output: The total distance between two samples, such that DD/rr are a distance 
        of 2, and Dr/rr are a distance of 1
        cx, ca, and cn will be filled in both structs as a byproduct of this function.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int get_distance(struct zygosity *z1, struct zygosity *z2)
{
  int dist = 0;
  char Hor;
  char S;
  char ch_dist;
  fill_zygosity(z1);
  fill_zygosity(z2);

  S = get_similarity_set(z1,z2);
  Hor = z1->ch | z2->ch;

  ch_dist = Hor | S;  // Force ones everywhere they are the same
  dist = get_zeros(S);  // Add one distance for every non-shared zygosity
  dist += get_zeros(ch_dist); // Add another one for every difference that has no heterozygotes

  return dist;
}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Counts the distance between two partially filled zygosity structs, using a custom
similarity set.
c1 and c2 must be filled in both structs prior to calling this function.

Input: A char representing the similarity set between two zygosity structs
       Two filled zygosity struct pointers representing the two sections to be compared.
Output: The total distance between two samples, such that DD/rr are a distance 
        of 2, and Dr/rr are a distance of 1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int get_distance_custom(char sim_set, struct zygosity *z1, struct zygosity *z2)
{
  int dist = 0;
  char Hor;
  char S;
  char ch_dist;

  S = sim_set;
  Hor = z1->ch | z2->ch;

  ch_dist = Hor | S;  // Force ones everywhere they are the same
  dist = get_zeros(S);  // Add one distance for every non-shared zygosity
  dist += get_zeros(ch_dist); // Add another one for every difference that has no heterozygotes

  return dist;
}
