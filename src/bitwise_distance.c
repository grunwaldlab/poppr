/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# This software was authored by Zhian N. Kamvar and Javier F. Tabima, graduate
# students at Oregon State University; and Dr. Nik Grünwald, an employee of
# USDA-ARS.
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
//  All SNPs are diploid.

struct zygosity
{
  char c1;  // A 32 bit fragment of one chromosome  
  char c2;  // The corresponding fragment from the outer chromosome

  char ch;  // Heterozygous sites are indicated by 1's
  char cd;  // Homozygous dominant sites are indicated by 1's
  char cr;  // Homozygous recessive sites are indicated by 1's
};

struct locus
{
  double d;  // Number of dominant alleles found at this locus across genotypes
  double r;  // Number of recessive alleles found at this locus across genotypes
  double h;  // Number of genotypes that are heterozygous at this locus
  double n;  // Number of genotypes that contributed data to this struct  

};

SEXP bitwise_distance_haploid(SEXP genlight, SEXP missing, SEXP requested_threads);
SEXP bitwise_distance_diploid(SEXP genlight, SEXP missing, SEXP differences_only, SEXP requested_threads);
SEXP association_index_haploid(SEXP genlight, SEXP missing, SEXP requested_threads, SEXP indices);
SEXP association_index_diploid(SEXP genlight, SEXP missing, SEXP differences_only, SEXP requested_threads, SEXP indices);
SEXP get_pgen_matrix_genind(SEXP genind, SEXP freqs, SEXP pops);
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

Input: A genlight object containing samples of happloids.
Output: A distance matrix representing the number of differences between each sample.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
SEXP bitwise_distance_haploid(SEXP genlight, SEXP missing, SEXP requested_threads)
{
  SEXP R_out;
  SEXP R_gen_symbol;
  SEXP R_chr_symbol;  
  SEXP R_nap_symbol;
  SEXP R_gen;
  int num_gens;
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

  
  PROTECT(R_gen_symbol = install("gen")); // Used for accessing the named elements of the genlight object
  PROTECT(R_chr_symbol = install("snp"));
  PROTECT(R_nap_symbol = install("NA.posi"));  

  // This will be a LIST of type LIST:RAW
  R_gen = getAttrib(genlight, R_gen_symbol);
  num_gens = XLENGTH(R_gen);

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
    R_chr1_1 = VECTOR_ELT(getAttrib(VECTOR_ELT(R_gen,i),R_chr_symbol),0); // Chromosome 1
    chr_length = XLENGTH(R_chr1_1);
    R_nap1 = getAttrib(VECTOR_ELT(R_gen,i),R_nap_symbol); // Vector of the indices of missing values 
    nap1_length = XLENGTH(R_nap1);
    // Loop through every other genotype
    #ifdef _OPENMP
    #pragma omp parallel for \
      private(j,cur_distance,R_chr2_1,R_nap2,next_missing_index_j,next_missing_j,next_missing_index_i,next_missing_i,\
              tmp_sim_set, k, mask, nap2_length) \
      shared(R_nap1, nap1_length, i, distance_matrix)
    #endif
    for(j = 0; j < i; j++)
    {
      cur_distance = 0;
      // These will be arrays of type RAW
      R_chr2_1 = VECTOR_ELT(getAttrib(VECTOR_ELT(R_gen,j),R_chr_symbol),0); // Chromosome 1
      R_nap2 = getAttrib(VECTOR_ELT(R_gen,j),R_nap_symbol); // Vector of the indices of missing values
      nap2_length = XLENGTH(R_nap2);
      if(nap2_length > 0)
      {
        next_missing_index_j = 0; // Next set of chromosomes start back at index 0
        next_missing_j = (int)INTEGER(R_nap2)[next_missing_index_j] - 1; //Compensate for Rs 1 based indexing
      }
      else
      {
        next_missing_index_j = -1; 
        next_missing_j = -1;
      }
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
      // Finally, loop through all the chunks of SNPs for these genotypes
      for(k = 0; k < chr_length; k++) 
      {
        // Find the sites where the two individuals differ. 1's for match, 0's for different
        tmp_sim_set = ~((char)RAW(R_chr1_1)[k] ^ (char)RAW(R_chr2_1)[k]);

        // Check for missing values and force them to match
        while(next_missing_index_i < nap1_length && next_missing_i < (k+1)*8 && next_missing_i >= (k*8) )
        {
          // Handle missing bit in both chromosomes and both samples with mask
          mask = 1 << (next_missing_i%8); 
          if(missing_match)
          {
            tmp_sim_set |= mask; // Force the missing bit to match
          }
          else
          {
            tmp_sim_set &= ~mask; // Force the missing bit to not match
          }
          next_missing_index_i++; 
          next_missing_i = (int)INTEGER(R_nap1)[next_missing_index_i] - 1;
        }       
        // Repeat for j
        while(next_missing_index_j < nap2_length && next_missing_j < (k+1)*8 && next_missing_j >= (k*8))
        {
          // Handle missing bit in both chromosomes and both samples with mask
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
        cur_distance += get_zeros(tmp_sim_set);
      }
      // Store the distance between these two genotypes in the distance matrix
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
zygosity.

Input: A genlight object containing samples of diploids.
Output: A distance matrix representing the number of differences in zygosity
        between each sample.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
SEXP bitwise_distance_diploid(SEXP genlight, SEXP missing, SEXP differences_only, SEXP requested_threads)
{
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

  
  PROTECT(R_gen_symbol = install("gen")); // Used for accessing the named elements of the genlight object
  PROTECT(R_chr_symbol = install("snp"));
  PROTECT(R_nap_symbol = install("NA.posi"));  

  // This will be a LIST of type LIST:RAW
  R_gen = getAttrib(genlight, R_gen_symbol);
  num_gens = XLENGTH(R_gen);

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
    R_chr1_1 = VECTOR_ELT(getAttrib(VECTOR_ELT(R_gen,i),R_chr_symbol),0); // Chromosome 1
    R_chr1_2 = VECTOR_ELT(getAttrib(VECTOR_ELT(R_gen,i),R_chr_symbol),1); // Chromosome 2
    chr_length = XLENGTH(R_chr1_1);
    R_nap1 = getAttrib(VECTOR_ELT(R_gen,i),R_nap_symbol); // Vector of the indices of missing values 
    nap1_length = XLENGTH(R_nap1);
    // Loop through every other genotype
    #ifdef _OPENMP
    #pragma omp parallel for \
      private(j,cur_distance,R_chr2_1,R_chr2_2,R_nap2,next_missing_index_j,next_missing_j,next_missing_index_i,next_missing_i,\
              set_1,set_2,tmp_sim_set, k, mask, nap2_length) \
      shared(R_nap1, nap1_length, i, distance_matrix)
    #endif
    for(j = 0; j < i; j++)
    {
      cur_distance = 0;
      // These will be arrays of type RAW
      R_chr2_1 = VECTOR_ELT(getAttrib(VECTOR_ELT(R_gen,j),R_chr_symbol),0); // Chromosome 1
      R_chr2_2 = VECTOR_ELT(getAttrib(VECTOR_ELT(R_gen,j),R_chr_symbol),1); // Chromosome 2
      R_nap2 = getAttrib(VECTOR_ELT(R_gen,j),R_nap_symbol); // Vector of the indices of missing values
      nap2_length = XLENGTH(R_nap2);
      if(nap2_length > 0)
      {
        next_missing_index_j = 0; // Next set of chromosomes start back at index 0
        next_missing_j = (int)INTEGER(R_nap2)[next_missing_index_j] - 1; //Compensate for Rs 1 based indexing
      }
      else
      {
        next_missing_index_j = -1;
        next_missing_j = -1;
      }
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
      // Finally, loop through all the chunks of SNPs for these genotypes
      for(k = 0; k < chr_length; k++) 
      {
        set_1.c1 = (char)RAW(R_chr1_1)[k];
        set_1.c2 = (char)RAW(R_chr1_2)[k];
        fill_zygosity(&set_1);

        set_2.c1 = (char)RAW(R_chr2_1)[k];
        set_2.c2 = (char)RAW(R_chr2_2)[k];
        fill_zygosity(&set_2);

        tmp_sim_set = get_similarity_set(&set_1,&set_2);

        // Check for missing values and force them to match
        while(next_missing_index_i < nap1_length && next_missing_i < (k+1)*8 && next_missing_i >= (k*8) )
        {
          // Handle missing bit in both chromosomes and both samples with mask
          mask = 1 << (next_missing_i%8); 
          if(missing_match)
          {
            tmp_sim_set |= mask; // Force the missing bit to match
            if((1&(set_1.c1>>(next_missing_i%8)))!=0 || (1&(set_1.c2>>(next_missing_i%8)))!=0)
            {
              Rprintf("\nLocus %d of genotype %d\tmod8:%d,actual:%d/%d\n",next_missing_i,i,next_missing_i%8,(1&(set_1.c1>>(next_missing_i%8))),(1&(set_1.c2>>(next_missing_i%8))));
            }
          }
          else
          {
            tmp_sim_set &= ~mask; // Force the missing bit to not match
          }
          next_missing_index_i++; 
          next_missing_i = (int)INTEGER(R_nap1)[next_missing_index_i] - 1;
        }       
        // Repeat for j
        while(next_missing_index_j < nap2_length && next_missing_j < (k+1)*8 && next_missing_j >= (k*8))
        {
          // Handle missing bit in both chromosomes and both samples with mask
          mask = 1 << (next_missing_j%8);
          if(missing_match)
          {
            tmp_sim_set |= mask; // Force the missing bit to match
            if((1&(set_2.c1>>(next_missing_j%8)))!=0 || (1&(set_2.c2>>(next_missing_j%8)))!=0)
            {
              Rprintf("\nLocus %d of genotype %d\tmod8:%d,actual:%d/%d\n",next_missing_j,j,next_missing_j%8,(1&(set_2.c1>>(next_missing_j%8))),(1&(set_2.c2>>(next_missing_j%8))));
            }
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
          cur_distance += get_zeros(tmp_sim_set);
        }
        else
        {
          cur_distance += get_distance_custom(tmp_sim_set,&set_1,&set_2);
        }
      }
      // Store the distance between these two genotypes in the distance matrix
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
       A boolean representing whether distances or differences should be counted.
       An integer representing the number of threads to be used.
       A vector of locus indices to be used in the calculations.
Output: The index of association for this genlight object over the specified loci
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
SEXP association_index_haploid(SEXP genlight, SEXP missing, SEXP requested_threads, SEXP indices)
{

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


  PROTECT(R_gen_symbol = install("gen")); // Used for accessing the named elements of the genlight object
  PROTECT(R_chr_symbol = install("snp"));
  PROTECT(R_nap_symbol = install("NA.posi"));  
  PROTECT(R_nloc_symbol = install("n.loc"));

  // This will be a LIST of type LIST:RAW
  R_gen = getAttrib(genlight, R_gen_symbol);
  num_gens = XLENGTH(R_gen);

  R_chr1_1 = VECTOR_ELT(getAttrib(VECTOR_ELT(R_gen,0),R_chr_symbol),0); // Chromosome 1
  num_chunks = XLENGTH(R_chr1_1);

  R_nloc = getAttrib(genlight,R_nloc_symbol);
  num_loci = INTEGER(R_nloc)[0];

  PROTECT(R_out = allocVector(REALSXP, 1));
  chunk_matrix = R_Calloc(num_gens,char*);
  for(i = 0; i < num_gens; i++)
  {
    chunk_matrix[i] = R_Calloc(num_chunks,char);
    R_chr1_1 = VECTOR_ELT(getAttrib(VECTOR_ELT(R_gen,i),R_chr_symbol),0); // Chromosome 1
    for(j = 0; j < num_chunks; j++)
    {
      chunk_matrix[i][j] = (char)RAW(R_chr1_1)[j];
    }
  }

  // This should be 8 for all chunks. If this assumption is wrong things will fail.
  chunk_length = 8;

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
  #pragma omp parallel for \
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
      R_nap1 = getAttrib(VECTOR_ELT(R_gen,j),R_nap_symbol); // Vector of the indices of missing values
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
          // Store it's inner-chunk location in the map
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
        Sn = (set_1.c1 ^ set_2.c1);

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
            // Store it's inner-chunk location in the map
            missing_mask_j |= (1<<(next_missing_j%8));
          }
          next_missing_index_j++;
        }
        // Use the missing data masks to correct for missing data
        if(missing_match)
        {
          Sn &= ~missing_mask_i; // Force the missing bit to 0 since Sn is matching-not
          Sn &= ~missing_mask_j; // Force the missing bit to 0 since Sn is matching-not
        }
        else
        {
          Sn |= missing_mask_i; // Force the missing bit to 1 since SN is matching-not
          Sn |= missing_mask_j; // Force the missing bit to 1 since SN is matching-not
        }
        // Loop through 0 through chunk_length, call this x
        for(x = 0; x < chunk_length; x++)
        {
          // Find the offset for the next bit/locus of interest 
          offset = x;
          // Create mask to retrieve just that bit/locus
          mask = 1<<offset;
          // If only_differences is set, subtract (mask&Hs)>>offset back off of val
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
  #pragma omp parallel for reduction(+ : D,D2) private(i,j) 
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
  //#pragma omp parallel for private(i) shared(Nc2,vars,M2,M,num_loci)
  for(i = 0; i < num_loci; i++)
  {
    vars[i] = (M2[i] - (M[i]*M[i])/Nc2) / Nc2;
  }
  // Calculate the expected variance
  Ve = 0;
  //#pragma omp parallel for reduction(+ : Ve) private(i)
  for(i = 0; i < num_loci; i++)
  {
    Ve += vars[i];
  }

  // Calculate the denominator for the index of association
  denom = 0;
  #ifdef _OPENMP
  #pragma omp parallel for reduction(+ : denom) private(i, j)
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
       A vector of locus indices to be used in the calculations.
Output: The index of association for this genlight object over the specified loci
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
SEXP association_index_diploid(SEXP genlight, SEXP missing, SEXP differences_only, SEXP requested_threads, SEXP indices)
{

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


  PROTECT(R_gen_symbol = install("gen")); // Used for accessing the named elements of the genlight object
  PROTECT(R_chr_symbol = install("snp"));
  PROTECT(R_nap_symbol = install("NA.posi"));  
  PROTECT(R_nloc_symbol = install("n.loc"));

  // This will be a LIST of type LIST:RAW
  R_gen = getAttrib(genlight, R_gen_symbol);
  num_gens = XLENGTH(R_gen);

  R_chr1_1 = VECTOR_ELT(getAttrib(VECTOR_ELT(R_gen,0),R_chr_symbol),0); // Chromosome 1
  num_chunks = XLENGTH(R_chr1_1);

  R_nloc = getAttrib(genlight,R_nloc_symbol);
  num_loci = INTEGER(R_nloc)[0];

  PROTECT(R_out = allocVector(REALSXP, 1));
  chunk_matrix = R_Calloc(num_gens*2,char*);
  for(i = 0; i < num_gens; i++)
  {
    chunk_matrix[i*2] = R_Calloc(num_chunks,char);
    chunk_matrix[i*2+1] = R_Calloc(num_chunks,char);
    R_chr1_1 = VECTOR_ELT(getAttrib(VECTOR_ELT(R_gen,i),R_chr_symbol),0); // Chromosome 1
    R_chr1_2 = VECTOR_ELT(getAttrib(VECTOR_ELT(R_gen,i),R_chr_symbol),1); // Chromosome 2
    for(j = 0; j < num_chunks; j++)
    {
      chunk_matrix[i*2][j] = (char)RAW(R_chr1_1)[j];
      chunk_matrix[i*2+1][j] = (char)RAW(R_chr1_2)[j];
    }
  }

  // This should be 8 for all chunks. If this assumption is wrong things will fail.
  chunk_length = 8;

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
  #pragma omp parallel for \
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
      R_nap1 = getAttrib(VECTOR_ELT(R_gen,j),R_nap_symbol); // Vector of the indices of missing values
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
          // Store it's inner-chunk location in the map
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
            // Store it's inner-chunk location in the map
            missing_mask_j |= (1<<(next_missing_j%8));
          }
          next_missing_index_j++;
        }
        // Use the missing data masks to correct for missing data
        if(missing_match)
        {
          Sn &= ~missing_mask_i; // Force the missing bit to 0 since Sn is matching-not
          Hs &= ~missing_mask_i; // Hs is differing homozygote sites, so missing should be 0
          Sn &= ~missing_mask_j; // Force the missing bit to 0 since Sn is matching-not
          Hs &= ~missing_mask_j; // Hs is differing homozygote sites, so missing should be 0
        }
        else
        {
          Sn |= missing_mask_i; // Force the missing bit to 1 since SN is matching-not
          Hs |= (~set_2.ch & missing_mask_i);// Hs should be 1 if set2.ch is 0 at this bit
                                  // This makes missing data a distance of 1 from a heterozygote site
                                  // and 2 from either homozygote site.
          Sn |= missing_mask_j; // Force the missing bit to 1 since SN is matching-not
          Hs |= (~set_1.ch & missing_mask_j);// Hs should be 1 if set2.ch is 0 at this bit
                                  // This makes missing data a distance of 1 from a heterozygote site
                                  // and 2 from either homozygote site.

        }
        // Loop through 0 through chunk_length, call this x
        for(x = 0; x < chunk_length; x++)
        {
          // Find the offset for the next bit/locus of interest 
          offset = x;
          // Create mask to retrieve just that bit/locus
          mask = 1<<offset;
    // val = (mask&Sn)>>offset + (mask&Hs)>>offset // Add one for not same, add one for opposite homozygotes
          // If only_differences is set, subtract (mask&Hs)>>offset back off of val
          val = (mask&Sn)>>offset;
          if(!only_differences)
          {
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
  #pragma omp parallel for reduction(+ : D,D2) private(i,j) 
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
  #pragma omp parallel for reduction(+ : denom) private(i, j)
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
Output: A matrix containing the Pgen value of each genotype at each locus.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
SEXP get_pgen_matrix_genind(SEXP genind, SEXP freqs, SEXP pops)
{
  SEXP R_out;
  SEXP R_tab_symbol;
  SEXP R_loc_symbol;
  SEXP R_pop_symbol;
  SEXP R_tab;
  PROTECT(R_tab_symbol = install("tab")); // Used for accessing the named elements of the genind object
  PROTECT(R_loc_symbol = install("loc.names"));
  PROTECT(R_pop_symbol = install("pop.names"));
  double *pgens;
  int *indices;
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
  num_pops = XLENGTH(getAttrib(genind, R_pop_symbol));
  missing_last = 0;
  size = num_gens*num_loci;
  pgens = R_Calloc(size, double);
  indices = R_Calloc(size*2, int);
  PROTECT(R_out = allocVector(REALSXP,size));

  // Fill indices with the indices inside freqs of alleles found in each genotype
  // Note that R_tab is column-major ordered due to being passed from R
  for(int i = 0; i < num_gens; i++)
  {
    index = 0;
    for(int j = 0; j < num_alleles; j++)
    {
      if(ISNA(REAL(R_tab)[i + j*num_gens]))
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
      else if(REAL(R_tab)[i + j*num_gens] > 0.75)
      {
        indices[i*num_loci*2 + index] = j;
        indices[i*num_loci*2 + index+1] = j;
        index += 2;
        missing_last = 0;
      }
      else if(REAL(R_tab)[i+ j*num_gens] < 0.75 && REAL(R_tab)[i + j*num_gens] > 0.25)
      {
        indices[i*num_loci*2 + index] = j;
        index += 1;
        missing_last = 0;
      }
    }
  } 

  for(int i = 0; i < num_gens; i++)
  {
    pop = INTEGER(pops)[i]-1;
    index = 0;
    for(int j = 0; j < num_loci; j++)
    {
      if(indices[i*num_loci*2 + index] == -1) // And therefore also [~ index+1] == -1
      {
        // Set the pgen value of this genotype at this locus to 0
        pgens[i*num_loci +j] = 1;
      }
      else
      {
        pgens[i*num_loci + j] = log(REAL(freqs)[pop + indices[i*num_loci*2 + index]*num_pops]) + log(REAL(freqs)[pop + indices[i*num_loci*2 + index+1]*num_pops]);
        // Account for both permutations of heterozygous loci
        if(indices[i*num_loci*2 + index] != indices[i*num_loci*2 + index+1])
        {
          pgens[i*num_loci + j] += log(2);
        }
      }
      index += 2;
    }
  }

  for(int i = 0; i < num_gens; i++)
  {
    for(int j = 0; j < num_loci; j++)
    {
      // Transpose the matrix into column-major ordering before returning to R
      REAL(R_out)[i + j*num_gens] = (pgens[i*num_loci + j]);
    }
  }

  R_Free(indices);
  R_Free(pgens);
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
