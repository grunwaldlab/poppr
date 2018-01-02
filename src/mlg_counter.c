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
#include <string.h>
#include <stdlib.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <R.h>

int mlg_round_robin_cmpr (const void *a, const void *b);
void SampleWithoutReplacement(int populationSize, int sampleSize, int* samples);
SEXP mlg_round_robin(SEXP mat);
SEXP genotype_curve_internal(SEXP mat, SEXP iter, SEXP maxloci, SEXP report);
// global variable indicating the size of array to use for comparison in memcmp.
int NLOCI = 0;


/*
* The mask struct will hold a single individual represented by an array of m-1
* columns, where m represents the number of loci in the data set. The integer i
* represents the initial position of the sample in the data set for back
* reference.
* 
* Note that there is a potential extra value for the struct: n
* This vaule would contain the comparison size needed for memcmp. If this code
* were to be parallelized, make sure this part of the struct is available.
*/
struct mask {
  int* ind;
  int i;
  // int n;
};

int mlg_round_robin_cmpr (const void *a, const void *b){
  struct mask *ia = (struct mask *)a;
  struct mask *ib = (struct mask *)b;
  // Note: this depends on the global variable NLOCI. This potentially means
  // that it's not possible in its current form to go parallel unless the
  // variable is moved into the struct.
  return memcmp( (const int *)(ia->ind), (const int *)(ib->ind), NLOCI );
}


// Adapted from http://stackoverflow.com/a/311716/2752888
// Algorithm 3.4.2S by Donald Knuth
void SampleWithoutReplacement(int populationSize, int sampleSize, int* samples)
{
    
    // Use Knuth's variable names
    int n = sampleSize;
    int N = populationSize;
    int t = 0; // total input records dealt with
    int m = 0; // number of items selected so far
    double u;
    
    GetRNGstate();
    while (m < n)
    {
        u = unif_rand(); // call a uniform(0,1) random number generator

        if ( (N - t)*u >= n - m )
        {
            t++;
        }
        else
        {
            samples[m] = t;
            t++; m++;
        }
    }
    PutRNGstate();
}


/*
* This will be a function to calculate round-robin multilocus genotypes using
* qsort. It currenlty takes in an integer matrix and spits out a vector that
* is the same length as the number of columns indicating the number of unique
* genotypes when masking that column.
*
* Process:
* Make an array of mask structs representing each sample and fill them with
* their respective genotypes except for the final locus.
*
* For each locus in loci:
*
* 	If this is the first iteration:
* 		set the value of mask_col to the last column of the genotype array
* 	Else:
* 		set the value of mask_col to the previous column of the genotype array.
* 	sort the mask array
*
* 	For each sample:
* 		If the index is not zero:
* 			compare the current sample with the previous sample.
* 			If the samples are not equal:
* 				increment the multilocus genotype number
* 		Else:
* 			Initialize the multilocus genotype number to 1.
*
* 		Replace the current locus in the previous sample in the mask matrix with
* 		mask_col in the previous sample from the genotype matrix.
*
* 		If this is the final sample:
* 			replace the current locus in the current sample in the mask matrix with
* 			mask_col in the current sample from the genotype matrix.
*
* 	Add the multilocus genotype count to the output array
* 	Zero out the multilocus genotype count.
*/
SEXP mlg_round_robin(SEXP mat)
{
  SEXP Rout;
  SEXP Rdim;
  int rows;
  int cols;
  int i;
  int j;
  int is_missing;
  int new_genotype;
  int mask_position;
  int nmlg;
  int* genotype_matrix;
  struct mask* mask_matrix;
  
  Rdim = getAttrib(mat, R_DimSymbol);
  rows = INTEGER(Rdim)[0];
  cols = INTEGER(Rdim)[1];
  PROTECT(Rout = allocMatrix(INTSXP, rows, cols));
  
  NLOCI = (cols - 1)*sizeof(int);
  genotype_matrix = INTEGER(mat);
  
  mask_matrix = R_Calloc(rows, struct mask);
  for (i = 0; i < rows; i++)
  {
    mask_matrix[i].ind = R_Calloc(cols, int);
    mask_matrix[i].i = i;
    // Initializing the pointers.
    // The mask matrix will start by masking the first column of the genotype.
    for (j = 1; j < cols; j++)
    {
      mask_matrix[i].ind[j - 1] = genotype_matrix[i + j*rows];
      if (j == cols - 1)
      {
        mask_matrix[i].ind[j] = 0;
      }
    }
  }
  
  // For each iteration, the genotypes will be calculated and then the current
  // locus in the mask_matrix (which is j + 1) will be replaced with the jth 
  // locus in the genotype_matrix, effectively masking the j + 1 locus. 
  for (j = 0; j < cols; j++)
  {
    R_CheckUserInterrupt();
    qsort(mask_matrix, rows, sizeof(struct mask), mlg_round_robin_cmpr);

    for (i = 0; i < rows; i++)
    {
      if (i != 0)
      {
        if (memcmp(mask_matrix[i].ind, mask_matrix[i - 1].ind, sizeof(int)*(cols - 1)) != 0)
        {
          nmlg++;
        }
        mask_position = mask_matrix[i - 1].i;
        is_missing = genotype_matrix[mask_position + j*rows] == NA_INTEGER;
        new_genotype = (is_missing) ? 0 : genotype_matrix[mask_position + j*rows];
        mask_matrix[i - 1].ind[j] = new_genotype;
      }
      else
      {
        nmlg = 1;
      }
      INTEGER(Rout)[mask_matrix[i].i + j*rows] = nmlg;
      if (i == (rows - 1))
      {
        mask_position = mask_matrix[i].i;
        is_missing = genotype_matrix[mask_position + j*rows] == NA_INTEGER;
        new_genotype = (is_missing) ? 0 : genotype_matrix[mask_position + j*rows];
        mask_matrix[i].ind[j] = new_genotype;
      }
    }
    nmlg = 0;
  }
  
  for (i = 0; i < rows; i++)
  {
    R_Free(mask_matrix[i].ind);
  }
  R_Free(mask_matrix);
  
  UNPROTECT(1);
  return(Rout);
}

/*
* This function will randomly sample with replacement 1..m-1 loci and calculate
* the number of multilocus genotypes to give a genotype accumulation curve.
* 
* Input:
*   - mat an n x m integer matrix where n = number of samples and m = number of 
*       loci. The integers represent different genotypes.
*   - iter an integer specifying the number of iterations per locus. This will
*       be the number of rows in the the output matrix.
*   - maxloci the maximum number of loci to be analyzed.
*   - report an integer specifying after how many steps you want the function
*       to report progess.
* Output:
*   - A matrix with iter rows and m - 1 columns filled with counts of the number
*       of multilocus genotypes for j loci. 
*/
SEXP genotype_curve_internal(SEXP mat, SEXP iter, SEXP maxloci, SEXP report)
{
  SEXP Rout;
  SEXP Rdim;
  int rows;
  int cols;
  int iteration = 0;
  int nloci = 1;
  int i;
  int j;
  int REPORT;
  int is_missing;
  int new_genotype;
  int mask_position;
  int nmlg;
  int nmax;
  int* genotype_matrix;
  int* sampled_loci;
  int selected_locus;
  struct mask* mask_matrix;
  
  Rdim = getAttrib(mat, R_DimSymbol);
  rows = INTEGER(Rdim)[0];
  cols = INTEGER(Rdim)[1];
  nmax = (INTEGER(maxloci)[0] < cols - 1) ? INTEGER(maxloci)[0] : cols - 1;
  REPORT = INTEGER(report)[0];
  PROTECT(Rout = allocMatrix(INTSXP, INTEGER(iter)[0], nmax));
  
  
  genotype_matrix = INTEGER(mat);
  sampled_loci = R_Calloc(nmax, int);
  mask_matrix = R_Calloc(rows, struct mask);
  for (i = 0; i < rows; i++)
  {
    mask_matrix[i].ind = R_Calloc(nmax, int);
    mask_matrix[i].i = i;
  }
  // Step 1: Loop over the number of loci
  while (nloci < nmax + 1)
  {
    R_CheckUserInterrupt();
    // Initialize the global variable to be the number of loci we want to 
    // compare. 
    NLOCI = nloci*sizeof(int);
      
    // iterate the number of times defined by the user. 
    while (iteration < INTEGER(iter)[0])
    {
      // sampled_loci here is an array of integers specifying the columns to
      // copy from the genotype_matrix.
      SampleWithoutReplacement(cols, nloci, sampled_loci);
      
      // If it's the first iteration, the matrix needs to be initialized.
      // We have to get three things for this:
      //  1. The position of the sample in the genotype matrix
      //  2. The locus in the genotype matrix
      //  3. If the value is missing.
      if (iteration == 0)
      {
        for (i = 0; i < rows; i++)
        {
          mask_position = mask_matrix[i].i; // position in genotype_matrix
          for (j = 0; j < nloci; j++)
          {
            selected_locus = sampled_loci[j]*rows; // locus in genotype_matrix
            is_missing = genotype_matrix[mask_position + selected_locus] == NA_INTEGER;
            new_genotype = (is_missing) ? 0 : genotype_matrix[mask_position + selected_locus];
            mask_matrix[i].ind[j] = new_genotype;
          }
        }
        // Since it's the first iteration, sample again to set up the next 
        // iteration.
        SampleWithoutReplacement(cols, nloci, sampled_loci);
      }
      if (REPORT > 0 && (iteration + 1) % REPORT == 0)
      {
        Rprintf("\rCalculating genotypes for %2d/%d loci. Completed iterations: %3.0f%%", nloci, nmax, (float)((iteration + 1)*100)/(INTEGER(iter)[0]));
      }
      // Here, we sort the mask_matrix and then iterate through, counting up the
      // number of times we see a change in genotype. We also fill the matrix 
      // with the values for the next iteration.
      qsort(mask_matrix, rows, sizeof(struct mask), mlg_round_robin_cmpr);
      for (i = 0; i < rows; i++)
      {
        // Here, we compare the current sample with the previous to determine if
        // the number of MLGs needs to go up. 
        if (i != 0)
        {
          if (memcmp(mask_matrix[i].ind, mask_matrix[i - 1].ind, NLOCI) != 0)
          {
            nmlg++;
          }
          // We don't need the previous sample any more, so we replace it with 
          // the genotype for the next iteration.
          mask_position = mask_matrix[i - 1].i;
          for (j = 0; j < nloci; j++)
          {
            selected_locus = sampled_loci[j]*rows;
            is_missing = genotype_matrix[mask_position + selected_locus] == NA_INTEGER;
            new_genotype = (is_missing) ? 0 : genotype_matrix[mask_position + selected_locus];
            mask_matrix[i - 1].ind[j] = new_genotype;
          }
        }
        else
        {
          nmlg = 1;
        }
        // If this is the last sample, we replace it with the genotype for the
        // next iteration.
        if (i == (rows - 1))
        {
          mask_position = mask_matrix[i].i;
          for (j = 0; j < nloci; j++)
          {
            selected_locus = sampled_loci[j]*rows;
            is_missing = genotype_matrix[mask_position + selected_locus] == NA_INTEGER;
            new_genotype = (is_missing) ? 0 : genotype_matrix[mask_position + selected_locus];
            mask_matrix[i].ind[j] = new_genotype;
          }
        }
      }
      // After going through the samples, we add the value into our output 
      // matrix, reset the MLG counter and increase the iterations.
      INTEGER(Rout)[iteration + (nloci - 1)*INTEGER(iter)[0]] = nmlg;
      nmlg = 0;
      iteration++;
    }
    // Once all the iterations are done, reset the counter and increase the
    // number of loci.
    iteration = 0;
    nloci++;
  }
  
  for (i = 0; i < rows; i++)
  {
    R_Free(mask_matrix[i].ind);
  }
  R_Free(mask_matrix);
  UNPROTECT(1);
  return(Rout);
}
