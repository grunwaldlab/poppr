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
#include <string.h>
#include <stdlib.h>
#include <Rinternals.h>
#include <R.h>

int mlg_round_robin_cmpr (const void *a, const void *b);
SEXP mlg_round_robin(SEXP mat);


/*
* The mask struct will hold a single individual represented by an array of m-1
* columns, where m represents the number of loci in the data set. The integer i
* represents the initial position of the sample in the data set for back
* reference.
*/
struct mask {
  int* ind;
  int i;
  int n;
};

int mlg_round_robin_cmpr (const void *a, const void *b){
  struct mask *ia = (struct mask *)a;
  struct mask *ib = (struct mask *)b;
  // Something here doesn't work :(
  return memcmp( (const int *)(ia->ind), (const int *)(ib->ind), ia->n );
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
  int k;
  int mask_col;
  int mask_position;
  int nmlg;
  int* genotype_matrix;
  struct mask* mask_matrix;
  
  Rdim = getAttrib(mat, R_DimSymbol);
  rows = INTEGER(Rdim)[0];
  cols = INTEGER(Rdim)[1];
  PROTECT(Rout = allocVector(INTSXP, cols));
  
  genotype_matrix = INTEGER(mat);
  
  mask_matrix = R_Calloc(rows, struct mask);
  for (i = 0; i < rows; i++)
  {
    mask_matrix[i].ind = R_Calloc(cols, int);
    mask_matrix[i].i = i;
    mask_matrix[i].n = (cols - 1)*sizeof(int);
    // Initializing the pointers.
    for (j = 0; j < cols - 1; j++)
    {
      mask_matrix[i].ind[j] = genotype_matrix[i + j*rows];
      if (j == cols - 2)
      {
        mask_matrix[i].ind[cols - 1] = 0;
      }
    }
  }
  
  for (j = 0; j < cols; j++)
  {
    mask_col = (j == 0) ? cols - 1 : j - 1;
    
    
    qsort(mask_matrix, rows, sizeof(struct mask), mlg_round_robin_cmpr);
    
    
    for (i = 0; i < rows; i++)
    {
      for (k = 0; k < cols; k++)
      {
        if (k == cols - 1)
        {
          Rprintf("%d\t", mask_matrix[i].i);
        }
        else
        {
          Rprintf("%d\t", mask_matrix[i].ind[k]);
        }
      }
      if (i != 0)
      {
        if (memcmp(mask_matrix[i].ind, mask_matrix[i - 1].ind, sizeof(int)*(cols - 1)) != 0)
        {
          nmlg++;
        }
        mask_position = mask_matrix[i - 1].i;
        mask_matrix[i - 1].ind[j] = genotype_matrix[mask_position + mask_col*rows];
      }
      else
      {
        nmlg = 1;
      }
      if (i == (rows - 1))
      {
        mask_position = mask_matrix[i].i;
        mask_matrix[i].ind[j] = genotype_matrix[mask_position + mask_col*rows];
      }
      Rprintf("\n");
    }
    Rprintf("\n");
    INTEGER(Rout)[mask_col] = nmlg;
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
