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


/*
* This will be a function to calculate round-robin multilocus genotypes using
* qsort.
*/
SEXP mlg_round_robin(SEXP mat)
{
  SEXP Rout;
  SEXP Rdim;
  int rows;
  int cols;
  int i;
  int j;
  int mask_position;
  int nmlg;
  int* genotype_matrix;
  int** mask_matrix;
  
  Rdim = getAttrib(mat, R_DimSymbol);
  rows = INTEGER(Rdim)[0];
  cols = INTEGER(Rdim)[1];
  PROTECT(Rout = allocVector(INTSXP, cols));
  
  genotype_matrix = INTEGER(mat);
  mask_matrix = R_Calloc(rows + 1, int*);
  for (i = 0; i < rows + 1; i++)
  {
    mask_matrix[i] = R_Calloc(cols + 1, int);
  }
  
  for (i = 0; i < rows; i++)
  {
    for (j = 0; j < cols; j++)
    {
      mask_matrix[i][j] = genotype_matrix[i + j*rows];
    }
    mask_matrix[i][cols] = i;
  }
  
  int mlg_round_robin_cmpr (const void *a, const void *b){
    return memcmp(*(const int**)a, *(const int**)b, sizeof(int)*cols);
  }
  
  for (j = 0; j < cols; j++)
  {
    for (i = 0; i < rows; i++)
    {
      mask_matrix[i][j] = 0;
    }
    qsort(mask_matrix, rows, sizeof(int*), mlg_round_robin_cmpr);
    
    for (i = 0; i < rows; i++)
    {
      if (i != 0)
      { 
        if (mlg_round_robin_cmpr(mask_matrix[i], mask_matrix[i - 1]) != 0)
        {
          nmlg++;
        }
        mask_position = mask_matrix[i - 1][cols];
        mask_matrix[i - 1][j] = genotype_matrix[mask_position + j*rows];
      }
      else
      {
        nmlg = 1;
      }
      if (i == (rows - 1))
      {
        mask_position = mask_matrix[i][cols];
        mask_matrix[i][j] = genotype_matrix[mask_position + j*rows];
      }
    }
    INTEGER(Rout)[j] = nmlg;
    nmlg = 0;
  }
  
  for (i = 0; i < rows; i++)
  {
    R_Free(mask_matrix[i]);
  }
  R_Free(mask_matrix);
  
  UNPROTECT(1);
  return(Rout);
  
}