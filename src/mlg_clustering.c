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

SEXP farthest_neighbor(SEXP dist, SEXP mlg, SEXP threshold)
{

  // TODO: More descriptive header
  // TODO: Make unit tests
  // TODO: Input validation

  // Assumptions:
  //  dist is an n by n matrix containing distances between individuals
  //  mlg is a vector of length n containing mlg assignments
  //  threshold is a real
  //  mlg is a vector of integers
  //  dist is a matrix of reals

  int num_individuals; // n above; total number of individuals
  int num_clusters; // Number of clusters; updated upon merge
  int num_mlgs; // Number of clusters already present
  int cur_mlg;  // Temporary number used for tracking mlgs
  double thresh; // Threshold for distance under which mlgs are clones
  double min_cluster_distance; // Used in finding the closest cluster pairing
  double tmp_max;     // Used in finding the closest cluster
  int closest_pair[2]; // Used in finding pair of clusters closest together
  int** cluster_matrix;
  double** cluster_distance_matrix;
  int* cluster_size; // Size of each cluster
  int* out_vector; // A copy of Rout for internal use

  SEXP Rout;
  SEXP Rdim;

  thresh = REAL(threshold)[0];
  Rdim = getAttrib(dist, R_DimSymbol);
  num_individuals = INTEGER(Rdim)[0]; // dist is a square matrix
  PROTECT(Rout = allocVector(INTSXP, num_individuals));  

  // TODO: Limit memory usage by counting initial mlgs from mlg
  num_mlgs = num_individuals;
  // Allocate empty matrix for storing clusters
  cluster_matrix = (int**)malloc(num_mlgs * sizeof(int*));
  cluster_distance_matrix = (double**)malloc(num_mlgs * sizeof(double*));
  if(cluster_matrix == NULL || cluster_distance_matrix == NULL)
  {
    // TODO: Throw error and exit
  }
  for(int i = 0; i < num_mlgs; i++)
  {
    cluster_matrix[i] = (int*)malloc(num_individuals * sizeof(int));
    memset(cluster_matrix[i], -1, sizeof(int)*num_individuals);
    cluster_distance_matrix[i] = (double*)malloc(num_individuals * sizeof(double));
    if(cluster_matrix[i] == NULL || cluster_distance_matrix[i] == NULL)
    {
      // TODO: Throw error and exit
    }
    for(int j = 0; j < num_individuals; j++)
    {
      // Using this instead of memset to preserve sentinel value
      cluster_distance_matrix[i][j] = -1.0;
    }
  }
  // Allocate memory for storing sizes of each cluster
  // calloc clears the memory allocated to 0s
  cluster_size = (int*)calloc(num_mlgs, sizeof(int));
  if(cluster_size == NULL)
  {
    // TODO: Throw error and exit
  }
  // Allocate memory for storing cluster assignments
  // calloc clears the memory allocated to 0s
  out_vector = (int*)calloc(num_individuals, sizeof(int));
  if(out_vector == NULL)
  {
    // TODO: Throw error and exit
  }

  // Fill initial clustering matrix via mlg
  // Steps through mlg.
  // Adds the index of each individual to cluster_matrix
  // Increments the cluster size for each individual added
  // For each new cluster added, increments num_clusters
  num_clusters = 0;
  for(int i = 0; i < num_individuals; i++)
  {
    cur_mlg = INTEGER(mlg)[i]; // Get the initial cluster of this individual / Casts as int implicitly
    // Insert index of individual in the result vector  
    out_vector[i] = cur_mlg-1;    
    // Then add this individual to the next open slot in its cluster
    cluster_matrix[cur_mlg-1][cluster_size[cur_mlg-1]] = i;
    cluster_size[cur_mlg-1]++;

    // If this is the first individual in this cluster, increment num_clusters
    if(cluster_size[cur_mlg-1] == 1)
    {
      num_clusters++;
    }
  }
  num_mlgs = num_clusters;

  // Main processing loop.
  // Finds the two closest cluster
  // then merges them together if they are within threshold of each other
  // Repeats until no clusters are withing threshold of each other
  //  or only one cluster remains
  min_cluster_distance = -1;
  while(min_cluster_distance < thresh && num_clusters > 1)
  {
    min_cluster_distance = -1;
    closest_pair[0] = -1;
    closest_pair[1] = -1;
    // Fill the distance matrix with the new distances between each cluster
    for(int i = 0; i < num_individuals; i++)
    {
      for(int j = i+1; j < num_individuals; j++)
      {
        if(out_vector[i] != out_vector[j])
        {
          if(REAL(dist)[(i)*num_individuals + (j)] > cluster_distance_matrix[out_vector[i]][out_vector[j]])
          {
            cluster_distance_matrix[out_vector[i]][out_vector[j]] = REAL(dist)[(i)*num_individuals + (j)];
          }
        } 
      }
    }
    for(int i = 0; i < num_mlgs; i++)
    {
      for(int j = 0; j < num_mlgs; j++)
      {
        if((cluster_distance_matrix[i][j] < min_cluster_distance && cluster_distance_matrix[i][j] > -0.5)
           || min_cluster_distance < -.05 && cluster_distance_matrix[i][j] > -0.5)
        {
          min_cluster_distance = cluster_distance_matrix[i][j];
          closest_pair[0] = i;
          closest_pair[1] = j;
        }
      }
    }
    // Merge the two closest clusters together into the lower numbered cluster
    //  if they are within the threshold distance from each other
    if(min_cluster_distance < 0 || closest_pair[0] < 0 || closest_pair[1] < 0)
    {
      printf("\nThe data resulted in a negative distance or cluster id. Please check your data.\n");
      num_clusters = 0;
    }
    else if(min_cluster_distance < thresh)
    {
      for(int i = 0; i < cluster_size[closest_pair[1]] && cluster_matrix[closest_pair[1]][i] > -1; i++)
      {
        // Change the assignment for this individual in the result vector
        out_vector[cluster_matrix[closest_pair[1]][i]] = closest_pair[0];
        // Then update the cluster_matrix and cluster_size given the moved individual
        cluster_matrix[closest_pair[0]][cluster_size[closest_pair[0]]] = cluster_matrix[closest_pair[1]][i];
        cluster_matrix[closest_pair[1]][i] = -1;
        cluster_size[closest_pair[0]]++;
      }
      cluster_size[closest_pair[1]] = 0;
      num_clusters--;
      // Erase the distance entries for all pairings containing the old cluster
      for(int i = 0; i < num_individuals; i++)
      {
        cluster_distance_matrix[closest_pair[1]][i] = -1.0;
        cluster_distance_matrix[i][closest_pair[1]] = -1.0;
      }
    }
  }

  // Fill return vector
  for(int i = 0; i < num_individuals; i++)
  {
    INTEGER(Rout)[i] = out_vector[i]+1;
  }

  // Free memory allocated for the various arrays and matrices
  for(int i = 0; i < num_mlgs; i++)
  { 
    free(cluster_matrix[i]);
    free(cluster_distance_matrix[i]);
  }
  free(cluster_matrix);
  free(cluster_distance_matrix);
  free(cluster_size);
  free(out_vector);
  UNPROTECT(1);
  
  return Rout;
}


