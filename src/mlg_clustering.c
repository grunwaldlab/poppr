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
#include <omp.h>

SEXP neighbor_clustering(SEXP dist, SEXP mlg, SEXP threshold, SEXP algorithm)
{
  // This function uses various clustering algorithms to
  // condense a set of multilocus genotypes into a potentially smaller
  // set. Clusters closer than threshold apart will be merged together,
  // and the process will repeat until no remaining clusters are close
  // enough to be merged. Farthest Neighbor clustering is the strictest,
  // (clustering only groups whose extremes are nearest) and is therefore 
  // used as the default. Nearest Neighbor is the least strict, and is
  // prone to linking vastly different individuals into the same cluster
  // if there is a path between them. Average Neighbor takes the average
  // distance between all pairs of points in two clusters, then merges
  // the pair of clusters with the smallest average distance.

  // Note: Input validation is handled in the R wrapper function.
  //  Take care if calling this function directly.

  // Assumptions:
  //  dist is an n by n matrix containing distances between individuals
  //  mlg is a vector of length n containing mlg assignments
  //  length(unique(mlg)) is at most n and at least 1
  //  threshold is a real
  //  mlg is a vector of integers with max value of n
  //  dist is a matrix of reals

  int num_individuals; // n above; total number of individuals
  int num_clusters; // Number of clusters; updated upon merge
  int num_mlgs; // Number of clusters already present
  int cur_mlg;  // Temporary number used for tracking mlgs
  double thresh; // Threshold for distance under which mlgs are clones
  double min_cluster_distance; // Used in finding the closest cluster pairing
  int closest_pair[2]; // Used in finding pair of clusters closest together
  int** cluster_matrix;
  double** cluster_distance_matrix;
  double*** private_distance_matrix;
  omp_lock_t** cluster_distance_locks;
  int* cluster_size; // Size of each cluster
  int* out_vector; // A copy of Rout for internal use
  double* dist_ij; // Variables to store distances inside loops
  double* dist_ji;
  char algo;  // Used for storing the first letter of algorithm

  SEXP Rout;
  SEXP Rdim;

  algo = *CHAR(STRING_ELT(algorithm,0));
  thresh = REAL(threshold)[0];
  Rdim = getAttrib(dist, R_DimSymbol);
  num_individuals = INTEGER(Rdim)[0]; // dist is a square matrix
  PROTECT(Rout = allocVector(INTSXP, num_individuals));  

  // Find the MLG with the highest index value
  // and use it as the number of MLGs in the data set 
  num_mlgs = 0;
  for(int i = 0; i < num_individuals; i++)
  {
    if(INTEGER(mlg)[i] > num_mlgs)
    {
      num_mlgs = INTEGER(mlg)[i];
    }
  }
  if(num_mlgs == 0)
  {
    num_mlgs = num_individuals;
  }
  
  // Allocate empty matrix for storing clusters
  cluster_matrix = R_Calloc(num_mlgs, int*);
  cluster_distance_matrix = R_Calloc(num_mlgs, double*);
  cluster_distance_locks = R_Calloc(num_mlgs, omp_lock_t*);

  for(int i = 0; i < num_mlgs; i++)
  {
    cluster_matrix[i] = R_Calloc(num_individuals, int);
    cluster_distance_matrix[i] = R_Calloc(num_individuals, double);
    cluster_distance_locks[i] = R_Calloc(num_individuals, omp_lock_t);
    for(int j = 0; j < num_individuals; j++)
    {
      // Using this instead of memset to preserve sentinel value
      cluster_distance_matrix[i][j] = -1.0;
      cluster_matrix[i][j] = -1;
      omp_init_lock(&cluster_distance_locks[i][j]);
    }
  }
  // Allocate memory for storing sizes of each cluster
  cluster_size = R_Calloc(num_mlgs, int);
  // Allocate memory for storing cluster assignments
  out_vector = R_Calloc(num_individuals, int);


  #pragma omp parallel \
      shared(num_individuals,num_mlgs)
  {
    #pragma omp single
    {
      private_distance_matrix = R_Calloc(omp_get_num_threads(), double**);
    }
    // All threads execute this loop
    private_distance_matrix[omp_get_thread_num()] = R_Calloc(num_mlgs, double*);
    for(int i = 0; i < num_mlgs; i++)
    {
      private_distance_matrix[omp_get_thread_num()][i] = R_Calloc(num_individuals, double);
      for(int j = 0; j < num_individuals; j++)
      {
        private_distance_matrix[omp_get_thread_num()][i][j] = -1.0;
      }
    }
  } // End parallel

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
    #pragma omp parallel private(dist_ij, dist_ji) \
      shared(out_vector,dist,thresh,algo,num_individuals,num_mlgs,cluster_size,cluster_distance_matrix)
    {
      for(int i = 0; i < num_individuals; i++)
      {
        // Parallel statistics for H3N2 data set
        //    No Locks:      50 seconds // New method. No locks, so unusable
        //    With Parallel: 60 seconds // No locks, so not usable
        //    With Critical: 65 seconds // Private matrices and critical section should make this safe
        //    Serial:       66 seconds
        //    With Locks:    90 seconds // New style of threading
        //    With Locks:    123 seconds
        //    Outer With Locks: 110 seconds
        //#pragma omp parallel for schedule(guided) \
        //  shared(dist,out_vector,private_distance_matrix,cluster_size,algo) \
        //  private(dist_ij,dist_ji)
        if(omp_get_thread_num() == out_vector[i] % omp_get_num_threads())
        {
          for(int j = i+1; j < num_individuals; j++)
          {
            if(out_vector[i] != out_vector[j])
            {
              dist_ij = &(private_distance_matrix[omp_get_thread_num()][out_vector[i]][out_vector[j]]);
              dist_ji = &(private_distance_matrix[omp_get_thread_num()][out_vector[j]][out_vector[i]]);
              //omp_set_lock(&cluster_distance_locks[out_vector[i]][out_vector[j]]);
              //omp_set_lock(&cluster_distance_locks[out_vector[j]][out_vector[i]]);
              if(ISNA(REAL(dist)[i*num_individuals + j]) || ISNAN(REAL(dist)[i*num_individuals + j]) 
                || !R_FINITE(REAL(dist)[i*num_individuals + j]))
              {
                error("Data set contains missing or invalid distances. Please check your data.\n");
              }
              else if(algo=='n' && ((REAL(dist)[(i)*num_individuals + (j)] < *dist_ij) || *dist_ij < -0.5))
              { // Nearest Neighbor clustering
                *dist_ij = REAL(dist)[(i)*num_individuals + (j)];
                *dist_ji = REAL(dist)[(i)*num_individuals + (j)];
              }
              else if(algo=='a')
              { // Average Neighbor clustering
                // The average distance will be sum(D(xi,yi))/(|x|*|y|)
                // Since |x| and |y| are constant for now, that term can be moved into the sum
                // Which lets us add the elements in one at a time divided by the product of cluster sizes
                if(*dist_ij < -0.5)
                { // This is the first pair to be considered between these two clusters
                  double portion = REAL(dist)[i*num_individuals+j] / (double)(cluster_size[out_vector[i]]*cluster_size[out_vector[j]]); 
                  *dist_ij = portion;
                  *dist_ji = portion;
                }
                else
                { 
                  double portion = REAL(dist)[i*num_individuals+j] / (double)(cluster_size[out_vector[i]]*cluster_size[out_vector[j]]); 
                  *dist_ij += portion;
                  *dist_ji += portion;
                }
              }
              else if(REAL(dist)[(i)*num_individuals + (j)] > *dist_ij)
              { // Farthest Neighbor clustering
                // This is the default, so it will execute even if the algorithm argument is invalid
                *dist_ij = REAL(dist)[(i)*num_individuals + (j)];
                *dist_ji = REAL(dist)[(i)*num_individuals + (j)];
                //printf("\n%f\n",private_distance_matrix[i][j]);
              }
              //omp_unset_lock(&cluster_distance_locks[out_vector[i]][out_vector[j]]);
              //omp_unset_lock(&cluster_distance_locks[out_vector[j]][out_vector[i]]);
            } 
          }
        }
      }
      // Merge the private distance matrices back into the global distance matrix 
      #pragma omp critical
      for(int i = 0; i < num_mlgs; i++)
      {
        for(int j = 0; j < num_mlgs; j++)
        {
          if(cluster_size[i] * cluster_size[j] != 0)
          {
            if(algo=='n')
            { // Min every element with this thread's distances
              if(private_distance_matrix[omp_get_thread_num()][i][j] > -0.5)
              {
                if(cluster_distance_matrix[i][j] < -0.5)
                {
                  cluster_distance_matrix[i][j] = private_distance_matrix[omp_get_thread_num()][i][j];
                }
                else
                { // Set the value to the min of the stored value and this thread's value
                  cluster_distance_matrix[i][j] = (private_distance_matrix[omp_get_thread_num()][i][j] < cluster_distance_matrix[i][j]) ? (private_distance_matrix[omp_get_thread_num()][i][j]) : (cluster_distance_matrix[i][j]);
                }
              }
            }
            else if(algo=='a')
            {
              if(private_distance_matrix[omp_get_thread_num()][i][j] > -0.5)
              {
                if(cluster_distance_matrix[i][j] < -0.5)
                {
                  cluster_distance_matrix[i][j] = private_distance_matrix[omp_get_thread_num()][i][j];
                }
                else
                {
                  cluster_distance_matrix[i][j] += private_distance_matrix[omp_get_thread_num()][i][j];
                }
              }
            }
            else // Or algo=='f'
            { // Max every element with this thread's distances
              cluster_distance_matrix[i][j] = (private_distance_matrix[omp_get_thread_num()][i][j] > cluster_distance_matrix[i][j]) ? (private_distance_matrix[omp_get_thread_num()][i][j]) : (cluster_distance_matrix[i][j]);
            }
          }
        }
      } // End critical section
    
    } // End parallel
    for(int i = 0; i < num_mlgs; i++)
    {
      for(int j = 0; j < num_mlgs; j++)
      {
        if((cluster_distance_matrix[i][j] < min_cluster_distance && cluster_distance_matrix[i][j] > -0.5)
           || (min_cluster_distance < -.05 && cluster_distance_matrix[i][j] > -0.5))
        {
          min_cluster_distance = cluster_distance_matrix[i][j];
          closest_pair[0] = i;
          closest_pair[1] = j;
        }
        // Erase this distance to prepare for the next loop
        cluster_distance_matrix[i][j] = -1.0;
      }
    }
    // Merge the two closest clusters together into the lower numbered cluster
    //  if they are within the threshold distance from each other
    if(min_cluster_distance < 0 || closest_pair[0] < 0 || closest_pair[1] < 0)
    {
      warning("\nThe data resulted in a negative or invalid distance or cluster id. Result vector is incomplete.\n");
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
      for(int i = 0; i < num_mlgs; i++)
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

  #pragma omp parallel
  {
    for(int i = 0; i < num_mlgs; i++)
    {
      R_Free(private_distance_matrix[omp_get_thread_num()][i]);
    }
    R_Free(private_distance_matrix[omp_get_thread_num()]);
  } // End parallel
  R_Free(private_distance_matrix);
  // Free memory allocated for the various arrays and matrices
  for(int i = 0; i < num_mlgs; i++)
  { 
    R_Free(cluster_matrix[i]);
    R_Free(cluster_distance_matrix[i]);

    for(int j = 0; j < num_individuals; j++)
    {
      omp_destroy_lock(&cluster_distance_locks[i][j]);
    }
    R_Free(cluster_distance_locks[i]);

  }
  R_Free(cluster_matrix);
  R_Free(cluster_distance_matrix);
  R_Free(cluster_distance_locks);
  R_Free(cluster_size);
  R_Free(out_vector);
  UNPROTECT(1);
  
  return Rout;
}


