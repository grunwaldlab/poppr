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
#include <R_ext/Utils.h>
#include <Rdefines.h>
#include <R.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>

// Thu Apr 13 08:42:12 2017 ------------------------------
// This code produces bugs when run on Fedora with multiple threads. Because of
// this, I am removing OMP functionality here and it will run serially. Details
// appear here: https://github.com/grunwaldlab/poppr/issues/138
// Include openMP if compiled with an openMP compatible compiler
// #ifdef _OPENMP
// #include <omp.h>
// #endif


SEXP neighbor_clustering(SEXP dist, SEXP mlg, SEXP threshold, SEXP algorithm, SEXP requested_threads);
void fill_distance_matrix(double** cluster_distance_martix, double*** private_distance_matrix, int* out_vector, int* cluster_size, SEXP dist, char algo, int num_individuals, int num_mlgs, int num_threads);

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Reassigns genotypes from mlg into new clusters based on a minimum genetic distance
between clusters using one of three clustering algorithms. 
The result will be a new list of mll assignments based on those rules.

Input: A square matrix of numeric distances between each pair of samples.
       A vector of initial mlg assignments for each sample.
       A real used to govern the minimum distance between clusters/mlls.
       A string determining the algorithm to use. Must start with a lowercase
        "n", "f", or "a". Representing "nearest neighbor", "farthest neighbor",
        and "average neighbor" (otherwise known as UPGMA) respectively.
       An integer representing the number of threads that should be used.
Output: A vector of mll assignments based on the algorithm and threshold used.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
SEXP neighbor_clustering(SEXP dist, SEXP mlg, SEXP threshold, SEXP algorithm, SEXP requested_threads)
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
  //  algorithm is a lowercase string starting with "n", "f", or "a"

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
  int* cluster_size; // Size of each cluster
  int* out_vector; // A copy of Rout for internal use
  int num_threads;
  char algo;  // Used for storing the first letter of algorithm

  SEXP Rout;
  SEXP Rout_vects;
  SEXP Rout_stats;
  SEXP Rout_dists;
  SEXP Rout_sizes;
  SEXP Rdim;

  // Convert the R object arguments into C data types
  algo = *CHAR(STRING_ELT(algorithm,0));
  thresh = REAL(threshold)[0];
  Rdim = getAttrib(dist, R_DimSymbol);
  num_individuals = INTEGER(Rdim)[0]; // dist is a square matrix

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
  
  // Generate and protect the output R object as well as its components
  PROTECT(Rout_vects = allocVector(INTSXP, num_individuals)); // MLG assignments 
  PROTECT(Rout_stats = allocVector(REALSXP, num_mlgs));       // Threshold for each merge
  PROTECT(Rout_dists = allocMatrix(REALSXP, num_mlgs, num_mlgs)); // Resulting distance matrix
  PROTECT(Rout_sizes = allocVector(INTSXP,  num_mlgs));           // Sizes of new clusters
  //PROTECT(Rout = CONS(Rout_vects, CONS(Rout_stats, CONS(Rout_dists, CONS(Rout_sizes, R_NilValue))))); 
  PROTECT(Rout = allocVector(VECSXP, 4));
  // Allocate empty matrix for storing clusters
  cluster_matrix = R_Calloc(num_mlgs, int*);
  cluster_distance_matrix = R_Calloc(num_mlgs, double*);

  // Initialize output and intermediate matrices with sentinel value -1
  for(int i = 0; i < num_mlgs; i++)
  {
    REAL(Rout_stats)[i] = -1;
    INTEGER(Rout_sizes)[i] = -1;
    cluster_matrix[i] = R_Calloc(num_individuals, int);
    for(int j = 0; j < num_individuals; j++)
    {
      // Using this instead of memset to preserve sentinel value
      cluster_matrix[i][j] = -1;
    }
    cluster_distance_matrix[i] = R_Calloc(num_mlgs, double);
    for(int j = 0; j < num_mlgs; j++)
    {
      // Using this instead of memset to preserve sentinel value
      cluster_distance_matrix[i][j] = -1.0;
    }
  }
  // Allocate memory for storing sizes of each cluster
  cluster_size = R_Calloc(num_mlgs, int);
  // Allocate memory for storing cluster assignments
  out_vector = R_Calloc(num_individuals, int);
  
  // Thu Apr 13 08:31:55 2017 ------------------------------
  // I am commenting out the following code because of an issue where this
  // code produces bugs on Fedora with 2 threads. Details of the failure
  // appear here: https://github.com/grunwaldlab/poppr/issues/138
  // 
  /*
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
  */
  num_threads = 1;
  // Initialize a 3D array so that each thread can work on its
  // own private distance matrix
  private_distance_matrix = R_Calloc(num_threads, double**);
  for(int i = 0; i < num_threads; i++)
  {
    private_distance_matrix[i] = R_Calloc(num_mlgs, double*);
    for(int j = 0; j < num_mlgs; j++)
    {
      private_distance_matrix[i][j] = R_Calloc(num_mlgs, double);
      for(int k = 0; k < num_mlgs; k++)
      {
        private_distance_matrix[i][j][k] = -1.0;
      }
    }
  }
  
  // Fill initial clustering matrix via mlg
  // Steps through mlg.
  // Adds the index of each individual to cluster_matrix
  // Increments the cluster size for each individual added
  // For each new cluster added, increments num_clusters
  num_clusters = 0;
  for(int i = 0; i < num_individuals; i++)
  {
    cur_mlg = INTEGER(mlg)[i]; // Get the initial cluster of this individual / Casts as int 
    // Insert index of individual in the result vector  
    out_vector[i] = cur_mlg-1;    
    // Then add this individual's index location to the next open slot in its cluster
    cluster_matrix[cur_mlg-1][cluster_size[cur_mlg-1]] = i;
    // And increase the size of this cluster
    cluster_size[cur_mlg-1]++;

    // If this is the first individual in this cluster, increment num_clusters
    if(cluster_size[cur_mlg-1] == 1)
    {
      num_clusters++;
    }
  }

  // Main processing loop.
  // Finds the two closest clusters
  // then merges them together if they are within threshold of each other
  // Repeats until no clusters are withing threshold of each other
  //  or only one cluster remains
  min_cluster_distance = -1;
  while(min_cluster_distance < thresh && num_clusters > 1)
  {
    R_CheckUserInterrupt();
    min_cluster_distance = -1;
    closest_pair[0] = -1;
    closest_pair[1] = -1;
    // Fill the distance matrix with the new distances between each cluster
    fill_distance_matrix(cluster_distance_matrix,private_distance_matrix,out_vector,cluster_size,dist,algo,num_individuals,num_mlgs,num_threads);
    // Loop through each pairing of MLGs to find the pair whose clusters are separated by the smallest distance
    for(int i = 0; i < num_mlgs; i++)
    {
      for(int j = 0; j < num_mlgs; j++)
      {
        // Check if this might be the minimum distance between two clusters
        if((cluster_distance_matrix[i][j] < min_cluster_distance && cluster_distance_matrix[i][j] > -0.5)
           || (min_cluster_distance < -0.5 && cluster_distance_matrix[i][j] > -0.5))
        {
          min_cluster_distance = cluster_distance_matrix[i][j];
          closest_pair[0] = i;
          closest_pair[1] = j;
        }
      }
    }
    // Merge the two closest clusters together based on which is closer to the other clusters
    //  if both clusters are within the threshold distance from each other
    if(min_cluster_distance < 0 || closest_pair[0] < 0 || closest_pair[1] < 0)
    {
      warning("\nThe data resulted in a negative or invalid distance or cluster id. Result vector is incomplete.\n");
      num_clusters = 0;
    }
    else if(min_cluster_distance < thresh)
    {
      // Store the distance at which this merge occurred in reverse order, since we
      // want them listed from largest distance to smallest
      REAL(Rout_stats)[num_mlgs-num_clusters] = min_cluster_distance;


      // Determine which cluster should "survive" the merger, and which should be consumed
      // Merge to the cluster that is closest to all other clusters
      //  closest_pair[0] will survive, and closest_pair[1] will move to closest_pair[0]
      double mean0 = 0.0;
      double mean1 = 0.0;
      for(int i = 0; i < num_clusters; i++)
      {
        if(cluster_distance_matrix[closest_pair[0]][i] > 0)
        {
          mean0 += cluster_distance_matrix[closest_pair[0]][i] / num_clusters;
        }
        if(cluster_distance_matrix[closest_pair[1]][i] > 0)
        {
          mean1 += cluster_distance_matrix[closest_pair[1]][i] / num_clusters;
        }
      }
      // If cluster 1 is closer to all others than cluster 0 is, switch the two. Otherwise leave 0 as the "host"
      if(mean1 < mean0)
      {
        int tmp;
        tmp = closest_pair[0];
        closest_pair[0] = closest_pair[1];
        closest_pair[1] = tmp;
      }

      // Now merge the two together, collapsing the (new) closest_pair[1] into closest_pair[0]
      // This is done by stepping through all the individuals in closest_pair[1] and assigning 
      //  them to closest_pair[0] instead.
      for(int i = 0; i < cluster_size[closest_pair[1]] && cluster_matrix[closest_pair[1]][i] > -1; i++)
      {
        // Change the assignment for this individual in the result vector
        out_vector[cluster_matrix[closest_pair[1]][i]] = closest_pair[0];
        // Then update the cluster_matrix and cluster_size given the moved individual
        cluster_matrix[closest_pair[0]][cluster_size[closest_pair[0]]] = cluster_matrix[closest_pair[1]][i];
        cluster_matrix[closest_pair[1]][i] = -1;
        cluster_size[closest_pair[0]]++;
      }
      // Now effectively erase the cluster that was merged into closest_pair[0]
      cluster_size[closest_pair[1]] = 0;
      num_clusters--;
      // Erase distance matrix to prepare for the next loop
      for(int i = 0; i < num_mlgs; i++)
      {
        for(int j = 0; j < num_mlgs; j++)
        {
          cluster_distance_matrix[i][j] = -1.0;
        }
      }
    }
  }

  // Fill return vector
  for(int i = 0; i < num_individuals; i++)
  {
    INTEGER(Rout_vects)[i] = out_vector[i]+1;
  }
  // Fill return distance matrix with updated cluster_distance_matrix
  fill_distance_matrix(cluster_distance_matrix,private_distance_matrix,out_vector,cluster_size,dist,algo,num_individuals,num_mlgs,num_threads);
  for(int i = 0; i < num_mlgs; i++)
  {
    // Fill return sizes
    INTEGER(Rout_sizes)[i] = cluster_size[i];
    //Fill return distance matrix
    for(int j = 0; j < num_mlgs; j++)
    {
      if(i == j)
      {
        REAL(Rout_dists)[i + j*num_mlgs] = 0;  
      }
      else if(cluster_distance_matrix[i][j] < -0.5)
      {
        REAL(Rout_dists)[i +j*num_mlgs] = NA_REAL;
      }
      else
      {
        REAL(Rout_dists)[i + j*num_mlgs] = cluster_distance_matrix[i][j];
      }
    }
  }
  for(int i = 0; i < num_threads; i++)
  {
    for(int j = 0; j < num_mlgs; j++)
    {
      R_Free(private_distance_matrix[i][j]);
    }
    R_Free(private_distance_matrix[i]);
  }
  R_Free(private_distance_matrix);
  
  // Free memory allocated for the various arrays and matrices
  for(int i = 0; i < num_mlgs; i++)
  { 
    R_Free(cluster_matrix[i]);
    R_Free(cluster_distance_matrix[i]);
  }
  R_Free(cluster_matrix);
  R_Free(cluster_distance_matrix);
  R_Free(cluster_size);
  R_Free(out_vector);
  
  SET_VECTOR_ELT(Rout, 0, Rout_vects);
  SET_VECTOR_ELT(Rout, 1, Rout_stats);
  SET_VECTOR_ELT(Rout, 2, Rout_dists);
  SET_VECTOR_ELT(Rout, 3, Rout_sizes);
  
  UNPROTECT(5);
  
  return Rout;
}

// Fill the distance matrix given the current cluster assignments
void fill_distance_matrix(double** cluster_distance_matrix, double*** private_distance_matrix, int* out_vector, int* cluster_size, SEXP dist, char algo, int num_individuals, int num_mlgs, int num_threads)
{
  double* dist_ij; // Variables to store distances inside loops
  double* dist_ji;
  int thread_id;

  // Thu Apr 13 08:37:40 2017 ------------------------------
  // Commenting this section out. See https://github.com/grunwaldlab/poppr/issues/138
  /*
  #ifdef _OPENMP
  {
    // Set the number of threads to be used in each omp parallel region
    omp_set_num_threads(num_threads);
  }
  #else
  {
    // Make sure it works in serial
    num_threads = 1;
  }
  #endif
  #ifdef _OPENMP
  #pragma omp parallel private(dist_ij, dist_ji, thread_id) shared(out_vector,dist,algo,num_individuals,num_mlgs,cluster_size,cluster_distance_matrix,num_threads)
  #endif
  */
  num_threads = 1;
  {
    /*
    #ifdef _OPENMP
    {
      thread_id = omp_get_thread_num();
    }
    #else
    {
      thread_id = 0;
    }
    #endif
    */
    thread_id = 0;
    // Clear the distance matrices before filling them
    for(int i = 0; i < num_mlgs; i++)
    {
      for(int j = 0; j < num_mlgs; j++)
      {
        private_distance_matrix[thread_id][i][j] = -1.0;
        if(thread_id==0)
        {
          cluster_distance_matrix[i][j] = -1.0;
        }
      }
    }
    for(int i = 0; i < num_individuals; i++)
    {
      // Divvy up the work via thread number and genotype assignment
      if(thread_id == out_vector[i] % num_threads)
      {
        for(int j = i+1; j < num_individuals; j++)
        {
          if(out_vector[i] != out_vector[j])
          {
            // Get the distance, as seen by this thread, between the clusters containing
            //  samples i and j.
            dist_ij = &(private_distance_matrix[thread_id][out_vector[i]][out_vector[j]]);
            dist_ji = &(private_distance_matrix[thread_id][out_vector[j]][out_vector[i]]);
            if(ISNA(REAL(dist)[(i) + (j)*num_individuals]) || ISNAN(REAL(dist)[(i) + (j)*num_individuals]) 
              || !R_FINITE(REAL(dist)[(i) + (j)*num_individuals]))
            {
              error("Data set contains missing or invalid distances. Please check your data.\n");
            }
            else if(algo=='n' && ((REAL(dist)[(i) + (j)*num_individuals] < *dist_ij) || *dist_ij < -0.5))
            { // Nearest Neighbor clustering
              // These will update and store the smallest distance between an individual in one cluster
              //  and an individual in another. Note that the pointer nature of dist_ij means this is
              //  actually updating the private distance matrix itself in order to keep track of the
              //  distances between every pairing of clusters.
              *dist_ij = REAL(dist)[(i) + (j)*num_individuals];
              *dist_ji = REAL(dist)[(i) + (j)*num_individuals];
            }
            else if(algo=='a')
            { // Average Neighbor clustering, otherwise known as UPGMA
              // The average distance will be sum(D(xi,yi))/(|x|*|y|)
              // Since |x| and |y| are constant for now, that term can be moved into the sum
              // Which lets us add the elements in one at a time divided by the product of cluster sizes
              if(*dist_ij < -0.5)
              { // This is the first pair to be considered between these two clusters
                double portion = REAL(dist)[i + j*num_individuals] / (double)(cluster_size[out_vector[i]]*cluster_size[out_vector[j]]); 
                *dist_ij = portion;
                *dist_ji = portion;
              }
              else
              { 
                // This is adding to the existing value in order to find the mean distance between all
                // individuals in cluster a with all individuals in cluster b, for all combinations of a and b.
                double portion = REAL(dist)[i + j*num_individuals] / (double)(cluster_size[out_vector[i]]*cluster_size[out_vector[j]]); 
                *dist_ij += portion;
                *dist_ji += portion;
              }
            }
            else if(algo=='f' && REAL(dist)[(i) + (j)*num_individuals] > *dist_ij)
            { // Farthest Neighbor clustering
              // This functions exactly like Nearest Neighbor, but using the maximum distance between
              // any individual in one cluster to any individual in another.
              *dist_ij = REAL(dist)[(i) + (j)*num_individuals];
              *dist_ji = REAL(dist)[(i) + (j)*num_individuals];
            }
          } 
        }
      }
    }
    // Merge the private distance matrices back into the global distance matrix 
    // Since each thread has its own private matrix, all matrices need to be merged
    //  to get the final distance matrix.
    // Thu Apr 13 08:39:24 2017 ------------------------------
    // see https://github.com/grunwaldlab/poppr/issues/138
    /*
    #ifdef _OPENMP
    #pragma omp barrier
    #pragma omp critical
    #endif
    */
    for(int i = 0; i < num_mlgs; i++)
    {
      for(int j = 0; j < num_mlgs; j++)
      {
        if(cluster_size[i] * cluster_size[j] != 0)
        {
          if(algo=='n')
          { // Min every element with this thread's distances, since the one we want in the 
            // final matrix is the minimum distance found by any thread.
            if(private_distance_matrix[thread_id][i][j] > -0.5)
            {
              if(cluster_distance_matrix[i][j] < -0.5)
              {
                cluster_distance_matrix[i][j] = private_distance_matrix[thread_id][i][j];
              }
              else
              { // Set the value to the min of the stored value and this thread's value
                // The ternary operator ?: functions like an in-line if else statement
                //  in this case updating cluster_distance_matrix with the smaller of the two operands
                cluster_distance_matrix[i][j] = (private_distance_matrix[thread_id][i][j] < cluster_distance_matrix[i][j]) ? (private_distance_matrix[thread_id][i][j]) : (cluster_distance_matrix[i][j]);
              }
            }
          }
          else if(algo=='a')
          {
            if(private_distance_matrix[thread_id][i][j] > -0.5)
            {
              if(cluster_distance_matrix[i][j] < -0.5)
              {
                cluster_distance_matrix[i][j] = private_distance_matrix[thread_id][i][j];
              }
              else
              {
                // Since average neighbor is the mean value, and each element has already been
                // divided by the total number, we just need to sum up all the values found
                // by the threads to get the full distance.
                cluster_distance_matrix[i][j] += private_distance_matrix[thread_id][i][j];
              }
            }
          }
          else if(algo=='f')
          { 
            if(private_distance_matrix[thread_id][i][j] > -0.5)
            {
              // Max every element with this thread's distances
              // This is exactly like the nearest neighbor version, but used to find the maximum.
              cluster_distance_matrix[i][j] = (private_distance_matrix[thread_id][i][j] > cluster_distance_matrix[i][j]) ? (private_distance_matrix[thread_id][i][j]) : (cluster_distance_matrix[i][j]);
            }
          }
        }
        // Reset private distance matrix for next loop
        private_distance_matrix[thread_id][i][j] = -1;
      }
    } // End critical section

  } // End parallel
}
