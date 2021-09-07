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
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <Rdefines.h>
#include <R.h>


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Finds any edges that could have been included in the given MST, then adds them
to a list of edges to be added in the future.

Input: A minimum spanning tree (as a matrix of edge weights) and the distance 
       matrix used to construct it. epsi is used as the epsilon value in determining
       how similar two branches must be to be considered equal. Note that this value
       should never be less than an reasonable integer multiple of the machine epsilon,
Output: A list of edges that lost the tiebreak while constructing the mst.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
SEXP msn_tied_edges(SEXP mst, SEXP bclone, SEXP epsi)
{
  // R code
  /*
  # Add any relevant edges that were cut from the mst while still being tied for the title of optimal edge
  # Loop through every vertex
  for(v1 in dimnames(mst[])[[1]])
  {
    # Store the minimum path out of this vertex, as represented in the mst
    mn <- min(mst[v1,][mst[v1,] > 0])
    # Loop through every possible path from this, as stored in the distance matrix
    for(v2 in names(which(as.matrix(bclone)[v1,] > 0)))
    {
      # If this path has the same weight (+ or - epsilon) as the minimum, add it to the tree unless it's already there
      if(isTRUE(all.equal(as.matrix(bclone)[v1,v2],mn,tolerance=.0005)) & !(mst[v1,v2] > 0))
      {
        mst <- mst + edge(c(v1,v2), weight=mn)
      }
    }
  }
  */
  // End R code

  SEXP R_out;
  int edges_size = 24;
  int num_edges = 0;
  int num_vertices = 0;
  double mn;
  double candidate_edge;
  double this_edge;
  double *edges = R_Calloc(edges_size,double);
  num_vertices = INTEGER(getAttrib(bclone,R_DimSymbol))[1];
  for(int i = 0; i < num_vertices; i++)
  {
    R_CheckUserInterrupt();
    // Find the shortest path out of this vertex
    mn = -1;
    for(int j = 0; j < num_vertices; j++)
    {
      this_edge = REAL(mst)[i + j*num_vertices];
      // Set mn to REAL(mst)[i][j] if it is smaller than mn, or it mn is still -1
      // ie, this loop finds the shortest path out of node i
      mn = ((mn < 0 || this_edge < mn) && this_edge > 0) ? this_edge : mn;
    }
    // Find all paths out of this vertex that are tied in length with the minimum
    for(int j = i+1; j < num_vertices; j++)
    {
      candidate_edge = REAL(bclone)[i + j*num_vertices];
      this_edge = REAL(mst)[i + j*num_vertices];
      // Check for matching edges that do not already exist in this graph
      // For each path out of node i, check to see if its length is no more than
      // epsi larger than the minimum
      if(fabs(candidate_edge - mn) < asReal(epsi) && !(this_edge > 0))
      {
        // Resize the edges array if it does not have room for another set of 
        // three values
        if(num_edges+2 >= edges_size)
        {
          edges = R_Realloc(edges, edges_size*2, double);
          edges_size *= 2;
        }
        // add i, j, and mn to a vector to return
        // This sets edges[x] and edges[x+1] to be end points of a path of size 
        // edges[x+2] = mn
        edges[num_edges] = (double)(i+1);
        edges[num_edges+1] = (double)(j+1);
        edges[num_edges+2] = candidate_edge;
        num_edges += 3;
      }
    }
  } 
  R_out = PROTECT(allocVector(REALSXP,num_edges));
  for(int i = 0; i < num_edges; i+=3)
  {
    REAL(R_out)[i] = edges[i];
    REAL(R_out)[i+1] = edges[i+1];
    REAL(R_out)[i+2] = edges[i+2];
  }   
  R_Free(edges);
  UNPROTECT(1);
  return R_out;
}

