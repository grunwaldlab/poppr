#include <Rinternals.h>
#include <Rdefines.h>
#include <R.h>

/*
#' Calculate correction for genetic distances
#'
#' @param nas a list of missing positions per sample
#' @param nloc the number of loci
#' @param mat a logical specifying whether or not a matrix should be returned
#'   (default: TRUE)
#'
#' @return an n x n matrix or a choose(n, 2) length vector of values that scale
#'   from 1 to the number of loci.
#' @noRd
 missing_correction <- function(nas, nloc, mat = TRUE){
  n          <- length(nas)
  correction <- matrix(0, nrow = n, ncol = n)
  for (i in seq(n - 1)) {
    for (j in seq(from = i + 1, to = n)) {
      res <- length(unique(unlist(nas[c(i, j)], use.names = FALSE)))
      correction[j, i] <- res -> correction[i, j]
    }
  }
  res <- nloc/(nloc - correction)
  if (mat) {
    return(res)
  } else {
    return(res[lower.tri(res)])
  }
}
 
*/

int count_unique(SEXP arr1, SEXP arr2);
SEXP adjust_missing(SEXP nas, SEXP nloc);
/*
 * Count all unique elements for the union of two arrays
 * 
 * Parameters:
 *  arr1 an integer array of any length
 *  arr2 an integer array of any length
 *
 * Return:
 *  an integer that is between the length of the smallest array and the length
 *  of the largest array. 
 */
int count_unique(SEXP arr1, SEXP arr2)
{
  int i = 0;
  int j = 0;
  int duplicates = 0;
  int len1 = length(arr1);
  int len2 = length(arr2);
  while (i < len1 && j < len2)
  {
    if (INTEGER(arr1)[i] < INTEGER(arr2)[j])
    {
      i++;
    }
    else if (INTEGER(arr1)[i] > INTEGER(arr2)[j])
    {
      j++;
    }
    else
    {
      i++;
      j++;
      duplicates++;
    }
  }
  return len1 + len2 - duplicates;
}
/*
 * Calculate adjustment for missing data in pairwise comparisons. This will
 * return a square matrix that is used to multiply the raw differences of a
 * distance matrix in order to scale the differences by the number of observed
 * loci. 
 * 
 * Parameters:
 *  nas a list where each element represents a sample containing an integer 
 *      vector representing positions of missing data for that individual
 *  nloc an integer specifying the number of loci observed in the entire set
 * 
 * Return:
 *  a square matrix
 */
SEXP adjust_missing(SEXP nas, SEXP nloc)
{
  int i;
  int j;
  int NLOC = asInteger(nloc);
  SEXP nai;
  SEXP naj;
  double u;
  int n    = length(nas);
  SEXP out = PROTECT(allocMatrix(REALSXP, n, n));
  for (i = 0; i < n - 1; i++)
  {
    // set diag to one
    REAL(out)[i + i*n] = 1.0;
    // GET NA list for i
    nai = VECTOR_ELT(nas, i);
    for (j = i + 1; j < n; j++)
    {
      // Get NA list for j
      naj = VECTOR_ELT(nas, j);
      // Scale by N/(N - M)
      u   = (double)NLOC/(double)(NLOC - count_unique(nai, naj));
      
      REAL(out)[i + j*n] = u;
      REAL(out)[i*n + j] = u;
    }
  }
  // fencepost for identity
  REAL(out)[(n*n) - 1] = 1.0;
  UNPROTECT(1);
  return(out);
}
