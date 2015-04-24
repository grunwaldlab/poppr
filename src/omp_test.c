#include <Rinternals.h>
#include <Rdefines.h>
#include <R.h>

// Simple function to check for openMP support
// Returns 1 if openMP is supported, 0 if it is not.
SEXP omp_test()
{
  SEXP Rout;
  PROTECT(Rout = allocVector(INTSXP,1));
  #ifdef _OPENMP
  INTEGER(Rout)[0] = 1;
  #else
  INTEGER(Rout)[0] = 0;
  #endif
  UNPROTECT(1);
  return(Rout);
}
