#==============================================================================#
# bootjack is a function that will calculate values of the index of association
# after sampling with or without replacement. The purpose of this function is
# to help determine distributions of I_A under non-random mating by creating 
# exact copies of samples by sampling with replacement and also to derive a
# distribution of the data by sampling a subset of the data (63% by default)
# without replacement. 
#
# Since the data itself is not being changed, we can use the distances observed.
# These distances are calculated per-locus and presented in a matrix (V) with the
# number of columns equal to the number of loci and the number of rows equal to
# choose(N, 2) where N is the number of samples in the data set. Sampling has
# been optimized with this strategy by first creating a lower triangle distance
# matrix containing indices from 1 to choose(N, 2). This is then converted to a
# square matrix (mat) and subset after sampling with or without replacement. The
# remaining indices in the distance matrix is used to subset the distance-per-
# locus matrix (V). Calculations are then performed on this matrix. It must be
# noted that for sampling with replacement, duplicated indices must be
# supplemented with rows of zeroes to indicate no distance.  
#==============================================================================#
#' @rdname ia
#' @param n an integer specifying the number of samples to be drawn. Defaults to
#'   \code{NULL}, which then uses the number of multilocus genotypes.
#' @param reps an integer specifying the number of replicates to perform. 
#'   Defaults to 999.
#' @param use_psex a logical. If \code{TRUE}, the samples will be weighted by the value 
#'   of psex. Defaults to \code{FALSE}.
#' @param ... arguments passed on to \code{\link{psex}}
#'
#' @return \subsection{resample.ia()}{a data frame with the index of association and standardized index of
#' association in columns. Number of rows represents the number of reps.}
#' @export
resample.ia <- function(gid, n = NULL, reps = 999, quiet = FALSE, use_psex = FALSE, ...){
  
  quiet   <- should_poppr_be_quiet(quiet)
  weights <- if (use_psex) psex(gid, ...) else NULL
  N       <- nInd(gid)
  numLoci <- nLoc(gid)
  if (is.null(n)){
    n <- suppressWarnings(nmll(gid))
  }
  if (n >= N){
    stop("n must be fewer than the number of observations", call. = FALSE)
  }
  if (n < 3){
    stop("n must be greater than 2 observations", call. = FALSE)
  }
  if (gid@type == "codom"){
    gid <- seploc(gid)
  }
  
  # Step 1: make a distance matrix defining the indices for the pairwise 
  # distance matrix of loci.
  np  <- choose(N, 2)
  dis <- seq.int(np)
  dis <- make_attributes(dis, N, seq(N), "dist", call("dist"))
  mat <- as.matrix(dis)
  # Step 2: calculate the pairwise distances for each locus. 
  V   <- pair_matrix(gid, numLoci, np)
  np  <- choose(n, 2)
  
  sample.data        <- matrix(numeric(reps*2), ncol = 2, nrow = reps)
  if(!quiet) progbar <- dplyr::progress_estimated(reps)
  for (i in seq(reps)){
    sample.data[i, ] <- run.jack(V, mat, N, n, np, replace = FALSE, weights = weights)
    if (!quiet) print(progbar$tick())
  }
  if(!quiet){
    print(progbar$stop())
    cat("\n")
  }
  colnames(sample.data) <- c("Ia", "rbarD")
  return(data.frame(sample.data))
}
#' Bootstrap the index of association
#' 
#' This function will perform the index of association on a bootstrapped data
#' set multiple times to create a distribution, showing the variation of the
#' index due to repeat observations.
#'
#' @param gid a genind or genclone object
#' @param how method of bootstrap. The default `how = "partial"` will include
#'   all the unique genotypes and sample with replacement from the unique
#'   genotypes until the total number of individuals has been reached. Using
#'   `how = "full"` will randomly sample with replacement from the data as it
#'   is. Using `how = "psex"` will sample from the full data set after first
#'   weighting the samples via the probability of encountering the nth occurence
#'   of a particular multilocus genotype. See [psex()] for details.
#' @param reps an integer specifying the number of replicates to perform. 
#'   Defaults to 999.
#' @param quiet a logical. If `FALSE`, a progress bar will be displayed. If
#'   `TRUE`, the progress bar is suppressed.
#' @param ... options passed on to [psex()]
#'   
#' @return a data frame with the index of association and standardized index of 
#'   association in columns. Number of rows represents the number of reps.
#' @note This function is experimental. Please do not use this unless you know
#'   what you are doing.
#' @export
#' @md
#' @seealso
#'   [ia()], 
#'   [pair.ia()],
#'   [psex()]
#'
#' @examples
#' data(Pinf)
#' boot.ia(Pinf, reps = 99)
boot.ia <- function(gid, how = "partial", reps = 999, quiet = FALSE, ...){
  
  METHOD  <- match.arg(how, c("partial", "full", "psex"))
  weights <- if (METHOD == "psex") psex(gid, ...) else NULL
  quiet   <- should_poppr_be_quiet(quiet)
  N       <- nInd(gid)
  numLoci <- nLoc(gid)
  if (METHOD == "partial") {
    n   <- suppressWarnings(nmll(gid))
    gid <- clonecorrect(gid, NA)
  } else {
    n <- N
  }
  if (n < 3) {
    stop("n must be greater than 2 observations", call. = FALSE)
  }
  if (gid@type == "codom") {
    gid <- seploc(gid)
  }
  
  # Step 1: make a distance matrix defining the indices for the pairwise 
  # distance matrix of loci.
  np  <- choose(n, 2)
  dis <- seq.int(np)
  dis <- make_attributes(dis, n, seq(n), "dist", call("dist"))
  mat <- as.matrix(dis)
  # Step 2: calculate the pairwise distances for each locus. 
  V   <- pair_matrix(gid, numLoci, np)
  np  <- choose(N, 2)
  
  sample.data         <- matrix(numeric(reps*2), ncol = 2, nrow = reps)
  if (!quiet) progbar <- dplyr::progress_estimated(reps)
  for (i in seq(reps)) {
    sample.data[i, ] <- run.jack(V, mat, N, n, np, replace = TRUE, 
                                 method = METHOD, weights = weights)
    if (!quiet) progbar$tick()$print()
  }
  if (!quiet) {
    cat("\n")
    progbar$stop()
  }
  colnames(sample.data) <- c("Ia", "rbarD")
  return(data.frame(sample.data))
}


bias.ia <- function(theta_hat, theta_star){
  vapply(theta_star, mean, numeric(1), na.rm = TRUE) - theta_hat
}

#' @rdname ia
#' @export
jack.ia <- function(gid, n = NULL, reps = 999, quiet = FALSE){
  msg <- paste0("jack.ia() is deprecated and will be removed in future versions.\n",
               "Please use resample.ia() instead.\n",
               "I am returning the results of resample.ia()")
  .Deprecated("resample.ia", msg = msg, package = "poppr")
  resample.ia(gid, n = n, reps = reps, quiet = quiet)
}

#' Helper function to calculate the index of association on reduced data.
#' 
#' Since the data itself is not being changed, we can use the distances
#' observed. These distances are calculated per-locus and presented in a matrix
#' (V) with the number of columns equal to the number of loci and the number of
#' rows equal to choose(N, 2) where N is the number of samples in the data set.
#' Sampling has been optimized with this strategy by first creating a lower
#' triangle distance matrix containing indices from 1 to choose(N, 2). This is
#' then converted to a square matrix (mat) and subset after sampling with or
#' without replacement. The remaining indices in the distance matrix is used to
#' subset the distance-per-locus matrix (V). Calculations are then performed on
#' this matrix. It must be noted that for sampling with replacement, duplicated
#' indices must be supplemented with rows of zeroes to indicate no distance. 
#' This implementation does not sample with replacement.
#'
#' @param V a matrix of distances for each locus in columns and observations in
#'   rows
#' @param mat a square matrix with the lower triangle filled with indices of the
#'   rows
#' @param N The number of observations in the original data
#' @param n The number of observations to sample
#' @param np the number of rows in V
#' @param replace logical whether or not to sample with replacement. Defaults to
#'   FALSE
#' @param method passed from boot.ia
#' @param weights a vector of weights given by [psex()]
#'
#' @return Estimates of the index of association
#' @noRd
#'
#' @examples
#' # No examples here
run.jack <- function(V, mat, N, n, np, replace = FALSE, method = "partial", weights = NULL){

  if (replace & method == "partial") {
    # For the partial boot method. In this case, the incoming data is clone
    # censored, so N is the desired number of individuals and n is the observed
    # number of individuals. Since we want to keep the number of MLG steady, we
    # are only resampling N - n individuals here.
    inds <- c(seq.int(n), sample(n, N - n, replace = TRUE))
  } else {
    inds <- sample(N, n, replace = replace, prob = weights)
  }
  newmat  <- mat[inds, inds]
  newInds <- newmat[lower.tri(newmat)]

  newV <- V[newInds, ]
  V    <- list(d.vector  = colSums(newV), 
               d2.vector = colSums(newV * newV), 
               D.vector  = rowSums(newV)
              )
  return(ia_from_d_and_D(V, np))
}

