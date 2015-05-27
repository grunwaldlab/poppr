#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
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
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

# Keeping this for compatability
new.read.genalex <- function(genalex, ploidy=2, geo=FALSE, region=FALSE){
  return(read.genalex(genalex, ploidy, geo, region))
}

##' Produce a table of diversity statistics
#' 
#' @param z a table of integers representing counts of MLGs (columns) per 
#' population (rows)
#' 
#' @return a numeric matrix with 4 columns giving the following statistics for 
#'   each population: \itemize{\item H - Shannon Diversity \item G Stoddart and
#'   Taylor's Diveristy (inverse Simpson's) \item unbiased Simpson Diversity \eqn{(N/(N-1))*(1 - D)} \item E.5 evenness}
#' @export
#' @author Zhian N. Kamvar
#' @examples
#' library(poppr)
#' data(Pinf)
#' tab <- mlg.table(Pinf, bar = FALSE)
#' get_stats(tab)
get_stats <- function(z, H = TRUE, G = TRUE, simp = TRUE, E5 = TRUE, ...){
  
  E.5 <- function(x){
    H <- vegan::diversity(x)
    G <- vegan::diversity(x, "inv")
    (G - 1)/(exp(H) - 1)
  }
  BASIC_FUNS <- list(H = vegan::diversity,
                     G = function(x) vegan::diversity(x, "inv"),
                     simp = function(x) vegan::diversity(x, "simp"),
                     E.5 = E.5)
  FUNS   <- c(BASIC_FUNS, list(...))
  STATS <- names(FUNS)
  STATS <- STATS[c(H, G, simp, E5, rep(TRUE, length(list(...))))]
  
  boot <- is.null(dim(z))
  
  nrows <- ifelse(boot, 1, nrow(z))
  if (boot){
    dims <- list(NULL, Index = STATS)
  } else {
    dims <- list(Pop = rownames(z), Index = STATS)
  }
  
  mat <- matrix(nrow = nrows, ncol = length(STATS), dimnames = dims)
  for (i in STATS){
    mat[, i] <- FUNS[[i]](z)
  }
  return(drop(mat))
}

#' @importFrom vegan diversity
boot_stats <- function(x, i, H = TRUE, G = TRUE, simp = TRUE, E5 = TRUE, ...){
  res <- get_stats(tabulate(x[i]), H, G, simp, E5, ...)
  return(res)
}

extract_samples <- function(x) rep(1:length(x), x)

#' Perform a bootstrap analysis on diversity statistics
#' 
#' @param tab a table produced from the \pkg{poppr} function \code{\link[poppr]{mlg.table}}. MLGs in columns and populations in rows
#' @param n an integer > 0 specifying the number of bootstrap replicates to perform (corresponds to \code{R} in the function \code{\link[boot]{boot}}.
#' @param ... other parameters passed on to \code{\link[boot]{boot}}.
#' 
#' @return a list of objects of class "boot". 
#' @seealso \code{\link{boot_ci}}
#' @export
#' @author Zhian N. Kamvar
#' @examples
#' library(poppr)
#' data(Pinf)
#' tab <- mlg.table(Pinf, bar = FALSE)
#' do_boot(tab, 10L)
#' \dontrun{
#' # This can be done in a parallel fasion (OSX uses "multicore", Windows uses "snow")
#' system.time(do_boot(tab, 10000L, parallel = "multicore", ncpus = 4L))
#' system.time(do_boot(tab, 10000L))
#' }
#' @importFrom boot boot
do_boot <- function(tab, n, H = TRUE, G = TRUE, simp = TRUE, E5 = TRUE, ...){
  res <- apply(tab, 1, function(x){
    boot::boot(extract_samples(x), boot_stats, n, H = H, G = G, simp = simp, E5 = E5, ...)
  })
  return(res)
}

get_ci <- function(x, lb, ub){
  res <- apply(x$t, 2, quantile, c(lb, ub), na.rm = TRUE)
  return(res)
}

get_all_ci <- function(res, ci = 95){
  lower_bound  <- (100 - ci)/200
  upper_bound  <- 1 - lower_bound
  funval       <- matrix(numeric(8), nrow = 2)
  CI           <- vapply(res, FUN = get_ci, FUN.VALUE = funval, 
                         lower_bound, upper_bound)
  dCI          <- dimnames(CI)
  dimnames(CI) <- list(CI    = dCI[[1]], 
                       Index = c("H", "G", "Hexp", "E.5"),
                       Pop   = dCI[[3]])
  
  return(CI)
}

#' Perform bootstrap statistics, calculate and plot confidence intervals.
#' 
#' @param tab a genind object OR a matrix produced from \code{\link[poppr]{mlg.table}}.
#' @param n an integer defining the number of bootstrap replicates (defaults to 1000).
#' @param ci the percent for confidence interval.
#' @param total argument to be passed on to \code{\link[poppr]{mlg.table}} if \code{tab} is a genind object.
#' @param ... parameters to be passed on to \code{\link[boot]{boot}}
#' 
#' @return an array of 3 dimensions giving the lower and upper bound, the index
#' measured, and the population. This also prints barplots for each population
#' faceted by each index using \pkg{ggplot2}. This plot can be retrieved by using
#' \code{p <- last_plot()}
#' 
#' @export
#' @author Zhian N. Kamvar
#' @examples
#' library(poppr)
#' data(Pinf)
#' boot_ci(Pinf, n = 100)
#' \dontrun{
#' # This can be done in a parallel fasion (OSX uses "multicore", Windows uses "snow")
#' system.time(boot_ci(tab, 10000L, parallel = "multicore", ncpus = 4L))
#' system.time(boot_ci(tab, 10000L))
#' }
#' 
boot_ci <- function(tab, n = 1000, ci = 95, total = TRUE, ...){
  if (!is.matrix(tab) & is.genind(tab)){
    tab <- mlg.table(tab, total = total, bar = FALSE)
  }
  res  <- do_boot(tab, n, ...)
  orig <- get_stats(tab)
  orig <- melt(orig)
  orig$Pop <- factor(orig$Pop)
  CI   <- get_all_ci(res, ci = ci)
  samp <- vapply(res, "[[", FUN.VALUE = res[[1]]$t, "t")
  dimnames(samp) <- list(NULL, 
                         Index = c("H", "G", "Hexp", "E.5"),
                         Pop = rownames(tab))
  sampmelt <- melt(samp)
  sampmelt$Pop <- factor(sampmelt$Pop)
  pl <- ggplot(sampmelt, aes_string(x = "Pop", y = "value", group = "Pop")) + 
    geom_boxplot() + 
    geom_point(aes_string(color = "Pop", x = "Pop", y = "value"), 
               size = 5, pch = 16, data = orig) +
    xlab("Population") + labs(color = "Observed") +
    facet_wrap(~Index, scales = "free_y") + myTheme
  print(pl)
  return(CI)
}

myTheme <- theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

read_allele_columns <- function(x, ploidy = 2){
  clm <- ncol(x)
  if (length(ploidy) == 1){
    if (clm %% ploidy == 0){
      ploidy <- rep(ploidy, clm/ploidy)
    } else {
      stop("ploidy specified does not match number of columns.")
    }
  }
  loci  <- first_allele_col(ploidy)
  x2    <- x[, loci]
  lapply(1:length(loci), function(a) x2[, a] <<-
           pop_combiner(x, hier = loci[a]:(loci[a]+ploidy[a]-1), sep = "/"))
  
  x2 <- data.frame(lapply(x2, add_zeroes,  ploidy = max(ploidy)))
  x2[apply(x2, 2, function(na) grepl("NA", na))] <- NA
  return(x2)
}

first_allele_col <- function(ploidies) cumsum(ploidies) - ploidies + 1

# Function that will add extra zeroes onto genotypes that are deficient in the 
# right number of alleles.
add_zeroes <- function(x, ploidy = 3){
  extras <- ploidy - vapply(strsplit(x, "/"), length, 69)
  vapply(1:length(extras), zero_adder, character(1), extras, x)
}

zero_adder <- function(index, extras, df){
  if(extras[index] > 0){
    return(paste(c(rep(0, extras[index]), df[index]), collapse = "/"))
  } else {
    return(df[index])
  }
}
#==============================================================================#
# Create a population hierarchy in a genind object if one is not there.
#==============================================================================#
genind_hierarchy <- function(x, df = NULL, dfname = "population_hierarchy"){
  if (is.null(df)){
    df <- data.frame(Pop = pop(x))
  } 
  other(x)[[dfname]] <- df
  return(x)
}

pair_ia <- function(pop){

  if(pop@type == "codom"){
    pop_loci <- seploc(pop)
    loci_pairs <- combn(pop@loc.names, 2)
    pair_ia_vector <- apply(loci_pairs, 2, function(x) .Ia.Rd(pop_loci[x]))
    colnames(pair_ia_vector) <- apply(loci_pairs, 2, paste, collapse = ":")
  } else {
    loci_pairs <- combn(1:nLoc(pop), 2)
    pair_ia_vector <- apply(loci_pairs, 2, function(x) .PA.Ia.Rd(pop[, x], missing = "ignore"))
    colnames(pair_ia_vector) <- apply(combn(pop@loc.names, 2), 2, paste, collapse = ":")
  }
  rownames(pair_ia_vector) <- c("Ia", "rbarD")    
  return(pair_ia_vector)
}


poppr_pair_ia <- function(pop){
  if(is.null(pop(pop))){
    return(pair_ia(pop))
  }
  pops       <- seppop(pop, drop = FALSE)
  loci_pairs <- choose(nLoc(pop), 2)
  res_mat    <- matrix(0.5, 2, loci_pairs)
  pops_array <- vapply(pops, pair_ia, res_mat)
  return(pops_array)
}

testing_funk <- function(){
  cat("This test worked...maybe.\n")
}

getinds <- function(x){
  # x is a vector indicating the size of loci in a data set (must be continuous).
  # The sum of x is equal to the number of observations.
  indices <- matrix(ncol = 2, nrow = length(x))
  to   <- cumsum(x)
  from <- c(1, to[-length(to)]+1)
  indices[] <- as.integer(c(from, to))
  return(indices)
}

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
bootjack <- function(gid, iterations = 999, half = NULL, progbar = NULL){
  if (!is.list(gid)){
    N <- nInd(gid)
    numLoci <- nLoc(gid)
  } else {
    N <- nInd(gid[[1]])
    numLoci <- length(gid)
  }
  # Step 1: make a distance matrix defining the indices for the pairwise 
  # distance matrix of loci.
  np  <- choose(N, 2)
  dis <- 1:np
  dis <- make_attributes(dis, N, 1:N, "dist", call("dist"))
  mat <- as.matrix(dis)
  # Step 2: calculate the pairwise distances for each locus. 
  V   <- jack.pair_diffs(gid, numLoci, np)
  if (!is.null(half)){
    np <- choose(half, 2)
    sample.data <- vapply(1:iterations, jackrun, numeric(2), V, mat, N, half, 
                          np, progbar, iterations)
  } else {
    sample.data <- vapply(1:iterations, bootrun, numeric(2), V, mat, N, np, 
                          progbar, iterations)
  }
  rownames(sample.data) <- c("Ia", "rbarD")
  return(data.frame(t(sample.data)))
}

#------------------------------------------------------------------------------#
# Variables for bootrun and jackrun
#
# count      = the current iteration
# V          = the matrix of pairwise distances per locus
# mat        = square matrix containing indices for V
# N          = number of observations in data set
# half       = the number of observations to be sampled without replacement
# np         = choose(N, 2)
# progbar    = a progress bar object
# iterations = the number of total samples
#------------------------------------------------------------------------------#
bootrun <- function(count, V, mat, N, np, progbar, iterations){
  if (!is.null(progbar)){
    setTxtProgressBar(progbar, count/iterations)
  }
  # Step 3: sample individuals with replacement and subset distance matrix
  # from step 1.
  inds    <- sample(N, replace = TRUE)
  newInds <- as.vector(as.dist(mat[inds, inds])) 
  
  # Step 4: Find the number of zeroes in from the distance matrix. This will be 
  # the number of duplicated genotypes.
  zeroes <- length(inds == 0)

  # Step 5: subset the pairwise locus matrix.
  newV <- V[newInds, ]
  V    <- list(d.vector  = colSums(newV), 
               d2.vector = colSums(newV * newV), 
               D.vector  = c(rowSums(newV), rep(0, zeroes))
              )
  # Step 6: Calculate the index of association.
  return(jack.calc(V, np))
}

jackrun <- function(count, V, mat, N, half, np, progbar, iterations){
  if (!is.null(progbar)){
    setTxtProgressBar(progbar, count/iterations)
  }

  inds    <- sample(N, half)
  newInds <- as.vector(as.dist(mat[inds, inds])) 

  newV <- V[newInds, ]
  V    <- list(d.vector  = colSums(newV), 
               d2.vector = colSums(newV * newV), 
               D.vector  = rowSums(newV)
              )
  return(jack.calc(V, np))
}

jack.calc <- function(V, np){
  varD <- ((sum(V$D.vector^2) - ((sum(V$D.vector))^2)/np))/np
  vard.vector <- ((V$d2.vector - ((V$d.vector^2)/np))/np)
  vardpair.vector <- .Call("pairwise_covar", vard.vector, PACKAGE = "poppr")
  sigVarj <- sum(vard.vector)
  rm(vard.vector)
  Ia <- (varD/sigVarj) - 1
  rbarD <- (varD - sigVarj)/(2 * sum(vardpair.vector))
  return(c(Ia, rbarD))
}

jack.pair_diffs <- function(pop, numLoci, np){
  temp.d.vector <- matrix(nrow = np, ncol = numLoci, data = as.numeric(NA))
  if(!is.list(pop)){
    ploid <- 2
    temp.d.vector <- vapply(seq(numLoci), 
                      function(x) as.vector(dist(pop@tab[,x])), 
                      temp.d.vector[,1])
    # checking for missing data and imputing the comparison to zero.
    if(any(is.na(temp.d.vector))){
      temp.d.vector[which(is.na(temp.d.vector))] <- 0
    }
  } else {
    temp.d.vector <- vapply(pop, function(x) .Call("pairdiffs", tab(x), PACKAGE = "poppr")/2, 
                            temp.d.vector[, 1])
    temp.d.vector <- ceiling(temp.d.vector)
  }
  return(temp.d.vector)
}


jackbootplot <- function(df, obs.df, index = c("rd", "ia", "both"), obsline = TRUE){
  ARGS                <- c("rd", "ia", "both")
  index               <- match.arg(index, ARGS)
  Indexfac            <- c("I[A]","bar(r)[d]")
  levels(df$variable) <- Indexfac
  obs.df$variable     <- Indexfac

  if (index == "rd"){
    df     <- df[df$variable == "bar(r)[d]", ]
    obs.df <- obs.df[obs.df$variable == "bar(r)[d]", ]
  } else if (index == "ia"){
    df     <- df[df$variable == "I[A]", ]
    obs.df <- obs.df[obs.df$variable == "I[A]", ]
  }

  distplot <- ggplot(df, aes_string(y = "value", x = "Distribution")) +
                geom_violin(aes_string(fill = "Distribution")) + 
                geom_boxplot(alpha = 0.25, width = 0.125)
  if (index == "both"){
    distplot <- distplot + facet_grid("variable ~ .", scales = "free", 
                                      labeller = label_parsed) 
  } 
  if (obsline){
    aesth    <- aes_string(yintercept = "value", group = "variable") 
    distplot <- distplot + geom_hline(aesth, data = obs.df, linetype = 2)
  }
  return(distplot)      
}

#' @importFrom reshape2 melt

jackbootcomp <- function(pop, sample = 999, quiet = FALSE, method = 1, 
                         jack = 0.63, plotindex = "rd"){
  ARGS <- c("rd", "ia", "both")
  plotindex <- match.arg(plotindex, ARGS)
  
  inds <- nInd(pop)
  message("Null Distribution...")

  null <- ia(pop, sample = sample, valuereturn = TRUE, quiet = quiet, 
             hist = FALSE, method = method)
  if(pop@type == "codom"){
    popx <- seploc(pop)
  } else {
    popx <- pop
  }

  message("Alternative Distribution (Bootstrap)...")
  if(!quiet){
    bootbar <- txtProgressBar(style = 3)
  } else {
    bootbar <- NULL
  }
  jackbar <- bootbar

  altboot <- bootjack(popx, sample, half = NULL, bootbar)

  message("Alternative Distribution (Jack Knife)...")
  half    <- round(inds*jack)
  altjack <- bootjack(popx, sample, half, jackbar)  
  datlist <- list(null = null$sample, bootstrap = altboot, `jack knife` = altjack)
  dat           <- melt(datlist, measure.vars = c("Ia", "rbarD"))
  names(dat)[3] <- "Distribution"
  message("Creating Plots")
  obs.df   <- data.frame(list(variable = names(null$samples), 
                              value = null$index[c(1,3)]))
  distplot <- jackbootplot(dat, obs.df, plotindex) + 
              ggtitle(paste("Data:", deparse(substitute(pop)), "\n", inds, 
                            "Individuals,", half, "Sampled for Jack knife")
                     )
  if (plotindex == "both"){
    distplot <- distplot + theme(strip.text = element_text(size = rel(2)))    
  } else if (plotindex == "rd"){
    distplot <- distplot + ylab(expression(bar(r)[d]))
  } else {
    distplot <- distplot + ylab(expression(I[A]))
  }
  return(list(observed = null$index, null_samples = null$samples, 
              alt_samples_boot = altboot, alt_samples_jack = altjack, 
              plot = distplot)
        )
} 



################################################################################
#################### Zhian's Functions above ###################################
################################################################################
################################################################################
################################################################################
#################### DO NOT CROSS THIS LINE!!!! ################################
################################################################################
################################################################################
################################################################################
#################### Javier's Functions below ##################################
################################################################################



javier<-function(x){
  cat ("http://www.youtube.com/watch?v=1-ctsxVXvO0")
}
