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
#==============================================================================#
#' Calculate a dendrogram with bootstrap support using any distance applicable 
#' to genind or genclone objects.
#' 
#' @param x a \linkS4class{genind}, \linkS4class{genpop},
#'   \linkS4class{genclone}, \linkS4class{genlight}, \linkS4class{snpclone}, or
#'   \link{matrix} object.
#'   
#' @param tree a text string or function that can calculate a tree from a 
#'   distance matrix. Defaults to "upgma". Note that you must load the package 
#'   with the function for it to work.
#'   
#' @param distance a character or function defining the distance to be applied 
#'   to x. Defaults to \code{\link{nei.dist}}.
#'   
#' @param sample An integer representing the number of bootstrap replicates 
#'   Default is 100.
#'   
#' @param cutoff An integer from 0 to 100 setting the cutoff value to return the
#'   bootstrap values on the nodes. Default is 0.
#'   
#' @param showtree If \code{TRUE} (Default), a dendrogram will be plotted. If 
#'   \code{FALSE}, nothing will be plotted.
#'   
#' @param missing any method to be used by \code{\link{missingno}}: "mean" 
#'   (default), "zero", "loci", "genotype", or "ignore".
#'   
#' @param mcutoff a value between 0 (default) and 1 defining the percentage of 
#'   tolerable missing data if the \code{missing} parameter is set to "loci" or 
#'   "genotype". This should only be set if the distance metric can handle 
#'   missing data.
#'   
#' @param quiet if \code{FALSE} (default), a progress bar will be printed to 
#'   screen.
#'   
#' @param root is the tree rooted? This is a parameter passed off to 
#'   \code{\link[ape]{boot.phylo}}. If the \code{tree} parameter returns a 
#'   rooted tree (like UPGMA), this should be \code{TRUE}, otherwise (like 
#'   neighbor-joining), it should be false. When set to \code{NULL} (default), 
#'   the tree will be considered unrooted unless it is a upgma tree.
#'   
#' @param ... any parameters to be passed off to the distance method.
#'   
#' @return an object of class \code{\link[ape]{phylo}}.
#'   
#' @details This function utilizes an internal class called 
#'   \code{\linkS4class{bootgen}} that allows bootstrapping of objects that 
#'   inherit the genind class. This is necessary due to the fact that columns in
#'   the genind matrix are defined as alleles and are thus interrelated. This 
#'   function will specifically bootstrap loci so that results are biologically 
#'   relevant. With this function, the user can also define a custom distance to
#'   be performed on the genind or genclone object.
#'   
#' @note \code{\link{provesti.dist}} and \code{\link{diss.dist}} are exactly the
#'   same, but \code{\link{diss.dist}} scales better for large numbers of 
#'   individuals (n > 125) at the cost of required memory. \subsection{missing 
#'   data}{Missing data is not allowed by many of the distances. Thus, one of 
#'   the first steps of this function is to treat missing data by setting it to 
#'   the average allele frequency in the data set. If you are using a distance 
#'   that can handle missing data (Provesti's distance), you can set 
#'   \code{missing = "ignore"} to allow the distance function to handle any 
#'   missing data. See \code{\link{missingno}} for details on missing 
#'   data.}\subsection{Bruvo's Distance}{While calculation of Bruvo's distance 
#'   is possible with this function, it is optimized in the function 
#'   \code{\link{bruvo.boot}}.}
#'   
#' @seealso \code{\link{nei.dist}} \code{\link{edwards.dist}} 
#'   \code{\link{rogers.dist}} \code{\link{reynolds.dist}} 
#'   \code{\link{provesti.dist}} \code{\link{diss.dist}} 
#'   \code{\link{bruvo.boot}} \code{\link[ape]{boot.phylo}} 
#'   \code{\link[adegenet]{dist.genpop}} \code{\link{dist}}
#'   
#' @export
#' @keywords bootstrap
#' @aliases bootstrap
#' @examples
#' 
#' data(nancycats)
#' nan9 <- popsub(nancycats, 9)
#' 
#' set.seed(9999)
#' # Generate a tree using nei's distance
#' neinan <- aboot(nan9, dist = nei.dist)
#' 
#' set.seed(9999)
#' # Generate a tree using custom distance
#' bindist <- function(x) dist(x$tab, method = "binary")
#' binnan <- aboot(nan9, dist = bindist)
#' 
#' \dontrun{
#' # AFLP data
#' data(Aeut)
#' 
#' # Nei's distance
#' anei <- aboot(Aeut, dist = nei.dist, sample = 1000, cutoff = 50)
#' 
#' # Rogers' distance
#' arog <- aboot(Aeut, dist = rogers.dist, sample = 1000, cutoff = 50)
#' 
#' # This can also be run on genpop objects
#' strata(Aeut) <- other(Aeut)$population_hierarchy[-1]
#' Aeut.gc <- as.genclone(Aeut)
#' setPop(Aeut.gc) <- ~Pop/Subpop
#' Aeut.pop <- genind2genpop(Aeut.gc)
#' set.seed(5000)
#' aboot(Aeut.pop, sample = 1000) # compare to Grunwald et al. 2006
#' 
#' # And genlight objects 
#' # From glSim:
#' ## 1,000 non structured SNPs, 100 structured SNPs
#' x <- glSim(100, 1e3, n.snp.struc=100, ploid=2)
#' aboot(x, distance = bitwise.dist)
#' 
#' # Utilizing other tree methods
#' 
#' library("ape")
#' 
#' aboot(Aeut.pop, tree = fastme.bal, sample = 1000)
#' 
#' # Utilizing options in other tree methods
#' 
#' myFastME <- function(x) fastme.bal(x, nni = TRUE, spr = FALSE, tbr = TRUE)
#' aboot(Aeut.pop, tree = myFastME, sample = 1000)
#' 
#' }
#==============================================================================#
aboot <- function(x, tree = "upgma", distance = "nei.dist", sample = 100,
                  cutoff = 0, showtree = TRUE, missing = "mean", mcutoff = 0,
                  quiet = FALSE, root = NULL, ...){
  if (!is(x, "genlight") && x@type == "PA"){
    xboot           <- x@tab
    colnames(xboot) <- locNames(x)
    if (is.genpop(x)){
      rownames(xboot) <- popNames(x)
    } else {
      rownames(xboot) <- indNames(x)
    }
  } else if (is(x, "gen")){
    if (is.genind(x)){
      if (missing %in% c("loci", "geno", "ignore")){
        x <- missingno(x, missing, quiet = quiet, cutoff = mcutoff)
        missing <- "asis"  
      } 
    }
    if (as.character(substitute(distance)) %in% "diss.dist"){
      xboot <- new("bootgen", x, na = missing, freq = FALSE)
      if (!is.integer(xboot@tab)){
        otab <- xboot@tab
        xboot@tab <- matrix(as.integer(otab), nrow = nrow(otab),
                            ncol = ncol(otab), dimnames = dimnames(otab))
      }
    } else {
      xboot <- new("bootgen", x, na = missing, freq = TRUE)
    }
  } else if (is(x, "genlight")){
    xboot <- x
    if (!as.character(substitute(distance)) %in% "bitwise.dist"){
      warning("distance from genlight objects can only be calculated by bitwise.dist.")
      distance <- bitwise.dist
    }
  } else {
    stop("x must be a genind, genpop, or genlight object.")
  }
  treefunk <- tree_generator(tree, distance, ...)
  xtree    <- treefunk(xboot)
  if (any(xtree$edge.len < 0)){
    xtree <- fix_negative_branch(xtree)
    warning(negative_branch_warning())
  }
  treechar <- paste(substitute(tree), collapse = "")
  if (is.null(root)){
    root <- grepl("upgma", treechar)
  }
  nodelabs <- boot.phylo(xtree, xboot, treefunk, B = sample, rooted = root, 
                         quiet = quiet)
  nodelabs <- (nodelabs/sample)*100
  nodelabs <- ifelse(nodelabs >= cutoff, nodelabs, NA)
  if (!is.genpop(x)){
    if (!is.null(indNames(x))){
      xtree$tip.label <- indNames(x)      
    }
  } else {
    xtree$tip.label <- popNames(x)
  }
  xtree$node.label <- nodelabs
  if (showtree){
    poppr.plot.phylo(xtree, treechar, root)
  }
  return(xtree)
}

#==============================================================================#
#' Produce a table of diversity statistics
#' 
#' @param z a table of integers representing counts of MLGs (columns) per 
#'   population (rows)
#'   
#' @param H logical whether or not to calculate Shannon's index
#' @param G logical whether or not to calculate Stoddart and Taylor's index
#' @param lambda logical whether or not to calculate Simpson's index
#' @param E5 logical whether or not to calculate Evenness
#' @param ... any functions that can be calculated on a vector or matrix of
#'   genotype counts.
#'   
#' @return a numeric matrix giving statistics (columns) for each population
#'   (rows).
#' @export
#' @seealso \code{\link{diversity_boot}} \code{\link{diversity_ci}} \code{\link{poppr}}
#' @author Zhian N. Kamvar
#' @examples
#' library(poppr)
#' data(Pinf)
#' tab <- mlg.table(Pinf, plot = FALSE)
#' diversity_stats(tab)
#' \dontrun{
#' # Example using the poweRlaw package to calculate the negative slope of the
#' # Pareto distribution.
#'
#' library("poweRlaw")
#' power_law_beta <- function(x){
#'   xpow <- displ(x[x > 0])                 # Generate the distribution
#'   xpow$setPars(estimate_pars(xpow))       # Estimate the parameters
#'   xdat <- plot(xpow, draw = FALSE)        # Extract the data
#'   xlm <- lm(log(y) ~ log(x), data = xdat) # Run log-log linear model for slope
#'   return(-coef(xlm)[2])
#' }
#' 
#' Beta <- function(x){
#'   x <- drop(as.matrix(x))
#'   if (length(dim(x)) > 1){
#'     res <- apply(x, 1, power_law_beta)
#'   } else {
#'     res <- power_law_beta(x)
#'   }
#'   return(res)
#' }
#'
#' diversity_stats(tab, B = Beta)
#' }
#==============================================================================#
diversity_stats <- function(z, H = TRUE, G = TRUE, lambda = TRUE, E5 = TRUE, ...){
  E.5 <- function(x){
    H <- vegan::diversity(x)
    G <- vegan::diversity(x, "inv")
    (G - 1)/(exp(H) - 1)
  }
  BASIC_FUNS <- list(H = vegan::diversity,
                     G = function(x) vegan::diversity(x, "inv"),
                     lambda = function(x) vegan::diversity(x, "simp"),
                     E.5 = E.5)
  FUNS   <- c(BASIC_FUNS, list(...))
  STATS <- names(FUNS)
  STATS <- STATS[c(H, G, lambda, E5, rep(TRUE, length(list(...))))]
  
  boot <- is.null(dim(z))
  
  nrows <- ifelse(boot, 1, nrow(z))
  if (boot){
    dims <- list(NULL, Index = STATS)
  } else {
    dims <- list(Pop = rownames(z), Index = STATS)
  }
  
  mat <- matrix(nrow = nrows, ncol = length(STATS), dimnames = dims)
  for (i in STATS){
    if (i == "E.5" && H && G){
      mat[, i] <- (mat[, "G"] - 1)/(exp(mat[, "H"]) - 1)
    } else {
      mat[, i] <- FUNS[[i]](z)
    }
  }
  return(drop(mat))
}

#==============================================================================#
#' Perform a bootstrap analysis on diversity statistics
#' 
#' @param tab a table produced from the \pkg{poppr} function 
#'   \code{\link[poppr]{mlg.table}}. MLGs in columns and populations in rows
#' @param n an integer > 0 specifying the number of bootstrap replicates to 
#'   perform (corresponds to \code{R} in the function \code{\link[boot]{boot}}.
#' @param n.rare a sample size at which all resamplings should be performed. 
#'   This should be no larger than the smallest sample size. Defaults to 
#'   \code{NULL}, indicating that each population will be sampled at its own 
#'   size.
#' @inheritParams diversity_stats
#' @param ... other parameters passed on to \code{\link[boot]{boot}} and 
#'   \code{\link{diversity_stats}}.
#'   
#' @return a list of objects of class "boot".
#' @seealso \code{\link{diversity_stats}} for basic statistic calculation, 
#'    \code{\link{diversity_ci}} for confidence intervals and plotting, and 
#'    \code{\link{poppr}}. For bootstrap sampling:
#'    \code{\link[stats]{rmultinom}} \code{\link[boot]{boot}}
#'    
#' @details Bootstrapping is performed in two ways depending on whether or not
#' rarefaction is being performed. 
#' \itemize{
#' \item if \code{n.rare = NULL}, bootstrapping is performed by sampling n
#' samples from a multinomial distribution weighted by the number of samples in
#' each MLG.
#' \item if \code{n.rare} is a number greater than zero, then bootstrapping is
#' performed by randomly sampling without replacement \emph{n.rare} samples from
#' the data.
#' }
#' @export
#' @author Zhian N. Kamvar
#' @examples
#' library(poppr)
#' data(Pinf)
#' tab <- mlg.table(Pinf, plot = FALSE)
#' diversity_boot(tab, 10L)
#' \dontrun{
#' # This can be done in a parallel fasion (OSX uses "multicore", Windows uses "snow")
#' system.time(diversity_boot(tab, 10000L, parallel = "multicore", ncpus = 4L))
#' system.time(diversity_boot(tab, 10000L))
#' }
#' @importFrom boot boot
#==============================================================================#
diversity_boot <- function(tab, n, n.rare = NULL, H = TRUE, G = TRUE, 
                           lambda = TRUE, E5 = TRUE, ...){
  if (!is.null(n.rare)){
    FUN <- rare_sim_boot
  } else {
    FUN <- multinom_boot
  }
  res <- apply(tab, 1, boot_per_pop, rg = FUN, n = n, 
               n.rare = n.rare, H = H, G = G, lambda = lambda, E5 = E5, ...)
  return(res)
}

boot_per_pop <- function(x, rg = multinom_boot, n, n.rare = NULL, H = TRUE, 
                         G = TRUE, lambda = TRUE, E5 = TRUE, ...){
  xi  <- extract_samples(x)
  res <- boot::boot(xi, boot_stats, R = n, sim = "parametric", ran.gen = rg, 
                    mle = n.rare, H = H, G = G, lambda = lambda, E5 = E5, ...)
  return(res)
}

multinom_boot <- function(x, mle = NULL){
  if (is.null(mle)){
    res <- rmultinom(1, length(x), prob = tabulate(x))
  } else {
    res <- rmultinom(1, mle, prob = tabulate(x))
  }
  extract_samples(res)
}

rare_sim_boot <- function(x, mle = 10){
  sample(x, mle)
}

boot_stats <- function(x, i, H = TRUE, G = TRUE, lambda = TRUE, E5 = TRUE, ...){
  xi  <- tabulate(x[i])
  res <- diversity_stats(xi, H, G, lambda, E5, ...)
  return(res)
}

extract_samples <- function(x) rep(1:length(x), x)



#' @importFrom boot boot.ci norm.ci
get_boot_ci <- function(index, x, type = "normal" , conf = 0.95, around_estimate = TRUE, ...){
  if (length(unique(x$t[, index])) == 1){
    return(c(NA_real_, NA_real_))
  } else if (around_estimate){
    res <- boot::norm.ci(x, conf = conf, index = index)[1, ]
  } else {
    res <- boot::boot.ci(x, conf, type, index, ...)[[type]][1, ]
  }
  return(tail(res, 2))
}

get_ci <- function(x, lb, ub, bci = TRUE, btype = "normal", center = TRUE){
  if (bci){
    lenout <- length(x$t0)
    conf   <- diff(c(lb, ub))
    res    <- vapply(1:lenout, get_boot_ci, numeric(2), x, btype, conf, center)
    if (!is.null(dim(res))) rownames(res) <- paste(c(lb, ub)*100, "%")
  } else {
    res <- apply(x$t, 2, quantile, c(lb, ub), na.rm = TRUE)
    if (all(apply(res, 2, function(i) i[1] == i[2]))){
      res[] <- NA_real_
    }
  }
  return(res)
}

get_all_ci <- function(res, ci = 95, index_names = c("H", "G", "Hexp", "E.5"),
                       center = TRUE, btype = "normal", bci = TRUE){
  lower_bound  <- (100 - ci)/200
  upper_bound  <- 1 - lower_bound
  n_indices    <- length(index_names)
  funval       <- matrix(numeric(n_indices*2), nrow = 2)
  CI           <- vapply(res, FUN = get_ci, FUN.VALUE = funval, 
                         lower_bound, upper_bound, bci = bci, btype = btype,
                         center = center)
  dCI          <- dimnames(CI)
  dimnames(CI) <- list(CI    = dCI[[1]], 
                       Index = index_names,
                       Pop   = dCI[[3]])
  
  return(CI)
}



sim_boot <- function(x, mle = 100){
  sample(x, mle, replace = TRUE)
}

#==============================================================================#
#' Perform bootstrap statistics, calculate and plot confidence intervals.
#' 
#' @param tab a genind object OR a matrix produced from 
#'   \code{\link[poppr]{mlg.table}}.
#' @param n an integer defining the number of bootstrap replicates (defaults to 
#'   1000).
#' @param ci the percent for confidence interval.
#' @param total argument to be passed on to \code{\link[poppr]{mlg.table}} if 
#'   \code{tab} is a genind object.
#' @param rarefy if \code{TRUE}, bootstrapping will be performed on the smallest
#'   population size. Defaults to \code{FALSE}, indicating that bootstrapping 
#'   will be performed respective to each population size.
#' @param plot If \code{TRUE} (default), boxplots will be produced for each 
#'   population, grouped by statistic. Colored dots will indicate the observed 
#'   value.This plot can be retrieved by using \code{p <- last_plot()} from the 
#'   \pkg{ggplot2} package.
#' @param raw if \code{TRUE} (default) a list containing three elements will be
#'   returned
#' @param ... parameters to be passed on to \code{\link[boot]{boot}} and 
#'   \code{\link{diversity_stats}}
#'   
#' @return \subsection{raw = TRUE}{
#' \itemize{
#' \item \strong{obs} - a matrix with observed statistics in columns, populations in rows
#' \item \strong{est} - a matrix with estimated statistics in columns, populations in rows
#' \item \strong{CI} - an array of 3 dimensions giving the lower and upper bound, the index measured, and the population.
#' \item \strong{boot} - a list containing the output of \code{\link[boot]{boot}} for each population.
#'   }
#'  }
#' \subsection{raw = FALSE}{ a data frame with the statistic observations,
#' estimates, and confidence intervals in columns, and populations in rows. Note
#' that the confidence intervals are converted to characters and rounded to
#' three decimal places. }
#' 
#' @details For details on the bootstrapping procedures, see
#'   \code{\link{diversity_boot}}. Confidence intervals are derived from the
#'   function \code{\link[boot]{norm.ci}}. The confidence interval is caluclated
#'   from the bootstrapped distribution and centered around the bias-corrected
#'   estimate. For rarefaction, the confidence interval is simply determined by
#'   caluclating the percentiles from the bootstrapped distribution. If you 
#'   want to calcuate your own confidence intervals, you can use the results of 
#'   the permutations stored in the \code{$boot} element of the output.
#'   
#' @note While it is possible to use custom functions with this, there are three
#'   important things to remember when using these functions:
#' \enumerate{
#' \item The function must return a single value.
#' \item The function must allow for both matrix and vector inputs
#' \item The function name cannot match or partially match any arguments from
#' \code{\link[boot]{boot}}
#' } Anonymous functions are okay \cr(e.g. \code{function(x) vegan::rarefy(t(as.matrix(x)), 10)}).
#'   
#' @export
#' @seealso \code{\link{diversity_boot}} \code{\link{diversity_stats}} \code{\link{poppr}} \code{\link[boot]{boot}} \code{\link[boot]{norm.ci}} \code{\link[boot]{boot.ci}}
#' @author Zhian N. Kamvar
#' @examples
#' library(poppr)
#' data(Pinf)
#' diversity_ci(Pinf, n = 100L)
#' \dontrun{
#' # With pretty results
#' diversity_ci(Pinf, n = 100L, raw = FALSE)
#' 
#' # This can be done in a parallel fasion (OSX uses "multicore", Windows uses "snow")
#' system.time(diversity_ci(Pinf, 10000L, parallel = "multicore", ncpus = 4L))
#' system.time(diversity_ci(Pinf, 10000L))
#' 
#' # We often get many requests for a clonal fraction statistic. As this is 
#' # simply the number of observed MLGs over the number of samples, we 
#' # recommeneded that people calculate it themselves. With this function, you
#' # can add it in:
#' 
#' CF <- function(x){
#'  x <- drop(as.matrix(x))
#'  if (length(dim(x)) > 1){
#'    res <- rowSums(x > 0)/rowSums(x)
#'  } else {
#'    res <- sum(x > 0)/sum(x)
#'  }
#'  return(res)
#' }
#' # Show pretty results
#' 
#' diversity_ci(Pinf, 1000L, CF = CF, raw = FALSE)
#' diversity_ci(Pinf, 1000L, CF = CF, rarefy = TRUE, raw = FALSE)
#' }
#' 
#==============================================================================#
diversity_ci <- function(tab, n = 1000, ci = 95, total = TRUE, rarefy = FALSE, 
                         plot = TRUE, raw = TRUE, ...){
  if (!is.matrix(tab) & is.genind(tab) | is(tab, "genlight")){
    tab <- mlg.table(tab, total = total, plot = FALSE)
  }
  rareval <- NULL
  if (rarefy){
    rareval <- min(rowSums(tab))
  }
  res  <- diversity_boot(tab, n, n.rare = rareval, ...)
  orig <- get_boot_stats(res)
  statnames <- colnames(orig)
  CI  <- get_all_ci(res, ci = ci, index_names = statnames,
                    center = ifelse(rarefy, FALSE, TRUE),
                    btype = ifelse(rarefy, "percent", "normal"),
                    bci = ifelse(rarefy, FALSE, TRUE))
  est <- get_boot_se(res, "mean")
  out <- list(obs = orig, est = est, CI = CI, boot = res)
  if (plot){
    if (rarefy){
      boot_plot(res, orig, statnames, rownames(tab), NULL)
    } else {
      boot_plot(res, orig, statnames, rownames(tab), CI)      
    }

  }
  if (!raw){
    out <- do.call("pretty_info", out)
  }
  return(out)
}

pretty_info <- function(obs, est, CI, boots = NULL){
  pretty_ci <- t(apply(round(CI, 3), 2:3, 
                       function(x){
                         if (all(is.na(x))) return(NA_character_)
                         paste0("(", paste(x, collapse = ", "), ")")
                      }))
  colnames(est) <- paste(colnames(est), "est", sep = ".")
  out <- vector(mode = "list", length = ncol(est)*3)
  colnames(pretty_ci) <- paste(colnames(pretty_ci), "ci", sep = ".")
  names(out)[(1:length(out)) %% 3 != 0] <- intersp(colnames(obs), colnames(est))
  names(out)[(1:length(out)) %% 3 == 0] <- colnames(pretty_ci)
  for (i in names(out)){
    out[[i]] <- obs[, 1]
  }
  out <- data.frame(out)
  out[colnames(obs)] <- obs
  out[colnames(est)] <- est
  out[colnames(pretty_ci)] <- pretty_ci
  class(out) <- c("popprtable", "data.frame")
  return(out)
}



boot_plot <- function(res, orig, statnames, popnames, CI){
  statnames <- colnames(orig)
  orig <- reshape2::melt(orig)
  orig$Pop <- factor(orig$Pop)
  if (!is.null(CI)){
    cidf <- reshape2::melt(CI)
    cidf <- reshape2::dcast(cidf, as.formula("Pop + Index ~ CI"))
    orig <- merge(orig, cidf)
    colnames(orig)[4:5] <- c("lb", "ub")    
  }
  samp <- vapply(res, "[[", FUN.VALUE = res[[1]]$t, "t")
  dimnames(samp) <- list(NULL, 
                         Index = statnames,
                         Pop = popnames)
  sampmelt <- melt(samp)
  sampmelt$Pop <- factor(sampmelt$Pop)
  pl <- ggplot(sampmelt, aes_string(x = "Pop", y = "value", group = "Pop")) + 
    geom_boxplot() + 
    geom_point(aes_string(color = "Pop", x = "Pop", y = "value"), 
               size = rel(4), pch = 16, data = orig) +
    xlab("Population") + labs(color = "Observed") +
    facet_wrap(~Index, scales = "free_y") + myTheme
  if (!is.null(CI)){
    pl <- pl +
      geom_errorbar(aes_string(color = "Pop", x = "Pop", ymin = "lb", ymax = "ub"),
                    data = orig)
  }
  print(pl)
}



get_boot_stats <- function(bootlist){
  npop   <- length(bootlist)
  bstats <- bootlist[[1]]$t0
  nstat  <- length(bstats)
  resmat <- matrix(nrow = npop, ncol = nstat,
                   dimnames = list(Pop = names(bootlist), Index = names(bstats)))
  resmat[] <- t(vapply(bootlist, FUN = "[[", FUN.VALUE = bstats, "t0"))
  return(resmat)
}

sd.boot <- function(x, na.rm = TRUE) apply(x$t, 2, sd, na.rm)
mean.boot <- function(x, ...) apply(x$t, 2, mean, na.rm = TRUE)

get_boot_se <- function(bootlist, res = "sd"){
  npop   <- length(bootlist)
  bstats <- bootlist[[1]]$t0
  nstat  <- length(bstats)
  THE_FUN <- match.fun(paste0(res, ".boot"))
  resmat <- matrix(nrow = npop, ncol = nstat,
                   dimnames = list(Pop = names(bootlist), Index = names(bstats)))
  resmat[] <- t(vapply(bootlist, FUN = THE_FUN, FUN.VALUE = bstats))
  if (res == "sd"){
    colnames(resmat) <- paste(colnames(resmat), res, sep = ".")  
  }
  return(resmat)
}

evens <- function(x){
  if (length(x) > 1){
    if (is.numeric(x)){
      res <- x %% 2 == 0
    } else {
      res <- 1:length(x) %% 2 == 0
    }
  } else {
    res <- seq(x) %% 2 == 0
  }
  return(res)
}

intersp <- function(v1, v2){
  v1l <- length(v1)
  v2l <- length(v2)
  the_evens <- evens(v1l + v2l)
  stopifnot(v1l == v2l)
  ov <- vector(length = v1l + v2l, mode = class(c(v1[1], v2[1])))
  ov[the_evens]  <- v2
  ov[!the_evens] <- v1
  return(ov)
}

boot_se_table <- function(tab, n = 1000, ci = 95, total = TRUE, rarefy = FALSE, 
                          rare.val = NULL, ...){
  if (!is.matrix(tab) & is.genind(tab)){
    tab <- mlg.table(tab, total = total, plot = FALSE)
  }
  rareval <- NULL
  
  if (rarefy){
    rareval <- ifelse(is.null(rare.val), min(rowSums(tab)), rare.val)
  }
  res <- diversity_boot(tab, n, n.rare = rareval, ...)
  if (rarefy){
    orig <- get_boot_se(res, "mean")
  } else {
    orig <- get_boot_stats(res)
  }
  out  <- matrix(nrow = nrow(orig), ncol = ncol(orig)*2, 
                 dimnames = list(rownames(orig), NULL))
  out[, !evens(1:ncol(out))]  <- orig
  se_out <- get_boot_se(res, "sd")
  out[, evens(1:ncol(out))] <- se_out
  colnames(out) <- intersp(colnames(orig), colnames(se_out))
  return(out)
}

rare_ia <- function(x, n = 1000, rare = 10, obs = FALSE){
  if (is.genind(x) || is.clone(x) || is(x, "genlight")){
    xloc <- seploc(x)
    est <- bootjack(xloc, n, rare, progbar = NULL)
    se  <- vapply(est, sd, numeric(1), na.rm = TRUE)
    if (obs){
      res <- matrix(nrow = 4, ncol = 1,
                    dimnames = list(c("Ia", "Ia.se", "rbarD", "rbarD.se"), NULL))
      est <- vapply(est, mean, numeric(1), na.rm = TRUE)
      res[c("Ia", "rbarD"), ] <- est
    } else {
      res <- matrix(nrow = 2, ncol = 1,
                    dimnames = list(c("Ia.se", "rbarD.se"), NULL))
    } 
    res[c("Ia.se", "rbarD.se"), ] <- se
  } else {
    if (obs){
      res <- matrix(nrow = 4, ncol = length(x),
                    dimnames = list(c("Ia", "Ia.se", "rbarD", "rbarD.se"), NULL))
    } else {
      res <- matrix(nrow = 2, ncol = length(x),
                    dimnames = list(c("Ia.se", "rbarD.se"), NULL))
    }
    for (i in seq(length(x))){
      if (nInd(x[[i]]) > rare){
        iloc <- seploc(x[[i]])
        est  <- bootjack(iloc, n , rare, progbar = NULL)
        se   <- vapply(est, sd, numeric(1), na.rm = TRUE)
        if (obs) est <- vapply(est, mean, numeric(1), na.rm = TRUE)
      } else {
        if (obs) est <- ia(x[[i]])
        se  <- c(0, 0)
      }
      if (obs) res[c("Ia", "rbarD"), i] <- est
      res[c("Ia.se", "rbarD.se"), i] <- se
    }
  }
  return(res)
}
