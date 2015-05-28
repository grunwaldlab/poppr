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
#' Aeut.gc <- as.genclone(Aeut, hierarchy=other(Aeut)$population_hierarchy[-1])
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
#' @param simp logical whether or not to calculate Simpson's index
#' @param E5 logical whether or not to calculate Evenness
#' @param ... any functions that can be calculated on a vector or matrix of
#'   genotype counts.
#'   
#' @return a numeric matrix giving statistics (columns) for each population
#'   (rows).
#' @export
#' @author Zhian N. Kamvar
#' @examples
#' library(poppr)
#' data(Pinf)
#' tab <- mlg.table(Pinf, plot = FALSE)
#' get_stats(tab)
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
#' get_stats(tab, B = Beta)
#' }
#==============================================================================#
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
    if (i == "E.5" && H && G){
      mat[, i] <- (mat[, "G"] - 1)/(exp(mat[, "H"]) - 1)
    } else {
      mat[, i] <- FUNS[[i]](z)
    }
  }
  return(drop(mat))
}

boot_stats <- function(x, i, H = TRUE, G = TRUE, simp = TRUE, E5 = TRUE, ...){
  res <- get_stats(tabulate(x[i]), H, G, simp, E5, ...)
  return(res)
}

extract_samples <- function(x) rep(1:length(x), x)

#==============================================================================#
#' Perform a bootstrap analysis on diversity statistics
#' 
#' @param tab a table produced from the \pkg{poppr} function
#'   \code{\link[poppr]{mlg.table}}. MLGs in columns and populations in rows
#' @param n an integer > 0 specifying the number of bootstrap replicates to
#'   perform (corresponds to \code{R} in the function \code{\link[boot]{boot}}.
#' @inheritParams get_stats
#' @param ... other parameters passed on to \code{\link[boot]{boot}} and
#'   \code{\link{get_stats}}.
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
#==============================================================================#
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

get_all_ci <- function(res, ci = 95, index_names = c("H", "G", "Hexp", "E.5")){
  lower_bound  <- (100 - ci)/200
  upper_bound  <- 1 - lower_bound
  n_indices    <- length(index_names)
  funval       <- matrix(numeric(n_indices*2), nrow = 2)
  CI           <- vapply(res, FUN = get_ci, FUN.VALUE = funval, 
                         lower_bound, upper_bound)
  dCI          <- dimnames(CI)
  dimnames(CI) <- list(CI    = dCI[[1]], 
                       Index = index_names,
                       Pop   = dCI[[3]])
  
  return(CI)
}

rare_sim_boot <- function(x, mle = 10){
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
#' @param ... parameters to be passed on to \code{\link[boot]{boot}} and
#'   \code{\link{get_stats}}
#'   
#' @return an array of 3 dimensions giving the lower and upper bound, the index 
#'   measured, and the population. This also prints barplots for each population
#'   faceted by each index using \pkg{ggplot2}. This plot can be retrieved by
#'   using \code{p <- last_plot()}
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
#==============================================================================#
boot_ci <- function(tab, n = 1000, ci = 95, total = TRUE, ...){
  if (!is.matrix(tab) & is.genind(tab)){
    tab <- mlg.table(tab, total = total, plot = FALSE)
  }
  res  <- do_boot(tab, n, ...)
  dotlist <- list(...)
  bootargs <- names(formals(boot))
  dotlist <- dotlist[!names(dotlist) %in% bootargs]
  get_stats_args <- formals(get_stats)
  get_stats_args$z <- tab
  matching_names <- names(dotlist) %in% names(get_stats_args)
  if (any(matching_names)){
    for (i in names(dotlist[matching_names])){
      if (!is.logical(dotlist[[i]])){
        warning("Extra functions cannot have the same name a defaults. Renaming")
        names(dotlist)[names(dotlist) == i] <- paste0("d", i)
      } else {
        get_stats_args[[i]] <- dotlist[[i]]
      }
    }
  }
  unique_names <- !names(dotlist) %in% names(get_stats_args)
  dotlist <- dotlist[unique_names]
  orig <- do.call("get_stats", c(get_stats_args[-6], dotlist))
  statnames <- colnames(orig)
  orig <- melt(orig)
  orig$Pop <- factor(orig$Pop)
  CI   <- get_all_ci(res, ci = ci, index_names = statnames)
  samp <- vapply(res, "[[", FUN.VALUE = res[[1]]$t, "t")
  dimnames(samp) <- list(NULL, 
                         Index = statnames,
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
