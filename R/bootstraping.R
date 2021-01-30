#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
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
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#==============================================================================#
#' Calculate a dendrogram with bootstrap support using any distance applicable 
#' to genind or genclone objects.
#' 
#' @param x a [genind-class][genind], [genpop-class][genpop],
#'   [genclone-class][genclone], [genlight-class], [snpclone-class] or
#'   \link{matrix} object.
#'   
#' @param strata a formula specifying the strata to be used to convert x to a
#'   genclone object if x is a genind object. Defaults to `NULL`. See details.
#' 
#' @param tree a text string or function that can calculate a tree from a 
#'   distance matrix. Defaults to "upgma". Note that you must load the package 
#'   with the function for it to work.
#'   
#' @param distance a character or function defining the distance to be applied 
#'   to x. Defaults to [nei.dist()].
#'   
#' @param sample An integer representing the number of bootstrap replicates 
#'   Default is 100.
#'   
#' @param cutoff An integer from 0 to 100 setting the cutoff value to return the
#'   bootstrap values on the nodes. Default is 0.
#'   
#' @param showtree If `TRUE` (Default), a dendrogram will be plotted. If 
#'   `FALSE`, nothing will be plotted.
#'   
#' @param missing any method to be used by [missingno()]: "mean" 
#'   (default), "zero", "loci", "genotype", or "ignore".
#'   
#' @param mcutoff a value between 0 (default) and 1 defining the percentage of 
#'   tolerable missing data if the `missing` parameter is set to "loci" or 
#'   "genotype". This should only be set if the distance metric can handle 
#'   missing data.
#'   
#' @param quiet if `FALSE` (default), a progress bar will be printed to 
#'   screen.
#'   
#' @param root is the tree rooted? This is a parameter passed off to 
#'   [ape::boot.phylo()]. If the `tree` parameter returns a 
#'   rooted tree (like UPGMA), this should be `TRUE`, otherwise (like 
#'   neighbor-joining), it should be false. When set to `NULL` (default), 
#'   the tree is considered rooted if [ape::is.ultrametric()] is true.
#'   
#' @param ... any parameters to be passed off to the distance method.
#'   
#' @return an object of class [ape::phylo()].
#'   
#' @details This function automates the process of bootstrapping genetic data to
#'   create a dendrogram with bootstrap support on the nodes. It will randomly
#'   sample with replacement the loci of a `gen` (genind/genpop) object or the columns of a 
#'   numeric matrix, **assuming that all loci/columns are independent**. The 
#'   process of randomly sampling `gen` objects with replacement is carried out
#'   through the use of an internal class called 
#'   [bootgen-class]. This is necessary due to the fact that columns in
#'   the genind matrix are defined as alleles and are thus interrelated. This 
#'   function will specifically bootstrap loci so that results are biologically 
#'   relevant. With this function, the user can also define a custom distance to
#'   be performed on the genind or genclone object. If you have a data frame-like
#'   object where all of the columns are independent or pairs of columns are
#'   independent, then it may be simpler to use [ape::boot.phylo()] to calculate
#'   your bootstrap support values. 
#'   
#'   \subsection{the strata argument}{
#'   There is an argument called `strata`. This argument is useful for when
#'   you want to bootstrap by populations from a [adegenet::genind()]
#'   object. When you specify strata, the genind object will be converted to
#'   [adegenet::genpop()] with the specified strata.
#'   }
#'   
#' @note [prevosti.dist()] and [diss.dist()] are exactly the
#'   same, but [diss.dist()] scales better for large numbers of 
#'   individuals (n > 125) at the cost of required memory. 
#'   \subsection{missing data}{
#'   Missing data is not allowed by many of the distances. Thus, one of 
#'   the first steps of this function is to treat missing data by setting it to 
#'   the average allele frequency in the data set. If you are using a distance 
#'   that can handle missing data (Prevosti's distance), you can set 
#'   `missing = "ignore"` to allow the distance function to handle any 
#'   missing data. See [missingno()] for details on missing 
#'   data.}
#'   \subsection{Bruvo's Distance}{
#'   While calculation of Bruvo's distance 
#'   is possible with this function, it is optimized in the function 
#'   [bruvo.boot()].}
#'   
#' @seealso [nei.dist()] [edwards.dist()] 
#'   [rogers.dist()] [reynolds.dist()] 
#'   [prevosti.dist()] [diss.dist()] 
#'   [bruvo.boot()] [ape::boot.phylo()] 
#'   [adegenet::dist.genpop()] [dist()]
#'   [bootgen2genind()] [bootgen-class]
#'   
#' @references
#' Kamvar ZN, Brooks JC and Grünwald NJ (2015) Novel R tools for analysis of 
#' genome-wide population genetic data with emphasis on clonality. 
#' Frontiers in Genetics 6:208. \doi{10.3389/fgene.2015.00208}
#' @export
#' @md
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
#' bindist <- function(x) dist(tab(x), method = "binary")
#' binnan <- aboot(nan9, dist = bindist)
#' 
#' \dontrun{
#' # Distances from other packages.
#' #
#' # Sometimes, distance functions from other packages will have the constraint
#' # that the incoming data MUST be genind. Internally, aboot uses the 
#' # bootgen class ( class?bootgen ) to shuffle loci, and will throw an error
#' # The function bootgen2genind helps fix that. Here's an example of a function
#' # that expects a genind class from above
#' bindist <- function(x){
#'   stopifnot(is.genind(x))
#'   dist(tab(x), method = "binary")
#' }
#' #
#' # Fails:
#' # aboot(nan9, dist = bindist)
#' ## Error: is.genind(x) is not TRUE
#' #
#' # Add bootgen2genind to get it working!
#' # Works:
#' aboot(nan9, dist = function(x) bootgen2genind(x) %>% bindist)
#' 
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
#' # You can also use the strata argument to convert to genpop inside the function.
#' set.seed(5000)
#' aboot(Aeut.gc, strata = ~Pop/Subpop, sample = 1000)
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
aboot <- function(x, strata = NULL, tree = "upgma", distance = "nei.dist", 
                  sample = 100, cutoff = 0, showtree = TRUE, missing = "mean", 
                  mcutoff = 0, quiet = FALSE, root = NULL, ...){
  if (!is.null(strata)){
    if (!is.genind(x)){
      warning("Sorry, the strata argument can only be used with genind objects.")
    } else {
      x <- genind2genpop(x, pop = strata, quiet = TRUE, process.other = FALSE)
    }
  }
  if (is.matrix(x)){
    if (is.null(rownames(x))) rownames(x) <- .genlab("", nrow(x))
    if (is.null(colnames(x))) colnames(x) <- .genlab("L", ncol(x))
    xboot <- x
  } else if (!is(x, "genlight") && x@type == "PA") {
    xboot           <- x@tab
    colnames(xboot) <- locNames(x)
    rownames(xboot) <- if (is.genpop(x)) popNames(x) else indNames(x)
  } else if (is(x, "gen")) {
    if (is.genind(x)) {
      if (missing %in% c("loci", "geno", "ignore")){
        x <- missingno(x, missing, quiet = quiet, cutoff = mcutoff)
        missing <- "asis"  
      } 
    }
    distname <- as.character(substitute(distance))
    if (length(distname) == 1 && distname %in% "diss.dist"){
      if (missing == "mean"){
        missing <- "asis"
        warning("Sorry, missing = 'mean' is incompatible with diss.dist(), setting missing to 'asis'.", 
                call. = FALSE)
      }
      xboot <- new("bootgen", x, na = missing, freq = FALSE)
    } else {
      xboot <- new("bootgen", x, na = missing, freq = TRUE)
    }
  } else if (is(x, "genlight")){
    xboot <- x
    my_dists <- c("diss.dist", "nei.dist", "prevosti.dist", "edwards.dist", 
                  "reynolds.dist", "rogers.dist", "provesti.dist")
    wrong_dist <- any(as.character(substitute(distance)) %in% my_dists)
    if (wrong_dist){
      warning("Sorry, distance from genlight objects can only be calculated by bitwise.dist.")
      distance <- bitwise.dist
    }
  } else {
    stop("Sorry, x must be a genind, genpop, or genlight object.")
  }
  treefunk <- tree_generator(tree, distance, ...)
  xtree    <- treefunk(xboot)
  if (any(xtree$edge.len < 0)){
    xtree <- fix_negative_branch(xtree)
    warning(negative_branch_warning())
  }
  treechar <- paste(substitute(tree), collapse = "")
  if (is.null(root)) {
    root <- ape::is.ultrametric(xtree)
  }
  nodelabs <- boot.phylo(xtree, xboot, treefunk, B = sample, rooted = root, 
                         quiet = quiet)
  nodelabs <- (nodelabs/sample)*100
  nodelabs <- ifelse(nodelabs >= cutoff, nodelabs, NA)
  if (!is.genpop(x)){
    if (is.matrix(x)){
      xtree$tip.label <- rownames(x)
    } else if (!is.null(indNames(x))) {
      xtree$tip.label <- indNames(x)
    }
  } else {
    xtree$tip.label <- popNames(x)
  }
  xtree$node.label <- nodelabs
  if (showtree) {
    poppr.plot.phylo(xtree, treechar, root)
  }
  return(xtree)
}

#==============================================================================#
#' Produce a table of diversity statistics
#' 
#' Calculate diversity statistics on a matrix containing counts of multilocus
#' genotypes per population.
#' 
#' @param z a table of integers representing counts of MLGs (columns) per 
#'   population (rows)
#'   
#' @param H logical whether or not to calculate Shannon's index
#' @param G logical whether or not to calculate Stoddart and Taylor's index (aka
#'   inverse Simpson's index).
#' @param lambda logical whether or not to calculate Simpson's index
#' @param E5 logical whether or not to calculate Evenness
#' @param ... any functions that can be calculated on a vector or matrix of 
#'   genotype counts.
#'   
#' @return a numeric matrix giving statistics (columns) for each population 
#'   (rows).
#'   
#' @details This function will calculate any diversity statistic for counts of 
#'   multilocus genotypes per population. This does not count allelic diversity.
#'   The calculations of H, G, and lambda are all performed by 
#'   [vegan::diversity()]. E5 is calculated as \deqn{E_{5} = 
#'   \frac{(1/\lambda) - 1}{e^{H} - 1}}{(G - 1)/(exp(H) - 1)}.
#'   
#' @export
#' @md
#' @seealso [diversity_boot()] [diversity_ci()]
#'   [poppr()]
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
  
  # Indicator to determine if the input is from diversity_boot
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
#' This function is intended to perform bootstrap statistics on a matrix of
#' multilocus genotype counts in different populations. Results from this 
#' function should be interpreted carefully as the default statistics are known
#' to have a downward bias. See the details for more information.
#' 
#' 
#' @param tab a table produced from the \pkg{poppr} function 
#'   [mlg.table()]. MLGs in columns and populations in rows
#' @param n an integer > 0 specifying the number of bootstrap replicates to 
#'   perform (corresponds to `R` in the function [boot::boot()].
#' @param n.boot an integer specifying the number of samples to be drawn in each
#'   bootstrap replicate. If `n.boot` < 2 (default), the number of samples 
#'   drawn for each bootstrap replicate will be equal to the number of samples in
#'   the data set. 
#' @param n.rare a sample size at which all resamplings should be performed. 
#'   This should be no larger than the smallest sample size. Defaults to 
#'   `NULL`, indicating that each population will be sampled at its own 
#'   size.
#' @inheritParams diversity_stats
#' @param ... other parameters passed on to [boot::boot()] and 
#'   [diversity_stats()].
#'   
#' @return a list of objects of class "boot".
#' @seealso [diversity_stats()] for basic statistic calculation, 
#'    [diversity_ci()] for confidence intervals and plotting, and 
#'    [poppr()]. For bootstrap sampling:
#'    [stats::rmultinom()] [boot::boot()]
#'    
#' @details 
#'   Bootstrapping is performed in three ways:
#'   \itemize{
#'     \item if `n.rare` is a number greater than zero, then bootstrapping
#'     is performed by randomly sampling without replacement *n.rare*
#'     samples from the data.
#'     
#'     \item if `n.boot` is greater than 1, bootstrapping is performed by
#'     sampling n.boot samples from a multinomial distribution weighted by the
#'     proportion of each MLG in the data.
#'     
#'     \item if `n.boot` is less than 2, bootstrapping is performed by 
#'     sampling N samples from a multinomial distribution weighted by the
#'     proportion of each MLG in the data.
#'   }
#'   
#'   \subsection{Downward Bias}{
#'     When sampling with replacement, the diversity statistics here present a 
#'     downward bias partially due to the small number of samples in the data. 
#'     The result is that the mean of the bootstrapped samples will often be 
#'     much lower than the observed value. Alternatively, you can increase the
#'     sample size of the bootstrap by increasing the size of `n.boot`. Both
#'     of these methods should be taken with caution in interpretation. There
#'     are several R packages freely available that will calculate and perform
#'     bootstrap estimates of Shannon and Simpson diversity metrics (eg.
#'     \pkg{entropart}, \pkg{entropy}, \pkg{simboot}, and
#'     \pkg{EntropyEstimation}. These packages also offer unbiased estimators of
#'     Shannon and Simpson diversity. Please take care when attempting to
#'     interpret the results of this function.
#'   }
#' @export
#' @md
#' @author Zhian N. Kamvar
#' @examples
#' library(poppr)
#' data(Pinf)
#' tab <- mlg.table(Pinf, plot = FALSE)
#' diversity_boot(tab, 10L)
#' \dontrun{
#' # This can be done in a parallel fashion (OSX uses "multicore", Windows uses "snow")
#' system.time(diversity_boot(tab, 10000L, parallel = "multicore", ncpus = 4L))
#' system.time(diversity_boot(tab, 10000L))
#' }
#' @importFrom boot boot boot.ci norm.ci
#==============================================================================#
diversity_boot <- function(tab, n, n.boot = 1L, n.rare = NULL, H = TRUE, 
                           G = TRUE, lambda = TRUE, E5 = TRUE, ...){
  if (!is.null(n.rare)){
    FUN <- rare_sim_boot
    mle <- n.rare
  } else {
    FUN <- multinom_boot
    mle <- n.boot
  }
  res <- apply(tab, 1, boot_per_pop, rg = FUN, n = n, mle = mle, H = H, G = G, 
               lambda = lambda, E5 = E5, ...)
  return(res)
}

#==============================================================================#
#' Perform bootstrap statistics, calculate, and plot confidence intervals.
#' 
#' This function is for calculating bootstrap statistics and their confidence 
#' intervals. It is important to note that the calculation of confidence 
#' intervals is not perfect (See Details). Please be cautious when interpreting 
#' the results.
#' 
#' @param tab a [genind()], [genclone()],
#'   [snpclone()], OR a matrix produced from 
#'   [poppr::mlg.table()].
#' @param n an integer defining the number of bootstrap replicates (defaults to 
#'   1000).
#' @param n.boot an integer specifying the number of samples to be drawn in each
#'   bootstrap replicate. If `n.boot` < 2 (default), the number of samples 
#'   drawn for each bootstrap replicate will be equal to the number of samples
#'   in the data set. See Details.
#' @param ci the percent for confidence interval.
#' @param total argument to be passed on to [poppr::mlg.table()] if 
#'   `tab` is a genind object.
#' @param rarefy if `TRUE`, bootstrapping will be performed on the smallest
#'   population size or the value of `n.rare`, whichever is larger. 
#'   Defaults to `FALSE`, indicating that bootstrapping will be performed 
#'   respective to each population size.
#' @param n.rare an integer specifying the smallest size at which to resample 
#'   data. This is only used if `rarefy = TRUE`.
#' @param plot If `TRUE` (default), boxplots will be produced for each 
#'   population, grouped by statistic. Colored dots will indicate the observed 
#'   value.This plot can be retrieved by using `p <- last_plot()` from the 
#'   \pkg{ggplot2} package.
#' @param raw if `TRUE` (default) a list containing three elements will be 
#'   returned
#' @param center if `TRUE` (default), the confidence interval will be 
#'   centered around the observed statistic. Otherwise, if `FALSE`, the 
#'   confidence interval will be bias-corrected normal CI as reported from 
#'   [boot::boot.ci()]
#' @param ... parameters to be passed on to [boot::boot()] and 
#'   [diversity_stats()]
#'   
#' @return \subsection{raw = TRUE}{
#' 
#'    - **obs** a matrix with observed statistics in columns,
#'   populations in rows
#'    - **est** a matrix with estimated statistics in columns,
#'   populations in rows
#'    - **CI** an array of 3 dimensions giving the lower and upper
#'   bound, the index measured, and the population.
#'    - **boot** a list containing the output of
#'   [boot::boot()] for each population.
#'   }
#' 
#' \subsection{raw = FALSE}{ a data frame with the statistic observations,
#' estimates, and confidence intervals in columns, and populations in rows. Note
#' that the confidence intervals are converted to characters and rounded to
#' three decimal places. }
#' 
#' @details 
#'   \subsection{Bootstrapping}{
#'   For details on the bootstrapping procedures, see 
#'   [diversity_boot()]. Default bootstrapping is performed by 
#'   sampling \strong{N} samples from a multinomial distribution weighted by the
#'   relative multilocus genotype abundance per population where \strong{N} is
#'   equal to the number of samples in the data set. If \strong{n.boot} > 2,
#'   then \strong{n.boot} samples are taken at each bootstrap replicate. When
#'   `rarefy = TRUE`, then samples are taken at the smallest population
#'   size without replacement. This will provide confidence intervals for all
#'   but the smallest population.
#'   }
#'   \subsection{Confidence intervals}{
#'   Confidence intervals are derived from the function 
#'   [boot::norm.ci()]. This function will attempt to correct for bias
#'   between the observed value and the bootstrapped estimate. When `center
#'   = TRUE` (default), the confidence interval is calculated from the
#'   bootstrapped distribution and centered around the bias-corrected estimate
#'   as prescribed in Marcon (2012). This method can lead to undesirable
#'   properties, such as the confidence interval lying outside of the maximum
#'   possible value. For rarefaction, the confidence interval is simply
#'   determined by calculating the percentiles from the bootstrapped
#'   distribution. If you want to calculate your own confidence intervals, you
#'   can use the results of the permutations stored in the `$boot` element
#'   of the output.
#'   }
#'   \subsection{Rarefaction}{
#'   Rarefaction in the sense of this function is simply sampling a subset of 
#'   the data at size **n.rare**. The estimates derived from this method
#'   have straightforward interpretations and allow you to compare diversity
#'   across populations since you are controlling for sample size.
#'   }
#'   \subsection{Plotting}{ Results are plotted as boxplots with point
#'   estimates. If there is no rarefaction applied, confidence intervals are
#'   displayed around the point estimates. The boxplots represent the actual
#'   values from the bootstrapping and will often appear below the estimates and
#'   confidence intervals. 
#'   }
#' @note 
#'   \subsection{Confidence interval calculation}{ Almost all of the statistics 
#'   supplied here have a maximum when all genotypes are equally represented. 
#'   This means that bootstrapping the samples will always be downwardly biased.
#'   In many cases, the confidence intervals from the bootstrapped distribution 
#'   will fall outside of the observed statistic. The reported confidence 
#'   intervals here are reported by assuming the variance of the bootstrapped 
#'   distribution is the same as the variance around the observed statistic. As
#'   different statistics have different properties, there will not always be
#'   one clear method for calculating confidence intervals. A suggestion for
#'   correction in Shannon's index is to center the CI around the observed
#'   statistic (Marcon, 2012), but there are theoretical limitations to this.
#'   For details, see <http://stats.stackexchange.com/q/156235/49413>.
#'   }
#' 
#'   \subsection{User-defined functions}{ 
#'   While it is possible to use custom functions with this, there are three
#'   important things to remember when using these functions:
#'   
#'     1. The function must return a single value. 
#'     2. The function must allow for both matrix and vector inputs 
#'     3. The function name cannot match or partially match any arguments 
#'     from [boot::boot()]
#'    
#'   Anonymous functions are okay \cr(e.g. `function(x)
#'   vegan::rarefy(t(as.matrix(x)), 10)`).
#'   }
#' @export
#' @seealso [diversity_boot()] [diversity_stats()]
#'   [poppr()] [boot::boot()] [boot::norm.ci()]
#'   [boot::boot.ci()]
#' @author Zhian N. Kamvar
#' @references
#' Marcon, E., Herault, B., Baraloto, C. and Lang, G. (2012). The Decomposition 
#' of Shannon’s Entropy and a Confidence Interval for Beta Diversity.
#' *Oikos* 121(4): 516-522.
#' 
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
#' # recommended that people calculate it themselves. With this function, you
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
#' diversity_ci(Pinf, 1000L, CF = CF, center = TRUE, raw = FALSE)
#' diversity_ci(Pinf, 1000L, CF = CF, rarefy = TRUE, raw = FALSE)
#' }
#' 
#==============================================================================#
diversity_ci <- function(tab, n = 1000, n.boot = 1L, ci = 95, total = TRUE, 
                         rarefy = FALSE, n.rare = 10, plot = TRUE, raw = TRUE, 
                         center = TRUE, ...){
  allowed_objects <- c("genind", "genclone", "snpclone")
  if (!is.matrix(tab) & inherits(tab, allowed_objects)){
    tab <- mlg.table(tab, total = total, plot = FALSE)
  }
  rareval <- NULL
  if (rarefy){
    rareval <- max(min(rowSums(tab)), n.rare)
    msg     <- paste("\nSamples for rarefaction:", rareval)
    message(msg)
  } else {
    if (center == TRUE){
      msg <- paste("\nConfidence Intervals have been centered around observed",
                   "statistic.\nPlease see ?diversity_ci for details.\n")
      message(msg)
    } else if (n.boot > 2){
      msg <- paste("\nTake caution in interpreting the results.\nWhile sampling",
                   n.boot, "samples will produce a confidence interval,\nit",
                   "might be unclear what the correct interpretation should be.\n")
      message(msg)
    }
    
  }
  res  <- diversity_boot(tab, n, n.rare = rareval, n.boot = n.boot, ...)
  orig <- get_boot_stats(res)
  statnames <- colnames(orig)
  CI  <- get_all_ci(res, ci = ci, index_names = statnames, 
                    center = if (rarefy) FALSE else center,
                    btype =  if (rarefy) "percent" else "normal",
                    bci =    if (rarefy) FALSE else TRUE)
  est <- get_boot_se(res, "mean")
  out <- list(obs = orig, est = est, CI = CI, boot = res)
  if (plot) {
    boot_plot(res = res, 
              orig = orig,  
              CI = if (rarefy) NULL else CI
              )
  }
  out <- if (raw) out else do.call("pretty_info", out)
  return(out)
}

