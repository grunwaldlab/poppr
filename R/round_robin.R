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
#' Round Robin Multilocus Genotypes
#' 
#' This function will mask each locus one by one and then calculate multilocus
#' genotypes from the remaining loci in a round-robin fashion. This is used for
#' calculating the round robin allele frequencies for pgen and psex.
#' 
#' @param gid a genind, genclone, or loci object.
#' 
#' @author Zhian N. Kamvar, Jonah Brooks, Stacy A. Krueger-Hadfield, Erik Sotka
#' 
#' @return a matrix of multilocus genotype assignments by masked locus. There 
#'   will be n rows and m columns where n = number of samples and m = number of
#'   loci.
#'
#' @export
#' @seealso \code{\link{rraf}}, \code{\link{pgen}}, \code{\link{psex}}
#' @references
#' 
#' Arnaud-Haond, S., Duarte, C. M., Alberto, F., & Serrão, E. A. 2007.
#' Standardizing methods to address clonality in population studies.
#' \emph{Molecular Ecology}, 16(24), 5115-5139.
#' 
#' Parks, J. C., & Werth, C. R. 1993. A study of spatial features of clones in a
#' population of bracken fern, \emph{Pteridium aquilinum} (Dennstaedtiaceae).
#' \emph{American Journal of Botany}, 537-544.
#' 
#' 
#' @examples
#' 
#' # Find out the round-robin multilocus genotype assignments for P. ramorum
#' data(Pram)
#' pmlg_rr <- rrmlg(Pram)
#' head(pmlg_rr)
#' \dontrun{
#' # You can find out how many unique genotypes are found without each locus:
#' 
#' colSums(!apply(pmlg_rr, 2, duplicated))
#' }
#==============================================================================#
rrmlg <- function(gid){
  if (class(gid)[1] %in% c("genind", "genclone")){
    gid <- pegas::as.loci(gid)
  }
  the_loci <- attr(gid, "locicol")
  res      <- integer(nrow(gid))
  suppressWarnings(gid <- vapply(gid[the_loci], as.integer, res))
  out <- .Call("mlg_round_robin", gid, PACKAGE = "poppr")
  dimnames(out) <- dimnames(gid)
  return(out)
}

#==============================================================================#
#' Round Robin Allele Frequencies
#' 
#' This function utilizes \code{\link{rrmlg}} to calculate multilocus genotypes 
#' and then subsets each locus by the resulting MLGs to calculate the 
#' round-robin allele frequencies used for pgen and psex.
#' 
#' @param gid a genind or genclone object
#' @param pop either a formula to set the population factor from the
#'   \code{\link{strata}} slot or a vector specifying the population factor for
#'   each sample. Defaults to NULL. 
#' @param res Either "list" (default), "vector", or "data.frame".
#' @param by_pop When this is \code{TRUE}, the calculation will be done by
#'   population. Defaults to \code{FALSE}
#' @param correction a logical indicating whether or not zero value allele 
#'   frequencies should be set to 1/n (Default: \code{TRUE})
#'   
#' @return a vector or list of allele frequencies
#' @note When \code{by_pop = TRUE}, the output will be a matrix of allele 
#'   frequencies. Additionally, when the argument \code{pop} is not \code{NULL},
#'   \code{by_pop} is automatically \code{TRUE}. When \code{correction = TRUE},
#'   \emph{n} refers to the population size.
#' 
#' @author Zhian N. Kamvar, Jonah Brooks, Stacy A. Krueger-Hadfield, Erik Sotka
#' @references
#' 
#' Arnaud-Haond, S., Duarte, C. M., Alberto, F., & Serrão, E. A. 2007.
#' Standardizing methods to address clonality in population studies.
#' \emph{Molecular Ecology}, 16(24), 5115-5139.
#' 
#' Parks, J. C., & Werth, C. R. 1993. A study of spatial features of clones in a
#' population of bracken fern, \emph{Pteridium aquilinum} (Dennstaedtiaceae).
#' \emph{American Journal of Botany}, 537-544.
#' 
#' @export
#' @seealso \code{\link{rrmlg}}, \code{\link{pgen}}, \code{\link{psex}}
#' @examples
#' 
#' data(Pram)
#' # Round robin allele frequencies.
#' rraf(Pram)
#' 
#' 
#' \dontrun{
#' 
#' # Compare to without round robin:
#' PrLoc <- seploc(Pram, res = "mat") # get locus by matrix
#' lapply(PrLoc, colMeans, na.rm = TRUE)
#' 
#' # Without round robin, clone corrected:
#' Pcc <- clonecorrect(Pram, strata = NA) # indiscriminantly clone correct
#' PccLoc <- seploc(Pcc, res = "mat")
#' lapply(PccLoc, colMeans, na.rm = TRUE)
#' 
#' # Get vector output.
#' rraf(Pram, res = "vector")
#' 
#' # Get frequencies per population (matrix only)
#' rraf(Pram, by_pop = TRUE, correction = FALSE)
#' 
#' # Get data frame output and plot.
#' (Prdf <- rraf(Pram, res = "data.frame"))
#' library("ggplot2")
#' ggplot(Prdf, aes(y = allele, x = frequency)) +
#'   geom_point() +
#'   facet_grid(locus ~ ., scale = "free_y", space = "free")
#' }
#==============================================================================#
rraf <- function(gid, pop = NULL, res = "list", by_pop = FALSE, correction = TRUE){
  RES     <- c("list", "vector", "data.frame")
  res     <- match.arg(res, RES)
  if (!is.null(pop)){
    if (!is.language(pop)){
      pop(gid) <- pop
    } else {
      setPop(gid) <- pop
    }
    by_pop <- TRUE
  }
  loclist <- seploc(gid)
  mlgs    <- rrmlg(gid)
  if (by_pop & !is.null(pop(gid))){
    out <- matrix(numeric(0), nrow = nPop(gid), ncol = ncol(tab(gid)))
    for (i in locNames(gid)){
      out[, locFac(gid) %in% i] <- rrccbp(i, loclist, mlgs, correction, popNames(gid))
    }
    rownames(out) <- popNames(gid)
    colnames(out) <- colnames(tab(gid))
    return(out)
  } else {
    out <- lapply(locNames(gid), rrcc, loclist, mlgs, correction)
  }
  names(out) <- locNames(gid)
  if (res == "vector"){
    out <- unlist(out, use.names = FALSE)
    names(out) <- colnames(tab(gid))
  } else if (res == "data.frame"){
    outdf <- reshape2::melt(out)
    names(outdf)        <- c("frequency", "locus")
    levels(outdf$locus) <- locNames(gid)
    outdf <- cbind(outdf, allele = unlist(lapply(out, names), use.names = FALSE))
    return(outdf)
  }
  return(out)
}
#==============================================================================#
#' Genotype Probability
#'
#' @inheritParams rraf
#' 
#' @param gid a genind or genclone object. 
#' 
#' @param by_pop When this is \code{TRUE} (default), the calculation will be
#'   done by population.
#'
#' @param log a \code{logical} if \code{log =TRUE} (default), the values
#'   returned will be log(Pgen). If \code{log = FALSE}, the values returned will
#'   be Pgen.
#'   
#' @param freq a vector or matrix of allele frequencies. This defaults to 
#'   \code{NULL}, indicating that the frequencies will be determined via 
#'   round-robin approach. \strong{If \code{by_pop = TRUE}, zero-value allele
#'   frequencies will automatically be corrected by 1/n.} See \code{\link{rraf}}
#'   for details.
#'   
#' @note For haploids, Pgen at a particular locus is the allele frequency. This
#'   function cannot handle polyploids. Additionally, when the argument
#'   \code{pop} is not \code{NULL}, \code{by_pop} is automatically \code{TRUE}.
#'
#' @return A vector containing Pgen values per locus for each genotype in the 
#'   object.
#'   
#' @details Pgen is the probability of a given genotype occuring in a population
#'   assuming HWE. Thus, the value for diploids is 
#'   
#'   \deqn{P_{gen} = \left(\prod_{i=1}^m p_i\right)2^h}{pgen = prod(p_i)*(2^h)} 
#'   
#'   where \eqn{p_i} are
#'   the allele frequencies and \emph{h} is the count of the number of
#'   heterozygous sites in the sample (Arnaud-Haond et al. 2007; Parks and
#'   Werth, 1993). The allele frequencies, by default, are calculated using a
#'   round-robin approach where allele frequencies at a particular locus are
#'   calculated on the clone-censored genotypes without that locus. 
#'   
#'   To avoid issues with numerical precision of small numbers, this function 
#'   calculates pgen per locus by adding up log-transformed values of allele 
#'   frequencies. These can easily be transformed to return the true value (see
#'   examples).
#'   
#' @author Zhian N. Kamvar, Jonah Brooks, Stacy A. Krueger-Hadfield, Erik Sotka
#' @seealso \code{\link{psex}}, \code{\link{rraf}}, \code{\link{rrmlg}}
#' @references
#' 
#' Arnaud-Haond, S., Duarte, C. M., Alberto, F., & Serrão, E. A. 2007.
#' Standardizing methods to address clonality in population studies.
#' \emph{Molecular Ecology}, 16(24), 5115-5139.
#' 
#' Parks, J. C., & Werth, C. R. 1993. A study of spatial features of clones in a
#' population of bracken fern, \emph{Pteridium aquilinum} (Dennstaedtiaceae).
#' \emph{American Journal of Botany}, 537-544.
#' 
#' @export
#' @examples
#' data(Pram)
#' head(pgen(Pram, log = FALSE))
#' 
#' \dontrun{
#' # You can get the Pgen values over all loci by summing over the logged results:
#' exp(rowSums(pgen(Pram, log = TRUE, na.rm = TRUE)))
#' 
#' # You can also take the product of the non-logged results:
#' apply(pgen(Pram, log = FALSE), 1, prod, na.rm = TRUE)
#' }
#==============================================================================#
pgen <- function(gid, pop = NULL, by_pop = TRUE, log = TRUE, freq = NULL){
  stopifnot(is.genind(gid))
  # Stop if the ploidy of the object is not diploid
  stopifnot(all(ploidy(gid) %in% 1:2)) 
  if (!is.null(pop)){
    if (!is.language(pop)){
      pop(gid) <- pop
    } else {
      setPop(gid) <- pop
    }
    by_pop <- TRUE
  }
  # Set single population if none exists OR if user explicitly said by_pop = FALSE
  if (is.null(pop(gid)) || !by_pop){
    pop(gid) <- rep(1, nInd(gid))
  }   
  pops <- pop(gid)
  if (is.null(freq)){
    freqs <- rraf(gid, by_pop = TRUE)
  } else if (is.matrix(freq)){
    if (nrow(freq) != nlevels(pops) || ncol(freq) != ncol(tab(gid))){
      stop("frequency matrix must have the same dimensions as the data.")
    }
    freqs <- freq
  } else if (is.numeric(freq) && length(freq) == ncol(tab(gid))){
    freqs <- freq
    pops <- factor(rep(1, nInd(gid)))
  } else {
    stop("frequencies must be the same length as the number of alleles.")
  }
  
  pgen_matrix <- .Call("get_pgen_matrix_genind", gid, freqs, pops, nlevels(pops), 
                       PACKAGE = "poppr")
  if (!log)
  {
    pgen_matrix <- exp(pgen_matrix)
  }
  dimnames(pgen_matrix) <- list(indNames(gid), locNames(gid))

  return(pgen_matrix);
}
#==============================================================================#
#' Probability of encountering a genotype more than once by chance
#' 
#' @inheritParams pgen
#' @param G an integer specifying the number of observed genets. If NULL, this 
#'   will be the number of original multilocus genotypes.
#' @param method which method of calculating psex should be used? Using 
#'   \code{method = "single"} (default) indicates that the calculation for psex 
#'   should reflect the probability of encountering a second genotype. Using 
#'   \code{method = "multiple"} gives the probability of encountering multiple 
#'   samples of the same genotype (see details).
#' @return a vector of Psex for each sample.
#' @note The values of Psex represent the value for each multilocus genotype. 
#'   Additionally, when the argument \code{pop} is not \code{NULL}, 
#'   \code{by_pop} is automatically \code{TRUE}.
#'   
#' @author Zhian N. Kamvar, Jonah Brooks, Stacy A. Krueger-Hadfield, Erik Sotka
#'   
#' @details Psex is the probability of encountering a given genotype more than 
#'   once by chance. The basic equation is
#'   
#'   \deqn{p_{sex} = 1 - (1 - p_{gen})^{G})}{psex = 1 - (1 - pgen)^G}
#'   
#'   where \emph{G} is the number of multilocus genotypes. See
#'   \code{\link{pgen}} for its calculation. For a given value of alpha (e.g.
#'   alpha = 0.05), genotypes with psex < alpha can be thought of as a single
#'   genet whereas genotypes with psex > alpha do not have strong evidence that
#'   members belong to the same genet (Parks and Werth, 1993).
#'   
#'   When \code{method = "multiple"}, the method from Arnaud-Haond et al. (1997)
#'   is used where the sum of the binomial density is taken:
#'   
#'   \deqn{p_{sex} = \sum_{i = 1}^N {N \choose i} \left(p_{gen}\right)^i\left(1
#'   - p_{gen}\right)^{N - i}}{psex = sum(dbinom(1:N, N, pgen))}
#'   
#'   where \emph{N} is the number of samples with the same genotype, \emph{i} is
#'   the ith sample, and \emph{pgen} is the value of pgen for that genotype.
#'   
#'   The function will automatically calculate the round-robin allele
#'   frequencies with \code{\link{rraf}} and \emph{G} with \code{\link{nmll}}.
#'   
#' @seealso \code{\link{pgen}}, \code{\link{rraf}}, \code{\link{rrmlg}}
#' @references
#' 
#' Arnaud-Haond, S., Duarte, C. M., Alberto, F., & Serrão, E. A. 2007. 
#' Standardizing methods to address clonality in population studies. 
#' \emph{Molecular Ecology}, 16(24), 5115-5139.
#' 
#' Parks, J. C., & Werth, C. R. 1993. A study of spatial features of clones in a
#' population of bracken fern, \emph{Pteridium aquilinum} (Dennstaedtiaceae). 
#' \emph{American Journal of Botany}, 537-544.
#' 
#' @export
#' @examples
#' 
#' data(Pram)
#' Pram_psex <- psex(Pram, by_pop = FALSE)
#' plot(Pram_psex, log = "y", col = ifelse(Pram_psex > 0.05, "red", "blue"))
#' abline(h = 0.05, lty = 2)
#' \dontrun{
#' 
#' # With multiple encounters
#' Pram_psex <- psex(Pram, by_pop = FALSE, method = "multiple")
#' plot(Pram_psex, log = "y", col = ifelse(Pram_psex > 0.05, "red", "blue"))
#' abline(h = 0.05, lty = 2)
#' 
#' # This can be also done assuming populations structure
#' Pram_psex <- psex(Pram, by_pop = TRUE)
#' plot(Pram_psex, log = "y", col = ifelse(Pram_psex > 0.05, "red", "blue"))
#' abline(h = 0.05, lty = 2)
#' 
#' ## An example of supplying previously calculated frequencies and G
#' # From Parks and Werth, 1993, using the first three genotypes.
#' 
#' # The row names indicate the number of samples found with that genotype
#' x <- "
#'  Hk Lap Mdh2 Pgm1 Pgm2 X6Pgd2
#' 54 12 12 12 23 22 11
#' 36 22 22 11 22 33 11
#' 10 23 22 11 33 13 13"
#' 
#' # Since we aren't representing the whole data set here, we are defining the
#' # allele frequencies before the analysis.
#' afreq <- c(Hk.1 = 0.167, Hk.2 = 0.795, Hk.3 = 0.038, 
#'            Lap.1 = 0.190, Lap.2 = 0.798, Lap.3 = 0.012,
#'            Mdh2.0 = 0.011, Mdh2.1 = 0.967, Mdh2.2 = 0.022,
#'            Pgm1.2 = 0.279, Pgm1.3 = 0.529, Pgm1.4 = 0.162, Pgm1.5 = 0.029,
#'            Pgm2.1 = 0.128, Pgm2.2 = 0.385, Pgm2.3 = 0.487,
#'            X6Pgd2.1 = 0.526, X6Pgd2.2 = 0.051, X6Pgd2.3 = 0.423)
#'
#' xtab <- read.table(text = x, header = TRUE, row.names = 1)
#' 
#' # Here we are expanding the number of samples to their observed values.
#' # Since we have already defined the allele frequencies, this step is actually
#' # not necessary. 
#' all_samples <- rep(rownames(xtab), as.integer(rownames(xtab)))
#' xgid        <- df2genind(xtab[all_samples, ], ncode = 1)
#' 
#' freqs <- afreq[colnames(tab(xgid))] # only used alleles in the sample
#' pSex  <- psex(xgid, by_pop = FALSE, freq = freqs, G = 45)
#' 
#' # Note, pgen returns log values for each locus, here we take the sum across
#' # all loci and take the exponent to give us the value of pgen for each sample
#' pGen <- exp(rowSums(pgen(xgid, by_pop = FALSE, freq = freqs)))
#' 
#' res  <- matrix(c(unique(pGen), unique(pSex)), ncol = 2)
#' colnames(res) <- c("Pgen", "Psex")
#' res <- cbind(xtab, nRamet = rownames(xtab), round(res, 5))
#' rownames(res) <- 1:3
#' res # Compare to the first three rows of Table 2 in Parks & Werth, 1993
#' }
#==============================================================================#
psex <- function(gid, pop = NULL, by_pop = TRUE, freq = NULL, G = NULL, 
                 method = c("single", "multiple")){
  stopifnot(is.genind(gid))
  if (!is.null(pop)){
    if (!is.language(pop)){
      pop(gid) <- pop
    } else {
      setPop(gid) <- pop
    }
    by_pop <- TRUE
  }
  if (!is.genclone(gid)){
    gid <- as.genclone(gid)
  }
  mll(gid) <- "original"
  METHOD <- c("single", "multiple")
  method <- match.arg(method, METHOD)
  xpgen  <- pgen(gid, by_pop = by_pop, freq = freq)
  xpgen  <- exp(rowSums(xpgen, na.rm = TRUE))

  if (method == "single"){
    # Only caculate for single encounter (Parks and Werth, 1993)
    G       <- ifelse(is.null(G), nmll(gid), G)
    pNotGen <- (1 - xpgen)^G
    return(1 - pNotGen)
  } else {
    # Calculate for n encounters (Arnaud-Haond et al, 1997)

    # Wed Sep  2 11:20:56 2015 ------------------------------
    # First step: get a vector of number of samples per mlg. Since pgen might
    # have been calculated using population-level frequencies, I'm not clone
    # correcting here, which results in a lot of redundant computations. I will
    # come up with a more clever way to handle this later.
    mlls       <- mll(gid, "original")
    mll_counts <- table(mlls)
    mll_counts <- mll_counts[match(as.character(mlls), names(mll_counts))]
    # Loop over each sample
    pSex <- vapply(seq(nInd(gid)), function(i){
      trials <- seq(mll_counts[i]) # number of samples in the MLG
      pgeni  <- xpgen[i]           # pgen for that MLG
      dens   <- dbinom(trials, mll_counts[i], pgeni) # vector of probabilites
      sum(dens, na.rm = TRUE) # return the sum probability.
    }, numeric(1))
    return(pSex)
  }
}