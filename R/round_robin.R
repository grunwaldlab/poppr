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
#' @author Zhian N. Kamvar, Jonah C. Brooks, Stacey Hadfield
#' 
#' @return a matrix of multilocus genotype assignments by masked locus. There 
#'   will be n rows and m columns where n = number of samples and m = number of
#'   loci.
#'
#' @export
#' @seealso \code{\link{rraf}}, \code{\link{pgen}}, \code{\link{psex}}
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
#' @param res Either "list" (default), "vector", or "data.frame".
#' @param by_pop When this is \code{TRUE}, the calculation will be done by
#'   population. Defaults to \code{FALSE}
#' @param correction a logical indicating whether or not zero value allele 
#'   frequencies should be set to 1/n (Default: \code{TRUE})
#'   
#' @return a vector or list of allele frequencies
#' @note When \code{by_pop = TRUE}, the output will be a matrix of allele 
#'   frequencies.
#' 
#' @author Zhian N. Kamvar, Jonah C. Brooks, Stacey Hadfield
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
#' rraf(Pram, type = "vector")
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
rraf <- function(gid, res = "list", by_pop = FALSE, correction = TRUE){
  RES     <- c("list", "vector", "data.frame")
  if ("loci" %in% class(gid)){
    gid <- pegas::loci2genind(gid)
  }
  res     <- match.arg(res, RES)
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
#' @param x a genind or genclone object. 
#'
#' @param log a \code{logical} if \code{log =TRUE} (default), the values
#'   returned will be log(Pgen). If \code{log = FALSE}, the values returned will
#'   be Pgen.
#'
#' @param by_pop a \code{logical} to determine whether allelic frequencies
#'   should be calculated per population (\code{TRUE}, default) or across all
#'   populations in the data (\code{FALSE}).
#'   
#' @param freq a vector or matrix of allele frequencies. This defaults to
#'   \code{NULL}, indicating that the frequencies will be determined via
#'   round-robin approach.
#'   
#' @note For haploids, Pgen at a particular locus is the allele frequency.
#'
#' @return A vector containing Pgen values per locus for each genotype in the 
#'   object.
#' @author Zhian N. Kamvar, Jonah Brooks, Stacy Hadfield
#' @seealso \code{\link{psex}}, \code{\link{rraf}}, \code{\link{rrmlg}}
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
pgen <- function(x, log = TRUE, by_pop = TRUE, freq = NULL){
  stopifnot(is.genind(x))
  # Stop if the ploidy of the object is not diploid
  stopifnot(all(ploidy(x) %in% 1:2)) 

  # Ensure there is a population assignment for every genotype
  if (is.null(pop(x)) || !by_pop){
    pop(x) <- rep(1, nInd(x))
  }   
  pops  <- pop(x)
  if (is.null(freq)){
    freqs <- rraf(x, by_pop = TRUE)
  } else if (is.matrix(freq)){
    if (nrow(freq) != nlevels(pops) || ncol(freq) != ncol(tab(x))){
      stop("frequency matrix must have the same dimensions as the data.")
    }
    freqs <- freq
  } else if (is.numeric(freq) && length(freq) == ncol(tab(x))){
    freqs <- freq
    pops <- factor(rep(1, nInd(x)))
  } else {
    stop("frequencies must be the same length as the number of alleles.")
  }
  
  pgen_matrix <- .Call("get_pgen_matrix_genind", x, freqs, pops, nlevels(pops), 
                       PACKAGE = "poppr")
  if (!log)
  {
    pgen_matrix <- exp(pgen_matrix)
  }
  dimnames(pgen_matrix) <- list(indNames(x), locNames(x))

  return(pgen_matrix);
}
#==============================================================================#
#' Probability of encountering a genotype more than once by chance
#' 
#' @param x a genind or genclone object.
#' @param by_pop a \code{logical} to determine whether allelic frequencies
#'   should be calculated per population (\code{TRUE}, default) or across all
#'   populations in the data (\code{FALSE}).
#' @param freq a vector or matrix of allele frequencies. This defaults to
#'   \code{NULL}, indicating that the frequencies will be determined via
#'   round-robin approach.
#' @param G an integer specifying the number of observed genets. If NULL, this
#'   will be the number of original multilocus genotypes.
#' @return a vector of Psex for each sample.
#' @note The values of Psex represent the value for each multilocus genotype.
#' 
#' @author Zhian N. Kamvar, Jonah Brooks, Stacy Hadfield
#' @seealso \code{\link{pgen}}, \code{\link{rraf}}, \code{\link{rrmlg}}
#' @export
#' @examples
#' 
#' data(Pram)
#' Pram_psex <- psex(Pram, by_pop = FALSE)
#' plot(Pram_psex, log = "y", col = ifelse(Pram_psex > 0.05, "red", "blue"))
#' abline(h = 0.05, lty = 2)
#' \dontrun{
#' # This can be also done assuming populations structure
#' Pram_psex <- psex(Pram, by_pop = TRUE)
#' plot(Pram_psex, log = "y", col = ifelse(Pram_psex > 0.05, "red", "blue"))
#' abline(h = 0.05, lty = 2)
#' 
#' ## An example of supplying previously calculated frequencies and G
#' # From Parks and Werth, 1993, using the first three genotypes.
#' x <- "
#'  Hk Lap Mdh2 Pgm1 Pgm2 X6Pgd2
#' 54 12 12 12 23 22 11
#' 36 22 22 11 22 33 11
#' 10 23 22 11 33 13 13"
#' 
#' xtab <- read.table(text = x, header = TRUE, row.names = 1)
#' xgid <- df2genind(xtab[rep(rownames(xtab), as.integer(rownames(xtab))), ], ncode = 1)
#' afreq <- c(Hk.1 = 0.167, Hk.2 = 0.795, Hk.3 = 0.038, 
#'            Lap.1 = 0.190, Lap.2 = 0.798, Lap.3 = 0.012,
#'            Mdh2.0 = 0.011, Mdh2.1 = 0.967, Mdh2.2 = 0.022,
#'            Pgm1.2 = 0.279, Pgm1.3 = 0.529, Pgm1.4 = 0.162, Pgm1.5 = 0.029,
#'            Pgm2.1 = 0.128, Pgm2.2 = 0.385, Pgm2.3 = 0.487,
#'            X6Pgd2.1 = 0.526, X6Pgd2.2 = 0.051, X6Pgd2.3 = 0.423)
#' freqs   <- afreq[colnames(tab(xgid))]
#' pNotGen <- psex(xgid, by_pop = FALSE, freq = freqs, G = 45)
#' pGen    <- exp(rowSums(pgen(xgid, by_pop = FALSE, freq = freqs)))
#' res     <- matrix(c(unique(pGen), unique(pNotGen)), ncol = 2)
#' colnames(res) <- c("Pgen", "Psex")
#' res
#' }
#==============================================================================#
psex <- function(x, by_pop = TRUE, freq = NULL, G = NULL){
  xpgen   <- pgen(x, by_pop = by_pop, freq = freq)
  G       <- ifelse(is.null(G), mlg(x, quiet = TRUE), G)
  xpgen   <- exp(rowSums(xpgen, na.rm = TRUE))
  pNotGen <- (1 - xpgen)^G
  return(1 - pNotGen)
}