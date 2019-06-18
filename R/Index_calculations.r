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
#' Produce a basic summary table for population genetic analyses.
#' 
#' @description
#' 
#' For the \pkg{poppr} package description, please see
#' \code{\link[=poppr-package]{package?poppr}}
#' 
#' This function allows the user to quickly view indices of heterozygosity, 
#' evenness, and linkage to aid in the decision of a path to further analyze 
#' a specified dataset. It natively takes \code{\linkS4class{genind}} and 
#' \code{\linkS4class{genclone}} objects, but can convert any raw data formats 
#' that adegenet can take (fstat, structure, genetix, and genpop) as well as 
#' genalex files exported into a csv format (see \code{\link{read.genalex}} for 
#' details).
#' 
#' 
#' @param dat a \code{\linkS4class{genind}} object OR a 
#'   \code{\linkS4class{genclone}} object OR any fstat, structure, genetix, 
#'   genpop, or genalex formatted file.
#'   
#' @param total When \code{TRUE} (default), indices will be calculated for the 
#'   pooled populations.
#'   
#' @param sublist a list of character strings or integers to indicate specific 
#'   population names (accessed via \code{\link[adegenet]{popNames}}). 
#'   Defaults to "ALL".
#'   
#' @param blacklist a list of character strings or integers to indicate specific
#'   populations to be removed from analysis. Defaults to NULL.
#'   
#' @param sample an integer indicating the number of permutations desired to 
#'   obtain p-values. Sampling will shuffle genotypes at each locus to simulate 
#'   a panmictic population using the observed genotypes. Calculating the 
#'   p-value includes the observed statistics, so set your sample number to one 
#'   off for a round p-value (eg. \code{sample = 999} will give you p = 0.001 
#'   and \code{sample = 1000} will give you p = 0.000999001).
#'   
#' @param method an integer from 1 to 4 indicating the method of sampling 
#'   desired. see \code{\link{shufflepop}} for details.
#'   
#' @param missing how should missing data be treated? \code{"zero"} and 
#'   \code{"mean"} will set the missing values to those documented in 
#'   \code{\link{tab}}. \code{"loci"} and \code{"geno"} will remove any loci or
#'   genotypes with missing data, respectively (see \code{\link{missingno}} for
#'   more information.
#'   
#' @param cutoff \code{numeric} a number from 0 to 1 indicating the percent 
#'   missing data allowed for analysis. This is to be used in conjunction with 
#'   the flag \code{missing} (see \code{\link{missingno}} for details)
#'   
#' @param quiet \code{FALSE} (default) will display a progress bar for each 
#'   population analyzed.
#'   
#' @param clonecorrect default \code{FALSE}. must be used with the \code{strata}
#'   parameter, or the user will potentially get undesired results. see
#'   \code{\link{clonecorrect}} for details.
#'   
#' @param strata a \code{formula} indicating the hierarchical levels to be used.
#'   The hierarchies should be present in the \code{strata} slot. See
#'   \code{\link{strata}} for details.
#'   
#' @param keep an \code{integer}. This indicates which strata you wish to keep 
#'   after clone correcting your data sets. To combine strata, just set keep 
#'   from 1 to the number of straifications set in strata. see 
#'   \code{\link{clonecorrect}} for details.
#'   
#' @param plot \code{logical} if \code{TRUE} (default) and \code{sampling > 0}, 
#'   a histogram will be produced for each population.
#'   
#' @param hist \code{logical} Deprecated. Use plot.
#'   
#' @param index \code{character} Either "Ia" or "rbarD". If \code{hist = TRUE}, 
#'   this will determine the index used for the visualization.
#'   
#' @param minsamp an \code{integer} indicating the minimum number of individuals
#'   to resample for rarefaction analysis. See \code{\link[vegan]{rarefy}} for 
#'   details.
#'   
#' @param legend \code{logical}. When this is set to \code{TRUE}, a legend 
#'   describing the resulting table columns will be printed. Defaults to 
#'   \code{FALSE}
#'   
#' @param ... arguments to be passed on to \code{\link{diversity_stats}}
#'   
#' @return A data frame with populations in rows and the following columns:
#'   \item{Pop}{A vector indicating the population factor} 
#'   \item{N}{An integer vector indicating the number of individuals/isolates in
#'   the specified population.}
#'   \item{MLG}{An integer vector indicating the number of multilocus genotypes 
#'   found in the specified population, (see: \code{\link{mlg}})}
#'   \item{eMLG}{The expected number of MLG at the lowest common sample size
#'   (set by the parameter \code{minsamp}).}
#'   \item{SE}{The standard error for the rarefaction analysis}
#'   \item{H}{Shannon-Weiner Diversity index}
#'   \item{G}{Stoddard and Taylor's Index} 
#'   \item{lambda}{Simpson's index} 
#'   \item{E.5}{Evenness} 
#'   \item{Hexp}{Nei's gene diversity (expected heterozygosity)}
#'   \item{Ia}{A numeric vector giving the value of the Index of Association for
#'   each population factor, (see \code{\link{ia}}).}
#'   \item{p.Ia}{A numeric vector indicating the p-value for Ia from the number
#'   of reshufflings indicated in \code{sample}. Lowest value is 1/n where n is
#'   the number of observed values.}
#'   \item{rbarD}{A numeric vector giving the value of the Standardized Index of
#'   Association for each population factor, (see \code{\link{ia}}).}
#'   \item{p.rD}{A numeric vector indicating the p-value for rbarD from the
#'   number of reshuffles indicated in \code{sample}. Lowest value is 1/n where
#'   n is the number of observed values.}
#'   \item{File}{A vector indicating the name of the original data file.}
#'   
#' @details This table is intended to be a first look into the dynamics of 
#'   mutlilocus genotype diversity. Many of the statistics (except for the the 
#'   index of association) are simply based on counts of multilocus genotypes 
#'   and do not take into account the actual allelic states.
#'   \strong{Descriptions of the statistics can be found in the Algorithms and
#'   Equations vignette}: \code{vignette("algo", package = "poppr")}.
#'   \subsection{sampling}{The sampling procedure is explicitly for testing the
#'   index of association. None of the other diversity statistics (H, G, lambda,
#'   E.5) are tested with this sampling due to the differing data types. To
#'   obtain confidence intervals for these statistics, please see 
#'   \code{\link{diversity_ci}}.}
#'   \subsection{rarefaction}{Rarefaction analysis is performed on the number of
#'   multilocus genotypes because it is relatively easy to estimate (Grünwald et
#'   al., 2003). To obtain rarefied estimates of diversity, it is possible to
#'   use \code{\link{diversity_ci}} with the argument \code{rarefy = TRUE}}
#'   \subsection{graphic}{This function outputs a \pkg{ggplot2} graphic of
#'   histograms. These can be manipulated to be visualized in another manner by
#'   retrieving the plot with the \code{\link{last_plot}} command from
#'   \pkg{ggplot2}. A useful manipulation would be to arrange the graphs into a
#'   single column so that the values of the statistic line up: \cr \code{p <-
#'   last_plot(); p + facet_wrap(~population, ncol = 1, scales = "free_y")}\cr
#'   The name for the groupings is "population" and the name for the x axis is
#'   "value".}
#'   
#' @note The calculation of \code{Hexp} has changed from \pkg{poppr} 1.x. It was
#'   previously calculated based on the diversity of multilocus genotypes, 
#'   resulting in a value of 1 for sexual populations. This was obviously not 
#'   Nei's 1978 expected heterozygosity. We have thus changed the statistic to 
#'   be the true value of Hexp by calculating \eqn{(\frac{n}{n-1}) 1 - \sum_{i =
#'   1}^k{p^{2}_{i}}}{(n/(n - 1))*(1 - sum(p^2))} where p is the allele
#'   frequencies at a given locus and n is the number of observed alleles (Nei,
#'   1978) in each locus and then returning the average. Caution should be 
#'   exercised in interpreting the results of Hexp with polyploid organisms with
#'   ambiguous ploidy. The lack of allelic dosage information will cause rare 
#'   alleles to be over-represented and artificially inflate the index. This is 
#'   especially true with small sample sizes.
#'   
#' @seealso \code{\link{clonecorrect}}, 
#'   \code{\link{poppr.all}}, 
#'   \code{\link{ia}}, 
#'   \code{\link{missingno}}, 
#'   \code{\link{mlg}}, 
#'   \code{\link{diversity_stats}},
#'   \code{\link{diversity_ci}}
#'   
#' @export
#' @author Zhian N. Kamvar
#' @references  Paul-Michael Agapow and Austin Burt. Indices of multilocus 
#'   linkage disequilibrium. \emph{Molecular Ecology Notes}, 1(1-2):101-102, 
#'   2001
#'   
#'   A.H.D. Brown, M.W. Feldman, and E. Nevo. Multilocus structure of natural 
#'   populations of \emph{Hordeum spontaneum}. \emph{Genetics}, 96(2):523-536,
#'   1980.
#'   
#'   Niklaus J. Gr\"unwald, Stephen B. Goodwin, Michael G. Milgroom, and William
#'   E. Fry. Analysis of genotypic diversity data for populations of 
#'   microorganisms. Phytopathology, 93(6):738-46, 2003
#'   
#'   Bernhard Haubold and Richard R. Hudson. Lian 3.0: detecting linkage 
#'   disequilibrium in multilocus data. Bioinformatics, 16(9):847-849, 2000.
#'   
#'   Kenneth L.Jr. Heck, Gerald van Belle, and Daniel Simberloff. Explicit 
#'   calculation of the rarefaction diversity measurement and the determination 
#'   of sufficient sample size. Ecology, 56(6):pp. 1459-1461, 1975
#'   
#'   Masatoshi Nei. Estimation of average heterozygosity and genetic distance 
#'   from a small number of individuals. Genetics, 89(3):583-590, 1978.
#'   
#'   S H Hurlbert. The nonconcept of species diversity: a critique and 
#'   alternative parameters. Ecology, 52(4):577-586, 1971.
#'   
#'   J.A. Ludwig and J.F. Reynolds. Statistical Ecology. A Primer on Methods and
#'   Computing. New York USA: John Wiley and Sons, 1988.
#'   
#'   Simpson, E. H. Measurement of diversity. Nature 163: 688, 1949 
#'   doi:10.1038/163688a0
#'   
#'   Good, I. J. (1953). On the Population Frequency of Species and the 
#'   Estimation of Population Parameters. \emph{Biometrika} 40(3/4): 237-264.
#'   
#'   Lande, R. (1996). Statistics and partitioning of species diversity, and 
#'   similarity among multiple communities. \emph{Oikos} 76: 5-13.
#'   
#'   Jari Oksanen, F. Guillaume Blanchet, Roeland Kindt, Pierre Legendre, Peter 
#'   R. Minchin, R. B. O'Hara, Gavin L. Simpson, Peter Solymos, M. Henry H. 
#'   Stevens, and Helene Wagner. vegan: Community Ecology Package, 2012. R 
#'   package version 2.0-5.
#'   
#'   E.C. Pielou. Ecological Diversity. Wiley, 1975.
#'   
#'   Claude Elwood Shannon. A mathematical theory of communication. Bell Systems
#'   Technical Journal, 27:379-423,623-656, 1948
#'   
#'   J M Smith, N H Smith, M O'Rourke, and B G Spratt. How clonal are bacteria? 
#'   Proceedings of the National Academy of Sciences, 90(10):4384-4388, 1993.
#'   
#'   J.A. Stoddart and J.F. Taylor. Genotypic diversity: estimation and 
#'   prediction in samples. Genetics, 118(4):705-11, 1988.
#'   
#'   
#' @examples
#' data(nancycats)
#' poppr(nancycats)
#' 
#' \dontrun{
#' # Sampling
#' poppr(nancycats, sample = 999, total = FALSE, plot = TRUE)
#' 
#' # Customizing the plot
#' library("ggplot2")
#' p <- last_plot()
#' p + facet_wrap(~population, scales = "free_y", ncol = 1)
#' 
#' # Turning off diversity statistics (see get_stats)
#' poppr(nancycats, total=FALSE, H = FALSE, G = FALSE, lambda = FALSE, E5 = FALSE)
#' 
#' # The previous version of poppr contained a definition of Hexp, which
#' # was calculated as (N/(N - 1))*lambda. It basically looks like an unbiased 
#' # Simpson's index. This statistic was originally included in poppr because it
#' # was originally included in the program multilocus. It was finally figured
#' # to be an unbiased Simpson's diversity metric (Lande, 1996; Good, 1953).
#' 
#' data(Aeut)
#' 
#' uSimp <- function(x){
#'   lambda <- vegan::diversity(x, "simpson")
#'   x <- drop(as.matrix(x))
#'   if (length(dim(x)) > 1){
#'     N <- rowSums(x)
#'   } else {
#'     N <- sum(x)
#'   }
#'   return((N/(N-1))*lambda)
#' }
#' poppr(Aeut, uSimp = uSimp)
#' 
#' 
#' # Demonstration with viral data
#' # Note: this is a larger data set that could take a couple of minutes to run
#' # on slower computers. 
#' data(H3N2)
#' strata(H3N2) <- data.frame(other(H3N2)$x)
#' setPop(H3N2) <- ~country
#' poppr(H3N2, total = FALSE, sublist=c("Austria", "China", "USA"), 
#' 				clonecorrect = TRUE, strata = ~country/year)
#' }
#==============================================================================#
#' @import adegenet ggplot2 vegan
poppr <- function(dat, total = TRUE, sublist = "ALL", blacklist = NULL, 
                  sample = 0, method = 1, missing = "ignore", cutoff = 0.05, 
                  quiet = FALSE, clonecorrect = FALSE, strata = 1, keep = 1, 
                  plot = TRUE, hist = TRUE, index = "rbarD", minsamp = 10, 
                  legend = FALSE, ...){

  if (inherits(dat, c("genlight", "snpclone"))){
    msg <- "The poppr function will not work with genlight or snpclone objects"
    msg <- paste0(msg, "\nIf you want to calculate genotypic diversity, use ",
                  "the function diversity_stats().")
    stop(msg)
  }
  quiet <- should_poppr_be_quiet(quiet)
  x <- process_file(dat, missing = missing, cutoff = cutoff, 
                    clonecorrect = clonecorrect, strata = strata,
                    keep = keep, quiet = TRUE)  
  # The namelist will contain information such as the filename and population
  # names so that they can easily be ported around.
  namelist <- NULL
  hist <- plot
  callpop  <- match.call()
  if (!is.na(grep("system.file", callpop)[1])){
    popsplt <- unlist(strsplit(dat, "/"))
    namelist$File <- popsplt[length(popsplt)]
  } else if (is.genind(dat)){
    namelist$File <- as.character(callpop[2])
  } else {
    namelist$File <- basename(x$X)
  }
  if (toupper(sublist[1]) == "TOTAL" & length(sublist) == 1){
    dat           <- x$GENIND
    pop(dat)      <- rep("Total", nInd(dat))
    poplist       <- NULL
    poplist$Total <- dat
  } else {
    dat <- popsub(x$GENIND, sublist = sublist, blacklist = blacklist)
    if (any(levels(pop(dat)) == "")) {
      levels(pop(dat))[levels(pop(dat)) == ""] <- "?"
      warning("missing population factor replaced with '?'")
    }
    pdrop   <- if (dat$type == "PA") FALSE else TRUE
    poplist <- if (is.null(pop(dat))) NULL else seppop(dat, drop = pdrop)
  }

  # Creating the genotype matrix for vegan's diversity analysis.
  pop.mat <- mlg.matrix(dat)
  if (total == TRUE & !is.null(poplist) & length(poplist) > 1){
    poplist$Total <- dat
    pop.mat       <- rbind(pop.mat, colSums(pop.mat))
  }
  sublist <- names(poplist)
  Iout    <- NULL
  total   <- toupper(total)
  missing <- toupper(missing)
  # For presence/absences markers, a different algorithm is applied. 
  if (legend) poppr_message()
  
  MLG.vec <- rowSums(ifelse(pop.mat > 0, 1, 0))
  N.vec   <- rowSums(pop.mat)
  datploid <- unique(ploidy(dat))
  Hexp_correction <- 1
  if (length(datploid) > 1 || any(datploid > 2)){
    datploid <- NULL
    Hexp_correction <- N.vec/(N.vec - 1)
  }
  divmat <- diversity_stats(pop.mat, ...)
  if (!is.matrix(divmat)){
    divmat <- matrix(divmat, nrow = 1, dimnames = list(NULL, names(divmat)))
  }
  
  if (!is.null(poplist)){
    # rarefaction giving the standard errors. This will use the minimum pop size
    # above a user-defined threshold.
    raremax <- ifelse(is.null(nrow(pop.mat)), sum(pop.mat), 
                      ifelse(min(rowSums(pop.mat)) > minsamp, 
                             min(rowSums(pop.mat)), minsamp))

    Hexp <- vapply(lapply(poplist, pegas::as.loci), FUN = get_hexp_from_loci, 
                   FUN.VALUE = numeric(1), ploidy = datploid, type = dat@type)

    Hexp   <- data.frame(Hexp = Hexp)
    N.rare <- suppressWarnings(vegan::rarefy(pop.mat, raremax, se = TRUE))
    IaList <- lapply(sublist, function(x){
      namelist <- list(file = namelist$File, population = x)
      .ia(poplist[[x]], 
          sample = sample, 
          method = method,
          quiet = quiet, 
          missing = missing, 
          hist = FALSE,
          namelist = namelist)
    })    
    names(IaList) <- sublist
    if (sample > 0){
      classtest <- summary(IaList)
      classless <- !classtest[, "Class"] %in% "ialist"
      if (any(classless)){
        no_class_pops <- paste(names(IaList[classless]), collapse = ", ")
        msg    <- paste0("values for ", no_class_pops, 
                         " could not be plotted.\n")
        IaList[classless] <- lapply(IaList[classless], function(x) list(index = x))
        warning(msg, call. = FALSE)
      }
      if (plot){
        try(print(poppr.plot(sample = IaList[!classless], file = namelist$File)))
      }
      IaList <- data.frame(t(vapply(IaList, "[[", numeric(4), "index")))
    } else {
      IaList <- t(as.data.frame(IaList))
    }
    Iout <- as.data.frame(
      list(
        Pop = sublist,
        N = N.vec,
        MLG = MLG.vec,
        eMLG = N.rare[1, ],
        SE = N.rare[2, ],
        divmat,
        Hexp,
        IaList,
        File = namelist$File
      )
    ) 
    rownames(Iout) <- NULL
  } else { 
    # rarefaction giving the standard errors. No population structure means that
    # the sample is equal to the number of individuals.
    N.rare <- rarefy(pop.mat, sum(pop.mat), se = TRUE)
    Hexp   <- get_hexp_from_loci(pegas::as.loci(dat), 
                                 ploidy = datploid, type = dat@type)
    Hexp   <- data.frame(Hexp = Hexp)
    IaList <-.ia(dat, 
                 sample = sample, 
                 method = method, 
                 quiet = quiet,
                 missing = missing, 
                 namelist = list(File = namelist$File, population = "Total"),
                 hist = plot
                )
    IaList <- if (sample > 0) IaList$index else IaList
    Iout <- as.data.frame(list(
      Pop = "Total",
      N = N.vec,
      MLG = MLG.vec,
      eMLG = N.rare[1, ],
      SE = N.rare[2, ],
      divmat,
      Hexp,
      as.data.frame(t(IaList)),
      File = namelist$File
    )) 
    rownames(Iout) <- NULL
  }
  class(Iout) <- c("popprtable", "data.frame")
  return(Iout) 
}

#==============================================================================#
#' Process a list of files with poppr
#'
#' poppr.all is a wrapper function that will loop through a list of files from
#' the working directory, execute \code{\link{poppr}}, and concatenate the
#' output into one data frame.
#'
#' @param filelist a list of files in the current working directory
#'
#' @param ... arguments passed on to poppr
#'
#' @return see \code{\link{poppr}}
#'
#' @seealso \code{\link{poppr}}, \code{\link{getfile}}
#' @export
#' @author Zhian N. Kamvar
#' @examples
#' \dontrun{
#' # Obtain a list of fstat files from a directory.
#' x <- getfile(multi=TRUE, pattern="^.+?dat$")
#'
#' # run the analysis on each file.
#' poppr.all(file.path(x$path, x$files))
#' }
#==============================================================================# 
poppr.all <- function(filelist, ...){
  result <- NULL
  for(a in seq(length(filelist))){
    cat(" \\    \n")
    input <- filelist[[a]]
    if (is.genind(input)){
      file <- names(filelist)[a]
      if (is.null(file)){
        file <- a
      }
      cat("  | Data: ")
    } else {
      file <- basename(input)
      cat("  | File: ")
    }
    cat(file, "\n /    \n")
    res      <- poppr(input, ...)
    res$File <- file
    result   <- rbind(result, res)
  }
  return(result)
}
#==============================================================================#
#' Index of Association
#' 
#' Calculate the Index of Association and Standardized Index of Association.
#' \itemize{
#'   \item \code{ia()} calculates the index of association over all loci in
#'   the data set.
#'   \item \code{pair.ia()} calculates the index of association in a pairwise
#'   manner among all loci.
#'   \item  \code{resample.ia()} calculates the index of association on a
#'   reduced data set multiple times to create a distribution, showing the
#'   variation of values observed at a given sample size (previously 
#'   \code{jack.ia}).
#' }
#' 
#' @param gid a \code{\link{genind}} or \code{\link{genclone}} object.
#'   
#' @param sample an integer indicating the number of permutations desired (eg 
#'   999).
#'   
#' @param method an integer from 1 to 4 indicating the sampling method desired. 
#'   see \code{\link{shufflepop}} for details.
#'   
#' @param quiet Should the function print anything to the screen while it is 
#'   performing calculations?
#'   
#'   \code{TRUE} prints nothing.
#'   
#'   \code{FALSE} (default) will print the population name and progress bar.
#'   
#' @param missing a character string. see \code{\link{missingno}} for details.
#'   
#' @param plot When \code{TRUE} (default), a heatmap of the values per locus
#'   pair will be plotted (for \code{pair.ia()}). When \code{sampling > 0}, 
#'   different things happen with \code{ia()} and \code{pair.ia()}. For 
#'   \code{ia()}, a histogram for the data set is plotted. For \code{pair.ia()},
#'   p-values are added as text on the heatmap. 
#'   
#' @param hist \code{logical} Deprecated. Use plot.
#' 
#' @param index \code{character} either "Ia" or "rbarD". If \code{hist = TRUE}, 
#'   this indicates which index you want represented in the plot (default:
#'   "rbarD").
#'   
#' @param valuereturn \code{logical} if \code{TRUE}, the index values from the 
#'   reshuffled data is returned. If \code{FALSE} (default), the index is 
#'   returned with associated p-values in a 4 element numeric vector.
#'   
#' @return 
#'   \subsection{for \code{pair.ia}}{
#'   A matrix with two columns and choose(nLoc(gid), 2) rows representing the
#'   values for Ia and rbarD per locus pair.
#'   }
#'   \subsection{If no sampling has occurred:}{ 
#'   A named number vector of length 2 giving the Index of Association, "Ia";
#'   and the Standardized Index of Association, "rbarD" 
#'   }
#'   \subsection{If there is sampling:}{ A a named number vector of length 4
#'   with the following values:
#'   \itemize{
#'     \item{Ia - }{numeric. The index of association.} 
#'     \item{p.Ia - }{A number indicating the p-value resulting from a
#'     one-sided permutation test based on the number of samples indicated in
#'     the original call.}
#'     \item{rbarD - }{numeric. The standardized index of association.}
#'     \item{p.rD - }{A factor indicating the p-value resulting from a
#'     one-sided permutation test based on the number of samples indicated in
#'     the original call.} 
#'   }
#'   }
#'   \subsection{If there is sampling and valureturn = TRUE}{ 
#'   A list with the following elements:
#'   \itemize{ 
#'     \item{index }{The above vector}
#'     \item{samples }{A data frame with s by 2 column data frame where s is the
#'     number of samples defined. The columns are for the values of Ia and
#'     rbarD, respectively.}
#'   }
#'   }
#'   
#' @note \code{jack.ia()} is deprecated as the name was misleading. Please use
#'   \code{resample.ia()}
#' @details The index of association was originally developed by A.H.D. Brown 
#'   analyzing population structure of wild barley (Brown, 1980). It has been widely 
#'   used as a tool to detect clonal reproduction within populations . 
#'   Populations whose members are undergoing sexual reproduction, whether it be
#'   selfing or out-crossing, will produce gametes via meiosis, and thus have a 
#'   chance to shuffle alleles in the next generation. Populations whose members
#'   are undergoing clonal reproduction, however, generally do so via mitosis. 
#'   This means that the most likely mechanism for a change in genotype is via 
#'   mutation. The rate of mutation varies from species to species, but it is 
#'   rarely sufficiently high to approximate a random shuffling of alleles. The 
#'   index of association is a calculation based on the ratio of the variance of
#'   the raw number of differences between individuals and the sum of those 
#'   variances over each locus . You can also think of it as the observed 
#'   variance over the expected variance. If they  are the same, then the index 
#'   is zero after subtracting one (from Maynard-Smith, 1993): \deqn{I_A = 
#'   \frac{V_O}{V_E}-1}{Ia = (Vo/Ve) - 1} Since the distance is more or less a binary 
#'   distance, any sort of marker can be used for this analysis. In the 
#'   calculation, phase is not considered, and any difference increases the 
#'   distance between two individuals. Remember that each column represents a 
#'   different allele and that each entry in the table represents the fraction 
#'   of the genotype made up by that allele at that locus. Notice also that the 
#'   sum of the rows all equal one. Poppr uses this to calculate distances by 
#'   simply taking the sum of the absolute values of the differences between 
#'   rows.
#'   
#'   The calculation for the distance between two individuals at a single locus 
#'   with \emph{a} allelic states and a ploidy of \emph{k} is as follows (except
#'   for Presence/Absence data): \deqn{ d = \displaystyle 
#'   \frac{k}{2}\sum_{i=1}^{a} \mid A_{i} - B_{i}\mid }{d(A,B) = (k/2)*sum(abs(Ai - Bi))} 
#'   To find the total number of differences 
#'   between two individuals over all loci, you just take \emph{d} over \emph{m}
#'   loci, a value we'll call \emph{D}:
#'   
#'   \deqn{D = \displaystyle \sum_{i=1}^{m} d_i }{D = sum(di)}
#'   
#'   These values are calculated over all possible combinations of individuals 
#'   in the data set, \eqn{{n \choose 2}}{choose(n, 2)} after which you end up 
#'   with \eqn{{n \choose 2}\cdot{}m}{choose(n, 2) * m} values of \emph{d} and 
#'   \eqn{{n \choose 2}}{choose(n, 2)} values of \emph{D}. Calculating the 
#'   observed variances is fairly straightforward (modified from Agapow and 
#'   Burt, 2001):
#'   
#'   \deqn{ V_O = \frac{\displaystyle \sum_{i=1}^{n \choose 2} D_{i}^2 - 
#'   \frac{(\displaystyle\sum_{i=1}^{n \choose 2} D_{i})^2}{{n \choose 2}}}{{n 
#'   \choose 2}}}{Vo = var(D)}
#'   
#'   Calculating the expected variance is the sum of each of the variances of 
#'   the individual loci. The calculation at a single locus, \emph{j} is the 
#'   same as the previous equation, substituting values of \emph{D} for 
#'   \emph{d}:
#'   
#'   \deqn{ var_j = \frac{\displaystyle \sum_{i=1}^{n \choose 2} d_{i}^2 - 
#'   \frac{(\displaystyle\sum_{i=1}^{n \choose 2} d_i)^2}{{n \choose 2}}}{{n 
#'   \choose 2}} }{Varj = var(dj)}
#'   
#'   The expected variance is then the sum of all the variances over all 
#'   \emph{m} loci:
#'   
#'   \deqn{ V_E = \displaystyle \sum_{j=1}^{m} var_j }{Ve = sum(var(dj))}
#'   
#'   Agapow and Burt showed that \eqn{I_A}{Ia} increases steadily with the
#'   number of loci, so they came up with an approximation that is widely used,
#'   \eqn{\bar r_d}{rbarD}. For the derivation, see the manual for
#'   \emph{multilocus}.
#'   
#'   \deqn{ \bar r_d = \frac{V_O - V_E} {2\displaystyle 
#'   \sum_{j=1}^{m}\displaystyle \sum_{k \neq j}^{m}\sqrt{var_j\cdot{}var_k}} 
#'   }{rbarD = (Vo - Ve)/(2*sum(sum(sqrt(var(dj)*var(dk))))}
#'   
#' @references  Paul-Michael Agapow and Austin Burt. Indices of multilocus 
#'   linkage disequilibrium. \emph{Molecular Ecology Notes}, 1(1-2):101-102, 
#'   2001
#'   
#'   A.H.D. Brown, M.W. Feldman, and E. Nevo. Multilocus structure of natural 
#'   populations of \emph{Hordeum spontaneum}. \emph{Genetics}, 96(2):523-536, 1980.
#'   
#'   J M Smith, N H Smith, M O'Rourke, and B G Spratt. How clonal are bacteria? 
#'   Proceedings of the National Academy of Sciences, 90(10):4384-4388, 1993.
#'   
#' @seealso \code{\link{poppr}}, \code{\link{missingno}}, 
#'   \code{\link{import2genind}}, \code{\link{read.genalex}}, 
#'   \code{\link{clonecorrect}}, \code{\link{win.ia}}, \code{\link{samp.ia}}
#'   
#' @export
#' @rdname ia
#' @author Zhian N. Kamvar
#' @examples
#' data(nancycats)
#' ia(nancycats)
#' 
#' # Pairwise over all loci:
#' data(partial_clone)
#' res <- pair.ia(partial_clone)
#' plot(res, low = "black", high = "green", index = "Ia")
#' 
#' # Resampling
#' data(Pinf)
#' resample.ia(Pinf, reps = 99)
#' 
#' \dontrun{
#' 
#' # Pairwise IA with p-values (this will take about a minute)
#' res <- pair.ia(partial_clone, sample = 999)
#' head(res)
#' 
#' # Plot the results of resampling rbarD. 
#' library("ggplot2")
#' Pinf.resamp <- resample.ia(Pinf, reps = 999)
#' ggplot(Pinf.resamp[2], aes(x = rbarD)) +
#'   geom_histogram() +
#'   geom_vline(xintercept = ia(Pinf)[2]) +
#'   geom_vline(xintercept = ia(clonecorrect(Pinf))[2], linetype = 2) +
#'   xlab(expression(bar(r)[d]))
#' 
#' # Get the indices back and plot the distributions.
#' nansamp <- ia(nancycats, sample = 999, valuereturn = TRUE)
#' 
#' plot(nansamp, index = "Ia")
#' plot(nansamp, index = "rbarD")
#' 
#' # You can also adjust the parameters for how large to display the text
#' # so that it's easier to export it for publication/presentations.
#' library("ggplot2")
#' plot(nansamp, labsize = 5, linesize = 2) +
#'   theme_bw() +                                      # adding a theme
#'   theme(text = element_text(size = rel(5))) +       # changing text size
#'   theme(plot.title = element_text(size = rel(4))) + # changing title size
#'   ggtitle("Index of Association of nancycats")      # adding a new title
#' 
#' # Get the index for each population.
#' lapply(seppop(nancycats), ia)
#' # With sampling
#' lapply(seppop(nancycats), ia, sample = 999)
#' 
#' # Plot pairwise ia for all populations in a grid with cowplot
#' # Set up the library and data
#' library("cowplot")
#' data(monpop)
#' splitStrata(monpop) <- ~Tree/Year/Symptom
#' setPop(monpop)      <- ~Tree
#' 
#' # Need to set up a list in which to store the plots.
#' plotlist        <- vector(mode = "list", length = nPop(monpop))
#' names(plotlist) <- popNames(monpop)
#' 
#' # Loop throgh the populations, calculate pairwise ia, plot, and then
#' # capture the plot in the list
#' for (i in popNames(monpop)){
#'   x <- pair.ia(monpop[pop = i], limits = c(-0.15, 1)) # subset, calculate, and plot
#'   plotlist[[i]] <- ggplot2::last_plot() # save the last plot
#' }
#' 
#' # Use the plot_grid function to plot.
#' plot_grid(plotlist = plotlist, labels = paste("Tree", popNames(monpop)))
#' 
#' }
#==============================================================================#
ia <- function(gid, sample = 0, method = 1, quiet = FALSE, missing = "ignore", 
               plot = TRUE, hist = TRUE, index = "rbarD", valuereturn = FALSE){
  namelist <- list(population = ifelse(nPop(gid) > 1 | is.null(gid@pop), 
                                       "Total", popNames(gid)),
                   File = as.character(match.call()[2])
                  )
  hist    <- plot
  popx    <- gid
  missing <- toupper(missing)
  type    <- gid@type
  quiet   <- should_poppr_be_quiet(quiet)
  if (type == "PA"){
    .Ia.Rd <- .PA.Ia.Rd
  } else {
    popx <- seploc(popx)
  }

  # if there are less than three individuals in the population, the calculation
  # does not proceed. 
  if (nInd(gid) < 3){
    IarD <- stats::setNames(as.numeric(c(NA, NA)), c("Ia", "rbarD"))
    if (sample == 0){
      return(IarD)
    } else {
      IarD <- stats::setNames(as.numeric(rep(NA, 4)), c("Ia","p.Ia","rbarD","p.rD"))
      return(IarD)
    }
  }
  
  IarD <- .Ia.Rd(popx, missing)
  names(IarD) <- c("Ia", "rbarD")
  # no sampling, it will simply return two named numbers.
  if (sample == 0){
    Iout   <- IarD
    result <- NULL
  } else {
  # sampling will perform the iterations and then return a data frame indicating
  # the population, index, observed value, and p-value. It will also produce a 
  # histogram.
    Iout     <- NULL 
    # idx      <- data.frame(Index = names(IarD))
    samp     <- .sampling(popx, sample, missing, quiet = quiet, type = type, 
                          method = method)
    p.val    <- sum(IarD[1] <= c(samp$Ia, IarD[1]))/(sample + 1)
    p.val[2] <- sum(IarD[2] <= c(samp$rbarD, IarD[2]))/(sample + 1)

    if (hist == TRUE){
      the_plot <- poppr.plot(samp, observed = IarD, pop = namelist$population, 
                             index = index, file = namelist$File, pval = p.val, 
                             N = nrow(gid@tab))
      print(the_plot)
    }
    result <- stats::setNames(vector(mode = "numeric", length = 4), 
                       c("Ia","p.Ia","rbarD","p.rD"))
    result[c(1, 3)] <- IarD
    result[c(2, 4)] <- p.val
    if (valuereturn == TRUE){
      iaobj        <- list(index = final(Iout, result), samples = samp)
      class(iaobj) <- "ialist"
      return(iaobj)
    } 
  }  
  return(final(Iout, result))
}

#==============================================================================#
#' @rdname ia
#' @param low (for pair.ia) a color to use for low values when \code{plot =
#'   TRUE}
#' @param high (for pair.ia) a color to use for low values when \code{plot =
#'   TRUE}
#' @param limits (for pair.ia) the limits to be used for the color scale. 
#'   Defaults to \code{NULL}. If you want to use a custom range, supply two
#'   numbers between -1 and 1, (e.g. \code{limits = c(-0.15, 1)})
#' @export
#==============================================================================#
pair.ia <- function(gid, sample = 0L, quiet = FALSE, plot = TRUE, low = "blue", 
                    high = "red", limits = NULL, index = "rbarD", method = 1L){
  N       <- nInd(gid)
  numLoci <- nLoc(gid)
  lnames  <- locNames(gid)
  np      <- choose(N, 2)
  nploci  <- choose(numLoci, 2)
  shuffle <- sample > 0L
  quiet   <- should_poppr_be_quiet(quiet)
  QUIET   <- if (shuffle) TRUE else quiet
  res     <- pair_ia_internal(gid, N, numLoci, lnames, np, nploci, QUIET, sample)
  if (shuffle) {
    # Initialize with 1 to account for the observed data.
    counts <- matrix(1L, nrow = nrow(res), ncol = ncol(res))
    if (!quiet) prog <- dplyr::progress_estimated(sample)
    for (i in seq_len(sample)) {
      tmp    <- shufflepop(gid, method = method)
      tmpres <- pair_ia_internal(tmp, N, numLoci, lnames, np, nploci, QUIET, sample)
      counts <- counts + as.integer(tmpres >= res)
      if (!quiet) print(prog$tick())
    }
    if (!quiet) print(prog$stop())
    if (!quiet) cat("\n")
    p   <- counts/(sample + 1)
    res <- cbind(Ia = res[, 1], 
                 p.Ia = p[, 1], 
                 rbarD = res[, 2], 
                 p.rD = p[, 2])
  }
  class(res) <- c("pairia", "matrix")
  if (plot) {
    tryCatch(plot(res, index = index, low = low, high = high, limits = limits),
             error = function(e) e)
  }
  res
}

pair_ia_internal <- function(gid, N, numLoci, lnames, np, nploci, quiet, sample) {
  if (gid@type == "codom") {
    V <- pair_matrix(seploc(gid), numLoci, np)
  } else { # P/A case
    V <- apply(tab(gid), 2, function(x) as.vector(dist(x)))
    # checking for missing data and imputing the comparison to zero.
    if (any(is.na(V))) {
      V[which(is.na(V))] <- 0
    }
  }
  colnames(V) <- lnames
  
  loci_pairs       <- matrix(NA, nrow = 3, ncol = nploci)
  loci_pairs[-3, ] <- combn(lnames, 2)
  loci_pairs[3, ]  <- as.character(1:nploci)
  prog             <- NULL
  if (!quiet) prog <- txtProgressBar(style = 3)
  pair_ia_vector   <- apply(loci_pairs, 2, ia_pair_loc, V, np, prog, nploci)
  if (!quiet) cat("\n")
  colnames(pair_ia_vector) <- apply(loci_pairs[-3, ], 2, paste, collapse = ":")
  rownames(pair_ia_vector) <- c("Ia", "rbarD")
  pair_ia_vector           <- t(pair_ia_vector)
  pair_ia_vector
}
#==============================================================================#
#' Create a table of summary statistics per locus. 
#' 
#' @param x a [genind-class] or [genclone-class]
#'   object.
#' 
#' @param index Which diversity index to use. Choices are \itemize{ \item
#'   `"simpson"` (Default) to give Simpson's index \item `"shannon"`
#'   to give the Shannon-Wiener index \item `"invsimpson"` to give the
#'   Inverse Simpson's index aka the Stoddard and Tayor index.}
#'   
#' @param lev At what level do you want to analyze diversity? Choices are
#'   `"allele"` (Default) or `"genotype"`.
#'   
#' @param population Select the populations to be analyzed. This is the
#'   parameter `sublist` passed on to the function [popsub()].
#'   Defaults to `"ALL"`.
#'   
#' @param information When `TRUE` (Default), this will print out a header
#'   of information to the R console.
#'   
#' @return a table with 4 columns indicating the Number of alleles/genotypes 
#'   observed, Diversity index chosen, Nei's 1978 gene diversity (expected
#'   heterozygosity), and Evenness.
#'   
#' @seealso [vegan::diversity()], [poppr()]
#' @md
#'   
#' @note The calculation of `Hexp` is \eqn{(\frac{n}{n-1}) 1 - \sum_{i = 
#'   1}^k{p^{2}_{i}}}{(n/(n - 1))*(1 - sum(p^2))} where p is the allele 
#'   frequencies at a given locus and n is the number of observed alleles (Nei, 
#'   1978) in each locus and then returning the average. Caution should be 
#'   exercised in interpreting the results of Hexp with polyploid organisms with
#'   ambiguous ploidy. The lack of allelic dosage information will cause rare 
#'   alleles to be over-represented and artificially inflate the index. This is 
#'   especially true with small sample sizes.
#'   
#'   If `lev = "genotype"`, then all statistics reflect **genotypic** diversity
#'   within each locus. This includes the calculation for `Hexp`, which turns
#'   into the unbiased Simpson's index.
#'   
#' @author Zhian N. Kamvar
#' 
#' @references
#'   Jari Oksanen, F. Guillaume Blanchet, Roeland Kindt, Pierre Legendre, Peter 
#'   R. Minchin, R. B. O'Hara, Gavin L. Simpson, Peter Solymos, M. Henry H. 
#'   Stevens, and Helene Wagner. vegan: Community Ecology Package, 2012. R 
#'   package version 2.0-5.
#' 
#'   Niklaus J. Gr\"unwald, Stephen B. Goodwin, Michael G. Milgroom, and William
#'   E. Fry. Analysis of genotypic diversity data for populations of 
#'   microorganisms. Phytopathology, 93(6):738-46, 2003
#'   
#'   J.A. Ludwig and J.F. Reynolds. Statistical Ecology. A Primer on Methods and
#'   Computing. New York USA: John Wiley and Sons, 1988.
#'   
#'   E.C. Pielou. Ecological Diversity. Wiley, 1975.
#'   
#'   J.A. Stoddart and J.F. Taylor. Genotypic diversity: estimation and
#'   prediction in samples. Genetics, 118(4):705-11, 1988.
#'   
#'   Masatoshi Nei. Estimation of average heterozygosity and genetic distance 
#'   from a small number of individuals. Genetics, 89(3):583-590, 1978.
#' 
#'   Claude Elwood Shannon. A mathematical theory of communication. Bell Systems
#'   Technical Journal, 27:379-423,623-656, 1948
#'   
#' @export
#' @examples
#' 
#' data(nancycats)
#' locus_table(nancycats[pop = 5])
#' \dontrun{
#' # Analyze locus statistics for the North American population of P. infestans.
#' # Note that due to the unknown dosage of alleles, many of these statistics
#' # will be artificially inflated for polyploids.
#' data(Pinf)
#' locus_table(Pinf, population = "North America")
#' }
#==============================================================================#
locus_table <- function(x, index = "simpson", lev = "allele", 
                        population = "ALL", information = TRUE){
  ploid   <- unique(ploidy(x))
  type    <- x@type
  INDICES <- c("shannon", "simpson", "invsimpson")
  index   <- match.arg(index, INDICES)
  x       <- popsub(x, population, drop = FALSE)
  x.loc   <- summary(as.loci(x))
  outmat  <- vapply(x.loc, locus_table_pegas, numeric(4), index, lev, ploid, type)
  loci    <- colnames(outmat)
  divs    <- rownames(outmat)
  res     <- matrix(0.0, nrow = ncol(outmat) + 1, ncol = nrow(outmat))
  dimlist <- list(`locus` = c(loci, "mean"), `summary` = divs)
  res[-nrow(res), ]     <- t(outmat)
  res[nrow(res), ]      <- colMeans(res[-nrow(res), ], na.rm = TRUE)
  attr(res, "dimnames") <- dimlist
  if (information){
    if (index == "simpson"){
      msg <- "Simpson index"
    } else if (index == "shannon"){
      msg <- "Shannon-Wiener index"
    } else {
      msg <- "Stoddard and Taylor index"
    }
    message("\n", divs[1], " = Number of observed ", paste0(divs[1], "s"), appendLF = FALSE)
    message("\n", divs[2], " = ", msg, appendLF = FALSE)
    message("\n", divs[3], " = Nei's 1978 gene diversity\n", appendLF = FALSE)
    message("------------------------------------------\n", appendLF = FALSE)
  }
  class(res) <- c("locustable", "matrix")
  return(res)
}

#==============================================================================#
#' Tabulate alleles the occur in only one population. 
#' 
#' @param gid a [genind-class] or [genclone-class]
#'   object.
#'   
#' @param form a [formula()] giving the levels of markers and 
#'   hierarchy to analyze. See Details.
#'   
#' @param report one of `"table", "vector",` or `"data.frame"`. Tables
#'   (Default) and data frame will report counts along with populations or 
#'   individuals. Vectors will simply report which populations or individuals 
#'   contain private alleles. Tables are matrices with populations or 
#'   individuals in rows and alleles in columns. Data frames are long form.
#'   
#' @param level one of `"population"` (Default) or `"individual"`.
#'   
#' @param count.alleles `logical`. If `TRUE` (Default), The report 
#'   will return the observed number of alleles private to each population. If 
#'   `FALSE`, each private allele will be counted once, regardless of 
#'   dosage.
#' 
#' @param drop `logical`. if `TRUE`, populations/individuals without 
#'   private alleles will be dropped from the result. Defaults to `FALSE`.
#'   
#' @return a matrix, data.frame, or vector defining the populations or
#'   individuals containing private alleles. If vector is chosen, alleles are
#'   not defined.
#'
#' @details the argument `form` allows for control over the strata at which
#'   private alleles should be computed. It takes a form where the left hand
#'   side of the formula can be either "allele", "locus", or "loci". The right
#'   hand of the equation, by default is ".". If you change it, it must
#'   correspond to strata located in the [adegenet::strata()] slot.
#'   Note, that the right hand side is disabled for genpop objects.
#' 
#' @export
#' @author Zhian N. Kamvar
#' @md
#' @examples
#' 
#' data(Pinf) # Load P. infestans data.
#' private_alleles(Pinf)
#' 
#' \dontrun{
#' # Analyze private alleles based on the country of interest:
#' private_alleles(Pinf, alleles ~ Country)
#' 
#' # Number of observed alleles per locus
#' private_alleles(Pinf, locus ~ Country, count.alleles = TRUE)
#' 
#' # Get raw number of private alleles per locus.
#' (pal <- private_alleles(Pinf, locus ~ Country, count.alleles = FALSE))
#' 
#' # Get percentages.
#' sweep(pal, 2, nAll(Pinf)[colnames(pal)], FUN = "/")
#' 
#' # An example of how these data can be displayed.
#' library("ggplot2")
#' Pinfpriv <- private_alleles(Pinf, report = "data.frame")
#' ggplot(Pinfpriv) + geom_tile(aes(x = population, y = allele, fill = count))
#' }
#==============================================================================#
private_alleles <- function(gid, form = alleles ~ ., report = "table", 
                            level = "population", count.alleles = TRUE,
                            drop = FALSE){
  REPORTARGS <- c("table", "vector", "data.frame")
  LEVELARGS  <- c("individual", "population")
  LHS_ARGS <- c("alleles", "locus", "loci")
  showform <- utils::capture.output(print(form))
  marker <- pmatch(as.character(form[[2]]), LHS_ARGS, nomatch = 0L, 
                   duplicates.ok = FALSE)
  if (all(marker == 0L)){
    stop("Left hand side of ", showform, " must be one of:\n ",
         paste(LHS_ARGS, collapse = " "))
  } else {
    marker <- LHS_ARGS[marker]
  }
  strataform <- form[c(1, 3)]
  the_strata <- all.vars(strataform[[2]])
  if (length(the_strata) > 1 || the_strata[1] != "."){
    if (!is.genpop(gid)){
      setPop(gid) <- strataform
    } else {
      warning("cannot set strata for a genpop object.")
    }
  } 
  report <- match.arg(report, REPORTARGS)
  level  <- match.arg(level, LEVELARGS)
  if (!is.genind(gid) & !is.genpop(gid)){
    stop(paste(gid, "is not a genind or genpop object."))
  }
  if (is.genind(gid) && !is.null(pop(gid)) | is.genpop(gid) && nPop(gid) > 1){
    if (is.genind(gid)){
      gid.pop <- tab(genind2genpop(gid, quiet = TRUE))
    } else {
      gid.pop <- tab(gid)
    }
    private_columns <- colSums(ifelse(gid.pop > 0, 1, 0), na.rm = TRUE) < 2
    privates <- gid.pop[, private_columns, drop = FALSE]
    if (level == "individual" & is.genind(gid)){
      gid.tab  <- tab(gid)
      privates <- gid.tab[, private_columns, drop = FALSE]
    } else if (!count.alleles){
      privates <- ifelse(privates > 0, 1, 0)
    }
    
    if (drop){
      privates <- privates[rowSums(privates, na.rm = TRUE) > 0, , drop = FALSE]
    }
    
    if (marker != "alleles"){
      private_fac <- locFac(gid)[private_columns]
      privates <- vapply(unique(private_fac), function(l){
        rowSums(privates[, private_fac == l, drop = FALSE], na.rm = TRUE)
      }, FUN.VALUE = numeric(nrow(privates))
      )
      colnames(privates) <- locNames(gid)[unique(private_fac)]
    }
    if (length(privates) == 0){
      privates <- NULL
      cat("No private alleles detected.")
      return(invisible(NULL))
    }
    if (report == "vector"){
      privates <- rownames(privates)
    } else if (report == "data.frame"){
      marker   <- if (marker == "alleles") "allele" else "locus"
      names(dimnames(privates)) <- c(level, marker)
      privates <- as.data.frame.table(privates, 
                                      responseName = "count",
                                      stringsAsFactors = FALSE)
    }
    return(privates)
  } else {
    stop("There are no populations detected")
  }
}
