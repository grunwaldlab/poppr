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
#' Produce a basic summary table for population genetic analyses.
#' 
#' This function allows the user to quickly view indicies of heterozygosity, 
#' evenness, and inbreeding to aid in the decision of a path to further analyze 
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
#'   population names (located in \code{$pop.names} within the 
#'   \code{\link{genind}} object) Defaults to "ALL".
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
#'   \code{\link{na.replace}}. \code{"loci"} and \code{"geno"} will remove any 
#'   loci or genotypes with missing data, respectively (see 
#'   \code{\link{missingno}} for more information.
#'   
#' @param cutoff \code{numeric} a number from 0 to 1 indicating the percent 
#'   missing data allowed for analysis. This is to be used in conjunction with 
#'   the flag \code{missing} (see \code{\link{missingno}} for details)
#'   
#' @param quiet \code{FALSE} (defualt) will display a progress bar for each
#'   population analyzed.
#'   
#' @param clonecorrect default \code{FALSE}. must be used with the \code{hier} 
#'   and \code{dfname} parameters, or the user will potentially get undesired 
#'   results. see \code{\link{clonecorrect}} for details.
#'   
#' @param hier \itemize{ \item \strong{for genclone objects} - a \code{formula} 
#'   indicating the hierarchical levels to be used. The hierarchies should be 
#'   present in the \code{hierarchy} slot. See \code{\link{sethierarchy}} for 
#'   details. \item \strong{for genind objects} - a \code{numeric or character} 
#'   vector OR a hierarchical formula. This is the list of columns within a data
#'   frame (specified in \code{dfname}) in the 'other' slot of the 
#'   \code{\link{genind}} object. The list should indicate the population 
#'   hierarchy to be used for clone correction.  }
#'   
#' @param dfname a \code{character string}. (Only for genind objects) This is
#'   the name of the data frame or heirarchy containing the vectors of the
#'   population hierarchy within the \code{other} slot of the
#'   \code{\link{genind}} object.
#'   
#' @param keep an \code{integer}. This indicates the levels of the population 
#'   hierarchy you wish to keep after clone correcting your data sets. To 
#'   combine the hierarchy, just set keep from 1 to the length of your 
#'   hierarchy. see \code{\link{clonecorrect}} for details.
#'   
#' @param hist \code{logical} if \code{TRUE} (default) and \code{sampling > 0}, 
#'   a histogram will be produced for each population.
#'   
#' @param minsamp an \code{integer} indicating the minimum number of individuals
#'   to resample for rarefaction analysis. See \code{\link[vegan]{rarefy}} for
#'   details.
#'   
#' @param legend \code{logical}. When this is set to \code{TRUE}, a legend 
#'   describing the resulting table columns will be printed. Defaults to 
#'   \code{FALSE}
#'   
#' @return \item{Pop}{A vector indicating the pouplation factor} \item{N}{An 
#'   integer vector indicating the number of individuals/isolates in the 
#'   specified population.} \item{MLG}{An integer vector indicating the number 
#'   of multilocus genotypes found in the specified poupulation, (see: 
#'   \code{\link{mlg}})} \item{eMLG}{The expected number of MLG at the lowest 
#'   common sample size (set by the parameter \code{minsamp}.} \item{SE}{The 
#'   standard error for the rarefaction analysis} \item{H}{Shannon-Weiner 
#'   Diversity index} \item{G}{Stoddard and Taylor's Index} \item{Hexp}{Expected
#'   heterozygosity or Nei's 1987 genotypic diversity corrected for sample 
#'   size.} \item{E.5}{Evenness} \item{Ia}{A numeric vector giving the value of 
#'   the Index of Association for each population factor, (see 
#'   \code{\link{ia}}).} \item{p.Ia}{A numeric vector indicating the p-value for
#'   Ia from the number of reshufflings indicated in \code{sample}. Lowest value
#'   is 1/n where n is the number of observed values.} \item{rbarD}{A numeric 
#'   vector giving the value of the Standardized Index of Association for each 
#'   population factor, (see \code{\link{ia}}).} \item{p.rD}{A numeric vector 
#'   indicating the p-value for rbarD from the number of reshufflings indicated 
#'   in \code{sample}. Lowest value is 1/n where n is the number of observed 
#'   values.} \item{File}{A vector indicating the name of the original data 
#'   file.}
#'   
#' @seealso \code{\link{clonecorrect}}, \code{\link{poppr.all}}, 
#'   \code{\link{ia}}, \code{\link{missingno}}, \code{\link{mlg}}
#'   
#' @export
#' @author Zhian N. Kamvar
#' @references  Paul-Michael Agapow and Austin Burt. Indices of multilocus 
#'   linkage disequilibrium. \emph{Molecular Ecology Notes}, 1(1-2):101-102, 
#'   2001
#'   
#'   A.H.D. Brown, M.W. Feldman, and E. Nevo. Multilocus structure of natural 
#'   populations of hordeum spontaneum. \emph{Genetics}, 96(2):523-536, 1980.
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
#'   S H Hurlbert. The nonconcept of species diversity: a critique and 
#'   alternative parameters. Ecology, 52(4):577-586, 1971.
#'   
#'   J.A. Ludwig and J.F. Reynolds. Statistical Ecology. A Primer on Methods and
#'   Computing. New York USA: John Wiley and Sons, 1988.
#'   
#'   Masatoshi Nei. Estimation of average heterozygosity and genetic distance 
#'   from a small number of individuals. Genetics, 89(3):583-590, 1978.
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
#' poppr(nancycats, sample=99, total=FALSE, quiet=FALSE)
#' 
#' # Note: this is a larger data set that could take a couple of minutes to run
#' # on slower computers. 
#' data(H3N2)
#' poppr(H3N2, total=FALSE, sublist=c("Austria", "China", "USA"), 
#' 				clonecorrect=TRUE, hier="country", dfname="x")
#' }
#==============================================================================#
#' @import adegenet ggplot2 vegan
poppr <- function(dat, total=TRUE, sublist="ALL", blacklist=NULL, sample=0,
                  method=1, missing="ignore", cutoff=0.05, quiet=FALSE,
                  clonecorrect=FALSE, hier=1, dfname="population_hierarchy", 
                  keep = 1, hist=TRUE, minsamp=10, legend=FALSE){
  METHODS = c("permute alleles", "parametric bootstrap",
              "non-parametric bootstrap", "multilocus")
  x <- process_file(dat, missing = missing, cutoff = cutoff, 
                  clonecorrect = clonecorrect, hier = hier, dfname = dfname, 
                  keep = keep, quiet = TRUE)  
  # The namelist will contain information such as the filename and population
  # names so that they can easily be ported around.
  namelist <- NULL
  callpop  <- match.call()
  if (!is.na(grep("system.file", callpop)[1])){
    popsplt <- unlist(strsplit(dat, "/"))
    namelist$File <- popsplt[length(popsplt)]
  } else if (is.genind(dat)){
    namelist$File <- as.character(callpop[2])
  } else {
    namelist$File <- basename(x$X)
  }
  #poplist <- x$POPLIST
  if(toupper(sublist[1]) == "TOTAL" & length(sublist) == 1){
    dat           <- x$GENIND
    pop(dat)      <- NULL
    poplist       <- NULL
    poplist$Total <- dat
  }
  else{
    dat <- popsub(x$GENIND, sublist=sublist, blacklist=blacklist)
    if (any(levels(pop(dat)) == "")){
      levels(pop(dat))[levels(pop(dat)) == ""] <- "?"
      warning("missing population factor replaced with '?'")
    }
    poplist <- .pop.divide(dat)
  }
  # Creating the genotype matrix for vegan's diversity analysis.
  pop.mat <- mlg.matrix(dat)
  if (total==TRUE & !is.null(poplist) & length(poplist) > 1){
    poplist$Total <- dat
    pop.mat <- rbind(pop.mat, colSums(pop.mat))
  }
  sublist <- names(poplist)
  Iout    <- NULL
  origpop <- x$GENIND
  total   <- toupper(total)
  missing <- toupper(missing)
  type    <- dat@type
  # For presence/absences markers, a different algorithm is applied. 
  if(type=="PA"){
    .Ia.Rd <- .PA.Ia.Rd
  }
  if (legend) poppr_message()
  
  MLG.vec <- rowSums(ifelse(pop.mat > 0, 1, 0))
  N.vec   <- rowSums(pop.mat)
  # Shannon-Weiner diversity index.
  H       <- vegan::diversity(pop.mat)
  # inverse Simpson's index aka Stoddard and Taylor: 1/lambda
  G       <- vegan::diversity(pop.mat, "inv")
  Hexp    <- (N.vec/(N.vec-1))*vegan::diversity(pop.mat, "simp")
  # E_5
  E.5     <- (G-1)/(exp(H)-1)
  
  if (!is.null(poplist)){
    # rarefaction giving the standard errors. This will use the minimum pop size
    # above a user-defined threshold.
    raremax <- ifelse(is.null(nrow(pop.mat)), sum(pop.mat), 
                      ifelse(min(rowSums(pop.mat)) > minsamp, 
                             min(rowSums(pop.mat)), minsamp))
    N.rare  <- suppressWarnings(rarefy(pop.mat, raremax, se=TRUE))
    IaList  <- NULL
    invisible(lapply(sublist, function(x){ 
      namelist <- list(File = namelist$File, population = x)
      IaList   <<- rbind(IaList, 
                       .ia(poplist[[x]], sample=sample, method=method, 
                        quiet=quiet, missing=missing, namelist=namelist, 
                        hist=hist) 
                       )
      })) 
    Iout <- as.data.frame(list(Pop=sublist, N=N.vec, MLG=MLG.vec, 
                               eMLG=N.rare[1, ], SE=N.rare[2, ], H=H, 
                               G=G, Hexp=Hexp, E.5=E.5, IaList, 
                               File=namelist$File)) 
    rownames(Iout) <- NULL
  } else { 
    # rarefaction giving the standard errors. No population structure means that
    # the sample is equal to the number of individuals.
    N.rare <- rarefy(pop.mat, sum(pop.mat), se=TRUE)
    IaList <- .ia(dat, sample=sample, method=method, quiet=quiet, missing=missing,
                  namelist=(list(File=namelist$File, population="Total")),
                  hist=hist)
    Iout <- as.data.frame(list(Pop="Total", N=N.vec, MLG=MLG.vec, 
                               eMLG=N.rare[1, ], SE=N.rare[2, ], H=H, G=G, 
                               Hexp=Hexp, E.5=E.5, as.data.frame(t(IaList)), 
                               File=namelist$File)) 
    rownames(Iout) <- NULL
  }
  class(Iout) <- c("popprtable", "data.frame")
  return(Iout) 
}

#==============================================================================#
#' Process a list of files with poppr
#'
#' poppr.all is a wrapper function that will loop through a list of files from
#' the workind directory, execute \code{\link{poppr}}, and concatenate the
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
#' # set the working directory to that directory.
#' setwd(x$path)
#'
#' # run the analysis on each file.
#' poppr.all(x$files)
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
#' Obtain p-values from one-sided permutation tests.
#' 
#' @param pop a \code{\link{genind}} object OR any fstat, structure, genetix, 
#'   genpop, or genalex formatted files.
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
#'   \code{FALSE} (defualt) will print the population name and progress bar.
#'   
#' @param missing a character string. see \code{\link{missingno}} for details.
#'   
#' @param hist \code{logical} if \code{TRUE}, a histogram will be printed for 
#'   each population if there is sampling.
#'   
#' @param valuereturn \code{logical} if \code{TRUE}, the index values from the 
#'   reshuffled data is returned. If \code{FALSE} (default), the index is 
#'   returned with associated p-values in a 4 element numeric vector.
#'   
#' @return \subsection{If no sampling has occured:}{ A named number vector of 
#'   length 2 giving the Index of Association, "Ia"; and the Standardized Index 
#'   of Association, "rbarD" } \subsection{If there is sampling:}{ A a named 
#'   number vector of length 4 with the following values: \itemize{\item{Ia -
#'   }{numeric. The index of association.} \item{p.Ia - }{A number indicating
#'   the p-value resulting from a one-sided permutation test based on the number
#'   of samples indicated in the original call.} \item{rbarD - }{numeric. The
#'   standardized index of association.} \item{p.rD - }{A factor indicating the
#'   p-value resutling from a one-sided permutation test based on the number of
#'   samples indicated in the original call.}} } \subsection{If there is
#'   sampling and valureturn = TRUE}{ A list with the following
#'   elements:\itemize{ \item{index}{The above vector} \item{samples}{A data
#'   frame with s by 2 column data frame where s is the number of samples
#'   defined. The columns are for the values of Ia and rbarD, respectively.}}}
#'   
#' @details The index of association was originally developed by A.H.D. Brown 
#'   analyzing population structure of wheat (Brown, 1980). It has been widely 
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
#'   \frac{V_O}{V_E}-1}{Ia = Vo/Ve} Since the distance is more or less a binary 
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
#'   with \eqn{{n \choose 2}\cdot{}m}{choose(n, 2) * m}. values of \emph{d} and 
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
#'   \deqn{ \bar{r_d} = \frac{V_O - V_E} {2\displaystyle 
#'   \sum_{j=1}^{m}\displaystyle \sum_{k \neq j}^{m}\sqrt{var_j\cdot{}var_k}} 
#'   }{rbarD = (Vo - Ve)/(2*sum(sum(sqrt(var(dj)*var(dk))))}
#'   
#' @references  Paul-Michael Agapow and Austin Burt. Indices of multilocus 
#'   linkage disequilibrium. \emph{Molecular Ecology Notes}, 1(1-2):101-102, 
#'   2001
#'   
#'   A.H.D. Brown, M.W. Feldman, and E. Nevo. Multilocus structure of natural 
#'   populations of hordeum spontaneum. \emph{Genetics}, 96(2):523-536, 1980.
#'   
#'   J M Smith, N H Smith, M O'Rourke, and B G Spratt. How clonal are bacteria? 
#'   Proceedings of the National Academy of Sciences, 90(10):4384-4388, 1993.
#'   
#' @seealso \code{\link{poppr}}, \code{\link{missingno}}, 
#'   \code{\link{import2genind}}, \code{\link{read.genalex}}, 
#'   \code{\link{clonecorrect}}
#'   
#' @export
#' @author Zhian N. Kamvar
#' @examples
#' data(nancycats)
#' ia(nancycats)
#' 
#' \dontrun{
#' # Get the indices back and plot them using base R graphics:
#' nansamp <- ia(nancycats, sample = 999, valuereturn = TRUE)
#' layout(matrix(c(1,1,2,2,), 2, 2, byrow = TRUE))
#' hist(nansamp$samples$Ia); abline(v = nansamp$index[1])
#' hist(nansamp$samples$rbarD); abline(v = nansamp$index[3])
#' 
#' # Get the index for each population.
#' lapply(seppop(nancycats), ia)
#' # With sampling
#' lapply(seppop(nancycats), ia, sample=999)
#' }
#==============================================================================#

ia <- function(pop, sample=0, method=1, quiet=FALSE, missing="ignore", 
                hist=TRUE, valuereturn = FALSE){
  METHODS = c("permute alleles", "parametric bootstrap",
              "non-parametric bootstrap", "multilocus")
  
  namelist            <- NULL
  namelist$population <- ifelse(length(levels(pop@pop)) > 1 | 
                                is.null(pop@pop), "Total", pop@pop.names)
  namelist$File       <- as.character(match.call()[2])
  
  popx    <- pop
  missing <- toupper(missing)
  type    <- pop@type
  
  if(type=="PA"){
    .Ia.Rd <- .PA.Ia.Rd
  } else {
    popx <- seploc(popx)
  }

  # if there are less than three individuals in the population, the calculation
  # does not proceed. 
  if (nInd(pop) < 3){
    IarD <- as.numeric(c(NA, NA))
    names(IarD) <- c("Ia", "rbarD")
    if (sample==0){
      return(IarD)
    } else {
      IarD <- as.numeric(rep(NA, 4))
      names(IarD) <- c("Ia","p.Ia","rbarD","p.rD")
      return(IarD)
    }
  }
  
  IarD <- .Ia.Rd(popx, missing)
  names(IarD) <- c("Ia", "rbarD")
  # no sampling, it will simply return two named numbers.
  if (sample==0){
    Iout   <- IarD
    result <- NULL
  } else {
  # sampling will perform the iterations and then return a data frame indicating
  # the population, index, observed value, and p-value. It will also produce a 
  # histogram.
    Iout     <- NULL 
    idx      <- as.data.frame(list(Index=names(IarD)))
    samp     <- .sampling(popx, sample, missing, quiet=quiet, type=type, method=method)
    p.val    <- sum(IarD[1] <= c(samp$Ia, IarD[1]))/(sample + 1)#ia.pval(index="Ia", samp2, IarD[1])
    p.val[2] <- sum(IarD[2] <= c(samp$rbarD, IarD[2]))/(sample + 1)#ia.pval(index="rbarD", samp2, IarD[2])
    if (hist == TRUE){
      poppr.plot(samp, observed=IarD, pop=namelist$population,
                        file=namelist$File, pval=p.val, N=nrow(pop@tab))
    }
    result         <- 1:4
    result[c(1,3)] <- IarD
    result[c(2,4)] <- p.val
    names(result)  <- c("Ia","p.Ia","rbarD","p.rD")
    if (valuereturn == TRUE){
      iaobj <- list(index = final(Iout, result), samples = samp)
      class(iaobj) <- "ialist"
      return(iaobj)
    } 
  }  
  return(final(Iout, result))
}



#==============================================================================#
#' Create a table of summary statistics per locus. 
#' 
#' @param x a \code{\linkS4class{genind}} or \code{\linkS4class{genclone}}
#'   object.
#' 
#' @param index Which diversity index to use. Choices are \itemize{ \item
#'   \code{"simpson"} (Default) to give Simpson's index \item \code{"shannon"}
#'   to give the Shannon-Wiener index \item \code{"invsimpson"} to give the
#'   Inverse Simpson's index aka the Stoddard and Tayor index.}
#'   
#' @param lev At what level do you want to analyze diversity? Choices are
#'   \code{"allele"} (Default) or \code{"genotype"}.
#'   
#' @param population Select the populations to be analyzed. This is the
#'   parameter \code{sublist} passed on to the function \code{\link{popsub}}.
#'   Defaults to \code{"ALL"}.
#'   
#' @param information When \code{TRUE} (Default), this will print out a header
#'   of information to the R console.
#'   
#' @return a table with 4 columns indicating the Number of alleles/genotypes
#'   observed, Diversity index chosen, Nei's 1978 expected heterozygosity, and
#'   Evenness.
#'   
#' @seealso \code{\link[vegan]{diversity}}, \code{\link{poppr}}
#'   
#' @note This will calculate statistics for polyploids as well by only counting
#'   observed allelic states.
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
#' # Analyze locus statistics for the North American population of P. infestans.
#' data(Pinf)
#' locus_table(Pinf, population = "North America")
#==============================================================================#
locus_table <- function(x, index = "simpson", lev = "allele", 
                        population = "ALL", information = TRUE){
  INDICES <- c("shannon", "simpson", "invsimpson")
  index   <- match.arg(index, INDICES)
  x       <- popsub(x, population, drop = FALSE)
  x.loc   <- summary(as.loci(x))
  outmat  <- vapply(x.loc, locus_table_pegas, numeric(4), index, lev, x@type)
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
    cat("\n", divs[1], "= Number of observed", paste0(divs[1], "s"))
    cat("\n", divs[2], "=", msg)
    cat("\n", divs[3], "= Nei's 1978 expected heterozygosity\n")
    cat("------------------------------------------\n")
  }
  class(res) <- c("locustable", "matrix")
  return(res)
}

#==============================================================================#
#' Tablulate alleles the occur in only one population. 
#' 
#' @param gid a \code{\linkS4class{genind}} or \code{\linkS4class{genclone}}
#'   object.
#'   
#' @param report one of \code{"table", "vector",} or \code{"data.frame"}. Tables
#'   (Default) and data frame will report counts along with populations or 
#'   individuals. Vectors will simply report which populations or individuals 
#'   contain private alleles. Tables are matrices with populations or 
#'   individuals in rows and alleles in columns. Data frames are long form.
#'   
#' @param level one of \code{"population"} (Default) or \code{"individual"}.
#'   
#' @return a matrix, data.frame, or vector defining the populations or
#'   individuals containing private alleles. If vector is chosen, alleles are
#'   not defined.
#'
#' @export
#' @examples
#' 
#' data(Pinf) # Load P. infestans data.
#' setpop(Pinf) <- ~Country # Set the population to be at the country level
#' private_alleles(Pinf)
#' \dontrun{
#' # An example of how this data can be displayed.
#' library(ggplot2)
#' Pinfpriv <- private_alleles(Pinf, report = "data.frame")
#' ggplot(Pinfpriv) + geom_tile(aes(x = population, y = allele, fill = count))
#' }
#==============================================================================#
private_alleles <- function(gid, report = "table", level = "population"){
  REPORTARGS <- c("table", "vector", "data.frame")
  LEVELARGS  <- c("individual", "population")
  report <- match.arg(report, REPORTARGS)
  level  <- match.arg(level, LEVELARGS)
  if (!is.genind(gid) & !is.genpop(gid)){
    stop(paste(gid, "is not a genind or genpop object."))
  }
  if (is.genind(gid) & !is.null(pop(gid)) | is.genpop(gid) & nrow(gid@tab) > 1){
    if (is.genind(gid)){
      gid.pop <- truenames(genind2genpop(gid, quiet = TRUE))
    } else {
      gid.pop <- truenames(gid)
    }
    privates <- gid.pop[, colSums(ifelse(gid.pop > 0, 1, 0), na.rm = TRUE) < 2]
    privates <- privates[rowSums(privates) > 0, ]
    if (level == "individual" & is.genind(gid)){
      gid.tab  <- truenames(gid)$tab
      privates <- gid.tab[, colnames(gid.pop) %in% colnames(privates)]
      privates <- privates[rowSums(privates, na.rm = TRUE) > 0, ]
    }
    if (length(privates) == 0){
      privates <- NULL
      cat("No private alleles detected.")
      return(invisible(NULL))
    }
    if (report == "vector"){
      privates <- rownames(privates)
    } else if (report == "data.frame"){
      privates <- melt(privates, varnames = c(level, "allele"), 
                       value.name = "count")
    }
    return(privates)
  } else {
    stop("There are no populations detected")
  }
}
