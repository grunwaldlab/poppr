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
#' Utilize all algorithms of mlg.filter
#' 
#' This function is a wrapper to mlg.filter. It will calculate all of the stats 
#' for mlg.filter utilizing all of the algorithms.
#' 
#' @param x a \code{\link{genind}}, \code{\link{genclone}}, \code{\link{genlight}}, or \code{\link{snpclone}} object
#' @param distance a distance function or matrix
#' @param threshold a threshold to be passed to \code{\link{mlg.filter}}
#'   (Default: 1e6)
#' @param stats what statistics should be calculated.
#' @param missing how to treat missing data with mlg.filter
#' @param plot If the threshold is a maximum threshold, should the statistics be
#'   plotted (Figure 2)
#' @param cols the colors to use for each algorithm (defaults to set1 of 
#'   \pkg{RColorBrewer}).
#' @param nclone the number of multilocus genotypes you expect for the data. 
#'   This will draw horizontal line on the graph at the value nclone and then 
#'   vertical lines showing the cutoff thresholds for each algorithm.
#' @param hist if you want a histogram to be plotted behind the statistics, 
#'   select a method here. Available methods are "sturges", "fd", or "scott"
#'   (default) as documented in \code{\link[graphics]{hist}}. If you don't want
#'   to plot the histogram, set \code{hist = NULL}.
#' @param ... extra parameters passed on to the distance function.
#'   
#' @return a list of results from mlg.filter from the three algorithms.
#' @export
#' @seealso \code{\link{mlg.filter}} 
#'   \code{\link{cutoff_predictor}} 
#'   \code{\link{bitwise.dist}} 
#'   \code{\link{diss.dist}}
#' @note This function originally appeared in 
#'   \href{http://dx.doi.org/10.5281/zenodo.17424}{DOI: 10.5281/zenodo.17424}
#' @references ZN Kamvar, JC Brooks, and NJ Grünwald. 2015. Supplementary
#' Material for Frontiers Plant Genetics and Genomics 'Novel R tools for
#' analysis of genome-wide population genetic data with emphasis on clonality'.
#' DOI: \href{http://dx.doi.org/10.5281/zenodo.17424}{10.5281/zenodo.17424}
#' 
#' Kamvar ZN, Brooks JC and Grünwald NJ (2015) Novel R tools for analysis of 
#' genome-wide population genetic data with emphasis on clonality. Front. Genet.
#' 6:208. doi:
#' \href{http://dx.doi.org/10.3389/fgene.2015.00208}{10.3389/fgene.2015.00208}
#' 
#' @author Zhian N. Kamvar, Jonah C. Brooks
#' @examples
#' \dontrun{
#' data(Pinf)
#' filter_stats(Pinf, distance = diss.dist, percent = TRUE, plot = TRUE)
#' }
#==============================================================================#
filter_stats <- function(x, distance = bitwise.dist, 
                         threshold = 1e6 + .Machine$double.eps^0.5, 
                         stats = "All", missing = "ignore", plot = FALSE, 
                         cols = NULL, nclone = NULL, hist = "Scott", ...){
  if (!inherits(distance, "dist")){
    DIST <- match.fun(distance)
    if (inherits(x, "genind")){
      x <- missingno(x, type = missing)
    }
    distmat <- DIST(x, ...)
  } else {
    distmat <- distance
  }
  f <- mlg.filter(x, threshold, missing, algorithm = "f", distance = distmat, 
                  stats = stats, ...)
  a <- mlg.filter(x, threshold, missing, algorithm = "a", distance = distmat, 
                  stats = stats, ...)
  n <- mlg.filter(x, threshold, missing, algorithm = "n", distance = distmat, 
                  stats = stats, ...)
  fanlist <- list(farthest = f, average = a, nearest = n)
  if (stats == "All"){
    if (plot){
      plot_filter_stats(x, fanlist, distmat, cols, nclone, hist)
    }
  }
  return(fanlist)
}
#==============================================================================#
#' Predict cutoff thresholds for use with mlg.filter
#' 
#' Given a series of thresholds for a data set that collapse it into one giant 
#' cluster, this will search the top fraction of threshold differences to find 
#' the largest difference. The average between the thresholds spanning that 
#' difference is the cutoff threshold defining the clonal lineage threshold.
#' 
#' @param thresholds a vector of numerics coming from mlg.filter where the 
#'   threshold has been set to the maximum threshold theoretically possible.
#' @param fraction the fraction of the data to seek the threshold.
#'   
#' @return a numeric value representing the threshold at which multilocus 
#'   lineages should be defined.
#' @seealso \code{\link{filter_stats}} \code{\link{mlg.filter}}
#' @note This function originally appeared in 
#'   \href{http://dx.doi.org/10.5281/zenodo.17424}{DOI: 10.5281/zenodo.17424}. 
#'   This is a bit of a blunt instrument.
#' @export
#' @references ZN Kamvar, JC Brooks, and NJ Grünwald. 2015. Supplementary
#' Material for Frontiers Plant Genetics and Genomics 'Novel R tools for
#' analysis of genome-wide population genetic data with emphasis on clonality'.
#' DOI: \href{http://dx.doi.org/10.5281/zenodo.17424}{10.5281/zenodo.17424}
#' 
#' Kamvar ZN, Brooks JC and Grünwald NJ (2015) Novel R tools for analysis of 
#' genome-wide population genetic data with emphasis on clonality. Front. Genet.
#' 6:208. doi:
#' \href{http://dx.doi.org/10.3389/fgene.2015.00208}{10.3389/fgene.2015.00208}
#' 
#' @author Zhian N. Kamvar
#' @examples
#' \dontrun{
#' data(Pinf)
#' pthresh <- mlg.filter(Pinf, distance = diss.dist, percent = TRUE, 
#'                       threshold = 1.1, stats = "THRESH")
#' cutoff_predictor(pthresh)
#' }
#==============================================================================#
cutoff_predictor <- function(thresholds, fraction = 0.5){
  frac <- 1:round(length(thresholds)*fraction)
  diffs <- diff(thresholds[frac])
  diffmax <- which.max(diffs)
  mean(thresholds[diffmax:(diffmax + 1)])
}

#==============================================================================#
#' Plot the results of filter_stats
#' 
#' @param x a genlight of genind object
#' @param fstats the list passed from \code{\link{filter_stats}}
#' @param distmat a distance matrix passed from \code{\link{filter_stats}}
#' @param cols colors to use for each algorithm (defaults to \pkg{RColorBrewer}
#'   set 1)
#' @param nclone see \code{\link{filter_stats}}
#'   
#' @return a plot depicting how many MLLs are collapsed as the genetic distance 
#'   increases for each algorithm.
#' @export
#' @seealso \code{\link{filter_stats}}
#' @note This function originally appeared in 
#'   \href{http://dx.doi.org/10.5281/zenodo.17424}{DOI: 10.5281/zenodo.17424}
#' @author Zhian N. Kamvar
#' @references ZN Kamvar, JC Brooks, and NJ Grünwald. 2015. Supplementary
#' Material for Frontiers Plant Genetics and Genomics 'Novel R tools for
#' analysis of genome-wide population genetic data with emphasis on clonality'.
#' DOI: \href{http://dx.doi.org/10.5281/zenodo.17424}{10.5281/zenodo.17424}
#' 
#' Kamvar ZN, Brooks JC and Grünwald NJ (2015) Novel R tools for analysis of 
#' genome-wide population genetic data with emphasis on clonality. Front. Genet.
#' 6:208. doi:
#' \href{http://dx.doi.org/10.3389/fgene.2015.00208}{10.3389/fgene.2015.00208}
#' 
#' @keywords internal
#==============================================================================#
plot_filter_stats <- function(x, fstats, distmat, cols = NULL, nclone = NULL, breaks = NULL){
  upper <- round(max(distmat), digits = 1)
  ylims <- c(ifelse(is.genind(x), mlg(x, quiet = TRUE), nInd(x)), 1)
  if (!is.null(breaks)){
    graphics::hist(distmat, breaks = breaks, xlab = "", ylab = "", axes = FALSE,
                   xlim = c(0, upper), main = "")
    par(new = TRUE)
  }
  plot(x = c(upper, 0), y = ylims, type = "n",
       ylab = "Number of Multilocus Lineages",
       xlab = "Genetic Distance Cutoff")
  a <- fstats$average$THRESHOLDS
  n <- fstats$nearest$THRESHOLDS
  f <- fstats$farthest$THRESHOLDS
  plotcols <- c("#E41A1C", "#377EB8", "#4DAF4A") # RColorBrewer::brewer.pal(3, "Set1")
  names(plotcols) <- c("f", "a", "n")
  points(x = rev(a), y = 1:length(a), col = plotcols["a"])
  points(x = rev(f), y = 1:length(f), col = plotcols["f"])
  points(x = rev(n), y = 1:length(n), col = plotcols["n"])
  if (!is.null(nclone)){
    abline(v = a[1 + length(a) - nclone] + .Machine$double.eps^0.5, lty = 2,
           col = plotcols["a"])
    abline(v = f[1 + length(f) - nclone] + .Machine$double.eps^0.5, lty = 2, 
           col = plotcols["f"])
    abline(v = n[1 + length(n) - nclone] + .Machine$double.eps^0.5, lty = 2, 
           col = plotcols["n"])
    abline(h = nclone)
    text(upper, nclone, labels = paste0("n = ", nclone), 
         adj = c(1, -0.5))
    legend("topright", 
           legend = c("Nearest Neighbor", "UPGMA", "Farthest Neighbor"), 
           col = plotcols[c("n", "a", "f")], 
           pch = 1, 
           lty = 2,
           title = "Clustering Method")
  } else {
    legend("topright", 
           legend = c("Nearest Neighbor", "UPGMA", "Farthest Neighbor"), 
           col = plotcols[c("n", "a", "f")], 
           pch = 1, 
           title = "Clustering Method")    
  }
}