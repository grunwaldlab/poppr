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
setClassUnion("mlgORnumeric", c("MLG", "numeric"))
#==============================================================================#
#' GENclone and SNPclone classes
#' 
#' @description \strong{GENclone} is an S4 class that extends the 
#'   \code{\linkS4class{genind}} object.\cr \strong{SNPclone} is an S4 class
#'   that extends the \code{\linkS4class{genlight}} object.\cr\cr They will have
#'   all of the same attributes as their parent classes, but they will contain
#'   one extra slot to store extra information about multilocus genotypes.
#'   
#' @section Extends: The \code{genclone} class extends class 
#'   \code{"\linkS4class{genind}"}, directly. \cr The \code{snpclone} class 
#'   extends class \code{"\linkS4class{genlight}"}, directly.
#'   
#' @details The genclone and snpclone classes will allow for more optimized 
#'   methods of clone correction.
#'   
#'   Previously for \linkS4class{genind} and \linkS4class{genlight} objects, 
#'   multilocus genotypes were not retained after a data set was subset by 
#'   population. The new \strong{\code{mlg}} slot allows us to assign the 
#'   multilocus genotypes and retain that information no matter how we subset 
#'   the data set. This new slot can either contain numeric values for 
#'   multilocus genotypes OR it can contain a special internal 
#'   \code{\linkS4class{MLG}} class that allows for custom multilocus genotype 
#'   definitions and filtering.
#'   
#' @name genclone-class
#' @rdname genclone-class
#' @aliases genclone
#' @export
#' @slot mlg a vector representing multilocus genotypes for the data set OR an 
#'   object of class \code{\linkS4class{MLG}}.
#' @author Zhian N. Kamvar
#' @seealso \code{\link{as.genclone}} \code{\link{as.snpclone}} 
#'   \code{\linkS4class{genind}} \code{\linkS4class{genlight}} 
#'   \code{\link[adegenet]{strata}} \code{\link[adegenet]{setPop}} 
#'   \code{\link{MLG}} \code{\link{mll}} \code{\link{mlg.filter}}
#' @import methods
#' @examples
#' \dontrun{
#' 
#' # genclone objects can be created from genind objects
#' #
#' data(partial_clone)
#' partial_clone
#' (pc <- as.genclone(partial_clone))
#' 
#' # snpclone objects can be created from genlight objects
#' #
#' set.seed(999)
#' (gl <- glSim(100, 0, n.snp.struc = 1e3, ploidy = 2, parallel = FALSE))
#' (sc <- as.snpclone(gl, parallel = FALSE))
#' # 
#' # Use mlg.filter to create a distance threshold to define multilocus genotypes.
#' mlg.filter(sc, threads = 1L) <- 0.25
#' sc # 82 mlgs
#' 
#' }
#==============================================================================#
setClass("genclone", 
         contains = "genind",
         representation = representation(mlg = "mlgORnumeric"),
         prototype(mlg = integer(0))
)

valid.genclone <- function(object){
  slots   <- slotNames(object)
  if (any(!"mlg" %in% slots)){
    return(FALSE)
  }
  inds    <- adegenet::nInd(object)
  mlgs    <- length(object@mlg)
  if (mlgs != inds){  
    message("Multilocus genotypes do not match the number of observations")
    return(FALSE)
  }
  return(TRUE)
}

setValidity("genclone", valid.genclone)

#==============================================================================#
#' @name snpclone-class
#' @rdname genclone-class
#' @aliases snpclone
#' @export
#==============================================================================#
setClass("snpclone",
         contains = "genlight",
         representation = representation(mlg = "mlgORnumeric"),
         prototype = prototype(mlg = integer(0))
)

setValidity("snpclone", valid.genclone)
#==============================================================================#
#' bruvomat object
#' 
#' An internal object used for bruvo's distance. 
#' Not intended for user interaction.
#' 
#' 
#' @name bruvomat-class
#' @rdname bruvomat-class
#' @export
#' @slot mat a matrix of genotypes with one allele per locus. Number of rows
#'   will be equal to (ploidy)*(number of loci)
#' @slot replen repeat length of microsatellite loci
#' @slot ploidy the ploidy of the data set
#' @slot ind.names names of individuals in matrix rows.
#' @keywords internal
#' @author Zhian N. Kamvar
#' @import methods
#==============================================================================#
setClass(
  Class = "bruvomat", 
  representation = representation(
    mat = "matrix", 
    replen = "numeric",
    ploidy = "numeric",
    ind.names = "character"
  ),
  prototype = prototype(
    mat = matrix(ncol = 0, nrow = 0),
    replen = integer(0),
    ploidy = integer(0),
    ind.names = character(0)
  )
)

#==============================================================================#
#' Bootgen object
#' 
#' An internal object used for bootstrapping. Not intended for user interaction.
#' 
#' @section Extends: 
#' Virtual Class \code{"\linkS4class{gen}"}.
#' 
#' @name bootgen-class
#' @rdname bootgen-class
#' @export
#' @slot type a character denoting Codominant ("codom") or Dominant data ("P/A")
#' @slot ploidy an integer denoting the ploidy of the data set. (>=1)
#' @slot alllist a list with numeric vectors, each representing a different
#'   locus where each element in the vector represents the index for a specific
#'   allele.
#' @slot names a vector containing names of the observed samples.
#' @keywords internal
#' @author Zhian N. Kamvar
#' @import methods
#==============================================================================#
setClass("bootgen", 
         contains = c("gen"),
         representation = representation(
                          type = "character",
                          ploidy = "integer",
                          names = "vector", 
                          alllist = "list"),
         prototype = prototype(
          type = character(0),
          ploidy = integer(0),
          names = character(0),
          alllist = list()
          )
)
