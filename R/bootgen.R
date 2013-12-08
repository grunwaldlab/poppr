#' Bootgen object
#' 
#' An internal object used for bootstrapping. Not intended for user interaction.
#' 
#' @section Extends: 
#' Class \code{"\linkS4class{genind}"}, directly.
#' 
#' @name bootgen-class
#' @rdname bootgen-class
#' @export
#' @slot replen repeat length of microsatellite loci
#' @author Zhian Kamvar
#' @import methods

setClass("bootgen", 
         contains = c("genind"),
         representation = representation(replen = "numeric"),
)


#' Methods used for the bootgen object. 
#' 
#' This is not designed for user interaction.
#' 
#' @rdname bootgen-methods
#' @param x a \code{"\linkS4class{bootgen}"} object
#' @param i vector of numerics indicating number of individuals desired
#' @param j a vector of numerics corresponding to the loci desired.
#' @param ... unused.
#' @param drop set to \code{FALSE}
setMethod(
  f = "[",
  signature(x = "bootgen"),
  definition = function(x, i, j, ..., drop = FALSE){
    if (missing(i)) i <- TRUE
    if (missing(j)) j <- TRUE
    loci <- levels(x@loc.fac)
    indices <- lapply(loci[j], function(locus){res <- which(x@loc.fac %in% locus)})
    names(indices) <- paste0(loci[1:length(indices)], ".")
    indices <- unlist(indices)
    res <- x@tab[i, indices, drop = drop]
    colnames(res) <- names(indices)
    return(new("bootgen", gen = genind(res), replen = x@replen[j]))
  }
)



#' @rdname bootgen-methods
setMethod(
  f = "dim",
  signature(x = "bootgen"),
  definition = function(x){
    return(c(nInd(x), nLoc(x)))
  }
)


#' @rdname bootgen-methods
#' @param .Object a character, "bootgen"
#' @param gen \code{"\linkS4class{genind}"} object
#' @param replen a vector of numbers indicating the repeat length for each microsatellite locus.
setMethod(
  f = "initialize",
  signature = "bootgen",
  definition = function(.Object, gen, replen){
    if (missing(gen)) gen <- new("genind")
    if (missing(replen)){
      replen <- vapply(gen@all.names, function(y) guesslengths(as.numeric(y)), 1)
    }
    lapply(names(gen), function(y) slot(.Object, y) <<- slot(gen, y))
    slot(.Object, "replen") <- replen
    return(.Object)
  })
