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

################################################################################
#------------------------------------------------------------------------------#
# BOOTGEN METHODS
#------------------------------------------------------------------------------#
################################################################################

#==============================================================================#
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
#==============================================================================#
setMethod(
  f = "[",
  signature(x = "bootgen"),
  definition = function(x, i, j, ..., drop = FALSE){
    if (missing(i)) i <- TRUE
    if (missing(j)) j <- TRUE
    if (length(j) > nLoc(x) | any(j > nLoc(x))){
      stop('subscript out of bounds')
    }
    
    # Taking Names
    locnall         <- x@loc.nall[j]
    allnames        <- x@all.names[j]
    names(allnames) <- names(x@all.names)[1:length(j)]
    names(locnall)  <- names(allnames)
    
    # Shuffling
    # matcols  <- apply(getinds(x@loc.nall), 1, function(ind) ind[1]:ind[2])[j]
    matcols  <- .Call("expand_indices", cumsum(x@loc.nall), nLoc(x))[j]
    indices  <- unlist(matcols)
    locnames <- rep(names(allnames), locnall)
    tabnames <- paste(locnames, unlist(allnames), sep = ".")
    res      <- slot(x, "tab")[i, indices, drop = drop]
    colnames(res) <- tabnames
    
    ## Resetting all factors that need to be set. 
    slot(x, "tab")       <- res
    slot(x, "loc.fac")   <- factor(locnames, names(allnames))
    slot(x, "loc.names") <- names(allnames)
    slot(x, "loc.nall")  <- locnall
    slot(x, "all.names") <- allnames
    slot(x, "replen")    <- slot(x, "replen")[j]
    return(x)
  }
)

#==============================================================================#
#' @rdname bootgen-methods
#==============================================================================#
setMethod(
  f = "dim",
  signature(x = "bootgen"),
  definition = function(x){
    return(c(nInd(x), nLoc(x)))
  }
)

#==============================================================================#
#' @rdname bootgen-methods
#' @param .Object a character, "bootgen"
#' @param gen \code{"\linkS4class{genind}"} object
#' @param replen a vector of numbers indicating the repeat length for each 
#' microsatellite locus.
#==============================================================================#
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

################################################################################
#------------------------------------------------------------------------------#
# BRUVOMAT METHODS
#------------------------------------------------------------------------------#
################################################################################

#==============================================================================#
#' @rdname bruvomat-methods
#' @param .Object a character, "bruvomat"
#' @param gen \code{"\linkS4class{genind}"} object
#' @param replen a vector of numbers indicating the repeat length for each 
#' microsatellite locus.
#==============================================================================#
setMethod(
  f = "initialize",
  signature = "bruvomat",
  definition = function(.Object, gen, replen){
    if (missing(gen)) gen <- new("genind")
    if (missing(replen)){
      replen <- vapply(gen@all.names, function(y) guesslengths(as.numeric(y)), 1)
    }
    ploid <- ploidy(gen)
    # This controlls for the user correcting missing data using "mean". 
    if (any(!gen@tab %in% c((0:ploid)/ploid, NA))){
      gen@tab[!gen@tab %in% c((0:ploid)/ploid, NA)] <- NA
    }
    # This will check for data that has missing scored as "zero".
    popcols <- ploid*nLoc(gen)
    if (!any(is.na(gen@tab)) & any(rowSums(gen@tab, na.rm=TRUE) < nLoc(gen))){
      mat1 <- as.matrix.data.frame(genind2df(gen, sep="/", usepop=FALSE))
      mat1[mat1 %in% c("", NA)] <- paste(rep(0, ploid), collapse="/")
      mat2 <- apply(mat1, 1, strsplit, "/")
      mat3 <- apply(as.matrix(t(sapply(mat2, unlist))), 2, as.numeric)
      vec1 <- suppressWarnings(as.numeric(unlist(mat3)))
      pop  <- matrix(vec1, nrow=nInd(gen), ncol=popcols)
    } else {
      popdf <- genind2df(gen, oneColPerAll=TRUE, usepop=FALSE)
      mat1  <- as.matrix.data.frame(popdf)
      pop   <- suppressWarnings(matrix(as.numeric(mat1), ncol=popcols))
    }
    slot(.Object, "mat")       <- pop
    slot(.Object, "replen")    <- replen
    slot(.Object, "ploidy")    <- ploid
    slot(.Object, "ind.names") <- indNames(gen)
    return(.Object)
  }
)

#==============================================================================#
#' @rdname bruvomat-methods
#==============================================================================#
setMethod(
  f = "dim",
  signature(x = "bruvomat"),
  definition = function(x){
    return(c(nrow(x@mat), ncol(x@mat)/x@ploidy))
  }
)

#==============================================================================#
#' Methods used for the bruvomat object. 
#' 
#' This is not designed for user interaction.
#' 
#' @rdname bruvomat-methods
#' @param x a \code{"\linkS4class{bruvomat}"} object
#' @param i vector of numerics indicating number of individuals desired
#' @param j a vector of numerics corresponding to the loci desired.
#' @param ... unused.
#' @param drop set to \code{FALSE}
#==============================================================================#
setMethod(
  f = "[",
  signature(x = "bruvomat"),
  definition = function(x, i, j, ..., drop = FALSE){
    if (missing(i)) i <- TRUE
    if (missing(j)) j <- TRUE
    x@replen    <- x@replen[j]
    x@ind.names <- x@ind.names[i]
    cols        <- rep(1:ncol(x), each = x@ploidy)
    replacement <- vapply(j, function(ind) which(cols == ind), 1:x@ploidy)
    x@mat       <- x@mat[i, as.vector(replacement), drop = FALSE]
    return(x)
  }
)

################################################################################
#------------------------------------------------------------------------------#
# GENCLONE METHODS
#------------------------------------------------------------------------------#
################################################################################
#==============================================================================#
#' Check for validity of a genclone object
#' 
#' In the hat
#' 
#' @export
#' @rdname is.genclone
#' @param x a genclone object
#' @examples
#' data(nancycats)
#' nanclone <- as.genclone(nancycats)
#' is.genclone(nanclone)
#==============================================================================#
is.genclone <- function(x){ 
  res <- (is(x, "genclone") & validObject(x))
  return(res)
}

#==============================================================================#
#' Methods used for the genclone object
#' 
#' words.
#' 
#' @rdname genclone-method
#' @param x a genclone object
#' @param i vector of numerics indicating number of individuals desired
#' @param j a vector of numerics corresponding to the loci desired.
#' @param ... passed on to the \code{\linkS4class{genind}} object.
#' @param drop set to \code{FALSE}
#==============================================================================#
setMethod(
  f = "[",
  signature(x = "genclone", i = "ANY", j = "ANY", drop = "ANY"),
  definition = function(x, i, j, ..., loc=NULL, treatOther=TRUE, quiet=TRUE, drop = FALSE){
    if (missing(i)) i <- TRUE
    if (missing(j)) j <- TRUE
    mlg       <- slot(x, "mlg")[i]
    hierarchy <- slot(x, "hierarchy")[i, , drop = FALSE]
    x <- callNextMethod()
    x <- new("genclone", x, hierarchy, mlg)
    return(x)
  }
)

#==============================================================================#
#' @rdname genclone-method
#' @param .Object a character, "genclone"
#' @param gen \code{"\linkS4class{genind}"} object
#' @param hierarchy a data frame where each row i represents the different
#' population assignments of individual i in the data set. If this is empty, the
#' hierarchy will be created from the population factor.
#' @param mlg a vector where each element assigns the multilocus genotype of
#' that individual in the data set. 
#==============================================================================#
setMethod(      
  f = "initialize",
  signature("genclone"),
  definition = function(.Object, gen, hierarchy, mlg){
    if (missing(gen)){
      gen <- new("genind")
      if (missing(mlg)) mlg <- 0
    } else {
      if (missing(mlg)) mlg <- mlg.vector(gen)
    }
    if (missing(hierarchy)){
      if (is.null(pop(gen))){
        hierarchy <- data.frame()
      } else {
        hierarchy <- data.frame(Pop = pop(gen))
      }
    }

    # No 'initialize' method for genind objects...
    lapply(names(gen), function(y) slot(.Object, y) <<- slot(gen, y))
    slot(.Object, "mlg")       <- mlg
    slot(.Object, "hierarchy") <- hierarchy
    return(.Object)
  }
)

#==============================================================================#
#' @rdname genclone-method
#' @param object a genclone object
#==============================================================================#
setMethod(
  f = "show",
  signature("genclone"),
  definition = function(object){
    ploid <- c("ha", "di", "tri", "tetra", "penta", "hexa", "hepta", "octa",
      "nona", "deca", "hendeca", "dodeca")
    ploid <- paste0(ploid[object@ploidy], "ploid")
    nind  <- nInd(object)
    type  <- ifelse(object@type == "PA", "dominant", "codominant")
    nmlg  <- length(unique(object@mlg))
    nloc  <- nLoc(object)
    npop  <- ifelse(is.null(object@pop), 0, length(object@pop.names))
    hier  <- length(object@hierarchy)
    chars <- nchar(c(nmlg, nind, nloc, hier, npop))
    ltab  <- max(chars) - chars
    ltab  <- vapply(ltab, function(x) substr("\t     ", 1, x + 1), character(1))
    pops  <- object@pop.names
    hiernames <- names(object@hierarchy)
    cat("\nThis is a genclone object\n")
    cat("-------------------------\n")
    cat("Genotype information:\n\n",
      nmlg, ltab[1], "multilocus genotypes\n",
      nind, ltab[2], ploid, "individuals\n", 
      nloc, ltab[3], type, "loci\n\n"
      )
    pophier <- ifelse(hier > 1, "hierarchies -", "hierarchy -")
    if (hier == 0) pophier <- "hierarchies."
    popdef  <- ifelse(npop > 0, "defined -", "defined.")
    cat("Population information:\n\n")
    cat("", hier, ltab[4], "population", pophier, hiernames, fill = TRUE)
    cat("", npop, ltab[5], "populations", popdef, pops, fill = TRUE)
    
  })

#==============================================================================#
#' Create a genclone object from a genind object.
#' 
#' Wrapper for genclone initializer. The hierarchy must have the same number of
#' rows as the number of observations in the genind object. If no hierarchy is
#' defined, the population will be set as the hierarchy under the label "Pop". 
#' 
#' @export 
#' @rdname coercion-methods
#' @aliases as.genclone,genind-method
#' @param x a \code{\linkS4class{genind}} or \code{\linkS4class{genclone}} object
#' @param hierarchy a data frame representing the population hierarchy.
#' @docType methods
#' @examples
#' data(Aeut)
#' Aeut
#' Aeut.gc <- as.genclone(Aeut)
#' Aeut.gc
#' Aeut.gc <- as.genclone(Aeut, other(Aeut)$population_hierarchy[-1])
#' Aeut.gc
#==============================================================================#
as.genclone <- function(x, hierarchy = NULL){
  standardGeneric("as.genclone")
}

#' @export
setGeneric("as.genclone")


setMethod(
  f = "as.genclone",
  signature(x = "genind"),
  definition = function(x, hierarchy){
    if (missing(hierarchy)){
      if ("population_hierarchy" %in% names(other(x))){
        hierarchy   <- other(x)[["population_hierarchy"]]
        newgenclone <- new("genclone", x, hierarchy)
      } else {
        newgenclone <- new("genclone", x)
      }
    } else {
      newgenclone   <- new("genclone", x, hierarchy)
    }
    return(newgenclone)
  })

#==============================================================================#
#' Access and manipulate the population hierarchy for genclone objects.
#' 
#' The following methods allow the user to quickly change the hierarchy or
#' population of a genclone object. 
#' 
#' @export 
#' @rdname hierarchy-methods
#' @aliases gethierarchy,genclone-method
#' @param x a genclone object
#' @param formula a nested formula indicating the order of the population
#' hierarchy.
#' @param combine if \code{TRUE}, the levels will be combined according to the
#' formula argument. If it is \code{FALSE}, the levels will not be combined.
#' @docType methods
#==============================================================================#
gethierarchy <- function(x, formula = NULL, combine = TRUE){
  standardGeneric("gethierarchy")
} 

#' @export
setGeneric("gethierarchy")

setMethod(
  f = "gethierarchy",
  signature(x = "genclone"),
  definition = function(x, formula = NULL, combine = TRUE){
    if (is.null(formula)) return(x@hierarchy)
    vars <- all.vars(formula)
    if (any(!vars %in% names(x@hierarchy))){
      stop(hier_incompatible_warning(vars, x@hierarchy))
    }
    if (combine){
      hier <- make_hierarchy(formula, x@hierarchy)
    } else {
      hier <- x@hierarchy[all.vars(formula)]
    }
    invisible(return(hier))
  })

#==============================================================================#
#' @export
#' @rdname hierarchy-methods
#' @aliases sethierarchy<-,genclone-method
#' @param value a data frame OR vector OR formula (see details).
#' @docType methods
#' 
#' @details
#' These functions allow the user to seamlessly assign the hierarchy of their
#' \code{\linkS4class{genclone}} object. Note that there are two ways of
#' performing \code{sethierarchy()} and \code{namehierarchy()}, they essentially
#' do the same thing except that the assignment method (the one with the "\code{<-}") 
#' will modify the object in place whereas the non-assignment method will not 
#' modify the original object. Due to convention, everything right of
#' the assignment is termed \code{value}. To avoid confusion, here is a guide:
#' \enumerate{
#'  \item \strong{sethierarchy()} This will be a \code{\link{data.frame}} that 
#'  defines the hierarchy for each individual in the rows.
#'  \item \strong{namehierarchy()} This will be either a \code{\link{vector}} or
#'  a \code{\link{formula}} that will define the names.
#' }
#' 
#' @examples
#' data(Aeut)
#' Aeut.gc <- as.genclone(Aeut)
#' head(gethierarchy(Aeut.gc))
#' popsub <- other(Aeut)$population_hierarchy[-1]
#' Aeut.gc <- sethierarchy(Aeut.gc, popsub)
#' sethierarchy(Aeut.gc) <- popsub
#' head(gethierarchy(Aeut.gc))
#' namehierarchy(Aeut.gc) <- ~Population/Subpopulation
#' Aeut.gc
#' Aeut.old <- namehierarchy(Aeut.gc, ~Pop/Subpop)
#' Aeut.old
#==============================================================================#
sethierarchy <- function(x, value){
  standardGeneric("sethierarchy")
} 

#' @export
setGeneric("sethierarchy")

setMethod(
  f = "sethierarchy",
  signature(x = "genclone"),
  definition = function(x, value){
    if (!inherits(value, "data.frame")){
      stop(paste(substitute(value), "is not a data frame"))
    }
    if (nrow(value) != nInd(x)){
      stop("Number of rows in data frame not equal to number of individuals in object.")
    }
    x@hierarchy <- value
    return(x)
  })

#==============================================================================#
#' @export 
#' @rdname hierarchy-methods
#' @aliases sethierarchy,genclone-method
#' @docType methods
#==============================================================================#
"sethierarchy<-" <- function(x, value){
  standardGeneric("sethierarchy<-")
}  

#' @export
setGeneric("sethierarchy<-")

setMethod(
  f = "sethierarchy<-",
  signature(x = "genclone"),
  definition = function(x, value){
    return(sethierarchy(x, value))
  })

#==============================================================================#
#' @export 
#' @rdname hierarchy-methods
#' @aliases namehierarchy,genclone-method
#' @docType methods
#==============================================================================#
"namehierarchy" <- function(x, value){
  standardGeneric("namehierarchy")
}  

#' @export
setGeneric("namehierarchy")

setMethod(
  f = "namehierarchy",
  signature(x = "genclone"),
  definition = function(x, value){
    if (is.language(value)){
      value <- all.vars(value)
    }
    if (!is.vector(value) | length(value) != length(x@hierarchy)){
      stop(paste("Hierarchy, needs a vector argument of length", length(x@hierarchy)))
    }
    names(x@hierarchy) <- value
    return(x)
  })

#==============================================================================#
#' @export 
#' @rdname hierarchy-methods
#' @aliases namehierarchy<-,genclone-method
#' @docType methods
#==============================================================================#
"namehierarchy<-" <- function(x, value){
  standardGeneric("namehierarchy<-")
}  

#' @export
setGeneric("namehierarchy<-")

setMethod(
  f = "namehierarchy<-",
  signature(x = "genclone"),
  definition = function(x, value){
    return(namehierarchy(x, value))
  })

#==============================================================================#
#' Manipulate the population factor of genclone objects.
#' 
#' The following methods allow the user to quickly change the population of a 
#' genclone object. 
#' 
#' @export 
#' @rdname population-methods
#' @param x a genclone object
#' @param formula a nested formula indicating the order of the population
#' hierarchy.
#' @param value same as formula
#' @aliases setpop,genclone-method
#' @docType methods
#' @examples
#' 
#' data(Aeut)
#' Aeut.gc <- new('genclone', Aeut, other(Aeut)$population_hierarchy)
#' 
#' # Notice that there are two hierarchies, Pop and Subpop
#' Aeut.gc 
#' 
#' # Currently set on just Pop
#' head(pop(Aeut.gc)) 
#' 
#' # setting the hierarchy to both Pop and Subpop
#' setpop(Aeut.gc) <- ~Pop/Subpop 
#' head(pop(Aeut.gc))
#' 
#' \dontrun{
#' 
#' # Can be used to create objects as well.
#' Aeut.old <- setpop(Aeut.gc, ~Pop) 
#' head(pop(Aeut.old))
#' }
#==============================================================================#
setpop <- function(x, formula = NULL) standardGeneric("setpop")

#' @export
setGeneric("setpop")

setMethod(
  f = "setpop",
  signature(x = "genclone"),
  definition = function(x, formula = NULL){
    if (is.null(formula) | !is.language(formula)){
      stop(paste(substitute(formula), "must be a valid formula object."))
    }
    vars <- all.vars(formula)
    if (!all(vars %in% names(x@hierarchy))){
      stop(hier_incompatible_warning(vars, x@hierarchy))
    }
    pop(x) <- make_hierarchy(formula, x@hierarchy)[[length(vars)]]
    return(x)
  })

#==============================================================================#
#' @export
#' @rdname population-methods
#' @aliases setpop<-,genclone-method
#' @docType methods
#==============================================================================#
"setpop<-" <- function(x, value) standardGeneric("setpop<-")

#' @export
setGeneric("setpop<-")

setMethod(
  f = "setpop<-",
  signature(x = "genclone"),
  definition = function(x, value){
    return(setpop(x, value))
  })
