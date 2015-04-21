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
#' @keywords internal
#' @author Zhian N. Kamvar
#==============================================================================#
setMethod(
  f = "[",
  signature(x = "bootgen"),
  definition = function(x, i, j, ..., drop = FALSE){
    if (missing(i)) i <- TRUE
    if (missing(j)) j <- TRUE
    loc <- dim(x)[2]
    if (length(j) > loc | any(j > loc)){
      stop('subscript out of bounds')
    }
    
    # Taking Names
    locnall         <- x@loc.nall[j]
    allnames        <- x@all.names[j]
    names(allnames) <- names(x@all.names)[1:length(j)]
    names(locnall)  <- names(allnames)
    alllist         <- slot(x, "alllist")[j]
    indices         <- unlist(alllist)
    locnames        <- rep(names(allnames), locnall)
    tabnames        <- paste(locnames, unlist(allnames), sep = ".")
    res             <- slot(x, "tab")[i, indices, drop = drop]
    colnames(res)   <- tabnames
    
    ## Resetting all factors that need to be set. 
    slot(x, "tab")       <- res
    slot(x, "loc.fac")   <- factor(locnames, names(allnames))
    slot(x, "loc.names") <- names(allnames)
    slot(x, "loc.nall")  <- locnall
    slot(x, "all.names") <- allnames
    slot(x, "alllist")   <- .Call("expand_indices", cumsum(locnall), length(j), PACKAGE = "poppr")
    slot(x, "names")     <- slot(x, "names")[i]
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
    return(c(length(slot(x, "names")), length(slot(x, "loc.names"))))
  }
)

#==============================================================================#
#' @rdname bootgen-methods
#==============================================================================#
setMethod(
  f = "$",
  signature(x = "bootgen"),
  definition = function(x, name){
    return(slot(x, name))
  }
)

#==============================================================================#
#' @rdname bootgen-methods
#' @param .Object a character, "bootgen"
#' @param gen a genind, genclone, or genpop object
#' @param na how missing data should be treated. Default is "mean", averaging 
#'   over other values in the matrix. Possible values are listed in 
#'   \code{\link[adegenet]{tab}}.
#' @param freq if \code{TRUE}, the matrix will be a genotype frequency matrix.
#'   If \code{FALSE}, the matrix will be allele counts.
#==============================================================================#
setMethod(
  f = "initialize",
  signature = "bootgen",
  definition = function(.Object, gen, na = "mean", freq = TRUE){
    if (missing(gen)){
      return(.Object)
    }
    if (is.genind(gen)){
      objnames <- slot(gen, "ind.names")
    } else if (is.genpop(gen)){
      objnames <- slot(gen, "pop.names")
    } else {
      stop("gen must be a valid genind or genpop object.")
    }
    num_alleles                <- slot(gen, "loc.nall")
    num_loci                   <- length(num_alleles)
    slot(.Object, "tab")       <- tab(gen, NA.method = na, freq = freq)     
    slot(.Object, "loc.fac")   <- slot(gen, "loc.fac")   
    slot(.Object, "loc.names") <- slot(gen, "loc.names") 
    slot(.Object, "loc.nall")  <- num_alleles  
    slot(.Object, "all.names") <- slot(gen, "all.names") 
    slot(.Object, "alllist")   <- .Call("expand_indices", cumsum(num_alleles), num_loci, PACKAGE = "poppr")
    slot(.Object, "names")     <- objnames
    slot(.Object, "type")      <- slot(gen, "type")
    slot(.Object, "ploidy")    <- as.integer(slot(gen, "ploidy"))
    return(.Object)
  })

setMethod(
  f = "tab",
  signature = "bootgen",
  definition = function(x, freq = TRUE){ # freq is a dummy argument
    return(x@tab)
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
#' @keywords internal
#' @author Zhian N. Kamvar
#==============================================================================#
setMethod(
  f = "initialize",
  signature = "bruvomat",
  definition = function(.Object, gen, replen){
    if (missing(gen)){
      return(.Object)
    }
    if (missing(replen)){
      replen <- vapply(gen@all.names, function(y) guesslengths(as.numeric(y)), 1)
    }
    ploid <- max(ploidy(gen))
    # This controlls for the user correcting missing data using "mean". 
    # if (any(!gen@tab %in% c((0:ploid)/ploid, NA))){
    #   gen@tab[!gen@tab %in% c((0:ploid)/ploid, NA)] <- NA
    # }
    # # This will check for data that has missing scored as "zero".
    # popcols <- ploid*nLoc(gen)
    # if (!any(is.na(gen@tab)) & any(rowSums(gen@tab, na.rm=TRUE) < nLoc(gen))){
    #   mat1 <- as.matrix.data.frame(genind2df(gen, sep="/", usepop=FALSE))
    #   mat1[mat1 %in% c("", NA)] <- paste(rep(0, ploid), collapse="/")
    #   mat2 <- apply(mat1, 1, strsplit, "/")
    #   mat3 <- apply(as.matrix(t(sapply(mat2, unlist))), 2, as.numeric)
    #   vec1 <- suppressWarnings(as.numeric(unlist(mat3)))
    #   pop  <- matrix(vec1, nrow=nInd(gen), ncol=popcols)
    # } else {
    #   popdf <- genind2df(gen, oneColPerAll=TRUE, usepop=FALSE)
    #   mat1  <- as.matrix.data.frame(popdf)
    #   pop   <- suppressWarnings(matrix(as.numeric(mat1), ncol=popcols))
    # }
    popdf <- genind2df(gen, sep = "/", usepop = FALSE)
    mat   <- generate_bruvo_mat(popdf, maxploid = ploid, sep = "/", mat = TRUE)
    mat[is.na(mat)] <- 0
    slot(.Object, "mat")       <- mat
    slot(.Object, "replen")    <- replen
    slot(.Object, "ploidy")    <- ploid
    slot(.Object, "ind.names") <- indNames(gen)
    return(.Object)
  }
)

#==============================================================================#
#' @rdname bruvomat-methods
#' @keywords internal
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
#' @keywords internal
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
# SNPCLONE METHODS
#------------------------------------------------------------------------------#
################################################################################
#==============================================================================#
#' Check for validity of a snpclone object
#' 
#' @note a \linkS4class{snpclone} object will always be a valid 
#' \linkS4class{genlight} object.
#' 
#' @export
#' @rdname is.snpclone
#' @param x a snpclone object 
#' @author Zhian N. Kamvar
#' @examples
#' (x <- as.snpclone(glSim(100, 1e3, ploid=2)))
#' is.snpclone(x)
#==============================================================================#
is.snpclone <- function(x){
  res <- (is(x, "snpclone"))
  return(res)
}

#==============================================================================#
#' Methods used for the snpclone object
#' 
#' Default methods for subsetting snpclone objects. 
#' 
#' @rdname snpclone-method
#' @param x a snpclone object
#' @param i vector of numerics indicating number of individuals desired
#' @param j a vector of numerics corresponding to the loci desired.
#' @param ... passed on to the \code{\linkS4class{genlight}} object.
#' @param drop set to \code{FALSE} 
#' @author Zhian N. Kamvar
#==============================================================================#
setMethod(
  f = "[",
  signature(x = "snpclone", i = "ANY", j = "ANY", drop = "ANY"),
  definition = function(x, i, j, ..., drop = FALSE){
    if (missing(i)) i <- TRUE
    ismlgclass <- "MLG" %in% class(x@mlg)

    if (ismlgclass){
      mlg <- x@mlg[i, all = TRUE]
    } else {
      mlg <- x@mlg[i]
    }
    x <- callNextMethod(x = x, i = i, j = j, ..., drop = drop)
    x@mlg <- mlg
    return(x)
  })

#==============================================================================#
#' @rdname snpclone-method
#' @param .Object a character, "snpclone"
#' @param mlg a vector where each element assigns the multilocus genotype of
#' that individual in the data set. 
#' @param mlgclass a logical value specifying whether or not to translate the
#' mlg object into an MLG class object. 
#' @author Zhian N. Kamvar
#' @keywords internal
#==============================================================================#
setMethod(      
  f = "initialize",
  signature("snpclone"),
  definition = function(.Object, ..., mlg, mlgclass = TRUE){

    .Object <- callNextMethod(.Object, ..., mlg)
    if (missing(mlg)){
      mlg <- mlg.vector(.Object)
    }
    if (mlgclass){
      mlg <- new("MLG", mlg)
    }
    slot(.Object, "mlg") <- mlg
    return(.Object)
  }
)

#==============================================================================#
#' @rdname snpclone-method
#' @param object a snpclone object
#==============================================================================#
setMethod(
  f = "show", 
  signature("snpclone"),
  definition = function(object){
    callNextMethod()
    cat(length(unique(object@mlg[])), "multilocus genotypes")
  }
)
#==============================================================================#
#' Create a snpclone object from a genlight object.
#' 
#' Wrapper for snpclone initializer.
#' 
#' @export
#' @rdname snpclone-coercion-methods
#' @aliases as.snpclone,genlight-method
#' @param x a \code{\linkS4class{genlight}} or \code{\linkS4class{snpclone}} 
#'   object
#' @param ... arguments to be passed on to the genlight constructor. These are
#'   not used if x is not missing.
#' @param parallel should the parallel package be used to construct the object?
#' @param n.cores how many cores should be utilized? See documentation for
#'   \code{\linkS4class{genlight}} for details.
#' @param mlg a vector of multilocus genotypes or an object of class MLG for the
#'   new snpclone object.
#' @param mlgclass if \code{TRUE} (default), the multilocus genotypes will be
#'   represented as an \code{\linkS4class{MLG}} object.
#'   
#' @docType methods
#'   
#' @author Zhian N. Kamvar
#' @examples
#' (x <- as.snpclone(glSim(100, 1e3, ploid=2)))
#==============================================================================#
as.snpclone <- function(x, ..., parallel = require('parallel'), n.cores = NULL, 
                        mlg, mlgclass = TRUE){
  standardGeneric("as.snpclone")
}

#' @export
setGeneric("as.snpclone")


setMethod(
  f = "as.snpclone",
  signature(x = "genlight"),
  definition = function(x, ..., parallel = require('parallel'), n.cores = NULL, 
                        mlg, mlgclass = TRUE){
    if (!missing(x) && is(x, "genlight")){
      if (missing(mlg)) mlg <- mlg.vector(x)
      res <- new("snpclone", 
                 gen = x@gen, 
                 n.loc = nLoc(x), 
                 ind.names = indNames(x), 
                 loc.names = locNames(x), 
                 loc.all = alleles(x), 
                 chromosome = chr(x), 
                 position = position(x), 
                 strata = strata(x), 
                 hierarchy = x@hierarchy, 
                 ploidy = ploidy(x), 
                 pop = pop(x), 
                 other = other(x), 
                 parallel = parallel,
                 n.cores = n.cores,
                 mlg = mlg,
                 mlgclass = mlgclass)
    } else {
      res <- new("snpclone", ..., mlg = mlg, mlgclass = mlgclass)
    }
    return(res)
  })
################################################################################
#------------------------------------------------------------------------------#
# GENCLONE METHODS
#------------------------------------------------------------------------------#
################################################################################
#==============================================================================#
#' Check for validity of a genclone object
#' 
#' @note a \linkS4class{genclone} object will always be a valid 
#' \linkS4class{genind} object.
#' 
#' @export
#' @rdname is.genclone
#' @param x a genclone object 
#' @author Zhian N. Kamvar
#' @examples
#' data(nancycats)
#' nanclone <- as.genclone(nancycats)
#' is.genclone(nanclone)
#==============================================================================#
is.genclone <- function(x){ 
  res <- (is(x, "genclone"))
  return(res)
}
#==============================================================================#
#' @rdname genclone-method
#' @param .Object a character, "genclone"
#' @param gen \code{"\linkS4class{genind}"} object
#' @param mlg a vector where each element assigns the multilocus genotype of 
#'   that individual in the data set.
#' @param mlgclass a logical value specifying whether or not to translate the 
#'   mlg object into an MLG class object.
#' @keywords internal
#==============================================================================#
setMethod(      
  f = "initialize",
  signature("genclone"),
  definition = function(.Object, ..., mlg, mlgclass = TRUE){

    .Object <- callNextMethod(.Object, ..., mlg = mlg)
    if (missing(mlg)){
      mlg <- mlg.vector(.Object)      
    } 
    if (mlgclass) {
      mlg <- new("MLG", mlg)
    }
    slot(.Object, "mlg") <- mlg
    return(.Object)
  }
)

#==============================================================================#
#' Methods used for the genclone object
#' 
#' Default methods for subsetting genclone objects. 
#' 
#' @rdname genclone-method
#' @param x a genclone object
#' @param i vector of numerics indicating number of individuals desired
#' @param j a vector of numerics corresponding to the loci desired.
#' @param ... passed on to the \code{\linkS4class{genind}} object.
#' @param drop set to \code{FALSE}
#' @param mlg.reset logical. Defaults to \code{FALSE}. If \code{TRUE}, the mlg
#'   vector will be reset
#' @param loc passed on to \code{\linkS4class{genind}} object.
#' @param treatOther passed on to \code{\linkS4class{genind}} object.
#' @param quiet passed on to \code{\linkS4class{genind}} object.
#' @author Zhian N. Kamvar
#==============================================================================#
setMethod(
  f = "[",
  signature(x = "genclone", i = "ANY", j = "ANY", drop = "ANY"),
  definition = function(x, i, j, ..., mlg.reset = FALSE, drop = FALSE){
    if (missing(i)) i <- TRUE
    x     <- callNextMethod(x = x, i = i, j = j, ..., drop = drop)
    if (!mlg.reset){
      x@mlg <- x@mlg[i]
    } else {
      if (class(x@mlg)[1] != "MLG"){
        x@mlg <- mlg.vector(x, reset = TRUE)        
      } else {
        x@mlg <- new("MLG", mlg.vector(x, reset = TRUE))
      }
    }
    return(x)
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
    ploid  <- c("ha", "di", "tri", "tetra", "penta", "hexa", "hepta", "octa",
      "nona", "deca", "hendeca", "dodeca")
    ploid  <- paste0(unique(ploid[object@ploidy]), "ploid")
    if (length(ploid) > 1){
      ploid  <- paste(ploid, paste0("(", table(object@ploidy), ")"))
      ploid1 <- paste0(ploid[-length(ploid)], collapse = ", ")
      sep    <- ifelse(length(ploid) > 2, ", ", " ")
      ploid  <- paste0(ploid1, sep, "and ", ploid[length(ploid)])
    }
    nind   <- nInd(object)     
    type   <- ifelse(object@type == "PA", "dominant", "codominant")
    nmlg   <- length(unique(object@mlg))
    nloc   <- nLoc(object)
    npop   <- ifelse(is.null(object@pop), 0, nPop(object))
    strata <- length(object@strata)
    chars  <- nchar(c(nmlg, nind, nloc, strata, 1, npop))
    ltab   <- max(chars) - chars
    ltab   <- vapply(ltab, function(x) substr("       ", 1, x+1), character(1))
    pops   <- popNames(object)
    poplen <- length(pops)
    if (poplen > 7) 
      pops <- c(pops[1:3], "...", pops[(poplen-2):poplen])
    stratanames <- names(object@strata)
    stratalen   <- length(stratanames)
    if (stratalen > 7) 
      stratanames <- c(stratanames[1:3], "...", stratanames[(stratalen-2):stratalen])
    cat("\nThis is a genclone object\n")
    cat("-------------------------\n")
    cat("Genotype information:\n\n",
      ltab[1], nmlg, "multilocus genotypes\n",
      ltab[2], nind, ploid, "individuals\n", 
      ltab[3], nloc, type, "loci\n\n"
      )
    popstrata <- ifelse(strata > 1, "strata -", "stratification -")
    if (strata == 0) popstrata <- "strata."
    popdef  <- ifelse(npop > 0, "defined -", "defined.")
    cat("Population information:\n\n")
    cat("", ltab[4], strata, popstrata, stratanames, fill = TRUE)
    if (!is.null(object@hierarchy)){
      cat("", ltab[5], "1", "hierarchy -", paste(object@hierarchy, collapse = ""), 
          fill = TRUE)
    }
    cat("", ltab[6], npop, "populations", popdef, pops, fill = TRUE)
    
  })

setGeneric("print")
#==============================================================================#
#' @rdname genclone-method
#' @export
#' @param x a genclone object
#' @param fullnames \code{logical}. If \code{TRUE}, then the full names of the
#'   populations will be printed. If \code{FALSE}, then only the first and last
#'   three population names are displayed.
#==============================================================================#
setMethod(
  f = "print",
  signature("genclone"),
  definition = function(x, ...){
    ploid  <- c("ha", "di", "tri", "tetra", "penta", "hexa", "hepta", "octa",
      "nona", "deca", "hendeca", "dodeca")
    ploid  <- paste0(unique(ploid[x@ploidy]), "ploid")
    if (length(ploid) > 1){
      ploid  <- paste(ploid, paste0("(", table(x@ploidy), ")"))
      ploid1 <- paste0(ploid[-length(ploid)], collapse = ", ")
      sep    <- ifelse(length(ploid) > 2, ", ", " ")
      ploid  <- paste0(ploid1, sep, "and ", ploid[length(ploid)])
    }
    nind   <- nInd(x)     
    type  <- ifelse(x@type == "PA", "dominant", "codominant")
    nmlg  <- length(unique(x@mlg))
    nloc  <- nLoc(x)
    npop  <- ifelse(is.null(x@pop), 0, nPop(x))
    strata  <- length(x@strata)
    chars <- nchar(c(nmlg, nind, nloc, strata, 1, npop))
    ltab  <- max(chars) - chars
    ltab  <- vapply(ltab, function(x) substr("       ", 1, x + 1), character(1))
    pops  <- popNames(x)
    stratanames <- names(x@strata)
    cat("\nThis is a genclone object\n")
    cat("-------------------------\n")
    cat("Genotype information:\n\n",
      ltab[1], nmlg, "multilocus genotypes\n",
      ltab[2], nind, ploid, "individuals\n", 
      ltab[3], nloc, type, "loci\n\n",
      fill = TRUE)
    popstrata <- ifelse(strata > 1, "strata -", "stratification -")
    if (strata == 0) popstrata <- "strata."
    popdef  <- ifelse(npop > 0, "defined -", "defined.")
    cat("Population information:\n\n")
    cat("", ltab[4], strata, popstrata, stratanames, fill = TRUE)
    if (!is.null(x@hierarchy)){
      cat("", ltab[5], "1", "hierarchy -", paste(x@hierarchy, collapse = ""), 
          fill = TRUE)
    }
    cat("", ltab[6], npop, "populations", popdef, pops, fill = TRUE)
  })

#==============================================================================#
#' Create a genclone object from a genind object.
#' 
#' Wrapper for genclone initializer.
#' 
#' @export
#' @rdname coercion-methods
#' @aliases as.genclone,genind-method
#' @param x a \code{\linkS4class{genind}} or \code{\linkS4class{genclone}} 
#'   object
#' @param ... arguments passed on to the \code{\linkS4class{genind}} constructor
#' @param mlg an optional vector of multilocus genotypes as integers
#' @param mlgclass should the mlg slot be of class MLG?
#' @docType methods
#'   
#' @seealso \code{\link{splitStrata}}, \code{\linkS4class{genclone}},
#'   \code{\link{read.genalex}}
#' @author Zhian N. Kamvar
#' @examples
#' data(Aeut)
#' Aeut
#' Aeut.gc <- as.genclone(Aeut)
#' Aeut.gc
#' Aeut.gc <- as.genclone(Aeut)
#' Aeut.gc
#==============================================================================#
as.genclone <- function(x, ..., mlg, mlgclass = TRUE){
  standardGeneric("as.genclone")
}

#' @export
setGeneric("as.genclone")


setMethod(
  f = "as.genclone",
  signature(x = "genind"),
  definition = function(x, ..., mlg, mlgclass = TRUE){
    theCall <- match.call()
    if (!missing(x) && is.genind(x)){
      theOther <- x@other
      if (missing(mlg)) mlg <- mlg.vector(x)
      res <- new("genclone", tab = tab(x), ploidy = ploidy(x), pop = pop(x), 
                 type = x@type, prevcall = theCall, strata = x@strata, 
                 hierarchy = x@hierarchy, mlg = mlg, mlgclass = mlgclass)
      res@other <- theOther
    } else {
      res <- new("genclone", ..., mlg = mlg, mlgclass = mlgclass)
    }
    return(res)
  })

#==============================================================================#
# Seploc method for genclone objects. 
# Note that this is not the most efficient, but it is difficult to work around
# the method for genind objects as they are forced to be genind objects later.
#==============================================================================#
setMethod(
  f = "seploc",
  signature(x = "genclone"),
  definition = function(x, truenames = TRUE, res.type = c("genind", "matrix")){
    ARGS <- c("genind", "matrix")
    res.type <- match.arg(res.type, ARGS)
    if (res.type == "matrix"){
      splitsville <- split(colnames(x@tab), x@loc.fac)
      listx       <- lapply(splitsville, function(i) x@tab[, i, drop = FALSE])
    } else {
      listx <- lapply(locNames(x), function(i) x[loc = i])
    }
    names(listx) <- locNames(x)
    return(listx)
  })



#==============================================================================#
#' Access and manipulate multilocus lineages.
#' 
#' The following methods allow the user to access and manipulate multilocus 
#' lineages in genclone or snpclone objects.
#' 
#' @export
#' @param x a \linkS4class{genclone} or \linkS4class{snpclone} object.
#' @param type a character specifying "original", "contracted", or "custom"
#'   defining they type of mlgs to return. Defaults to what is set in the
#'   object.
#' @param value a character specifying which mlg type is visible in the object.
#'   See details.
#'   
#' @return an object of the same type as x.
#' 
#' @details \linkS4class{genclone} and \linkS4class{snpclone} objects have a
#'   slot for an internal class of object called \linkS4class{MLG}. This class
#'   allows the storage of flexible mll definitions: \itemize{ \item "original"
#'   - naive mlgs defined by string comparison. This is default. \item
#'   "contracted" - mlgs defined by a genetic distance threshold. \item "custom"
#'   - user-defined MLGs }
#'   
#' @rdname mll-method
#' @aliases mll,genclone-method mll,snpclone-method
#' @docType methods
#' @author Zhian N. Kamvar
#' @seealso \code{\link{mll.custom}} \code{\link{mlg.table}}
#' @examples
#' 
#' data(partial_clone)
#' pc <- as.genclone(partial_clone)
#' mll(pc)
#' mll(pc) <- "custom"
#' mll(pc)
#' mll.levels(pc) <- LETTERS
#' mll(pc)
#==============================================================================#
mll <- function(x, type = NULL) standardGeneric("mll")

#' @export
setGeneric("mll")

mll.internal <- function(x, type = NULL){
  mlg <- x@mlg
  if (!"MLG" %in% class(mlg)){
    return(mlg)
  }
  if (!is.null(type)){
    TYPES <- c("original", "expanded", "contracted", "custom")
    type <- match.arg(type, TYPES)
  } else {
    type <- mlg@visible
  }
  return(mlg[, type])
}

setMethod(
  f = "mll",
  signature(x = "genclone"),
  definition = function(x, type = NULL){
    mll.internal(x, type)
  })

setMethod(
  f = "mll",
  signature(x = "snpclone"),
  definition = function(x, type = NULL){
    mll.internal(x, type)
  })

setMethod(
  f = "mll",
  signature(x = "snpclone"),
  definition = function(x, type = NULL){
    mlg <- x@mlg
    if (!"MLG" %in% class(mlg)){
      return(mlg)
    }
    if (!is.null(type)){
      TYPES <- c("original", "expanded", "contracted", "custom")
      type <- match.arg(type, TYPES)
    } else {
      type <- mlg@visible
    }
    return(mlg[, type])
  })

#==============================================================================#
#' @export
#' @rdname mll-method
#' @aliases mll<-,genclone-method mll<-,snpclone-method
#' @docType methods
#==============================================================================#
"mll<-" <- function(x, value) standardGeneric("mll<-")

#' @export
setGeneric("mll<-")

setMethod(
  f = "mll<-",
  signature(x = "genclone"),
  definition = function(x, value){
    TYPES <- c("original", "expanded", "contracted", "custom")
    value <- match.arg(value, TYPES)
    x@mlg@visible <- value
    return(x)
  })

setMethod(
  f = "mll<-",
  signature(x = "snpclone"),
  definition = function(x, value){
    TYPES <- c("original", "expanded", "contracted", "custom")
    value <- match.arg(value, TYPES)
    x@mlg@visible <- value
    return(x)
  })

#==============================================================================#
#' Define custom multilocus lineages
#' 
#' This function will allow you to define custom multilocus lineages for your
#' data set.
#' 
#' @export
#' @param x a \linkS4class{genclone} or \linkS4class{snpclone} object.
#' @param set logical. If \code{TRUE} (default), the visible mlls will be set to
#' 'custom'. 
#' @param value a vector that defines the multilocus lineages for your data.
#' This can be a vector of ANYTHING that can be turned into a factor. 
#' 
#' @return an object of the same type as x
#' @rdname mll.custom
#' @aliases mll.custom,genclone-method mll.custom,snpclone-method
#' @docType methods
#' @author Zhian N. Kamvar
#' @seealso \code{\link{mll}} \code{\link{mlg.table}}
#' @examples 
#' data(partial_clone)
#' pc <- as.genclone(partial_clone)
#' mll.custom(pc) <- LETTERS[mll(pc)]
#' mll(pc)
#' 
#' # Let's say we had a mistake and the A mlg was actually B. 
#' mll.levels(pc)[mll.levels(pc) == "A"] <- "B"
#' mll(pc)
#' 
#' # Set the MLL back to the original definition.
#' mll(pc) <- "original"
#' mll(pc)
#==============================================================================#
mll.custom <- function(x, set = TRUE, value) standardGeneric("mll.custom")

#' @export
setGeneric("mll.custom")

setMethod(
  f = "mll.custom",
  signature(x = "genclone"),
  definition = function(x, set = TRUE, value){
    mll.custom.internal(x, set, value)
  })

setMethod(
  f = "mll.custom",
  signature(x = "snpclone"),
  definition = function(x, set = TRUE, value){
    mll.custom.internal(x, set, value)
  })

mll.custom.internal <- function(x, set = TRUE, value){
  if (!is(x@mlg, "MLG")){
    x@mlg <- new("MLG", x@mlg)
  }
  if (missing(value)){
    return(x@mlg[, "custom"])
  }
  mlgs <- x@mlg
  if (length(value) != length(mlgs)){
    stop("value must be the same length as the mlls")
  }
  if (!is.factor(value)){
    value <- factor(value)
  }
  if (set){
    mlgs@visible <- "custom"
  }
  mlgs@mlg[, "custom"] <- value
  x@mlg <- mlgs
  return(x)
}

#' @export
#' @rdname mll.custom
#' @aliases mll.custom<-,genclone-method mll.custom<-,snpclone-method
#' @docType methods
"mll.custom<-" <- function(x, set = TRUE, value) standardGeneric("mll.custom<-")

#' @export
setGeneric("mll.custom<-")

setMethod(
  f = "mll.custom<-",
  signature(x = "genclone"),
  definition = function(x, set = TRUE, value){
    mll.custom.internal(x, set, value)
  })

setMethod(
  f = "mll.custom<-",
  signature(x = "snpclone"),
  definition = function(x, set = TRUE, value){
    mll.custom.internal(x, set, value)
  })

#' @export
#' @rdname mll.custom
#' @aliases mll.levels,genclone-method mll.levels,snpclone-method
#' @docType methods
mll.levels <- function(x, set = TRUE, value) standardGeneric("mll.levels")

#' @export
setGeneric("mll.levels")

setMethod(
  f = "mll.levels",
  signature(x = "genclone"),
  definition = function(x, set = TRUE, value){
    mll.levels.internal(x, set, value)
  }
)

setMethod(
  f = "mll.levels",
  signature(x = "snpclone"),
  definition = function(x, set = TRUE, value){
    mll.levels.internal(x, set, value)
  }
)

mll.levels.internal <- function(x, set = TRUE, value){
  if (!is(x@mlg, "MLG")){
    x@mlg <- new("MLG", x)
  }
  mlgs <- x@mlg
  if (missing(value)){
    return(levels(mlgs))
  }
  if (length(value) != nlevels(mlgs@mlg[, "custom"])){
    stop("value length should match the number of values in mlg")
  }
  if (set){
    vis <- "custom"
  } else {
    vis <- mlgs@visible
  }
  mlgs@visible <- "custom"
  levels(mlgs) <- value
  x@mlg        <- mlgs
  mll(x)       <- vis
  return(x)
}

#' @export
#' @rdname mll.custom
#' @aliases mll.levels<-,genclone-method mll.levels<-,snpclone-method
#' @docType methods
"mll.levels<-" <- function(x, set = TRUE, value) standardGeneric("mll.levels<-")

#' @export
setGeneric("mll.levels<-")

setMethod(
  f = "mll.levels<-",
  signature(x = "genclone"),
  definition = function(x, set = TRUE, value){
    mll.levels.internal(x, set, value)
  }
)

setMethod(
  f = "mll.levels<-",
  signature(x = "snpclone"),
  definition = function(x, set = TRUE, value){
    mll.levels.internal(x, set, value)
  }
)
