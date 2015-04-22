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
    # Shuffling
    # matcols  <- apply(getinds(x@loc.nall), 1, function(ind) ind[1]:ind[2])[j]
    # matcols  <- .Call("expand_indices", cumsum(x@loc.nall), nLoc(x))[j]
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
#==============================================================================#
setMethod(
  f = "initialize",
  signature = "bootgen",
  definition = function(.Object, gen){
    if (missing(gen)){
      stop("gen must be specified.")
    }
    if (is.genind(gen)){
      objnames <- slot(gen, "ind.names")
      tab      <- slot(gen, "tab")
    } else if (is.genpop(gen)){
      objnames <- slot(gen, "pop.names")
      tab      <- makefreq(gen, missing = "mean", quiet = TRUE)$tab
    } else {
      stop("gen must be a valid gen object.")
    }
    num_alleles                <- slot(gen, "loc.nall")
    num_loci                   <- length(num_alleles)
    slot(.Object, "tab")       <- tab       
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
    if (length(x@hierarchy) > 0){
      hierarchy <- x@hierarchy[i, , drop = FALSE]      
    } else {
      hierarchy <- x@hierarchy
    }


    x <- callNextMethod(x = x, i = i, j = j, ..., drop = drop)
    if (!"snpclone" %in% class(x)){
      x <- new("snpclone", x, hierarchy, mlg, mlgclass = ismlgclass)
    } else {
      x@mlg <- mlg
      x@hierarchy <- hierarchy      
    }

    return(x)
  })

#==============================================================================#
#' @rdname snpclone-method
#' @param .Object a character, "snpclone"
#' @param gen \code{"\linkS4class{genlight}"} object
#' @param hierarchy a data frame where each row i represents the different
#' population assignments of individual i in the data set. If this is empty, the
#' hierarchy will be created from the population factor.
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
  definition = function(.Object, gen, hierarchy, mlg, mlgclass = TRUE){
    if (missing(gen)){
      gen <- new("genlight")
    }
    .Object <- callNextMethod(.Object)
    invisible(lapply(slotNames(gen), function(x) slot(.Object, x) <<- slot(gen, x)))

    if (missing(mlg)){
      mlg <- mlg.vector(gen)      
    }
    if (missing(hierarchy)){
      hierarchy <- data.frame()
    }
    if (mlgclass){
      mlg <- new("MLG", mlg)
    }
    slot(.Object, "mlg")       <- mlg
    slot(.Object, "hierarchy") <- hierarchy
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
    cat(" --- snpclone contents ---\n")
    if (length(x@hierarchy) > 0){
      hiernames <- names(x@hierarchy)
      hierlen <- length(hiernames)
      nameshow <- ifelse(hierlen > 6, 
                         c(head(hiernames, 3), "...", tail(hiernames, 3)), 
                         hiernames)
      nameshow <- paste(nameshow, collapse = ", ")
      cat(" @hierarchy:", "a data frame with", hierlen, 
          "levels: (", nameshow, ")\n")
    }
    mlgtype <- ifelse(is(object@mlg, "MLG"), paste0(object@mlg@visible, " "), "")
    mlgtype <- paste0(mlgtype, "multilocus genotypes")
    cat(" @mlg:", length(unique(object@mlg[])), mlgtype)
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
#' @param hierarchy a data frame representing the population hierarchy.
#' @docType methods
#'   
#' @note The hierarchy must have the same number of rows as the number of 
#'   observations in the genlight object. If no hierarchy is defined, the function
#'   will search for a data frame in the \code{\link{other}} slot called 
#'   "population_hierarchy" and set that as the hieararchy. If none is defined,
#'   the population will be set as the hierarchy under the label "Pop". Use the 
#'   function \code{\link{splithierarchy}} to split up any population 
#'   hierarchies that might be combined in the population factor.
#'   
#' @author Zhian N. Kamvar
#' @examples
#' (x <- as.snpclone(glSim(100, 1e3, ploid=2)))
#==============================================================================#
as.snpclone <- function(x, hierarchy = NULL){
  standardGeneric("as.snpclone")
}

#' @export
setGeneric("as.snpclone")


setMethod(
  f = "as.snpclone",
  signature(x = "genlight"),
  definition = function(x, hierarchy){
    if (missing(hierarchy)){
      if ("population_hierarchy" %in% names(other(x))){
        hierarchy   <- other(x)[["population_hierarchy"]]
        newsnpclone <- new("snpclone", x, hierarchy)
      } else {
        newsnpclone <- new("snpclone", x)
      }
    } else {
      newsnpclone   <- new("snpclone", x, hierarchy)
    }
    return(newsnpclone)
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
#' @param loc passed on to \code{\linkS4class{genind}} object.
#' @param treatOther passed on to \code{\linkS4class{genind}} object.
#' @param quiet passed on to \code{\linkS4class{genind}} object. 
#' @author Zhian N. Kamvar
#==============================================================================#
setMethod(
  f = "[",
  signature(x = "genclone", i = "ANY", j = "ANY", drop = "ANY"),
  definition = function(x, i, j, ..., loc=NULL, treatOther=TRUE, quiet=TRUE, drop = FALSE){
    if (missing(i)) i <- TRUE
    if (missing(j)) j <- TRUE
    if (class(slot(x, "mlg")) %in% "MLG"){
      mlg <- slot(x, "mlg")[i, all = TRUE]
    } else {
      mlg <- slot(x, "mlg")[i]  
    }
    hierarchy <- slot(x, "hierarchy")[i, , drop = FALSE]
    ## The following is lifted directly from the adegenet source code as
    ## callNextMethod() was throwing the error:
    ##
    ## Error in callNextMethod() : bad object found as method (class "function")
    ##
    ## The reason for this is the fact that the `genind` function calls
    ## `new("genind")` within the function. Because of this, using
    ## `callNextMethod()` returns a genind object instead of modifying the
    ## genclone object as it should. 
    pop <- NULL
    if (is.null(x@pop)) { 
      tab <- truenames(x) 
    } else {
      temp <- truenames(x)
      tab  <- temp$tab
      pop  <- temp$pop
      pop  <- factor(pop[i])
    }
    nrowx     <- nrow(x@tab)
    old.other <- other(x)

    ## handle loc argument
    if(!is.null(loc)){
      loc  <- as.character(loc)
      temp <- !loc %in% x@loc.fac
      if (any(temp)) { # si mauvais loci
        noloc <- paste(loc[temp], collapse = " ")
        warning(paste("the following specified loci do not exist:", noloc))
      } 
      j    <- x$loc.fac %in% loc
    } # end loc argument

    prevcall <- x@call
    tab      <- tab[i, j, ..., drop=FALSE]

    if(drop){
      allNb  <- apply(tab, 2, sum, na.rm=TRUE) # allele absolute frequencies
      toKeep <- (allNb > 1e-10)
      tab    <- tab[ , toKeep, drop=FALSE]
    }

    res <- genind(tab, pop=pop, prevcall=prevcall, ploidy=x@ploidy, type=x@type)
    res <- new("genclone", res, hierarchy, mlg, mlgclass = FALSE)
  
    ## handle 'other' slot
    nOther     <- length(x@other)
    namesOther <- names(x@other)
    counter    <- 0
    if (treatOther){
      f1 <- function(obj, n = nrowx){
        counter <<- counter + 1
        if (!is.null(dim(obj)) && nrow(obj) == n){ # if the element is a matrix-like obj
          obj <- obj[i , , drop=FALSE]
        } else if (length(obj) == n){ # if the element is not a matrix but has a length == n
          obj <- obj[i]
          if (is.factor(obj)){
            obj <- factor(obj)
          }
        } else {
          if (!quiet){
            warning(paste("cannot treat the object", namesOther[counter]))
          } 
        }
        return(obj)
      } # end f1

      res@other <- lapply(x@other, f1) # treat all elements
    } else {
      res@other <- old.other
    } # end treatOther

    return(res)
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
#' @param mlgclass a logical value specifying whether or not to translate the
#' mlg object into an MLG class object. 
#' @keywords internal
#==============================================================================#
setMethod(      
  f = "initialize",
  signature("genclone"),
  definition = function(.Object, gen, hierarchy, mlg, mlgclass = TRUE){
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
    } else {
      hierarchy <- data.frame(lapply(hierarchy, function(f) factor(f, unique(f))))
    }

    # No 'initialize' method for genind objects...
    lapply(names(gen), function(y) slot(.Object, y) <<- slot(gen, y))

    if (mlgclass) mlg <- new("MLG", mlg)

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
    ploid  <- c("ha", "di", "tri", "tetra", "penta", "hexa", "hepta", "octa",
      "nona", "deca", "hendeca", "dodeca")
    ploid  <- paste0(ploid[object@ploidy], "ploid")
    nind   <- nInd(object)
    type   <- ifelse(object@type == "PA", "dominant", "codominant")
    nmlg   <- length(unique(object@mlg))
    nloc   <- nLoc(object)
    npop   <- ifelse(is.null(object@pop), 0, length(object@pop.names))
    hier   <- length(object@hierarchy)
    chars  <- nchar(c(nmlg, nind, nloc, hier, npop))
    ltab   <- max(chars) - chars
    ltab   <- vapply(ltab, function(x) substr("       ", 1, x+1), character(1))
    pops   <- object@pop.names
    poplen <- length(pops)
    if (poplen > 7) 
      pops <- c(pops[1:3], "...", pops[(poplen-2):poplen])
    hiernames <- names(object@hierarchy)
    hierlen   <- length(hiernames)
    if (hierlen > 7) 
      hiernames <- c(hiernames[1:3], "...", hiernames[(hierlen-2):hierlen])
    cat("\nThis is a genclone object\n")
    cat("-------------------------\n")
    cat("Genotype information:\n\n",
      ltab[1], nmlg, "multilocus genotypes\n",
      ltab[2], nind, ploid, "individuals\n", 
      ltab[3], nloc, type, "loci\n\n"
      )
    pophier <- ifelse(hier > 1, "levels -", "level -")
    if (hier == 0) pophier <- "levels."
    popdef  <- ifelse(npop > 0, "defined -", "defined.")
    cat("Population information:\n\n")
    cat("", ltab[4], hier, "hierarchical", pophier, hiernames, fill = TRUE)
    cat("", ltab[5], npop, "populations", popdef, pops, fill = TRUE)
    
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
    ploid <- c("ha", "di", "tri", "tetra", "penta", "hexa", "hepta", "octa",
      "nona", "deca", "hendeca", "dodeca")
    ploid <- paste0(ploid[x@ploidy], "ploid")
    nind  <- nInd(x)
    type  <- ifelse(x@type == "PA", "dominant", "codominant")
    nmlg  <- length(unique(x@mlg))
    nloc  <- nLoc(x)
    npop  <- ifelse(is.null(x@pop), 0, length(x@pop.names))
    hier  <- length(x@hierarchy)
    chars <- nchar(c(nmlg, nind, nloc, hier, npop))
    ltab  <- max(chars) - chars
    ltab  <- vapply(ltab, function(x) substr("       ", 1, x + 1), character(1))
    pops  <- x@pop.names
    hiernames <- names(x@hierarchy)
    cat("\nThis is a genclone object\n")
    cat("-------------------------\n")
    cat("Genotype information:\n\n",
      ltab[1], nmlg, "multilocus genotypes\n",
      ltab[2], nind, ploid, "individuals\n", 
      ltab[3], nloc, type, "loci\n\n"
      )
    pophier <- ifelse(hier > 1, "levels -", "level -")
    if (hier == 0) pophier <- "levels."
    popdef  <- ifelse(npop > 0, "defined -", "defined.")
    cat("Population information:\n\n")
    cat("", ltab[4], hier, "hierarchical", pophier, hiernames, fill = TRUE)
    cat("", ltab[5], npop, "populations", popdef, pops, fill = TRUE)
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
#' @param hierarchy a data frame representing the population hierarchy.
#' @docType methods
#'   
#' @note The hierarchy must have the same number of rows as the number of 
#'   observations in the genind object. If no hierarchy is defined, the function
#'   will search for a data frame in the \code{\link{other}} slot called 
#'   "population_hierarchy" and set that as the hieararchy. If none is defined,
#'   the population will be set as the hierarchy under the label "Pop". Use the 
#'   function \code{\link{splithierarchy}} to split up any population 
#'   hierarchies that might be combined in the population factor.
#'   
#' @seealso \code{\link{splithierarchy}}, \code{\linkS4class{genclone}},
#'   \code{\link{read.genalex}}
#' @author Zhian N. Kamvar
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
# Seploc method for genclone objects. 
#==============================================================================#
setMethod(
  f = "seploc",
  signature(x = "genclone"),
  definition = function(x, ...){
    mlg       <- x@mlg
    hierarchy <- x@hierarchy
    listx     <- callNextMethod()
    if (is.genind(listx[[1]])){
      listx <- lapply(listx, function(gid) new("genclone", gid, hierarchy, mlg))
    }
    return(listx)
  })
#==============================================================================#
#' Access and manipulate the population hierarchy for genclone objects.
#' 
#' The following methods allow the user to quickly change the hierarchy or
#' population of a genclone object. 
#' 
#' @export 
#' @rdname hierarchy-methods
#' @aliases gethierarchy,genclone-method gethierarchy,snpclone-method
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
    gethierarchy.internal(x, formula, combine)
  })

setMethod(
  f = "gethierarchy",
  signature(x = "snpclone"),
  definition = function(x, formula = NULL, combine = TRUE){
    gethierarchy.internal(x, formula, combine)
  })

#==============================================================================#
#' @export
#' @rdname hierarchy-methods
#' @aliases sethierarchy<-,genclone-method sethierarchy<-,snpclone-method
#' @param value a data frame OR vector OR formula (see details).
#' @docType methods
#'   
#' @details \subsection{Function Specifics}{ \itemize{ \item
#' \strong{gethierarchy()} - This will retrieve the data from the
#' \emph{hierarchy} slot in the \linkS4class{genclone} object. You have the
#' option to choose specific heirarchical levels using a formula (see below) and
#' you can choose to combine the hierarchical levels (default) \item
#' \strong{sethierarchy()} - Set or reset the hierarchical levels in your
#' \linkS4class{genclone} object. \item \strong{namehierarchy()} - Rename the
#' hierarchical levels. \item \strong{splithierarchy()} - This is conceptually
#' similar to the default method of \code{\link{splitcombine}}. It is often
#' difficult to import files with several levels of hierarchy as most data
#' formats do not allow unlimited population levels. This is circumvented by
#' collapsing all hierarchical levels into a single population factor with a
#' common separator for each observation. This function will then split those
#' hierarchies for you, but it works best on a hierarchy that only has a single
#' column in it. See the rootrot example below. \item \strong{addhierarchy()} -
#' Add levels to your population hierarchy. If you have extra hierarchical
#' levels you want to add to your population hierarchy, you can use this method
#' to do so. You can input a data frame or a vector, but if you put in a vector,
#' you have the option to name it . }}
#' 
#' \subsection{Argument Specifics}{
#' 
#' These functions allow the user to seamlessly assign the hierarchical levels
#' of their \code{\linkS4class{genclone}} object. Note that there are two ways
#' of performing all methods (except for \code{gethierarchy()}). They
#' essentially do the same thing except that the assignment method (the one with
#' the "\code{<-}") will modify the object in place whereas the non-assignment 
#' method will not modify the original object. Due to convention, everything 
#' right of the assignment is termed \code{value}. To avoid confusion, here is a
#' guide to the inputs: \itemize{ \item \strong{sethierarchy()} This will be a 
#' \code{\link{data.frame}} that defines the hierarchy for each individual in 
#' the rows. \item \strong{namehierarchy()} This will be either a 
#' \code{\link{vector}} or a \code{\link{formula}} that will define the names. 
#' \item \strong{splithierarchy()} This will be a \code{\link{formula}} argument
#' with the same number of levels as the hierarchy you wish to split. \item 
#' \strong{addhierarchy()} This will be a \code{\link{vector}} or 
#' \code{\link{data.frame}} with the same length as the number of individuals in
#' your data. }}
#' 
#' \subsection{Details on Formulas}{
#' 
#' The preferred use of these functions is with a \code{\link{formula}} object. 
#' Specifically, a hierarchical formula argument is used to assign the levels of
#' the hierarchy. An example of a hierarchical formula would be:\cr 
#' \code{~Country/City/Neighborhood}\cr or \cr \code{~Country + Country:City + 
#' Country:City:Neighborhood}\cr of course, the first method is slightly easier 
#' to read. It is important to use hiearchical formulas when specifying 
#' hierarchies as other types of formulas (eg. 
#' \code{~Country*City*Neighborhood}) might give spurious results.}
#' 
#' @seealso \code{\link{setpop}} \code{\link{genclone}}
#'   \code{\link{as.genclone}}
#'   
#' @author Zhian N. Kamvar
#' @examples
#' # let's look at the microbov data set:
#' data(microbov)
#' microgc <- as.genclone(microbov)
#' microgc
#' 
#' # We see that we have three vectors of different names here. 
#' ?microbov
#' # These are Country, Breed, and Species
#' names(other(microgc))
#' 
#' # Let's set the hierarchy
#' sethierarchy(microgc) <- data.frame(other(microgc))
#' microgc
#' 
#' # And change the names so we know what they are
#' namehierarchy(microgc) <- ~Country/Breed/Species
#' 
#' # let's see what the hierarchy looks like by Species and Breed:
#' head(gethierarchy(microgc, ~Breed/Species))
#' 
#' \dontrun{
#' # Load our data set and convert it to a genclone object.
#' Aeut.gc <- read.genalex(system.file("files/rootrot.csv", package = "poppr"))
#' 
#' # we can see the hierarchy is set to Population_Subpopulation.
#' head(gethierarchy(Aeut.gc))
#' 
#' # We can use splithierarchy() to split them.
#' splithierarchy(Aeut.gc) <- ~Pop/Subpop
#' Aeut.gc
#' head(gethierarchy(Aeut.gc))
#' 
#' # We can also use gethierarchy to combine the hierarchy.
#' head(gethierarchy(Aeut.gc, ~Pop/Subpop))
#' 
#' # We can also give it a more descriptive name. 
#' namehierarchy(Aeut.gc) <- ~Population/Subpopulation
#' Aeut.gc
#' Aeut.gc <- namehierarchy(Aeut.gc, ~Pop/Subpop)
#' Aeut.gc
#' }
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
    sethierarchy.internal(x, value)
  })

setMethod(
  f = "sethierarchy",
  signature(x = "snpclone"),
  definition = function(x, value){
    sethierarchy.internal(x, value)
  })

#==============================================================================#
#' @export 
#' @rdname hierarchy-methods
#' @aliases sethierarchy,genclone-method sethierarchy,snpclone-method
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
    return(sethierarchy.internal(x, value))
  })

setMethod(
  f = "sethierarchy<-",
  signature(x = "snpclone"),
  definition = function(x, value){
    return(sethierarchy.internal(x, value))
  })

#==============================================================================#
#' @export 
#' @rdname hierarchy-methods
#' @aliases namehierarchy,genclone-method namehierarchy,snpclone-method
#' @docType methods
#==============================================================================#
namehierarchy <- function(x, value){
  standardGeneric("namehierarchy")
}  

#' @export
setGeneric("namehierarchy")

setMethod(
  f = "namehierarchy",
  signature(x = "genclone"),
  definition = function(x, value){
    namehierarchy.internal(x, value)
  })

setMethod(
  f = "namehierarchy",
  signature(x = "snpclone"),
  definition = function(x, value){
    namehierarchy.internal(x, value)
  })

#==============================================================================#
#' @export 
#' @rdname hierarchy-methods
#' @aliases namehierarchy<-,genclone-method namehierarchy<-,snpclone-method
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

setMethod(
  f = "namehierarchy<-",
  signature(x = "snpclone"),
  definition = function(x, value){
    return(namehierarchy(x, value))
  })

#==============================================================================#
#' @export 
#' @rdname hierarchy-methods
#' @aliases splithierarchy,genclone-method splithierarchy,snpclone-method
#' @docType methods
#' @param sep a \code{character} indicating the character used to separate
#' hierarchical levels. This defaults to "_".
#' @importFrom reshape2 colsplit
#==============================================================================#
splithierarchy <- function(x, value, sep = "_"){
  standardGeneric("splithierarchy")
}  

#' @export
setGeneric("splithierarchy")

setMethod(
  f = "splithierarchy",
  signature(x = "genclone"),
  definition = function(x, value, sep = "_"){
    splithierarchy.internal(x, value, sep) 
  })

setMethod(
  f = "splithierarchy",
  signature(x = "snpclone"),
  definition = function(x, value, sep = "_"){
    splithierarchy.internal(x, value, sep) 
  })

#==============================================================================#
#' @export 
#' @rdname hierarchy-methods
#' @aliases splithierarchy<-,genclone-method splithierarchy<-,snpclone-method
#' @docType methods
#==============================================================================#
"splithierarchy<-" <- function(x, sep = "_", value){
  standardGeneric("splithierarchy<-")
}  

#' @export
setGeneric("splithierarchy<-")

setMethod(
  f = "splithierarchy<-",
  signature(x = "genclone"),
  definition = function(x, sep = "_", value){
    return(splithierarchy.internal(x, value, sep))
  })


setMethod(
  f = "splithierarchy<-",
  signature(x = "snpclone"),
  definition = function(x, sep = "_", value){
    return(splithierarchy.internal(x, value, sep))
  })

#==============================================================================#
#' @export 
#' @rdname hierarchy-methods
#' @aliases addhierarchy,genclone-method addhierarchy,snpclone-method
#' @param name an optional name argument for use with addhierarchy if supplying
#'   a vector. Defaults to "NEW".
#' @docType methods
#==============================================================================#
addhierarchy <- function(x, value, name = "NEW"){
  standardGeneric("addhierarchy")
}  

#' @export
setGeneric("addhierarchy")

setMethod(
  f = "addhierarchy",
  signature(x = "genclone"),
  definition = function(x, value, name = "NEW"){
    addhierarchy.internal(x, value, name)
  })

setMethod(
  f = "addhierarchy",
  signature(x = "snpclone"),
  definition = function(x, value, name = "NEW"){
    addhierarchy.internal(x, value, name)
  })

#==============================================================================#
#' @export 
#' @rdname hierarchy-methods
#' @aliases addhierarchy<-,genclone-method addhierarchy<-,snpclone-method
#' @docType methods
#==============================================================================#
"addhierarchy<-" <- function(x, name = "NEW", value){
  standardGeneric("addhierarchy<-")
}  

#' @export
setGeneric("addhierarchy<-")

setMethod(
  f = "addhierarchy<-",
  signature(x = "genclone"),
  definition = function(x, name = "NEW", value){
    return(addhierarchy.internal(x, value, name))
  })

setMethod(
  f = "addhierarchy<-",
  signature(x = "snpclone"),
  definition = function(x, name = "NEW", value){
    return(addhierarchy.internal(x, value, name))
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
#' @aliases setpop,genclone-method setpop,snpclone-method
#' @docType methods 
#' @author Zhian N. Kamvar
#' @examples
#' 
#' data(Aeut)
#' Aeut.gc <- as.genclone(Aeut)
#' 
#' # Notice that there are two hierarchical levels, Pop and Subpop
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
    setpop.internal(x, formula)
  })

setMethod(
  f = "setpop",
  signature(x = "snpclone"),
  definition = function(x, formula = NULL){
    setpop.internal(x, formula)
  })
#==============================================================================#
#' @export
#' @rdname population-methods
#' @aliases setpop<-,genclone-method setpop<-,snpclone-method
#' @docType methods
#==============================================================================#
"setpop<-" <- function(x, value) standardGeneric("setpop<-")

#' @export
setGeneric("setpop<-")

setMethod(
  f = "setpop<-",
  signature(x = "genclone"),
  definition = function(x, value){
    return(setpop.internal(x, value))
  })

setMethod(
  f = "setpop<-",
  signature(x = "snpclone"),
  definition = function(x, value){
    return(setpop.internal(x, value))
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

setMethod(
  f = "mll",
  signature(x = "genclone"),
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

#==============================================================================#
#' @export
#' @rdname mll.custom
#' @aliases mll.custom<-,genclone-method mll.custom<-,snpclone-method
#' @docType methods
#==============================================================================#
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

#==============================================================================#
#' @export
#' @rdname mll.custom
#' @aliases mll.levels,genclone-method mll.levels,snpclone-method
#' @docType methods
#==============================================================================#
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


#==============================================================================#
#' @export
#' @rdname mll.custom
#' @aliases mll.levels<-,genclone-method mll.levels<-,snpclone-method
#' @docType methods
#==============================================================================#
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