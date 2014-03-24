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
    if (length(j) > nLoc(x) | any(j > nLoc(x))){
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
    slot(x, "alllist")   <- alllist
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
#==============================================================================#
setMethod(
  f = "initialize",
  signature = "bootgen",
  definition = function(.Object, gen){
    if (missing(gen)) gen <- new("genclone")
    if (!is.genclone(gen)) gen <- as.genclone(gen)
    lapply(names(gen), function(y) slot(.Object, y) <<- slot(gen, y))
    slot(.Object, "alllist") <- .Call("expand_indices", cumsum(.Object@loc.nall), 
                                      nLoc(.Object))
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
# GENCLONE METHODS
#------------------------------------------------------------------------------#
################################################################################
#==============================================================================#
#' Check for validity of a genclone object
#' 
#' Note that a \linkS4class{genclone} object will always be a valid 
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
    mlg       <- slot(x, "mlg")[i]
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
    res <- new("genclone", res, hierarchy, mlg)
  
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
#' @keywords internal
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
    ltab  <- vapply(ltab, function(x) substr("       ", 1, x + 1), character(1))
    pops  <- object@pop.names
    poplen <- length(pops)
    if (poplen > 7) 
      pops <- c(pops[1:3], "...", pops[(poplen-3):poplen])
    hiernames <- names(object@hierarchy)
    hierlen <- length(hiernames)
    if (hierlen > 7) 
      hiernames <- c(hiernames[1:3], "...", hiernames[(hierlen-3):hierlen])
    cat("\nThis is a genclone object\n")
    cat("-------------------------\n")
    cat("Genotype information:\n\n",
      ltab[1], nmlg, "multilocus genotypes\n",
      ltab[2], nind, ploid, "individuals\n", 
      ltab[3], nloc, type, "loci\n\n"
      )
    pophier <- ifelse(hier > 1, "hierarchies -", "hierarchy -")
    if (hier == 0) pophier <- "hierarchies."
    popdef  <- ifelse(npop > 0, "defined -", "defined.")
    cat("Population information:\n\n")
    cat("", ltab[4], hier, "population", pophier, hiernames, fill = TRUE)
    cat("", ltab[5], npop, "populations", popdef, pops, fill = TRUE)
    
  })

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
    pophier <- ifelse(hier > 1, "hierarchies -", "hierarchy -")
    if (hier == 0) pophier <- "hierarchies."
    popdef  <- ifelse(npop > 0, "defined -", "defined.")
    cat("Population information:\n\n")
    cat("", ltab[4], hier, "population", pophier, hiernames, fill = TRUE)
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
#'   "population\_hierarchy" and set that as the hieararchy. If none is defined,
#'   the population will be set as the hierarchy under the label "Pop". Use the
#'   function \code{\link{splithierarchy}} to split up any population
#'   hierarchies that might be combined in the population factor.
#' 
#' @seealso \code{\link{splithierarchy}}
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
#' 
#' \strong{Function Specifics}
#' \itemize{
#' \item \strong{gethierarchy()} - This will retrieve the data from the
#' \emph{hierarchy} slot in the \linkS4class{genclone} object. You have the
#' option to choose specific parts of the hierarchy using a formula (see below)
#' and you can choose to combine the hierarchy (default)
#' \item \strong{sethierarchy()} - Set or reset a population hierarchy
#' in your \linkS4class{genclone} object.
#' \item \strong{namehierarchy()} - Rename the levels of your population 
#' hierarchy.
#' \item \strong{splithierarchy()} - This is conceptually similar to the default
#' method of \code{\link{splitcombine}}. It is often difficult to import files
#' with several levels of hierarchy as most data formats do not allow unlimited
#' population levels. This is circumvented by collapsing all hierarchical levels
#' into a single population factor with a common separator for each observation.
#' This function will then split those hierarchies for you, but it works best
#' on a hierarchy that only has a single column in it. See the rootrot example
#' below.
#' \item \strong{addhierarchy()} - Add levels to your population hierarchy. If
#' you have extra hierarchical levels you want to add to your population
#' hierarchy, you can use this method to do so. You can input a data frame or a
#' vector, but if you put in a vector, you have the option to name it (if you
#' are using the functional version and not the assignment version).
#' }
#' 
#' \strong{Argument Specifics}
#' 
#' These functions allow the user to seamlessly assign the hierarchy of their 
#' \code{\linkS4class{genclone}} object. Note that there are two ways of 
#' performing all methods (except for \code{gethierarchy()}). They essentially 
#' do the same thing except that the assignment method (the one with the
#' "\code{<-}") will modify the object in place whereas the non-assignment
#' method will not modify the original object. Due to convention, everything
#' right of the assignment is termed \code{value}. To avoid confusion, here is a
#' guide:
#' \enumerate{
#'  \item \strong{sethierarchy()} This will be a \code{\link{data.frame}} that 
#'  defines the hierarchy for each individual in the rows.
#'  \item \strong{namehierarchy()} This will be either a \code{\link{vector}} or
#'  a \code{\link{formula}} that will define the names.
#'  \item \strong{splithierarchy()} This will be a \code{\link{formula}}
#'  argument with the same number of levels as the hierarchy you wish to split.
#'  \item \strong{addhierarchy()} This will be a \code{\link{vector}} or
#'  \code{\link{data.frame}} with the same length as the number of individuals
#'  in your data.
#' }
#' 
#' \strong{Details on Formulas}
#' 
#' The preferred use of these functions is with a \code{\link{formula}} object.
#' Specifically, a hierarchical formula argument is used to assign the levels of
#' the hierarchy. An example of a hierarchical formula would be:\cr
#'    \code{~Country/City/Neighborhood}\cr
#' or \cr
#'    \code{~Country + Country:City + Country:City:Neighborhood}\cr
#' of course, the first method is slightly easier to read. It is important to
#' use hiearchical formulas when specifying hierarchies as other types of formulas
#' (eg. \code{~Country*City*Neighborhood}) might give spurious results. 
#' 
#' @seealso \code{\link{setpop}} \code{\link{genclone}} \code{\link{as.genclone}}
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
namehierarchy <- function(x, value){
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
#' @export 
#' @rdname hierarchy-methods
#' @aliases splithierarchy,genclone-method
#' @docType methods
#' @param sep a \code{character} indicating the character used to separate
#' population hierarchies. This defaults to "_".
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
    if (is.language(value)){
      # valterms <- attr(terms(value), "term.labels")
      # valterms <- valterms[length(valterms)]
      # valterms <- gsub(":", sep, valterms)
      value    <- all.vars(value)
    } else {
      stop("value must be a formula.")
    }
    if (length(value) < 1){
      stop("value must have more than one hierarchical level.")
    }
    hierarchy  <- x@hierarchy
    if (length(hierarchy) > 1){
      warning("Hierarchy must be length 1. Taking the first column.")
      hierarchy <- hierarchy[1]
    }
    seps     <- gregexpr(sep, hierarchy[[1]])
    sepmatch <- vapply(seps, function(val) all(as.integer(val) > 0), logical(1))
    seps     <- vapply(seps, length, numeric(1))
    all_seps_match <- all(sepmatch)
    given_seps     <- length(value) - 1 
    if (!all_seps_match | all(seps != given_seps)){
      seps <- ifelse(all_seps_match, seps[1], 0) + 1
      msg1 <- paste("\n  Data has", seps, ifelse(seps == 1, "level", "levels"),
                    "of hierarchy with the separator", sep, ".")
      msg2 <- paste("Here is the fist column of the data:", hierarchy[1, ])
      stop(paste(msg1, "\n ", msg2))
    }
    hierarchy        <- pop_splitter(hierarchy, sep = sep)[-1]
    names(hierarchy) <- value
    x@hierarchy      <- hierarchy
    return(x) 
  })

#==============================================================================#
#' @export 
#' @rdname hierarchy-methods
#' @aliases splithierarchy<-,genclone-method
#' @docType methods
#==============================================================================#
"splithierarchy<-" <- function(x, value){
  standardGeneric("splithierarchy<-")
}  

#' @export
setGeneric("splithierarchy<-")

setMethod(
  f = "splithierarchy<-",
  signature(x = "genclone"),
  definition = function(x, value){
    return(splithierarchy(x, value))
  })

#==============================================================================#
#' @export 
#' @rdname hierarchy-methods
#' @aliases addhierarchy,genclone-method
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
    
    hierarchy  <- x@hierarchy
    if ((is.vector(value) | is.factor(value)) & length(value) == nrow(hierarchy)){
      NEW <- data.frame(value)
      names(NEW) <- name
      hierarchy <- cbind(hierarchy, NEW)
    } else if (is.data.frame(value) && nrow(value) == nrow(hierarchy)){
      hierarchy <- cbind(hierarchy, value)
    } else {
      stop("value must be a vector or data frame.")
    }
    x@hierarchy <- hierarchy
    return(x) 
  })

#==============================================================================#
#' @export 
#' @rdname hierarchy-methods
#' @aliases addhierarchy<-,genclone-method
#' @docType methods
#==============================================================================#
"addhierarchy<-" <- function(x, value){
  standardGeneric("addhierarchy<-")
}  

#' @export
setGeneric("addhierarchy<-")

setMethod(
  f = "addhierarchy<-",
  signature(x = "genclone"),
  definition = function(x, value){
    return(addhierarchy(x, value))
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
#' @author Zhian N. Kamvar
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
