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
# The point of this file serves as a way of parsing the internal functions that
# are specific to methods since there are now two classes that share common 
# slots. Having the methods specific to these internal slots ensures that the 
# code is a) DRY and b) debuggable (by calling the internal function).
#==============================================================================#
#==============================================================================#
# Internal functions to deal with the MLG class.
#==============================================================================#

#==============================================================================#
# internal function for mll method of clone type objects
#
# Public functions utilizing this function:
# ## mll mll<- (genclone and snpclone)
#
# Internal functions utilizing this function:
# ## none
#==============================================================================#
mll.internal <- function(x, type = NULL, the_call = match.call()){
  mlg <- x@mlg
  if (!"MLG" %in% class(mlg)){
    the_obj <- as.character(the_call[["x"]])
    the_type <- as.character(the_call[["type"]])
    if (length(the_type) == 0) the_type <- "original"
    msg <- paste("\n The @mlg slot does not contain an MLG class object.\n",
                 "Returning the original mlgs. Please use:\n\n",
                 paste0('mll(', the_obj, ') <- "', the_type, '"\n'),
                 "\n to convert your object.")
    warning(msg, call. = FALSE)
    return(mlg)
  }
  if (!is.null(type)){
    TYPES <- c("original", "expanded", "contracted", "custom")
    type <- match.arg(type, TYPES)
  } else {
    type <- visible(mlg)
  }
  return(mlg[, type])
}

#==============================================================================#
# wrapper for mlg.vector
#
# This serves as a wrapper due to the fact that mll should be available for
# genind and genlight objects to avoid unnecessary checking of type.
# 
# Public functions utilizing this function:
# ## mll (genind and genlight)
#
# Internal functions utilizing this function:
# ## none
#==============================================================================#
mll.gen.internal <- function(x, type = NULL){
  if (!is.null(type)){
    msg <- paste("The object you are using is a genind object and does not",
                 "contain an mlg slot. Returning the results of mlg.vector().")
    warning(msg)
  }
  mlg.vector(x)
}

#==============================================================================#
# internal function to reset multilocus genotypes of clone type objects
#
# Public functions utilizing this function:
# ## mll.reset mll.reset<-
#
# Internal functions utilizing this function:
# ## none
#==============================================================================#
mll.reset.internal <- function(x, value){
  if (missing(value)){
    stop("please specify a value to reset MLGs")
  }
  true_value <- is.logical(value) && length(value) == 1 && value == TRUE
  if (!is(x@mlg, "MLG") | true_value){
    x@mlg <- new("MLG", mlg.vector(x, reset = TRUE))
    return(x)
  }
  new_mlg <- x@mlg@mlg
  TYPES   <- c("original", "expanded", "contracted", "custom")
  types   <- match.arg(value, TYPES, several.ok = TRUE)
  
  if ("original" %in% types){
    newvec <- mlg.vector(x, reset = TRUE)
  } else {
    newvec <- mll(x, "original")
  }
  for (i in types){
    if (i == "custom"){
      new_mlg[i] <- as.factor(newvec)
    } else {
      new_mlg[i] <- newvec
    }
  }
  x@mlg@mlg <- new_mlg
  if ("contracted" %in% types){
    cutoff(x@mlg)["contracted"] <- 0
    distname(x@mlg) <- if (inherits(x, "genclone")) "diss.dist" else "bitwise.dist"
    distargs(x@mlg) <- list()
  }
  return(x)
}

#==============================================================================#
# internal function to get and set custom values for clone type objects
#
# Public functions utilizing this function:
# ## mll.custom mll.custom<-
#
# Internal functions utilizing this function:
# ## none
#==============================================================================#
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

#==============================================================================#
# internal function to get and set levels for custom mlls
#
# Public functions utilizing this function:
# ## mll.levels mll.levels<-
#
# Internal functions utilizing this function:
# ## none
#==============================================================================#
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

#==============================================================================#
# internal function for mlg.filter
#
# Public functions utilizing this function:
# ## mlg.filter mlg.filter<-
#
# Internal functions utilizing this function:
# ## none
#==============================================================================#
mlg.filter.internal <- function(gid, threshold = 0.0, missing = "asis", 
                                memory = FALSE, algorithm = "farthest_neighbor", 
                                distance = "diss.dist", threads = 0, 
                                stats = "MLGs", the_call = match.call(), ...){

  # This will return a vector indicating the multilocus genotypes after applying
  # a minimum required distance threshold between multilocus genotypes.
  if (is.character(distance) || is.function(distance)) {
    if (memory==TRUE && identical(c(gid, distance, ...), .last.value.param$get())){
      dis <- .last.value.dist$get()
    } else {
      if (is.genind(gid)) {
        the_dist <- as.character(the_call[["distance"]])
        call_len <- length(the_dist)
        is_diss_dist <- the_dist %in% "diss.dist"
        any_dist <- the_dist %in% c("diss.dist", "nei.dist", "prevosti.dist",
                                    "edwards.dist", "reynolds.dist", 
                                    "rogers.dist", "provesti.dist")
        if (missing == "mean" && call_len == 1 && is_diss_dist){
          disswarn <- paste("Cannot use function diss.dist and correct for", 
                            "mean values.", "diss.dist will automatically",
                            "ignore missing data.") 
          warning(disswarn, call. = FALSE)
          mpop <- gid
        } else if (call_len == 1 && any_dist) {
          mpop <- new("bootgen", gid, na = missing, 
                      freq = ifelse(is_diss_dist, FALSE, TRUE))
        } else {
          mpop <- missingno(gid, type = missing, quiet = TRUE)
        }
      } else {
        mpop <- gid
      }
      DISTFUN <- match.fun(distance)
      dis <- DISTFUN(mpop, ...)
      dis <- as.matrix(dis)
      if (memory == TRUE)
      {
        .last.value.param$set(c(gid, distance, ...))
        .last.value.dist$set(dis)
      }
    }
  } else {
    # Treating distance as a distance table 
    # Warning: Missing data in distance matrix or data uncorrelated with gid may
    # produce unexpected results.
    dis <- as.matrix(eval(distance))
  }
  
  if (!is.clone(gid)) {
    if (is(gid, "genlight")){
      gid <- as.snpclone(gid, mlg = seq_len(nInd(gid)))
    } else {
      gid <- as.genclone(gid)
    }
  } else {
    if (!is(gid@mlg, "MLG")){
      gid@mlg <- new("MLG", gid@mlg)
    }
    mll(gid) <- "original"
  }
  basemlg <- mlg.vector(gid)
  
  
  
  # Input validation --------------------------------------------------------
  # 
  if (any(is.na(dis))){
    msg <- paste("The resulting distance matrix contains missing data.\n",
                 "Please treat your missing data by using the missing",
                 "argument.\n")
    stop(msg, call. = FALSE)
  }
  
  if (!is.numeric(dis) && !is.integer(dis)){
    stop("Distance matrix must be a square matrix of numeric or integer values.",
         call. = FALSE)
  }
  if (nrow(dis) != ncol(dis)){
    stop("The distance matrix must be a square matrix", call. = FALSE)
  }
  
  if (nrow(dis) != nInd(gid)){
    msg <- paste0("The number of observations in the distance matrix (",
                  nrow(dis), ") are not equal to the number of observations in",
                  " the data (", nInd(gid), ").")
    if (nrow(dis) == nPop(gid)){
      msg <- paste(msg, "\n\nPlease check your distance function to make sure",
                   "it calculates distances between samples, not populations.")
    }
    stop(msg, call. = FALSE)
  }
    # Threshold must be something that can cast to numeric
  if (!is.numeric(threshold) && !is.integer(threshold)){
    stop("Threshold must be a numeric or integer value", .call = FALSE)
  } 
    # Threads must be something that can cast to integer
  if (!is.numeric(threads) && !is.integer(threads) && threads >= 0){
    stop("Threads must be a non-negative numeric or integer value", .call  = FALSE)
  } 
    # Stats must be logical
  STATARGS <- c("MLGS", "THRESHOLDS", "DISTANCES", "SIZES", "ALL")
  stats <- match.arg(toupper(stats), STATARGS, several.ok = TRUE)

  # Cast parameters to proper types before passing them to C
  dis_dim   <- dim(dis)
  dis       <- as.numeric(dis)
  dim(dis)  <- dis_dim # Turn it back into a matrix
  threshold <- as.numeric(threshold)
  algo      <- tolower(as.character(algorithm))
  threads   <- as.integer(threads)
  
  if (!isTRUE(all.equal(basemlg, as.integer(basemlg))))
  {
    warning("MLG contains non-integer values. MLGs differing only in decimal values will be merged.")
  }
  basemlg <- as.integer(basemlg)
  
  result_list <- .Call("neighbor_clustering", dis, basemlg, threshold, algo, threads) 
  
  # Cut out empty values from result_list[[2]]
  result_list[[2]] <- result_list[[2]][result_list[[2]] > -0.05]
  # Format result_list[[3]]
  mlgs  <- unique(result_list[[1]])
  dists <- result_list[[3]]
  dists <- dists[mlgs,mlgs]
  if (length(mlgs) > 1){
    rownames(dists) <- mlgs
    colnames(dists) <- mlgs
  }
  result_list[[3]] <- dists
  names(result_list) <- c("MLGS", "THRESHOLDS", "DISTANCES", "SIZES")
  if (length(stats) == 1){
    if (toupper(stats) == "ALL"){
      return(result_list)
    } else {
      return(result_list[[stats]])
    } 
  } else {
    return(result_list[stats])
  }
}
