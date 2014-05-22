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
#
# The clone correct function will need a parameter for the lowest population
# level in order to keep at least one individual represented in each population.
# It takes a popper object and will return a poppr object.
#
#' Remove potential bias caused by cloned genotypes in genind object.
#' 
#' This function removes any duplicated multi locus genotypes from any specified
#' population hierarchy.
#'
#' @param pop a \code{\link{genind}} object
#'
#' @param hier a \code{numeric or character list}. This is the list of vectors
#' within a data frame (specified in \code{dfname}) in the 'other' slot of the
#' \code{\link{genind}} object. The list should indicate the population
#' hierarchy to be used for clone correction.
#'
#' @param dfname a \code{character string}. This is the name of the data frame
#' or list containing the vectors of the population hierarchy within the
#' \code{other} slot of the \code{\link{genind}} object. 
#'
#' @param combine \code{logical}. When set to TRUE, the heirarchy will be
#' combined to create a new population for the genind object.
#'
#' @param keep \code{integer}. When \code{combine} is set to \code{FALSE}, you
#' can use this flag to choose the levels of your population hierarchy. For
#' example: if your clone correction hierarchy is set to "Pop", "Subpop", and 
#' "Year", and you want to analyze your populations with respect to year, you 
#' can set \code{keep = c(1,3)}.
#' 
#' @return a clone corrected \code{\link{genind}} object. 
#' 
#' @note 
#' This function will clone correct to the population level indicated in
#' the \code{pop} slot of the \code{\link{genind}} object if there is no data
#' frame specified in dfname. If there is no population structure and there is
#' no specified data frame, it will clone correct the entire
#' \code{\link{genind}} object. 
#' 
#'
#' @export
#' @author Zhian N. Kamvar
#' @examples
#' # LOAD A. euteiches data set
#' data(Aeut)
#' 
#' # Check the number of multilocus genotypes
#' mlg(Aeut)
#' Aeut$pop.names
#' 
#' # Clone correct at the population level.
#' Aeut.pop <- clonecorrect(Aeut, hier="Pop")
#' mlg(Aeut.pop)
#' Aeut.pop$pop.names
#' 
#' \dontrun{
#' # Clone correct at the subpopulation level with respect to population and
#' # combine.
#' Aeut.subpop <- clonecorrect(Aeut, hier=c("Pop", "Subpop"), combine=TRUE)
#' mlg(Aeut.subpop)
#' Aeut.subpop$pop.names
#' 
#' # Do the same, but set to the population level.
#' Aeut.subpop2 <- clonecorrect(Aeut, hier=c("Pop", "Subpop"), keep=1)
#' mlg(Aeut.subpop2)
#' Aeut.subpop2$pop.names
#' 
#' # LOAD H3N2 dataset
#' data(H3N2)
#'
#' # Extract only the individuals located in China
#' country <- clonecorrect(H3N2, hier=c("country"), dfname="x")
#'
#' # How many isolates did we have from China before clone correction?
#' length(which(other(H3N2)$x$country=="China")) # 155
#'
#' # How many unique isolates from China after clone correction?
#' length(which(other(country)$x$country=="China")) # 79
#' 
#' # Something a little more complicated. (This could take a few minutes on
#' # slower computers)
#'
#' # setting the hierarchy to be Country > Year > Month  
#' c.y.m <- clonecorrect(H3N2, hier=c("year","month","country"), dfname="x")
#'
#' # How many isolates in the original data set?
#' length(other(H3N2)$x$country) # 1903
#'
#' # How many after we clone corrected for country, year, and month?
#' length(other(c.y.m)$x$country) # 1190
#' }
#==============================================================================#

clonecorrect <- function(pop, hier=1, dfname="population_hierarchy", 
                         combine = FALSE, keep = 1){
  clonecall <- match.call()$pop
  if(!is.genind(pop)){
    stop(paste(paste(substitute(pop), collapse=""), "is not a genind object.\n"))
  }
  if (is.language(hier)){
    hierformula <- hier
    hier <- all.vars(hier)
  }
  popcall <- pop@call
  if (is.na(hier[1])){
    return(pop[.clonecorrector(pop), ])
  }
  if (is.genclone(pop)){
    if (is.numeric(hier)){
      hier <- names(gethierarchy(pop))[hier]
      hierformula <- as.formula(paste0("~", paste(hier, collapse = "/")))
    }
    if (!all(hier %in% names(gethierarchy(pop)))){
      stop(hier_incompatible_warning(hier, gethierarchy(pop)))
    }
    setpop(pop) <- hierformula
  } else if (is.null(other(pop)[[dfname]])) {
    if(length(hier) == 1 & hier[1] == 1){
      if(length(levels(pop(pop))) == 1 | is.null(pop(pop))){
        pop <- pop[.clonecorrector(pop), ]
        return(pop)
      }
      else if(length(levels(pop(pop))) > 1){
        other(pop)[[dfname]] <- as.data.frame(list(Pop = as.character(pop(pop))))
        warning(paste("There is no data frame in ",
                      paste(substitute(clonecall), collapse=""), 
                      "@other called ",dfname,
                      ".\nOne is being created from the population factor.", sep=""))
      }
    }
    else{
      stop(paste("There is no data frame in ",paste(substitute(clonecall), collapse=""), 
                 "@other called ",dfname,".\n", sep=""))
    }
  }
  # Corrects the individual names of the object. This is fo the fact that the
  # clone corrector relies on the unique individual names for it to work.
  if(all(pop@ind.names == "")){
    pop@ind.names <- as.character(1:nInd(pop))
  }
  
  # Combining the population factor by the hierarchy
  if(!is.genclone(pop)){
    pop <- splitcombine(pop, method=2, dfname=dfname, hier=hier)
  }
  cpop <- length(pop$pop.names)
  
  # Steps for correction:
  # Subset by population factor.
  # Run subsetted population by the .clonecorrector
  # Profit!
  corWrecked <- function(x, pop){
    subbed <- popsub(pop, x) # population to be...corrected.
    subbed <- subbed[.clonecorrector(subbed), ] 
    # Return the indices based off of the individual names.
    return(which(pop@ind.names %in% subbed@ind.names))
  }
  
  ccpop <- unlist(lapply(1:cpop, corWrecked, pop))
  pop <- pop[ccpop, ]
  
  if(!combine){
    # When the combine flag is not true, the default is to keep the first level
    # of the hierarchy. The keep flag is a numeric vector corresponding to the
    # hier flag indicating which levels the user wants to keep.
    if (is.genclone(pop)){
      hier <- hier[keep]
      newformula <- as.formula(paste0("~", paste(hier, collapse = "/")))
      setpop(pop) <- newformula
    } else {
      if(length(keep) > 1){
        pop <- splitcombine(pop, hier=hier[keep], method=2, dfname=dfname)
      }
      else{
        pop(pop) <- pop$other[[dfname]][[hier[keep]]]
      }
      names(pop$pop.names) <- levels(pop$pop)
    }
  }
  pop@call <- popcall
  return(pop)
}


#==============================================================================#
# subset a population with a combination of sublists and blacklists. Either one
# is optional, and the default is to do nothing. The structure will allow the
# user to select a range of populations and exclude a small number of them
# without having to use the total. 
# eg pop <- popsub(gid, sublist=1:50, blacklist=c(17, 33))
# 
#' Subset a \code{\link{genind}} object by population
#' 
#' Create a new dataset with specified populations or exclude specified
#' populations from the dataset.
#' 
#' @param gid a \code{\link{genind}} object.
#' 
#' @param sublist a \code{vector} of population names or indexes that the user
#' wishes to keep. Default to "ALL".
#'
#' @param blacklist a \code{vector} of population names or indexes that the user
#' wishes to discard. Default to \code{NULL}
#'
#' @param mat a \code{matrix} object produced by \code{\link{mlg.table}} to be
#' subsetted. If this is present, the subsetted matrix will be returned instead
#' of the genind object 
#'
#' @param drop \code{logical}. If \code{TRUE}, unvariate alleles will be dropped
#' from the population.
#' 
#' @return A \code{genind} object or a matrix.
#'
#' @examples
#' # Load the dataset microbov.
#' data(microbov)
#' 
#' # Analyze only the populations with exactly 50 individuals
#' mic.50 <- popsub(microbov, sublist=c(1:6, 11:15), blacklist=c(3,4,13,14))
#'
#' \dontrun{
#' # Analyze the first 10 populations, except for "Bazadais"
#' mic.10 <- popsub(microbov, sublist=1:10, blacklist="Bazadais")
#' 
#' # Take out the two smallest populations
#' micbig <- popsub(microbov, blacklist=c("NDama", "Montbeliard"))
#' 
#' # Analyze the two largest populations
#' miclrg <- popsub(microbov, sublist=c("BlondeAquitaine", "Charolais"))
#' }
#' @export
#==============================================================================#

popsub <- function(gid, sublist="ALL", blacklist=NULL, mat=NULL, drop=TRUE){

  if (!is.genind(gid)){
    stop("popsub requires a genind object\n")
  }
  if (is.null(pop(gid))){
    if(sublist[1] != "ALL")
      warning("No population structure. Subsetting not taking place.")
    return(gid)
  }
  orig_list <- sublist 
  popnames <- gid@pop.names
  if (toupper(sublist[1]) == "ALL"){
    if (is.null(blacklist)){
      return(gid)
    } else {
      # filling the sublist with all of the population names.
      sublist <- popnames 
    }
  }

  # Checking if there are names for the population names. 
  # If there are none, it will give them names. 
  if (is.null(names(popnames))){
    if (length(popnames) == length(levels(gid@pop))){
      names(popnames) <- levels(gid@pop)
    }
    else{
      stop("Population names do not match population factors.")
    }
  }

  # Treating anything present in blacklist.
  if (!is.null(blacklist)){

    # If both the sublist and blacklist are numeric or character.
    if (is.numeric(sublist) & is.numeric(blacklist) | class(sublist) == class(blacklist)){
      sublist <- sublist[!sublist %in% blacklist]
    }
    
    # if the sublist is numeric and blacklist is a character. eg s=1:10, b="USA"
    else if (is.numeric(sublist) & class(blacklist) == "character"){
      sublist <- sublist[sublist %in% which(!popnames %in% blacklist)]
    } else {
      # no sublist specified. Ideal situation
      if(all(popnames %in% sublist)){
        sublist <- sublist[-blacklist]
      }
      # weird situation where the user will specify a certain sublist, yet index
      # the blacklist numerically. Interpreted as an index of populations in the
      # whole data set as opposed to the sublist.
      else{
        warning("Blacklist is numeric. Interpreting blacklist as the index of the population in the total data set.")
        sublist <- sublist[!sublist %in% popnames[blacklist]]
      }
    }
  }
  if (!is.null(mat)){
    mat <- mat[sublist, , drop=FALSE]
    return(mat[, which(colSums(mat) > 0), drop=FALSE])
  } else {
    # subsetting the population. 
    if (is.numeric(sublist)){
      sublist <- popnames[sublist]
    } else {
      sublist <- popnames[popnames %in% sublist]
    }
    sublist <- pop(gid) %in% sublist
    if (!any(sublist)){
      if (!is.numeric(orig_list) & !any(gid@pop.names %in% orig_list)){
        stop(unmatched_pops_warning(gid@pop.names, orig_list))
      } else {
        nothing_warn <- paste("Nothing present in the sublist.\n",
                            "Perhaps the sublist and blacklist arguments have",
                            "duplicate entries?\n",
                            "Subsetting not taking place.")
        warning(nothing_warn)
        return(gid)
      }
    }
    gid <- gid[sublist, , drop = drop]
    gid@call <- match.call()
    return(gid)
  }
}

#==============================================================================#
# missigno simply applies one of four methods to deal with missing data.
# default is to remove missing loci. 
#' How to deal with missing data in a genind object.
#' 
#' missingno gives the user four options to deal with missing data.
#'
#' @param pop a \code{\link{genind}} object.
#'
#' @param type a \code{character} string: can be "zero", "mean", "loci", or "geno"
#' (see \code{Details} for definitions).]
#' 
#' @param cutoff \code{numeric}. A number from 0 to 1 indicating the allowable
#' rate of missing data in either genotypes or loci. This will be ignored for
#' \code{type} values of \code{"mean"} or \code{"zero"}.
#' 
#' @param quiet if \code{TRUE}, it will print to the screen the action performed.
#'
#' @section Details: The default way that functions in \code{poppr} deal with
#' missing data is to simply ignore it. These methods provide a way to deal with
#' systematic missing data and to give a wrapper for \code{adegenet}'s \code{
#' \link{na.replace}} function. ALL OF THESE ARE TO BE USED WITH CAUTION.
#'
#' \strong{\code{"loci"}} - removes all loci containing missing data in the entire data
#' set. 
#'
#' \strong{\code{"genotype"}} - removes any genotypes/isolates/individuals with missing data.
#'
#' \strong{\code{"mean"}} - replaces all NA's with the mean of the alleles for the entire
#' data set.
#'
#' \strong{\code{"zero"}} or \strong{\code{"0"}} - replaces all NA's with "0". 
#' Introduces more diversity.
#'
#' @return a \code{\link{genind}} object.
#'
#' @note
#' \emph{"wild missingno appeared!"}
#'
#' @seealso \code{\link{na.replace}}, \code{\link{poppr}}
#'
#' @export
#' @author Zhian N. Kamvar
#' @examples
#'
#' data(nancycats)
#' 
#' nancy.locina <- missingno(nancycats, type = "loci")
#' 
#' ## Found 617 missing values.
#' ## 2 loci contained missing values greater than 5%.
#' ## Removing 2 loci : fca8 fca45 
#' 
#' nancy.genona <- missingno(nancycats, type = "geno")
#'
#' ## Found 617 missing values.
#' ## 38 genotypes contained missing values greater than 5%.
#' ## Removing 38 genotypes : N215 N216 N188 N189 N190 N191 N192 N302 N304 N310 
#' ## N195 N197 N198 N199 N200 N201 N206 N182 N184 N186 N298 N299 N300 N301 N303 
#' ## N282 N283 N288 N291 N292 N293 N294 N295 N296 N297 N281 N289 N290  
#'
#' # Replacing all NA with "0" (see na.replace in the adegenet package).
#' nancy.0 <- missingno(nancycats, type = "0")
#'
#' ## Replaced 617 missing values 
#' 
#' # Replacing all NA with the mean of each column (see na.replace in the
#' # adegenet package).
#' nancy.mean <- missingno(nancycats, type = "mean")
#' 
#' ## Replaced 617 missing values 
#==============================================================================#
  #####
  #####
  #####
  #####
#######
#######
#######
#######
#######
#######
missingno <- function(pop, type = "loci", cutoff = 0.05, quiet=FALSE){
  if(sum(is.na(pop@tab)) > 0){
    # removes any loci (columns) with missing values.
    MISSINGOPTS <- c("loci", "genotypes", "mean", "zero", "0", "ignore")
    type        <- match.arg(tolower(type), MISSINGOPTS)
    if (type == "ignore"){
      return(pop)
    }
    navals      <- percent_missing(pop, type = type, cutoff = cutoff)
    if (type == "loci"){
      if(quiet != TRUE){
        # if(all(naloci < 0)){
        #   remloc <- pop@loc.names[which(cumsum(pop@loc.nall) %in% -naloci)]
        #   cat("\n Found", sum(is.na(pop@tab)),"missing values.")
        #   loci <- paste(length(remloc), ifelse(length(remloc) == 1, "locus", "loci"))
        #   cat("\n",loci,"contained missing values greater than",paste(cutoff*100,"%.",sep=""))
        #   cat("\n Removing",loci,":", remloc,"\n", fill = 80)
        # }
        # else{
        #   cat("\n No loci with missing values above",paste(cutoff*100,"%",sep=""),"found.\n")
        # }
        if (length(navals) == ncol(pop@tab)){
          cat("\n No loci with missing values above",
              paste0(cutoff*100,"%"),"found.\n")
        } else {
          remloc <- pop@loc.names[!cumsum(pop@loc.nall) %in% navals]
          cat("\n Found", sum(is.na(pop@tab)),"missing values.")
          loci   <- paste(length(remloc), ifelse(length(remloc) == 1, "locus", 
                          "loci"))
          cat("\n", loci, "contained missing values greater than",
              paste0(cutoff*100,"%."))
          cat("\n Removing", loci, ":", remloc,"\n", fill = 80)
        }
      }
      pop <- pop[, navals]
    }  
    # removes any genotypes (rows) with missing values.
    else if (type == "genotypes"){
      # nageno <- percent_missing(pop, type=type, cutoff=cutoff)
      if(quiet != TRUE){
        # if(all(nageno < 0)){
        #   remgeno <- pop@ind.names[-nageno]
        #   cat("\n Found", sum(is.na(pop@tab)),"missing values.")
        #   genotypes <- paste(length(remgeno), ifelse(length(remgeno) == 1, 
        #                                              "genotype", "genotypes"))
        #   cat("\n",genotypes,"contained missing values greater than",
        #       paste(cutoff*100,"%.",sep=""))
        #   cat("\n Removing",genotypes,":",remgeno,"\n", fill = 80)
        # }
        # else{
        #   cat("\n No genotypes with missing values above",
        #       paste(cutoff*100,"%",sep=""),"found.\n")
        # }
        if (length(navals) == nInd(pop)){
          cat("\n No genotypes with missing values above",
              paste0(cutoff*100, "%"),"found.\n")
        } else {
          remgeno <- indNames(pop)[-navals]
          cat("\n Found", sum(is.na(pop@tab)),"missing values.")
          geno    <- paste(length(remgeno), ifelse(length(remgeno) == 1, 
                           "genotype", "genotypes"))
          cat("\n", geno, "contained missing values greater than",
              paste0(cutoff*100,"%."))
          cat("\n Removing", geno, ":", remgeno,"\n", fill = 80)
        }
      }
      pop <- pop[navals, ]
    }
    # changes all NA's to the mean of the column. NOT RECOMMENDED
    else if (type == "mean"){
      pop <- na.replace(pop, "mean", quiet=quiet)
    }
    # changes all NA's to 0. NOT RECOMMENDED. INTRODUCES MORE DIVERSITY.
    else if (type %in% c("zero","0")){
      pop <- na.replace(pop, "0", quiet=quiet)
    }
  }
  else{
    if(quiet == FALSE){
      cat("\n No missing values detected.\n")
    }
  }
  return(pop)
}

#==============================================================================#
#' Split a or combine items within a data frame in \code{\link{genind}} objects (DEPRECATED).
#'
#' Often, one way a lot of file formats fail is that they do not allow multiple
#' population hierarchies. This can be circumvented, however, by coding all of
#' the hierarchies in one string in the input file with a common separator (eg.
#' "_"). \code{splitcombine} will be able to recognise those separators and
#' create a data frame of all the population structures for whatever subsetting
#' you might need. 
#'
#' @param pop a \code{\link{genind}} object.
#'
#' @param method an \code{integer}, 1 for splitting, 2 for combining.
#'
#' @param dfname the name of the data frame containing the population structure.
#' for the splitting method, the combined population structure must be in the
#' first column. 
#'
#' @param sep The separator used for separating or combining the data. See note.
#'
#' @param hier a \code{vector} containing the population hierarchy you wish to
#' split or combine. 
#'
#' @param setpopulation \code{logical}. if \code{TRUE}, the population of the
#' resulting genind object will be that of the highest population structure
#' (split method) or the combined populations (combine method).
#'
#' @param fixed \code{logical}. An argument to be passed onto
#' \code{\link{strsplit}}. If \code{TRUE}, \code{sep} must match exactly to the
#' populations for the split method. 
#'
#' @return a \code{\link{genind}} object with a modified data frame in the
#' \code{\link{other}} slot.
#'
#' @note 
#' This function has been deprecated and replaced by functions like
#' \code{\link{splithierarchy}}. Please consider using the
#' \code{\linkS4class{genclone}} object for storing hierarchies. 
#'
#' The separator
#' field is sensitive to regular expressions. If you do not know what those are,
#' please use the default underscore to separate your populations. Use \code{fixed
#' = TRUE} to ignore regular expressions. If you do not set the \code{hier} flag
#' for the split method, your new data frame will have the names "comb", "h1", "h2"
#' and so on; for the combine method, your data frame will return the first column
#' of your data frame.
#'
#' @export
#' @author Zhian N. Kamvar
#' @examples
#' \dontrun{
#' # Method 1: Splitting.
#' Aeut <- read.genalex(system.file("files/rootrot.csv", package="poppr"))
#' 
#' # We have 19 different "populations", but really, there is a hierarchy.
#' Aeut$pop.names
#' 
#' # Let's split them up. The default data frame from read.genalex is the same
#' # as the default for this function. 
#' Aeut <- splitcombine(Aeut, hier=c("Pop", "Subpop"))
#' 
#' # Much better!
#' Aeut$pop.names
#' }
#'
#' # Method 2: Combining.
#' 
#' data(H3N2)
#' # Create a new data set combining the population factors of year and country
#' H.comb <- splitcombine(H3N2, method=2, dfname="x", hier=c("year", "country"))
#' 
#' # Checking to make sure they were actually combined.
#' head(H.comb$other$x$year_country)
#' \dontrun{
#' # Creating new data frame in the object to mess around with. 
#' H.comb$other$year_country <- data.frame(H.comb$other$x$year_country)
#' 
#' # Splitting those factors into their original components and setting the
#' # population to year.
#' H.comb <- splitcombine(H.comb, method=1, dfname="year_country", hier=c("year", "country"))
#' }
#==============================================================================#
splitcombine <- function(pop, method=1, dfname="population_hierarchy", sep="_", hier=1, setpopulation=TRUE, fixed=TRUE){
  if (!is.genind(pop)){
    stop(paste(paste(substitute(pop), collapse=""), "is not a genind object.\n"))
  }
  if (!is.data.frame(pop$other[[dfname]])){
    stop(paste("There is no data frame in ",paste(substitute(pop), collapse=""), "@other called ",dfname,".\n", sep=""))
  }
  METHODS = c("Split", "Combine")
  if (all((1:2)!=method)) {
    cat("1 = Split\n")
    cat("2 = Combine\n")
    cat("Select an integer (1 or 2): ")
    method <- as.integer(readLines(n = 1))
  }
  if (all((1:2)!=method)) (stop ("Non convenient method number"))
  # Splitting !
  if(method == 1){
    # Remember, this acts on the first column of the data frame.
    # It will split the population factors in that data frame by sep. 
    df <- pop_splitter(pop$other[[dfname]], sep=sep)
    
    # This makes sure that the length of the given hierarchy is the same as the
    # hierarchy in the original population factor. It's renaming the original
    # factor the population hierarchy separated by sep.
    # eg. df: Pop_Subpop_Year, Pop, Subpop, Year
    if(length(df) - 1 == length(hier)){
      names(df) <- c(paste(hier, collapse=sep), hier)
    }
    # In the case that the number of columns added to the new df is the same
    # as the number of hierarchies specified.
    else if(length(df) - length(hier) == length(pop$other[[dfname]]) ){
      names(df)[(length(df) - length(hier) + 1):length(df)] <- hier
    }

    # Checking to see if there was only one column in the original data frame
    # This is necessary to avoid overwriting data.
    if(length(pop$other[[dfname]] == 1)){
      pop$other[[dfname]] <- df
    }
    
    # If there are other columns in the data frame, we check to see if any of
    # the names in the new data frame are the same. 
    else if(any(names(pop$other[[dfname]]) %in% names(df[-1]))){
      # Removing the combined column. It will remain the same.
      df <- df[-1]
      # Gathering the index of names that are the same as the split hierarchy.
      dfcols <- which(names(pop$other[[dfname]]) %in% names(df))
      
      # All the names are equal, so we replace those columns with the newly
      # subsetted columns.
      if(length(names(df)) == length(dfcols))
        pop$other[[dfname]][dfcols] <- df
      
      # Not all columns match, so you replace the ones that do and then bind
      # the new ones.
      else{
        popothernames <- which(names(df) %in% names(pop$other[[dfname]][dfcols]))
        pop$other[[dfname]][dfcols] <- df[popothernames]
        pop$other[[dfname]] <- cbind(pop$other[[dfname]], df[-dfcols])
      }
    }
    
    # If there are no names in the new data frame that match the original,
    # simply tack the new data frame on to the end. These will have the names
    # h1, h2, h3, etc.
    else{  
      pop$other[[dfname]] <- cbind(pop$other[[dfname]], df[-1])
    }
    
    # Set the population to the highest level of the hierarchy.
    if(setpopulation == TRUE){
      hier <- ifelse(is.numeric(hier), as.character(hier), hier)
      pop(pop) <- pop$other[[dfname]][[hier[1]]]
      names(pop$pop.names) <- levels(pop$pop)
    }
    return(pop)
  }
  
  # Combining !
  else if(method == 2){
    newdf <- pop_combiner(pop$other[[dfname]], hier=hier, sep=sep)
    if(all(is.na(newdf))){
      stop(paste("\n\nError in combining population factors.\nCheck your hier flag and make sure that the columns",paste(hier, collapse=" and "),
                 "exist within the data frame called",dfname,
                 "in the @other slot of your genind object.\n\n"))
    }
    # Combining the hierarchy. This will no longer give you numbers as names,
    # rather it will return the names.
    pop$other[[dfname]][[paste(names(pop$other[[dfname]][hier]), collapse=sep)]] <- newdf
    if(setpopulation == TRUE){
      pop(pop) <- newdf
      names(pop$pop.names) <- levels(pop$pop)
    }
    return(pop)
  }
}
#==============================================================================#
#' Remove all non-phylogentically informative loci
#' 
#' This function will facilitate in removing phylogenetically uninformative loci
#' from a \code{\link{genind}} object. The user can specify what is meant by
#' phylogenetically uninformative with a specification of the cutoff percentage.
#' Any loci under the cutoff will be removed. For convenience's sake, the
#' default cutoff is set to 2 individuals.
#' 
#' @param pop a \code{\link{genind}} object.
#' 
#' @param cutoff \code{numeric}. This is a number from 0 to 1 representing the
#' minimum percentage of differentiating individuals. Defaults is 2 individuals.
#'
#' @param quiet \code{logical}. When \code{quiet = TRUE}, messages indicating
#' the loci removed will be printed to screen. When \code{quiet = FALSE}, 
#' nothing will be printed to screen.
#' 
#' @return A \code{genind} object with user-defined informative loci.
#'
#' @note This will have a few side effects that affect certain analyses. First,
#' the number of multilocus genotypes might be reduced due to the reduced number
#' of markers. Second, if you plan on using this data for analysis of the index
#' of association, be sure to use the standardized version (rbarD) that corrects
#' for the number of observed loci. 
#'
#' @examples
#' # Load the data set H3N2
#' data(H3N2)
#' pop(H3N2) <- H3N2$other$x$country
#' Nepal <- popsub(H3N2, "Nepal")
#'
#' # Using the default 2 individuals.
#' N.inform <- informloci(Nepal)
#'
#' # 5 individuals.
#' N.informfive <- informloci(Nepal, cutoff = 5/nInd(Nepal))
#'
#' # 10 individuals. Too many. Gives warning.
#' N.informten <- informloci(Nepal, cutoff = 10/nInd(Nepal))
#'
#' # Decimate (10%)
#' N.informdecimated <- informloci(Nepal, cutoff = 0.1)
#' @export
#==============================================================================#
#' @importFrom pegas as.loci
informloci <- function(pop, cutoff = 2/nInd(pop), quiet = FALSE){
  if(!is.genind(pop)){
    stop("This function only works on genind objects.")
  }
  MLG <- mlg(pop, quiet = TRUE)
  if(MLG < 3){
    if(!isTRUE(quiet)){
      cat("Not enough multilocus genotypes to be meaningful.\n")
    }
    return(pop)
  }
  cutoff <- ifelse(cutoff > 0.5, 1 - cutoff, cutoff)
  min_ind = round(cutoff*nInd(pop))
  if(!isTRUE(quiet)){
    cat("cutoff value:", cutoff*100, "percent (",min_ind,"individuals ).\n")
  }
  if(pop@type == "PA"){
    # cutoff applies to too many or too few typed individuals in AFLP cases.
    locivals <- apply(pop@tab, 2, sum) %in% min_ind:(nInd(pop) - min_ind)
    if(!isTRUE(quiet)){
      if(all(locivals == TRUE)){
        cat("No sites found with fewer than", min_ind, 
            "different individuals.\n", fill = 80)
      }
      else{
        cat(sum(!locivals), "uninformative", 
            ifelse(sum(!locivals) > 1, "loci", "locus"), "found:", 
            pop@loc.names[!locivals],"\n", fill = 80)
      }
    }
    return(pop[, locivals])
  }
  else{
    # as.loci will put the population factor first when creating the data frame.
    if(is.null(pop@pop)){
      locivals <- apply(as.loci(pop), 2, test_table, min_ind, nInd(pop))
    }
    else{
      locivals <- apply(as.loci(pop)[-1], 2, test_table, min_ind, nInd(pop))
    }
    if(!isTRUE(quiet)){
      if(all(locivals == TRUE)){
        cat("No sites found with fewer than", min_ind, 
            "different individuals.\n", fill = 80)
      }
      else if(sum(locivals) < 2){
        cat("Fewer than 2 loci found informative. Perhaps you should choose a",
             "lower cutoff value?\nReturning with no changes.\n")
            return(pop)
      }
      else{
        cat(sum(!locivals), "uninformative", 
            ifelse(sum(!locivals) > 1, "loci", "locus"), "found:", 
            pop@loc.names[!locivals],"\n", fill = 80)
      }
    }
    return(pop[, loc = names(pop@loc.names[locivals])])
  }
}


