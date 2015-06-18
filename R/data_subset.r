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
#' Remove potential bias caused by cloned genotypes in genind or genclone 
#' object.
#' 
#' This function removes any duplicated multilocus genotypes from any specified 
#' population strata.
#' 
#' @param pop a \code{\link{genind}} object
#'   
#' @param strata a hierarchical formula or numeric vector. This will define the columns of the 
#'   data frame in the strata slot to use. 
#'   
#' @param combine \code{logical}. When set to TRUE, the strata will be 
#'   combined to create a new population for the clone-corrected genind or 
#'   genclone object.
#'   
#' @param keep \code{integer}. When \code{combine} is set to \code{FALSE}, you 
#'   can use this flag to choose the levels of your population strata. For 
#'   example: if your clone correction strata is set to "Pop", "Subpop", and 
#'   "Year", and you want to analyze your populations with respect to year, you 
#'   can set \code{keep = c(1,3)}.
#'   
#' @return a clone corrected \code{\linkS4class{genclone}} or
#'   \code{\linkS4class{genind}} object.
#'   
#' @details This function will clone correct based on the stratification 
#'   provided. To clone correct indiscriminately of population structure, set 
#'   \code{strata = NA}. 
#'   
#'   
#' @export
#' @author Zhian N. Kamvar
#' @examples
#' # LOAD A. euteiches data set
#' data(Aeut)
#' 
#' # Redefine it as a genclone object
#' Aeut <- as.genclone(Aeut)
#' strata(Aeut) <- other(Aeut)$population_hierarchy[-1]
#' 
#' # Check the number of multilocus genotypes
#' mlg(Aeut)
#' popNames(Aeut)
#' 
#' # Clone correct at the population level.
#' Aeut.pop <- clonecorrect(Aeut, strata =  ~Pop)
#' mlg(Aeut.pop)
#' popNames(Aeut.pop)
#' 
#' \dontrun{
#' # Clone correct at the subpopulation level with respect to population and
#' # combine.
#' Aeut.subpop <- clonecorrect(Aeut, strata = ~Pop/Subpop, combine=TRUE)
#' mlg(Aeut.subpop)
#' popNames(Aeut.subpop)
#' 
#' # Do the same, but set to the population level.
#' Aeut.subpop2 <- clonecorrect(Aeut, strata = ~Pop/Subpop, keep=1)
#' mlg(Aeut.subpop2)
#' popNames(Aeut.subpop2)
#' 
#' # LOAD H3N2 dataset
#' data(H3N2)
#' 
#' strata(H3N2) <- other(H3N2)$x
#'
#' # Extract only the individuals located in China
#' country <- clonecorrect(H3N2, strata =  ~country)
#' 
#' # How many isolates did we have from China before clone correction?
#' sum(strata(H3N2, ~country) == "China") # 155
#' 
#' # How many unique isolates from China after clone correction?
#' sum(strata(country, ~country) == "China") # 79
#' 
#' # Something a little more complicated. (This could take a few minutes on
#' # slower computers)
#' 
#' # setting the hierarchy to be Country > Year > Month  
#' c.y.m <- clonecorrect(H3N2, strata =  ~year/month/country)
#' 
#' # How many isolates in the original data set?
#' nInd(H3N2) # 1903
#' 
#' # How many after we clone corrected for country, year, and month?
#' nInd(c.y.m) # 1190
#' }
#==============================================================================#

clonecorrect <- function(pop, strata = 1, combine = FALSE, keep = 1){
  clonecall <- match.call()$pop
  if(!is.genind(pop) & !is(pop, "snpclone")){
    stop(paste(paste(substitute(pop), collapse=""), "is not a genind or snpclone object.\n"))
  }
  if (is.language(strata)){
    strataformula <- strata
    strata        <- all.vars(strata)
  }
  if (is.genind(pop)) popcall <- match.call()
  if (is.na(strata[1])){
    return(pop[.clonecorrector(pop), ])
  }

  if (is.numeric(strata)){
    strata        <- names(strata(pop))[strata]
    strataformula <- as.formula(paste0("~", paste(strata, collapse = "/")))
  }
  if (!all(strata %in% names(strata(pop)))){
    stop(hier_incompatible_warning(strata, strata(pop)))
  }
  setPop(pop) <- strataformula
  # Corrects the individual names of the object. This is fo the fact that the
  # clone corrector relies on the unique individual names for it to work.
  if(all(indNames(pop) == "")){
    indNames(pop) <- as.character(1:nInd(pop))
  }

  cpop <- nPop(pop)  

  # Steps for correction:
  # Subset by population factor.
  # Run subset population by the .clonecorrector
  # Profit!
  corWrecked <- function(x, pop){
    subbed <- popsub(pop, x) # population to be...corrected.
    subbed <- subbed[.clonecorrector(subbed), ] 
    # Return the indices based off of the individual names.
    return(which(indNames(pop) %in% indNames(subbed)))
  }
  
  ccpop <- unlist(lapply(1:cpop, corWrecked, pop))
  pop <- pop[ccpop, ]
  
  if(!combine){
    # When the combine flag is not true, the default is to keep the first level
    # of the strata. The keep flag is a numeric vector corresponding to the
    # strata flag indicating which levels the user wants to keep.
    strata <- strata[keep]
    newformula <- as.formula(paste0("~", paste(strata, collapse = "/")))
    setPop(pop) <- newformula
  }
  if (is.genind(pop)){
    pop@call <- popcall    
  }

  return(pop)
}


#==============================================================================#
#' Subset a \code{\linkS4class{genclone}} or \code{\linkS4class{genind}} object by population
#' 
#' Create a new dataset with specified populations or exclude specified
#' populations from the dataset.
#' 
#' @param gid a \code{\linkS4class{genclone}} or \code{\linkS4class{genind}} object.
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
#' @author Zhian N. Kamvar
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

  if (!is.genind(gid) & !is(gid, "genlight")){
    stop("popsub requires a genind or genlight object\n")
  }
  if (is.null(pop(gid))){
    if(sublist[1] != "ALL")
      warning("No population structure. Subsetting not taking place.")
    return(gid)
  }
  orig_list <- sublist 
  popnames  <- popNames(gid)
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
    if (length(popnames) == nPop(gid)){
      names(popnames) <- popNames(gid)
    } else {
      stop("Population names do not match population factors.")
    }
  }

  # Treating anything present in blacklist.
  if (!is.null(blacklist)){

    # If both the sublist and blacklist are numeric or character.
    if (is.numeric(sublist) & is.numeric(blacklist) | class(sublist) == class(blacklist)){
      sublist <- sublist[!sublist %in% blacklist]
    } else if (is.numeric(sublist) & class(blacklist) == "character"){
      # if the sublist is numeric and blacklist is a character. eg s=1:10, b="USA"
      sublist <- sublist[sublist %in% which(!popnames %in% blacklist)]
    } else {
      # no sublist specified. Ideal situation
      if(all(popnames %in% sublist)){
        sublist <- sublist[-blacklist]
      } else {
        # weird situation where the user will specify a certain sublist, yet
        # index the blacklist numerically. Interpreted as an index of
        # populations in the whole data set as opposed to the sublist.
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
      if (!is.numeric(orig_list) & !any(popNames(gid) %in% orig_list)){
        stop(unmatched_pops_warning(popNames(gid), orig_list))
      } else {
        nothing_warn <- paste("Nothing present in the sublist.\n",
                            "Perhaps the sublist and blacklist arguments have",
                            "duplicate entries?\n",
                            "Subset not taking place.")
        warning(nothing_warn)
        return(gid)
      }
    }
    gid <- gid[sublist, drop = drop]
    if (is.genind(gid)){
      gid@call <- match.call()      
    }
    return(gid)
  }
}

#==============================================================================#
#'How to deal with missing data in a genind object.
#'
#'missingno gives the user four options to deal with missing data.
#'
#'@param pop a \code{\linkS4class{genclone}} or \code{\linkS4class{genind}} 
#'  object.
#'  
#'@param type a \code{character} string: can be "ignore", "zero", "mean", 
#'  "loci", or "geno" (see \code{Details} for definitions).
#'  
#'@param cutoff \code{numeric}. A number from 0 to 1 indicating the allowable 
#'  rate of missing data in either genotypes or loci. This will be ignored for 
#'  \code{type} values of \code{"mean"} or \code{"zero"}.
#'  
#'@param quiet if \code{TRUE}, it will print to the screen the action performed.
#'  
#'@param freq defaults to \code{FALSE}. This option is passed on to the
#'  \code{\link[adegenet]{tab}} function. If \code{TRUE}, the matrix in the
#'  genind object will be replaced by a numeric matrix (as opposed to integer).
#'  THIS IS NOT RECOMMENDED. USE THE FUNCTION \code{\link[adegenet]{tab}}
#'  instead.
#'  
#'@details These methods provide a way to deal with systematic missing data and 
#'  to give a wrapper for \code{adegenet}'s \code{ \link{tab}} function. 
#'  ALL OF THESE ARE TO BE USED WITH CAUTION.
#'  
#'  \subsection{Treatment types}{ \itemize{ \item{\code{"ignore"} - does not 
#'  remove or replace missing data.} \item{\code{"loci"} - removes all loci 
#'  containing missing data in the entire data set. } \item{\code{"genotype"} - 
#'  removes any genotypes/isolates/individuals with missing data.} 
#'  \item{\code{"mean"} - replaces all NA's with the mean of the alleles for the
#'  entire data set.} \item{\code{"zero"} or \code{"0"} - replaces all NA's with
#'  "0". Introduces more diversity.}}}
#'@return a \code{\linkS4class{genclone}} or \code{\linkS4class{genind}} object.
#'  
#'@note \emph{"wild missingno appeared!"}
#'  
#'@seealso \code{\link{tab}}, \code{\link{poppr}}, \code{\link{poppr.amova}},
#'  \code{\link{nei.dist}}, \code{\link{aboot}}
#'  
#'@export
#'@author Zhian N. Kamvar
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
#' # Replacing all NA with "0" (see tab in the adegenet package).
#' nancy.0 <- missingno(nancycats, type = "0")
#'
#' ## Replaced 617 missing values 
#' 
#' # Replacing all NA with the mean of each column (see tab in the
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
missingno <- function(pop, type = "loci", cutoff = 0.05, quiet=FALSE, freq = FALSE){
  if(sum(is.na(pop@tab)) > 0){
    # removes any loci (columns) with missing values.
    MISSINGOPTS <- c("loci", "genotypes", "mean", "zero", "0", "ignore")
    freq_warning <- paste("Objects of class 'genind' must have integers in the",
      "genotype matrix. Setting freq = TRUE will force the matrix to be numeric.",
      "please see help('tab') for alternatives.")
    type        <- match.arg(tolower(type), MISSINGOPTS)
    if (type == "ignore"){
      return(pop)
    }
    navals <- percent_missing(pop, type = type, cutoff = cutoff)
    if (type == "loci"){
      if(quiet != TRUE){
        if (length(navals) == ncol(tab(pop))){
          message("\n No loci with missing values above ",
              paste0(cutoff*100,"%")," found.\n")
        } else {
          remloc <- locNames(pop)[!cumsum(pop@loc.nall) %in% navals]
          message("\n Found ", sum(is.na(tab(pop)))," missing values.")
          loci   <- paste(length(remloc), ifelse(length(remloc) == 1, "locus", 
                          "loci"))
          cat("\n", loci, "contained missing values greater than",
              paste0(cutoff*100,"%."))
          cat("\n Removing", loci, ":", remloc,"\n", fill = 80)
        }
      }
      pop <- pop[, navals]
    } else if (type == "genotypes"){
      if(quiet != TRUE){
        if (length(navals) == nInd(pop)){
          message("\n No genotypes with missing values above ",
              paste0(cutoff*100, "%")," found.\n")
        } else {
          remgeno <- indNames(pop)[-navals]
          message("\n Found ", sum(is.na(tab(pop)))," missing values.")
          geno    <- paste(length(remgeno), ifelse(length(remgeno) == 1, 
                           "genotype", "genotypes"))
          cat("\n", geno, "contained missing values greater than",
              paste0(cutoff*100,"%."))
          cat("\n Removing", geno, ":", remgeno,"\n", fill = 80)
        }
      }
      pop <- pop[navals, ]
    } else if (type == "mean"){
      if (!quiet){
        message("\n Replaced ", sum(is.na(tab(pop)))," missing values.")
      }
      pop@tab <- tab(pop, freq = freq, NA.method = "mean", quiet=quiet)
      if (freq){
        warning(freq_warning)
      }
    }
    # changes all NA's to 0. NOT RECOMMENDED. INTRODUCES MORE DIVERSITY.
    else if (type %in% c("zero","0")){
      if (!quiet){
        message("\n Replaced ", sum(is.na(tab(pop)))," missing values.")
      }
      pop@tab <- tab(pop, freq = freq, NA.method = "zero", quiet=quiet)

      if (freq){
        warning(freq_warning)
      }
    } else {
      message("\nInvalid type.\n")
    } 
  } else {
    if(quiet == FALSE){
      message("\n No missing values detected.\n")
    }
  }
  return(pop)
}

# #==============================================================================#
# #~ Split a or combine items within a data frame in \code{\link{genind}} objects (DEPRECATED).
# #~
# #~ Often, one way a lot of file formats fail is that they do not allow multiple
# #~ population hierarchies. This can be circumvented, however, by coding all of
# #~ the hierarchies in one string in the input file with a common separator (eg.
# #~ "_"). \code{splitcombine} will be able to recognise those separators and
# #~ create a data frame of all the population structures for whatever sub-setting
# #~ you might need. 
# #~
# #~ @param pop a \code{\link{genind}} object.
# #~
# #~ @param method an \code{integer}, 1 for splitting, 2 for combining.
# #~
# #~ @param dfname the name of the data frame containing the population structure.
# #~ for the splitting method, the combined population structure must be in the
# #~ first column. 
# #~
# #~ @param sep The separator used for separating or combining the data. See note.
# #~
# #~ @param hier a \code{vector} containing the population hierarchy you wish to
# #~ split or combine. 
# #~
# #~ @param setPopulation \code{logical}. if \code{TRUE}, the population of the
# #~ resulting genind object will be that of the highest population structure
# #~ (split method) or the combined populations (combine method).
# #~
# #~ @param fixed \code{logical}. An argument to be passed onto
# #~ \code{\link{strsplit}}. If \code{TRUE}, \code{sep} must match exactly to the
# #~ populations for the split method. 
# #~
# #~ @return a \code{\link{genind}} object with a modified data frame in the
# #~ \code{\link{other}} slot.
# #~
# #~ @note 
# #~ This function has been deprecated and replaced by functions like
# #~ \code{\link{splitStrata}}. Please consider using the
# #~ \code{\linkS4class{genclone}} object for storing hierarchies. 
# #~
# #~ The separator
# #~ field is sensitive to regular expressions. If you do not know what those are,
# #~ please use the default underscore to separate your populations. Use \code{fixed
# #~ = TRUE} to ignore regular expressions. If you do not set the \code{hier} flag
# #~ for the split method, your new data frame will have the names "comb", "h1", "h2"
# #~ and so on; for the combine method, your data frame will return the first column
# #~ of your data frame.
# #~
# #~ @export
# #~ @author Zhian N. Kamvar
# #~ @examples
# #~ \dontrun{
# #~ # Method 1: Splitting.
# #~ Aeut <- read.genalex(system.file("files/rootrot.csv", package="poppr"))
# #~ 
# #~ # We have 19 different "populations", but really, there is a hierarchy.
# #~ popNames(Aeut)
# #~ 
# #~ # Let's split them up. The default data frame from read.genalex is the same
# #~ # as the default for this function. 
# #~ Aeut <- splitcombine(Aeut, hier=c("Pop", "Subpop"))
# #~ 
# #~ # Much better!
# #~ popNames(Aeut)
# #~
# #~ # Method 2: Combining.
# #~ 
# #~ data(H3N2)
# #~ # Create a new data set combining the population factors of year and country
# #~ H.comb <- splitcombine(H3N2, method=2, dfname="x", hier=c("year", "country"))
# #~ 
# #~ # Checking to make sure they were actually combined.
# #~ head(H.comb$other$x$year_country)
# #~ # Creating new data frame in the object to mess around with. 
# #~ H.comb$other$year_country <- data.frame(H.comb$other$x$year_country)
# #~ 
# #~ # Splitting those factors into their original components and setting the
# #~ # population to year.
# #~ H.comb <- splitcombine(H.comb, method=1, dfname="year_country", hier=c("year", "country"))
# #~ }
# #==============================================================================#
# splitcombine <- function(pop, method=1, dfname="population_hierarchy", sep="_", hier=1, setPopulation=TRUE, fixed=TRUE){
  
#   # As this function is deprecated, the following lines will attempt to give the
#   # user equivalent code in poppr >v.1.1.0 
#   ##########
#   sccall <- substitute(pop)
#   depmsg <- paste0("Use:\n'", sccall, " <- as.genclone(", sccall, ", hier = ")
#   if (method == 1){
#     newhier <- paste0("other(", sccall, ")$", dfname, "[1]")
#     if (!is.numeric(hier)){
#       newform <- paste0("~", paste(hier, collapse = "/"))
#       sugg <- paste0(newhier, "); splitStrata(", sccall, ") <- ", newform, "'\n")
#     } else {
#       sugg <- paste0(newhier, ")' and 'splitStrata'\n")
#     }
#   } else {
#     newhier <- paste0("other(", sccall, ")$", dfname)
#     if (!is.numeric(hier)){
#       newform <- paste0("~", paste(hier, collapse = "/"))
#       sugg <- paste0(newhier, "); setPop(", sccall, ") <- ", newform, "'\n")
#     } else {
#       sugg <- paste0(newhier, ")' and 'setPop'\n") 
#     }
#   }
#   dmsg <- c("'splitcombine' is deprecated.\n\n",
#             paste0(depmsg, sugg), "\n",
#             "See help(\"genclone\"), help(\"setPop\"), ",
#             "and help(\"splitStrata\") for details.")
#   .Deprecated("genclone", package = "poppr", dmsg)
#   ###########
  
#   if (!is.genind(pop)){
#     stop(paste(paste(substitute(pop), collapse=""), "is not a genind object.\n"))
#   }
#   if (!is.data.frame(pop$other[[dfname]])){
#     stop(paste("There is no data frame in ",paste(substitute(pop), collapse=""), "@other called ",dfname,".\n", sep=""))
#   }
#   METHODS = c("Split", "Combine")
#   if (all((1:2)!=method)) {
#     cat("1 = Split\n")
#     cat("2 = Combine\n")
#     cat("Select an integer (1 or 2): ")
#     method <- as.integer(readLines(n = 1))
#   }
#   if (all((1:2)!=method)) (stop ("Non convenient method number"))
#   # Splitting !
#   if(method == 1){
#     # Remember, this acts on the first column of the data frame.
#     # It will split the population factors in that data frame by sep. 
#     df <- pop_splitter(pop$other[[dfname]], sep=sep)
    
#     # This makes sure that the length of the given hierarchy is the same as the
#     # hierarchy in the original population factor. It's renaming the original
#     # factor the population hierarchy separated by sep.
#     # eg. df: Pop_Subpop_Year, Pop, Subpop, Year
#     if(length(df) - 1 == length(hier)){
#       names(df) <- c(paste(hier, collapse=sep), hier)
#     }
#     # In the case that the number of columns added to the new df is the same
#     # as the number of hierarchies specified.
#     else if(length(df) - length(hier) == length(pop$other[[dfname]]) ){
#       names(df)[(length(df) - length(hier) + 1):length(df)] <- hier
#     }

#     # Checking to see if there was only one column in the original data frame
#     # This is necessary to avoid overwriting data.
#     if(length(pop$other[[dfname]] == 1)){
#       pop$other[[dfname]] <- df
#     }
    
#     # If there are other columns in the data frame, we check to see if any of
#     # the names in the new data frame are the same. 
#     else if(any(names(pop$other[[dfname]]) %in% names(df[-1]))){
#       # Removing the combined column. It will remain the same.
#       df <- df[-1]
#       # Gathering the index of names that are the same as the split hierarchy.
#       dfcols <- which(names(pop$other[[dfname]]) %in% names(df))
      
#       # All the names are equal, so we replace those columns with the newly
#       # subsetted columns.
#       if(length(names(df)) == length(dfcols))
#         pop$other[[dfname]][dfcols] <- df
      
#       # Not all columns match, so you replace the ones that do and then bind
#       # the new ones.
#       else{
#         popothernames <- which(names(df) %in% names(pop$other[[dfname]][dfcols]))
#         pop$other[[dfname]][dfcols] <- df[popothernames]
#         pop$other[[dfname]] <- cbind(pop$other[[dfname]], df[-dfcols])
#       }
#     }
    
#     # If there are no names in the new data frame that match the original,
#     # simply tack the new data frame on to the end. These will have the names
#     # h1, h2, h3, etc.
#     else{  
#       pop$other[[dfname]] <- cbind(pop$other[[dfname]], df[-1])
#     }
    
#     # Set the population to the highest level of the hierarchy.
#     if(setPopulation == TRUE){
#       hier <- ifelse(is.numeric(hier), as.character(hier), hier)
#       pop(pop) <- pop$other[[dfname]][[hier[1]]]
#       names(popNames(pop)) <- levels(pop$pop)
#     }
#     return(pop)
#   }
  
#   # Combining !
#   else if(method == 2){
#     newdf <- pop_combiner(pop$other[[dfname]], hier=hier, sep=sep)
#     if(all(is.na(newdf))){
#       stop(paste("\n\nError in combining population factors.\nCheck your hier flag and make sure that the columns",paste(hier, collapse=" and "),
#                  "exist within the data frame called",dfname,
#                  "in the @other slot of your genind object.\n\n"))
#     }
#     # Combining the hierarchy. This will no longer give you numbers as names,
#     # rather it will return the names.
#     pop$other[[dfname]][[paste(names(pop$other[[dfname]][hier]), collapse=sep)]] <- newdf
#     if(setPopulation == TRUE){
#       pop(pop) <- newdf
#       names(popNames(pop)) <- levels(pop$pop)
#     }
#     return(pop)
#   }
# }
#==============================================================================#
#' Remove all non-phylogentically informative loci
#' 
#' This function will facilitate in removing phylogenetically uninformative loci
#' from a \code{\linkS4class{genclone}} or \code{\linkS4class{genind}} object. 
#' The user has the ability to define what uninformative means by setting a 
#' cutoff value for either percentage of differentiating genotypes or minor 
#' allele frequency.
#' 
#' @param pop a \code{\linkS4class{genclone}} or \code{\linkS4class{genind}} 
#'   object.
#'   
#' @param cutoff \code{numeric}. A number from 0 to 1 defining the minimum 
#'   number of differentiating samples.
#'   
#' @param MAF \code{numeric}. A number from 0 to 1 defining the minimum minor
#'   allele frequency. This is passed as the \code{thresh} parameter of
#'   \code{\link[adegenet]{isPoly}}.
#'   
#' @param quiet \code{logical}. When \code{quiet = TRUE} (default), messages 
#'   indicating the loci removed will be printed to screen. When \code{quiet = 
#'   FALSE}, nothing will be printed to screen.
#'   
#' @return A \code{genind} object with user-defined informative loci.
#'   
#' @details This function will remove uninformative loci using a traditional MAF
#'   cutoff (using \code{\link[adegenet]{isPoly}} from \pkg{adegenet}) as well
#'   as analyzing the number of observed genotypes in a locus. This is important
#'   for clonal organisms that can have fixed heterozygous sites not detected by
#'   MAF methods.
#'   
#' @note This will have a few side effects that affect certain analyses. First, 
#'   the number of multilocus genotypes might be reduced due to the reduced 
#'   number of markers (if you are only using a genind object). Second, if you 
#'   plan on using this data for analysis of the index of association, be sure 
#'   to use the standardized version (rbarD) that corrects for the number of 
#'   observed loci.
#'   
#' @author Zhian N. Kamvar
#' @examples
#' # We will use a dummy data set to demonstrate how this detects uninformative
#' # loci using both MAF and a cutoff.
#' 
#' genos <- c("A/A", "A/B", "A/C", "B/B", "B/C", "C/C")
#' 
#' v <- sample(genos, 100, replace = TRUE)
#' w <- c(rep(genos[2], 99), genos[3])           # found by cutoff
#' x <- c(rep(genos[1], 98), genos[3], genos[2]) # found by MAF
#' y <- c(rep(genos[1], 99), genos[2])           # found by both
#' z <- sample(genos, 100, replace = TRUE)
#' dat <- df2genind(data.frame(v = v, w = w, x = x, y = y, x = z), sep = "/")
#' 
#' informloci(dat)
#' 
#' \dontrun{
#' # Ignore MAF
#' informloci(dat, MAF = 0)
#' 
#' # Ignore cutoff
#' informloci(dat, cutoff = 0)
#' 
#' # Real data
#' data(H3N2)
#' informloci(H3N2)
#' 
#' }
#' @export
#==============================================================================#
#' @importFrom pegas as.loci
informloci <- function(pop, cutoff = 2/nInd(pop), MAF = 0.01, quiet = FALSE){
  if (!is.genind(pop)){
    stop("This function only works on genind objects.")
  }
  MLG <- mlg(pop, quiet = TRUE)
  if (MLG < 3){
    if(!isTRUE(quiet)){
      cat("Not enough multilocus genotypes to be meaningful.\n")
    }
    return(pop)
  }
  cutoff <- ifelse(cutoff > 0.5, 1 - cutoff, cutoff)
  MAF    <- ifelse(MAF > 0.5, 1 - MAF, MAF)
  min_ind = round(cutoff * nInd(pop))
  if (!isTRUE(quiet)){
    ind <- ifelse(min_ind == 1, "sample", "samples")
    message("cutoff value: ", cutoff*100, " % ( ",min_ind, " ", ind," ).")
    message("MAF         : ", MAF)
  }
  if (pop@type == "PA"){
    # cutoff applies to too many or too few typed individuals in AFLP cases.
    locivals <- apply(pop@tab, 2, sum) %in% min_ind:(nInd(pop) - min_ind)
    if(!isTRUE(quiet)){
      if(all(locivals == TRUE)){
        message("No sites found with fewer than", min_ind, 
                "different samples.")
      }
      else{
        message(sum(!locivals), "uninformative", 
            ifelse(sum(!locivals) > 1, "loci", "locus"), "found:", 
            paste(locNames(pop)[!locivals], collapse = ", "))
      }
    }
    return(pop[, locivals])
  } 
  # as.loci will put the population factor first when creating the data frame.
  if (is.null(pop(pop))){
    glocivals <- apply(as.loci(pop), 2, test_table, min_ind, nInd(pop))
  } else {
    glocivals <- apply(as.loci(pop)[-1], 2, test_table, min_ind, nInd(pop))
  }

  alocivals <- isPoly(pop, "locus", thres = MAF)
  
  locivals  <- alocivals & glocivals
  
  if (all(locivals == TRUE)){
    msg <- paste("\nAll sites polymorphic")
  } else if (sum(locivals) < 2){
    msg <- paste0("\nFewer than 2 loci found informative.",
                  "\nPerhaps you should choose a ",
                  "lower cutoff value?\nReturning with no changes.")
    locivals <- rep(TRUE, nLoc(pop))
  } else {
    glocsum <- sum(!glocivals)
    alocsum <- sum(!alocivals)
    locsum  <- sum(!locivals)
    fmsg <- paste("Found", locsum, "uninformative", 
                  ifelse(locsum != 1, "loci", "locus"), "\n",
                  "============================")
    gmsg <- paste(glocsum, 
                  ifelse(glocsum != 1, "loci", "locus"), "found with",
                  "a cutoff of", min_ind, ind, 
                  ifelse(glocsum == 0, "", ":\n"),
                  paste(locNames(pop)[!glocivals], collapse = ", "))
    amsg <- paste(alocsum, 
                  ifelse(alocsum != 1, "loci", "locus"), "found with",
                  "found with MAF <", signif(MAF, 3), 
                  ifelse(alocsum == 0, "", ":\n"),
                  paste(locNames(pop)[!alocivals], collapse = ", "))
    msg <- paste("\n", fmsg, "\n", gmsg, "\n", amsg)
  }

  if (!isTRUE(quiet)){
    message(msg)
  }

  return(pop[, loc = locNames(pop)[locivals]])
}
#==============================================================================#
#' Recode polyploid microsatellite data for use in frequency based statistics.
#' 
#' As the genind object requires ploidy to be consistent across loci, a 
#' workaround to importing polyploid data was to code missing alleles as "0" 
#' (for microsatellite data sets). The advantage of this is that users would be 
#' able to calculate Bruvo's distance, the index of association, and genotypic 
#' diversity statistics. The tradeoff was the fact that this broke all other 
#' analyses as they relied on allele frequencies and the missing alleles are 
#' treated as extra alleles. This function removes those alleles and returns a 
#' \code{\linkS4class{genclone}} or \code{\linkS4class{genind}} object where 
#' allele frequencies are coded based on the number of alleles observed at a 
#' single locus per individual. See the examples for more details.
#' 
#' @param poly a \code{\linkS4class{genclone}}, \code{\linkS4class{genind}}, or
#'   \code{\linkS4class{genpop}} object that has a ploidy of > 2
#' @param newploidy for genind or genclone objects: if \code{FALSE} (default), the user-defined ploidy will stay
#'   contstant. if \code{TRUE}, the ploidy for each sample will be determined by
#'   the maximum ploidy observed for each genotype.
#'   
#' @param addzero add zeroes onto genind or genclone objects with uneven ploidy?
#' if \code{TRUE}, objects with uneven ploidies will have zeroes appended to all
#' loci to allow conversion to genpop objects. Defaults to \code{FALSE}.
#'   
#' @return a \code{\linkS4class{genclone}}, \code{\linkS4class{genind}}, or
#'   \code{\linkS4class{genpop}} object.
#'   
#' @details The genind object has two caveats that make it difficult to work 
#'   with polyploid data sets: \enumerate{\item ploidy must be constant
#'   throughout the data set \item missing data is treated as "all-or-none"} In
#'   an ideal world, polyploid genotypes would be just as unambigouous as
#'   diploid or haploid genotypes. Unfortunately, the world we live in is far
#'   from ideal and a genotype of AB in a tetraploid organism could be AAAB,
#'   AABB, or ABBB. In order to get polyploid data in to \pkg{adegenet} or
#'   \pkg{poppr}, we must code all loci to have the same number of allelic
#'   states as the ploidy or largest observed heterozygote (if ploidy is
#'   unknown). The way to do this is to insert zeroes to pad the alleles. So, to
#'   import two genotypes of:
#'
#' \tabular{rrrr}{
#' NA \tab 20 \tab 23 \tab 24\cr
#' 20 \tab 24 \tab 26 \tab 43
#' } 
#' they should be coded as: 
#' \tabular{rrrr}{
#'  0 \tab 20 \tab 23 \tab 24\cr
#' 20 \tab 24 \tab 26 \tab 43
#' } 
#' This zero is treated as an extra allele and is represented in the genind object as so:
#' \tabular{rrrrrr}{
#' \strong{0} \tab \strong{20} \tab \strong{23} \tab \strong{24} \tab \strong{26} \tab \strong{43}\cr
#' 1 \tab 1 \tab 1 \tab 1 \tab 0 \tab 0\cr
#' 0 \tab 1 \tab 0 \tab 1 \tab 1 \tab 1
#' }
#' 
#'   This function remedies this problem by removing the zero column.
#'   The above table would become:
#'
#' \tabular{rrrrr}{
#' \strong{20} \tab \strong{23} \tab \strong{24} \tab \strong{26} \tab \strong{43}\cr
#' 1 \tab 1 \tab 1 \tab 0 \tab 0\cr
#' 1 \tab 0 \tab 1 \tab 1 \tab 1
#' }
#' 
#' With this, the user is able to calculate frequency based statistics on the
#' data set.
#' 
#' @note This is an approximation, and a bad one at that. \pkg{Poppr} was not
#' originally intended for polyploids, but with the inclusion of Bruvo's
#' distance, it only made sense to attempt something beyond single use.
#'
#' @author Zhian N. Kamvar
#' @export
#' @examples
#' data(Pinf)
#' iPinf <- recode_polyploids(Pinf)
#' 
#' # Note that the difference between the number of alleles.
#' Pinf@@loc.nall
#' iPinf@@loc.nall
#' 
#' \dontrun{
#' library("ape")
#' 
#' # Removing missing data. 
#' setPop(Pinf) <- ~Country
#' 
#' # Calculating Rogers' distance. 
#' rog <- rogers.dist(genind2genpop(Pinf))
#' irog <- rogers.dist(recode_polyploids(genind2genpop(Pinf)))
#' 
#' # We will now plot neighbor joining trees. Note the decreased distance in the
#' # original data.
#' plot(nj(rog), type = "unrooted")
#' add.scale.bar(lcol = "red", length = 0.02)
#' plot(nj(irog), type = "unrooted")
#' add.scale.bar(lcol = "red", length = 0.02)
#' }
#==============================================================================#
recode_polyploids <- function(poly, newploidy = FALSE, addzero = FALSE){
  if (!is.genind(poly) && !is.genpop(poly)){
    stop("input must be a genind object.")
  } else if (!test_zeroes(poly) && addzero == FALSE){
    warning("Input is not a polyploid data set, returning original.")
    return(poly)
  }
  if (addzero && is.genind(poly)){
    
    if (test_zeroes(poly)){
      warning("addzero = TRUE, but data already has zeroes. Returning original.")
      return(poly)
    }
    polydf   <- genind2df(poly, sep = "/", usepop = FALSE)
    maxploid <- max(ploidy(poly))
    polydf   <- generate_bruvo_mat(polydf, maxploid, sep = "/")
    res      <- df2genind(polydf, sep = "/", ploidy = maxploid,
                          pop = pop(poly), strata = strata(poly), 
                          hierarchy = poly@hierarchy)
    other(res) <- other(poly)
    if (is.genclone(poly)){
      res <- as.genclone(res, mlg = poly@mlg)
    }
    return(res)
  }
  MAT <- tab(poly)
  fac <- poly@loc.fac

  non_zero_cols_list   <- lapply(poly@all.names, function(x) as.numeric(x) > 0)
  non_zero_cols_vector <- unlist(non_zero_cols_list, use.names = FALSE)

  poly@loc.fac   <- fac[non_zero_cols_vector]
  poly@loc.nall  <- setNames(tabulate(poly@loc.fac), locNames(poly))
  poly@tab       <- MAT[, non_zero_cols_vector, drop = FALSE]
  poly@all.names <- mapply("[", poly@all.names, non_zero_cols_list,
                           SIMPLIFY = FALSE)

  if (newploidy && is.genind(poly)){
    ploc         <- seploc(poly, res.type = "matrix")
    ploid_mat    <- vapply(ploc, FUN = apply, FUN.VALUE = integer(nInd(poly)), 
                           MARGIN = 1, sum)
    ploidy(poly) <- apply(ploid_mat, MARGIN = 1, max, na.rm = TRUE)
  }
  return(poly)
}



