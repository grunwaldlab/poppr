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
#' Oomycete root rot pathogen \emph{Aphanomyces euteiches} AFLP data
#' 
#' @name Aeut
#' @docType data
#' @usage data(Aeut)
#' @description The Aeut dataset consists of 187 isolates of the Oomycete root 
#'   rot pathogen, \emph{Aphanomyces euteiches} collected from two different 
#'   fields in NW Oregon and W Washington, USA.
#' @format a \code{\link{genind}} object with two populations containing a data
#'   frame in the \code{other} slot called \code{population_hierarchy}. This
#'   data frame gives indices of the populations and subpopulations for the data
#'   set.
#' @references Grunwald, NJ and Hoheisel, G.A. 2006. Hierarchical Analysis of
#'   Diversity, Selfing, and Genetic Differentiation in Populations of the 
#'   Oomycete \emph{Aphanomyces euteiches}. Phytopathology 96:1134-1141
#==============================================================================#
NULL
#==============================================================================#
#' Simulated data illustrating a Minimum Spanning Network based on Bruvo's
#' Distance
#' 
#' @name partial_clone
#' @docType data
#' @usage data(partial_clone)
#' @description These data were simulated using SimuPOP version 1.0.8 with 
#'   99.9\% clonal reproduction over 10,000 generations. Populations were
#'   assigned post-hoc and are simply present for the purposes of demonstrating
#'   a minimum spanning network with Bruvo's distance.
#' @format a \code{\link{genind}} object with 50 individuals, 10 loci, and four 
#'   populations.
#' @references Bo Peng and Christopher Amos (2008) Forward-time simulations of 
#'   nonrandom mating populations using simuPOP. \emph{bioinformatics}, 24 (11):
#'   1408-1409.
#==============================================================================#
NULL
#==============================================================================#
#' Phytophthora infestans data from Mexico and South America.
#' 
#' @name Pinf
#' @docType data
#' @usage data(Pinf)
#' @description The Pinf data set contains 86 isolates genotyped over 11 
#'   microsatellite loci collected from Mexico, Peru, Columbia, and Ecuador.
#'   This is a subset of the data used for the reference below.
#' @format a \code{\linkS4class{genclone}} object with 2 hierarchical levels
#'   called "Continent" and "Country" that contain 2 and 4 populations,
#'   respectively.
#' @references Goss, Erica M., Javier F. Tabima, David EL Cooke, Silvia
#'   Restrepo, William E. Fry, Gregory A. Forbes, Valerie J. Fieland, Martha
#'   Cardenas, and Niklaus J. Grünwald. "The Irish potato famine pathogen
#'   \emph{Phytophthora infestans} originated in central Mexico rather than the Andes."
#'   Proceedings of the National Academy of Sciences 111:8791-8796.
#==============================================================================#
NULL
#==============================================================================#
#' Peach brown rot pathogen \emph{Monilinia fructicola}
#' 
#' @name monpop
#' @docType data
#' @usage data(monpop)
#' @description This is microsatellite data for a population of the haploid 
#'   plant pathogen \emph{Monilinia fructicola} that causes disease within peach
#'   tree canopies (Everhart & Scherm, 2014). Entire populations within trees
#'   were sampled across 3 years (2009, 2010, and 2011) in a total of four
#'   trees, where one tree was sampled in all three years, for a total of 6
#'   within-tree populations. Within each year, samples in the spring were taken
#'   from affected blossoms (termed "BB" for blossom blight) and in late summer
#'   from affected fruits (termed "FR" for fruit rot). There are a total of 694 
#'   isolates with 65 to 173 isolates within each canopy population that were 
#'   characterized using a set of 13 microsatellite markers.
#' @format a \code{\linkS4class{genclone}} object with 3 hierarchical levels 
#'   coded into one population factor. These are named "Tree", "Year", and 
#'   "Symptom"
#' @references SE Everhart, H Scherm, (2014) Fine-scale genetic structure of 
#'   \emph{Monilinia fructicola} during brown rot epidemics within individual peach 
#'   tree canopies. Phytopathology 105:542-549 doi:10.1094/PHYTO-03-14-0088-R
#' @examples
#' data(monpop)
#' splithierarchy(monpop) <- ~Tree/Year/Symptom
#' setpop(monpop) <- ~Symptom/Year
#' monpop
#==============================================================================#
NULL
#==============================================================================#
# This function will extract all relevant information from the files
# 
# Public functions utilizing this function:
# # none
#
# Internal functions utilizing this function:
# # new.poppr (in testing)
#==============================================================================#
extract.info <- function(x) {
	if (length(grep("^clone.+?dat$", x$File)) != 0){
		# Rate of clonal reproduction
		x$Clone <- as.numeric(sub("^clone.(\\d{3}.\\d{2}).+?dat$","\\1", x$File))
		# Rate of Sexual reproduction
		x$Sex.Rate <- (100-x$Clone)/100
	}
	if (length(grep(".+?rep.+?dat$", x$File)) != 0){
		# Replicate indicators
		x$Replicate <- sub(".+?rep.(\\d{2}).+?dat$", "\\1", x$File)
	}
	if (length(grep(".+?pop.+?dat$", x$File)) != 0){
		# Population size indicators
		x$Pop.Size <- sub(".+?pop.(\\d+?).+?dat$","\\1", x$File)
	}
	if (length(grep(".+?sam.+?dat$", x$File)) != 0){
		# Sample size
		x$Samp.Size <- sub(".+?sam.+?(\\d{2,3}).+?dat$", "\\1", x$File)
	}
	return(x)
}
#==============================================================================#
# This function will attempt to convert external files of the following types:
#
# Structure (*.str, *.stru)
# Fstat (*.dat)
# Gentix (*.gtx)
# Genpop (*.gen)
# Genalex (*.csv)
#
# The output is a poppr object and the original filename that was imported. 
# Missing data is handled via the function missingno and the ability for clone
# correction is also possible.
# If quiet is set to false and you are importing non-genind or poppr objects,
# you will see many warnings. 
#
# Public functions utilizing this function:
# # poppr
#
# Internal functions utilizing this function:
# # new.poppr (in testing)
#==============================================================================#
process_file <- function(input, quiet=TRUE, missing="ignore", cutoff=0.05, keep=1,
                            clonecorrect=FALSE, hier=c(1), dfname="hier"){
  if (!is.genind(input)){
    x <- input
    if (toupper(.readExt(x)) == "CSV"){
      try(input <- read.genalex(x), silent=quiet)
      try(input <- read.genalex(x, region=TRUE), silent=quiet)
      try(input <- read.genalex(x, geo=TRUE), silent=quiet)
      try(input <- read.genalex(x, geo=TRUE, region=TRUE), silent=quiet)
    } else {
      try(input <- import2genind(x, quiet=quiet), silent=quiet)
    }
    stopifnot(is.genind(input))
    input@call[2] <- x
    popcall       <- input@call
    input         <- missingno(input, type=missing, cutoff=cutoff, quiet=quiet)
    input@call    <- popcall
    if (clonecorrect == TRUE){
      poplist    <- clonecorrect(input, hier=hier, dfname=dfname, keep=keep)
      input      <- poplist
      input@call <- popcall
    }
  } else if (is.genind(input)) {
    x         <- as.character(match.call()[2])
    popcall   <- input@call
    input     <- missingno(input, type=missing, cutoff=cutoff, quiet=quiet)
    if (clonecorrect == TRUE){
      poplist    <- clonecorrect(input, hier=hier, dfname=dfname, keep=keep)
      input      <- poplist
      input@call <- popcall
    }
  }
  return(list(X=x, GENIND=input))
}

#==============================================================================#
# .clonecorrector will simply give a list of individuals (rows) that are
# duplicated within a genind object. This can be used for clone correcting a
# single genind object.
#
# Public functions utilizing this function:
# # clonecorrect, bruvo.msn
#
# Internal functions utilizing this function:
# # none
#==============================================================================#

.clonecorrector <- function(x){
  if (is.genclone(x)){
    is_duplicated <- duplicated(x@mlg)
  } else {
    is_duplicated <- duplicated(x@tab[, 1:ncol(x@tab)])
  }
  res <- -which(is_duplicated)
  # conditional for the case that all individuals are unique.
  if(is.na(res[1])){
    res <- which(!is_duplicated)
  }
  return(res)
}

#==============================================================================#
# geno.na will find the genotypes in the population that contain na's and 
# remove them.
#
# Public functions utilizing this function:
# # none
#
# Internal functions utilizing this function:
# # percent_missing
#
# DEPRECATED
#==============================================================================#

geno.na <- function(pop){
  pop2 <- -unique(which(is.na(pop@tab), arr.ind=T)[,1])  
  if(is.na(pop2[1])){
    return(unique(which(!is.na(pop@tab), arr.ind=T)[,1]))
  } else {
    return(pop2)
  } 
}

#==============================================================================#
# loci.na will find the loci in the population that contain na's and remove
# them.
# Public functions utilizing this function:
# # none
#
# Internal functions utilizing this function:
# # percent_missing
#
# DEPRECATED
#==============================================================================#

loci.na <- function(pop) {
  pop2 <- -unique(which(is.na(pop@tab), arr.ind=T)[,2])  
  if(is.na(pop2[1])){
    return(unique(which(!is.na(pop@tab), arr.ind=T)[,2]))
  } else {
    return(pop2)
  } 
}

#==============================================================================#
# This will remove either loci or genotypes containing missing values above the
# cutoff percent.
# 
# Public functions utilizing this function:
# # missingno
#
# Internal functions utilizing this function:
# # none.
#==============================================================================#

percent_missing <- function(pop, type="loci", cutoff=0.05){
  if (toupper(type) == "LOCI"){
    missing_loci        <- 1 - propTyped(pop, "loc")
    names(missing_loci) <- levels(pop@loc.fac)
    missing_loci        <- missing_loci[missing_loci > cutoff]
    misslist            <- 1:ncol(pop@tab)
    filter              <- !pop@loc.fac %in% names(missing_loci)
  } else {
    missing_geno <- 1 - propTyped(pop, "ind")
    misslist     <- 1:nInd(pop)
    filter       <- missing_geno <= cutoff
  }
  return(misslist[filter])
}

#==============================================================================#
# This implements rounding against the IEEE standard and rounds 0.5 up
# Public functions utilizing this function:
# # none
#
# Internal functions utilizing this function:
# # .PA.pairwise.differences, .pairwise.differences
#==============================================================================#

round.poppr <- function(x){
  if (x - as.integer(x) == 0.5 & as.integer(x)%%2 == 0)
    x <- round(x) + 1
  else if(-x + as.integer(x) == 0.5 & as.integer(x)%%2 == 0)  
    x <- round(x) - 1
  else
    x <- round(x)
  return(x)
}


#==============================================================================#
# This will caluclulate p-values for permutation tests. 
# Public functions utilizing this function:
# # ia
#
# Internal functions utilizing this function:
# # .ia
# DEPRECATED
#==============================================================================#

ia.pval <- function(index="index", sampled, observed){
  if (all(is.nan(sampled[[index]]))){
    return(NA)
  }
  pval <- mean(ifelse(!is.na(sampled[[index]]) & sampled[[index]] >= observed,1,0))
  return(pval)
}

#==============================================================================#
# This will be used to split heirarchical population vectors that are separated
# by a given separator (normally "_"). It's useful for maintaining the
# population structure after clone correction. The input data is a data frame
# where the first column is a character vector of the combined population
# heirarchy. 
# Public functions utilizing this function:
# # splitcombine
#
# Internal functions utilizing this function:
# # none
#==============================================================================#

pop_splitter <- function(df, sep="_"){
  if(is.vector(df))
    df <- as.data.frame(list(comb=df), stringsAsFactors=FALSE)
  if(is.factor(df[[1]]))
    df[[1]] <- as.character(df[[1]])
  # iterating through the number of items separated by the given separator.
  for(x in seq(length(strsplit(df[[1]], sep)[[1]]))){
    # creating a column on the data frame called h# for the heirarchical level.
    # These levels are arbitrary and labeled as they are arranged in the
    # original vector. 
    df[[paste0("h",x)]] <- "NA"
    df[[paste0("h",x)]] <- vapply(strsplit(df[[1]],sep), 
                                          function(y) y[x], "1")
  }
  return(df)
}

#==============================================================================#
# This will be used to join heirarchical population vectors for the purposes of
# maintaining hierarchy. 
# Public functions utilizing this function:
# # splitcombine
#
# Internal functions utilizing this function:
# # none
#==============================================================================#

pop_combiner <- function(df, hier=c(1), sep="_"){
  if(!is.list(df)){
    warning("df must be a data frame or a list")
    return(df)
  }
  else{
    if(length(hier)==1){
      return(df[[hier]])
    }
    else{
      comb <- vector(length=length(df[[hier[1]]]))
      comb <- df[[hier[1]]]
      lapply(hier[-1], function(x) comb <<- paste(comb, df[[x]], sep=sep))
      return(comb)
    }
  }
}
#==============================================================================#
# Subsetting the population. 
# Public functions utilizing this function:
# # mlg.crosspop
#
# Internal functions utilizing this function:
# # none
#==============================================================================#
sub_index <- function(pop, sublist="ALL", blacklist=NULL){
  if (!is.genind(pop)){
    stop("pop.subset requires a genind object\n")
  }
  if (is.null(pop(pop))){
    warning("No population structure. Subsetting not taking place.")
    return(1:length(pop@ind.names))
  }
  if(toupper(sublist[1]) == "ALL"){
    if (is.null(blacklist)){
      return(1:length(pop@ind.names))
    }
    else {
      # filling the sublist with all of the population names.
      sublist <- pop@pop.names 
    }
  }

  # Checking if there are names for the population names. 
  # If there are none, it will give them names. 
  if (is.null(names(pop@pop.names))){
    if (length(pop@pop.names) == length(levels(pop@pop))){
      names(pop@pop.names) <- levels(pop@pop)
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
    } else if (is.numeric(sublist) & class(blacklist) == "character"){
    # if the sublist is numeric and blacklist is a character. eg s=1:10, b="USA"
      sublist <- sublist[sublist %in% which(!pop@pop.names %in% blacklist)]
    } else {
      # no sublist specified. Ideal situation
      if(all(pop@pop.names %in% sublist)){
        sublist <- sublist[-blacklist]
      } else {
      # weird situation where the user will specify a certain sublist, yet index
      # the blacklist numerically. Interpreted as an index of populations in the
      # whole data set as opposed to the sublist.
        warning("Blacklist is numeric. Interpreting blacklist as the index of the population in the total data set.")
        sublist <- sublist[!sublist %in% pop@pop.names[blacklist]]
      }
    }
  }

  # subsetting the population. 
  if (is.numeric(sublist))
    sublist <- names(pop@pop.names[sublist])
  else
    sublist <- names(pop@pop.names[pop@pop.names %in% sublist])
  sublist <- (1:length(pop@pop))[pop@pop %in% sublist]
  if(is.na(sublist[1])){
    warning("All items present in Sublist are also present in the Blacklist.\nSubsetting not taking place.")
    return(1:length(pop@ind.names))
  }
  #cat("Sublist:\n",sublist,"\n")
  return(sublist)
}


#==============================================================================#
# Internal function to create mlg.table.
# Public functions utilizing this function:
# # mlg.table, mlg.crosspop
#
# Internal functions utilizing this function:
# # none
#==============================================================================#
mlg.matrix <- function(x){
  if (is.genclone(x)){
    mlgvec <- x@mlg
  } else {
    mlgvec <- mlg.vector(x)
  }
  mlgs   <- length(unique(mlgvec))
  if (!is.null(pop(x))){
    mlg.mat <- table(pop(x), mlgvec)
  } else {
    mlg.mat <- matrix(table(mlgvec), nrow = 1)
    rownames(mlg.mat) <- "Total"
  }
  names(attr(mlg.mat, "dimnames")) <- NULL
  if (is.null(colnames(mlg.mat))){
    colnames(mlg.mat) <- 1:mlgs
  }
  colnames(mlg.mat) <- paste("MLG", colnames(mlg.mat), sep=".")
  return(mlg.mat)
}
#==============================================================================#
# DEPRECATED
#==============================================================================#
old.mlg.matrix <- function(x){
  mlgvec <- mlg.vector(x)
  mlgs   <- length(unique(mlgvec))
  
  if (!is.null(x@pop)){
    # creating a new population matrix. Rows are the population indicator and 
    # columns are the genotype indicator.
    mlg.mat <- matrix(ncol=mlgs, nrow=length(levels(x@pop)), data=0L)
    # populating (no, pun intended.) the matrix with genotype counts.
    lapply(levels(x@pop),function(z){
                           # This first part gets the index for the row names. 
                           count <- as.numeric(substr(z, 2, nchar(z)))
                           sapply(mlgvec[which(x@pop==z)], 
                                  function(a) mlg.mat[count, a] <<-
                                              mlg.mat[count, a] + 1L)
                         })
    rownames(mlg.mat) <-x@pop.names
  } else {
    # if there are no populations to speak of.
    mlg.mat <- t(as.matrix(vector(length=mlgs, mode="numeric")))
    sapply(mlgvec, function(a) mlg.mat[a] <<- mlg.mat[a] + 1)
    rownames(mlg.mat) <- "Total"
  }
  colnames(mlg.mat) <- paste("MLG", 1:mlgs, sep=".")
  return(mlg.mat)
}
#==============================================================================#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
# 
# The reason for this section of code is for the fact that Presence/Absence
# markers are dealt with in a different way for adegenet (to save memory) and
# so the calculations must be different as implemented in these mostly identical
# functions.
#
# Public functions utilizing this function:
# # ia
#
# Internal functions utilizing this function:
# # .ia
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
#==============================================================================#

.PA.Ia.Rd <- function(pop, missing=NULL){
	vard.vector <- NULL
	numLoci     <- ncol(pop@tab)
	numIsolates <- nrow(pop@tab)
	# Creating this number is necessary because it is how the variance is
	# calculated.
	np <- choose(numIsolates, 2)
  if(np < 2){
    return(as.numeric(c(NaN, NaN)))
  }  
  #''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''#
  # Starting the actual calculations. 
  #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,#
	V <- .PA.pairwise.differences(pop,numLoci, np, missing=missing)
	# First, set the variance of D	
	varD <- ((sum(V$D.vector^2)-((sum(V$D.vector))^2)/np))/np
	# Next is to create a vector containing all of the variances of d (there
	# will be one for each locus)
	vard.vector <- ((V$d2.vector-((V$d.vector^2)/np))/np)
	vardpair.vector <- .Call("pairwise_covar", vard.vector, PACKAGE = "poppr")
	# The sum of the variances necessary for the calculation of Ia is calculated
	sigVarj <- sum(vard.vector)
	rm(vard.vector)
	# Finally, the Index of Association and the standardized Index of associati-
	# on are calculated.
	Ia    <- (varD/sigVarj)-1
	rbarD <- (varD - sigVarj)/(2*sum(vardpair.vector))
	return(c(Ia, rbarD))
}

#==============================================================================#
# .PA.pairwise.differences will calculate three vectors that will be used for the
# calculation of the Index of Association and standardized Index of Association
# Later.
# pop = genind object 
# numLoci = should read numLoci. This will be fixed later.
# temp.d.vector = temporary vector to store the differences
# d.vector = a vector of the sum of the differences at each locus. The length
# 			 of this vector will be the same as the number of loci.
# d2.vector = the same as d.vector, except it's the sum of the squares
# D.vector = a vector of the the pairwise distances over all loci. The length
#			 of this vector will be the same as n(n-1)/2, where n is number of
# 			isolates.
# Public functions utilizing this function:
# # ia
#
# Internal functions utilizing this function:
# # .ia
#
#==============================================================================#

.PA.pairwise.differences <- function(pop,numLoci,np, missing){  
  temp.d.vector <- matrix(nrow=np, ncol=numLoci, data=as.numeric(NA))
  if( missing == "MEAN" ){
    # this will round all of the values if the missing indicator is "mean"  
    temp.d.vector <- vapply(seq(numLoci), 
                            function(x) as.vector(dist(pop@tab[,x])), 
                            temp.d.vector[,1])
    # since the replacement was with "mean", the missing data will not produce
    # a binary distance. The way we will handle this is to replace numbers that
    # are not one or zero with a rounded value. 
    tempz <- !temp.d.vector %in% 0:1
    temp.d.vector[tempz] <- vapply(temp.d.vector[tempz], round.poppr, 1)
  } else {    
    temp.d.vector <- vapply(seq(numLoci), 
                          function(x) as.vector(dist(pop@tab[,x])), 
                          temp.d.vector[,1])
    # checking for missing data and imputing the comparison to zero.
    if(any(is.na(temp.d.vector))){
      temp.d.vector[which(is.na(temp.d.vector))] <- 0
    }
  }
  if (ploidy(pop) > 1){
    # multiplying by two is the proper way to evaluate P/A diploid data because
    # one cannot detect heterozygous loci (eg, a difference of 1).
    temp.d.vector <- temp.d.vector*ploidy(pop)
    d.vector  <- as.vector(colSums(temp.d.vector))
    d2.vector <- as.vector(colSums(temp.d.vector^2))
    D.vector  <- as.vector(rowSums(temp.d.vector))
  } else {
    d.vector  <- as.vector(colSums(temp.d.vector))
    d2.vector <- d.vector
    D.vector  <- as.vector(rowSums(temp.d.vector))
  }
  vectors <- list(d.vector=d.vector, d2.vector=d2.vector, D.vector=D.vector)
  return(vectors)
}

#==============================================================================#
# Function for parsing output of poppr function.
#
# Public functions utilizing this function:
# # poppr, poppr.all
#
# Internal functions utilizing this function:
# # none
#
#==============================================================================#

final <- function(Iout, result){
  if (is.null(result)){
    return(Iout)
  } else {
    return(result)
  }
}

#==============================================================================#
# The internal version of ia. 
# Public functions utilizing this function:
# # ia, poppr
#
# Internal functions utilizing this function:
# # none
# 
#==============================================================================#

.ia <- function(pop, sample=0, method=1, quiet=FALSE, namelist=NULL, 
                missing="ignore", hist=TRUE){
  METHODS = c("permute alleles", "parametric bootstrap",
              "non-parametric bootstrap", "multilocus")
  if(pop@type!="PA"){
    type <- pop@type
    popx <- seploc(pop)
  }
  else {
    type   <- pop@type
    popx   <- pop
    .Ia.Rd <- .PA.Ia.Rd
  }
  # if there are less than three individuals in the population, the calculation
  # does not proceed. 
  if (nInd(pop) < 3){
    IarD <- as.numeric(c(NA,NA))
    names(IarD) <- c("Ia", "rbarD")
    if(sample==0){
      return(IarD)
    }
    else{
      IarD <- as.numeric(rep(NA,4))
      names(IarD) <- c("Ia","p.Ia","rbarD","p.rD")
      return(IarD)
    }
  }
  IarD <- .Ia.Rd(popx, missing)
  # data vomit options.
  if (!quiet){
    cat("|", namelist$population ,"\n")
  }
  names(IarD) <- c("Ia", "rbarD")
  # no sampling, it will simply return two named numbers.
  if (sample==0){
    Iout   <- IarD
    result <- NULL
  }
  # sampling will perform the iterations and then return a data frame indicating
  # the population, index, observed value, and p-value. It will also produce a 
  # histogram.
  else{
    Iout     <- NULL 
    idx      <- as.data.frame(list(Index=names(IarD)))
    samp     <- .sampling(popx, sample, missing, quiet=quiet, type=type, method=method)
    p.val    <- sum(IarD[1] <= c(samp$Ia, IarD[1]))/(sample + 1)#ia.pval(index="Ia", samp2, IarD[1])
    p.val[2] <- sum(IarD[2] <= c(samp$rbarD, IarD[2]))/(sample + 1)#ia.pval(index="rbarD", samp2, IarD[2])
    if(hist == TRUE){
      poppr.plot(samp, observed=IarD, pop=namelist$population,
                        file=namelist$File, pval=p.val, N=nrow(pop@tab))
    }
    result <- 1:4
    result[c(1,3)] <- IarD
    result[c(2,4)] <- p.val
    names(result)  <- c("Ia","p.Ia","rbarD","p.rD")
  } 
  return(final(Iout, result))
}
#==============================================================================#
#==============================================================================#
#=====================Index of Association Calculations========================#
#==============================================================================#
#==============================================================================#
# .pairwise.differences will calculate three vectors that will be used for the
# calculation of the Index of Association and standardized Index of Association
# Later. Note that all NA's must be treated or removed before this step.
# pop = genind object 
# numLoci = should read numLoci. This will be fixed later.
# temp.d.vector = temporary vector to store the differences
# d.vector = a vector of the sum of the differences at each locus. The length
# 			 of this vector will be the same as the number of loci.
# d2.vector = the same as d.vector, except it's the sum of the squares
# D.vector = a vector of the the pairwise distances over all loci. The length
#			 of this vector will be the same as n(n-1)/2, where n is number of
# 			isolates.
#
#
# DEPRECATED
#==============================================================================#
.pairwise.differences <- function(pop,numLoci,np, missing){  
  temp.d.vector <- matrix(nrow=np, ncol=numLoci, data=as.numeric(NA))
  if( missing == "MEAN" )
    temp.d.vector <- matrix(nrow=np, ncol=numLoci,
                            data=vapply(vapply(pop, pairwisematrix, 
                                        temp.d.vector[,1], np),round.poppr,1))  
  else    
    temp.d.vector <- vapply(pop, pairwisematrix, temp.d.vector[,1], np)
  d.vector  <- as.vector(colSums(temp.d.vector))
  d2.vector <- as.vector(colSums(temp.d.vector^2))
  D.vector  <- as.vector(rowSums(temp.d.vector))
  vectors   <- list(d.vector=d.vector, d2.vector=d2.vector, D.vector=D.vector)
  return(vectors)
}
#==============================================================================#
# pairwisematrix performs a pairwise comparison over all individuals per locus
# and returns a vector that will make its way into the final matrix for d.
# the conditional for this is that each locus must not be completely fixed for
# one allele. In that case, the resulting pairwise differences will all be zero.
#
#
# DEPRECATED 
#==============================================================================#
pairwisematrix <- function(pop, np){
  temp.d.vector <- vector(mode="numeric", length=np)
  if ( ncol(pop@tab) != 1 )
    temp.d.vector <- as.numeric(colSums(.pairwise.diffs(t(pop@tab)), na.rm=TRUE))
  return(temp.d.vector)
}
#==============================================================================#
# The original function pairwise.diffs can be found here
# https://stat.ethz.ch/pipermail/r-help/2004-August/055324.html
#
#
# DEPRECATED
#==============================================================================#
.pairwise.diffs <- function(x){
  stopifnot(is.matrix(x))

  # create column combination pairs
  prs <- cbind(rep(1:ncol(x), each = ncol(x)), 1:ncol(x))
  col.diffs <- prs[prs[, 1] < prs[, 2], , drop = FALSE]

  # do pairwise differences 
  result <- abs(x[, col.diffs[, 1]] - x[, col.diffs[, 2], drop = FALSE])

  return(result)
}
#==============================================================================#
# To calculate rbarD, the pairwise variances for each locus needs to be
# caluclated. 
#
#
# DEPRECATED
#==============================================================================#
.pairwise.variances <- function(vard.vector, pair.alleles){  
  # Here the roots of the products of the variances are being produced and
  # the sum of those values is taken. 
  vardpair.vector <- vector(length=pair.alleles)
  vardpair.vector <- sqrt(combn(vard.vector, 2, prod))
  return(vardpair.vector)
}
#==============================================================================#
# The actual calculation of Ia and rbarD. This allows for multiple populations
# to be calculated.
# pop: A list of genind objects consisting of one locus each over a population.
# Public functions utilizing this function:
# # none
#
# Internal functions utilizing this function:
# # .ia
#
#==============================================================================#
.Ia.Rd <- function (pop, missing = NULL) 
{
  vard.vector <- NULL
  numLoci     <- length(pop)
  numIsolates <- length(pop[[1]]@ind.names)
  np          <- choose(numIsolates, 2)
  if (np < 2) {
    return(as.numeric(c(NaN, NaN)))
  }
  V               <- pair_diffs(pop, numLoci, np)
  varD            <- ((sum(V$D.vector^2) - ((sum(V$D.vector))^2)/np))/np
  vard.vector     <- ((V$d2.vector - ((V$d.vector^2)/np))/np)
  vardpair.vector <- .Call("pairwise_covar", vard.vector, PACKAGE = "poppr")
  sigVarj         <- sum(vard.vector)
  rm(vard.vector)
  Ia              <- (varD/sigVarj) - 1
  rbarD           <- (varD - sigVarj)/(2 * sum(vardpair.vector))
  return(c(Ia, rbarD))
}

#==============================================================================#
# This creates a pairwise difference matrix via the C function pairdiffs in
# src/poppr_distance.c
# 
# Public functions utilizing this function:
# # none
#
# Internal functions utilizing this function:
# # .Ia.Rd
#
#==============================================================================#
pair_diffs <- function(pop, numLoci, np)
{
  ploid <- ploidy(pop[[1]])
  temp.d.vector <- matrix(nrow = np, ncol = numLoci, data = as.numeric(NA))
  temp.d.vector <- vapply(pop, function(x) .Call("pairdiffs", x@tab, PACKAGE = "poppr")*(ploid/2), 
                          temp.d.vector[, 1])
  d.vector  <- colSums(temp.d.vector)
  d2.vector <- colSums(temp.d.vector^2)
  D.vector  <- rowSums(temp.d.vector)
  return(list(d.vector = d.vector, d2.vector = d2.vector, D.vector = D.vector))
}

#==============================================================================#
# Internal counter...probably DEPRECATED.
#==============================================================================#
.new_counter <- function() {
  i <- 0
  function() {
    i <<- i + 1
    i
  }
}

#==============================================================================#
# Bruvo's distance calculation that takes in an SSR matrix. Note the conditions
# below.
#
# Public functions utilizing this function:
# # bruvo.boot
#
# Internal functions utilizing this function:
# # none
#
# DEPRECATED
#==============================================================================#
phylo.bruvo.dist <- function(ssr.matrix, replen=c(2), ploid=2, add = TRUE, loss = TRUE){
  # Preceeding functions should take care of this:
  # ssr.matrix <- genind2df(pop, sep="/", usepop=FALSE)
  # ssr.matrix[is.na(ssr.matrix)] <- paste(rep(0, ploid), collapse="/")
  # Bruvo's distance needs a matrix with the number of columns equal to the
  # number of loci multiplied by the polidy. 
  indnames <- rownames(ssr.matrix)
  ssr.matrix <- apply(ssr.matrix, 1, strsplit, "/")
  # Getting the values into numeric form.
  ssr.matrix <- apply(as.matrix(t(sapply(ssr.matrix, unlist))), 2, as.numeric)
  # Dividing each column by the repeat length and changing the values to integers.
  ssr.matrix <- apply(ssr.matrix / rep(replen, each=ploid*nrow(ssr.matrix)), 2, round)
  ssr.matrix <- apply(ssr.matrix, 2, as.integer)
  perms <- .Call("permuto", ploid, PACKAGE = "poppr")
  distmat <- .Call("bruvo_distance", ssr.matrix, perms, ploid, add, loss, PACKAGE = "poppr")
  distmat[distmat == 100] <- NA
  avg.dist.vec <- apply(distmat, 1, mean, na.rm=TRUE)
  # presenting the information in a lower triangle distance matrix.
  dist.mat <- matrix(ncol=nrow(ssr.matrix), nrow=nrow(ssr.matrix))
  dist.mat[which(lower.tri(dist.mat)==TRUE)] <- avg.dist.vec
  dist.mat <- as.dist(dist.mat)
  attr(dist.mat, "labels") <- indnames
  return(dist.mat)
}

#==============================================================================#
# This will transform the data to be in the range of [0, 1]
#
# Public functions utilizing this function:
# poppr.msn
# 
# Internal functions utilizing this function:
# # adjustcurve
#==============================================================================#

rerange <- function(x){
  minx <- min(x, na.rm = TRUE)
  maxx <- max(x, na.rm = TRUE)
  if (minx < 0)
    x <- x + abs(minx)
    maxx <- maxx + abs(minx)
  if (maxx > 1)
    x <- x/maxx
  return(x)
}

#==============================================================================#
# This will scale the edge widths
#
# Public functions utilizing this function:
# none
# 
# Internal functions utilizing this function:
# # update_edge_scales
#==============================================================================#

make_edge_width <- function(mst){
  edgewidth <- rerange(E(mst)$weight)
  if (any(edgewidth < 0.08)){
    edgewidth <- edgewidth + 0.08
  } 
  return(1/edgewidth)
}

#==============================================================================#
# This will scale the edge widths and edge color for a graph
#
# Public functions utilizing this function:
# poppr.msn bruvo.msn plot_poppr_msn
# 
# Internal functions utilizing this function:
# # singlepop_msn
#==============================================================================#
update_edge_scales <- function(mst, wscale = TRUE, gscale = TRUE, glim, gadj){
  if(gscale == TRUE){
    E(mst)$color <- gray(adjustcurve(E(mst)$weight, glim=glim, correction=gadj, 
                                     show=FALSE))
  } else {
    E(mst)$color <- rep("black", length(E(mst)$weight))
  }
  E(mst)$width <- 2
  if (wscale==TRUE){
    E(mst)$width <- make_edge_width(mst)
  }
  return(mst)
}

#==============================================================================#
# This will adjust the grey scale with respect to the edge weights for igraph.
# This is needed because the length of the edges do not correspond to weights.
# If show is set to TRUE, it will show a graph giving the equation used for con-
# version from the original scale to the grey scale, the grey scale itself in 
# the background, and the curve.
#
# Public functions utilizing this function:
# # bruvo.msn (soon)
#
# Internal functions utilizing this function:
# # new.bruvo.msn
# # new.poppr.msn
#
#==============================================================================#

adjustcurve <- function(weights, glim = c(0,0.8), correction = 3, show=FALSE, 
  scalebar = FALSE, smooth = TRUE){
  w    <- weights
  w    <- rerange(w)
  maxg <- max(glim)
  ming <- 1-(min(glim)/maxg)
  if (correction < 0){
    adj <- (w^abs(correction))/(1/ming) 
    adj <- (adj + 1-ming) / ((1 / maxg))
  } else {
    adj <- (1 - (((1-w)^abs(correction))/(1/ming)) )
    adj <- adj / (1/maxg)
  }
  if (!show){
    return(adj)
  } else if (!scalebar){
    with_quantiles <- sort(weights)
    wq_raster      <- t(as.raster(as.matrix(gray(sort(adj)), nrow = 1)))
    xlims <- c(min(weights), max(weights))
    plot(xlims, 0:1, type = "n", ylim = 0:1, xlim = xlims, xlab = "", ylab = "")
    rasterImage(wq_raster, xlims[1], 0, xlims[2], 1)
    points(x = sort(weights), y = sort(adj), col=grey(rev(sort(adj))), pch=20)
    title(xlab="Observed Value", ylab="Grey Adjusted", 
          main=paste("Grey adjustment\n min:", 
                     min(glim), 
                     "max:", max(glim), 
                     "adjust:",abs(correction)))
    if (correction < 0){
      text(bquote(frac(bgroup("(",frac(scriptstyle(x)^.(abs(correction)),
                                       .(ming)^-1),")") + .(1-ming), 
                       .(maxg)^-1)) , 
           x = min(weights) + (0.25*max(weights)), y=0.75, col="red")
    } else {
      text(bquote(frac(1-bgroup("(",frac((1-scriptstyle(x))^.(abs(correction)),
                                         .(ming)^-1),")"), 
                       .(maxg)^-1)) , 
           x= min(weights) + (0.15*max(weights)), y=0.75, col="red")
    }
    lines(x=xlims, y=c(min(glim),min(glim)), col="yellow")
    lines(x=xlims, y=c(max(glim),max(glim)), col="yellow")    
  } else {
    with_quantiles <- sort(weights)
    wq_raster      <- t(as.raster(as.matrix(gray(sort(adj)), nrow = 1)))
    no_quantiles   <- seq(min(weights), max(weights), length = 1000)
    nq_raster      <- adjustcurve(no_quantiles, glim, correction, show = FALSE)
    nq_raster      <- t(as.raster(as.matrix(gray(nq_raster), nrow = 1)))
    layout(matrix(1:2, nrow = 2))
    plot.new()
    rasterImage(wq_raster, 0, 0.5, 1, 1)
    polygon(c(0, 1, 1), c(0.5, 0.5, 0.8), col = "white", border = "white", lwd = 2)
    axis(3, at = c(0, 0.25, 0.5, 0.75, 1), labels = round(quantile(with_quantiles), 3))
    text(0.5, 0, labels = "Quantiles From Data", font = 2, cex = 1.5, adj = c(0.5, 0))
    plot.new()
    rasterImage(nq_raster, 0, 0.5, 1, 1)
    polygon(c(0, 1, 1), c(0.5, 0.5, 0.8), col = "white", border = "white", lwd = 2)
    axis(3, at = c(0, 0.25, 0.5, 0.75, 1), labels = round(quantile(no_quantiles), 3))
    text(0.5, 0, labels = "Quantiles From Smoothing", font = 2, cex = 1.5, adj = c(0.5, 0))
    # Return top level plot to defau lts.
    layout(matrix(c(1), ncol=1, byrow=T))
    par(mar=c(5,4,4,2) + 0.1) # number of lines of margin specified.
    par(oma=c(0,0,0,0)) # Figure margins
  }
}

#==============================================================================#
# This will guess the repeat lengths of the microsatellites for Bruvo's distance
# 
# Public functions utilizing this function:
# # bruvo.boot bruvo.dist
#
# Internal functions utilizing this function:
# # none
#
#==============================================================================#

guesslengths <- function(vec){
  if (length(vec) > 1){
    lens <- vapply(2:length(vec), function(x) abs(vec[x] - vec[x - 1]), 1)
    if (all(lens == 1)){
      return(1)
    } else {
      return(min(lens[lens > 1]))
    }
  } else {
    return(1)
  }
}
#==============================================================================#
# Specifically used to find loci with an excessive number of the same genotype.
#
# Public functions utilizing this function:
# # informloci
#
# Internal functions utilizing this function:
# # none
#==============================================================================#
test_table <- function(loc, min_ind, n){
  tab <- table(loc)
  return(ifelse(any(tab > n - min_ind), FALSE, TRUE))
}

#==============================================================================#
# Normalize negative branch lenght by converting the negative branch to zero
# and adding the negative value to the sibling branch.
#
# Public functions utilizing this function:
# # bruvo.boot
#
# Internal functions utilizing this function:
# # none
#==============================================================================#

fix_negative_branch <- function(tre){
  # Creating a matrix from the tree information: Tree edges and edge length
  all.lengths <- matrix(c(tre$edge, tre$edge.length), ncol = 3,
                        dimnames = list(1:length(tre$edge.length), 
                                        c("parent", "child", "length")
                        ))
  # Looking at the edges that are zero.
  zero.edges  <- all.lengths[tre$edge.length < 0, , drop = FALSE]
  # Checking which negative edges are included in all the edges
  all.edges   <- all.lengths[all.lengths[, "parent"] %in% zero.edges[, "parent"], , drop = FALSE]
  # Ordering all the edges
  index.table <- all.edges[order(all.edges[, "parent"]), , drop = FALSE]
  # Loop to change the NJ branch length
  for (i in (unique(index.table[, "parent"]))){
    the_parents <- index.table[, "parent"] == i
    fork <- index.table[, "length"][the_parents]
    index.table[, "length"][the_parents] <- abs(fork) + min(fork)      
  }
  # replacing the branch length for each negative value in the total table
  name_match <- match(rownames(index.table), rownames(all.lengths))
  all.lengths[, "length"][name_match] <- index.table[, "length"]
  # replacing the branch lengths to the original tree
  tre$edge.length <- all.lengths[, "length"]
  return(tre)
}

#==============================================================================#
# Plot Minimum Spanning Networks for data sets with no or single population
# factors.
#
# Public functions utilizing this function:
# # bruvo.msn, poppr.msn
#
# Internal functions utilizing this function:
# # none
#==============================================================================#


singlepop_msn <- function(pop, vertex.label, replen = NULL, distmat = NULL, gscale = TRUE, 
                      glim = c(0, 0.8), gadj = 3, wscale = TRUE, palette = topo.colors, showplot = TRUE, ...){
  # First, clone correct and get the number of individuals per MLG in order.
  cpop <- pop[.clonecorrector(pop), ]
  if (is.genclone(pop)){
    mlgs <- pop$mlg
    cmlg <- cpop$mlg
    mlg.number <- table(mlgs)[rank(cmlg)]
  } else {
    mlgs <- pop$other$mlg.vec
    cmlg <- cpop$other$mlg.vec
    mlg.number <- table(mlgs)[rank(cmlg)]
  }
  
  # Calculate distance matrix if not supplied (Bruvo's distance)
  if (is.null(distmat) & !is.null(replen)){
    distmat <- as.matrix(bruvo.dist(cpop, replen=replen))
  }
  
  # Create the graphs.
  g   <- graph.adjacency(distmat, weighted=TRUE, mode="undirected")
  mst <- minimum.spanning.tree(g, algorithm="prim", weights=E(g)$weight)
  
  # Create the vertex labels
  if (!is.na(vertex.label[1]) & length(vertex.label) == 1){
    if (toupper(vertex.label) == "MLG"){
      vertex.label <- paste0("MLG.", cmlg)
    } else if(toupper(vertex.label) == "INDS") {
      vertex.label <- cpop$ind.names
    }
  } 
  mst <- update_edge_scales(mst, wscale, gscale, glim, gadj)
  populations <- ifelse(is.null(pop(pop)), NA, pop$pop.names)
  
  # Plot everything
  if (showplot){
    plot.igraph(mst, edge.width = E(mst)$width, edge.color = E(mst)$color,  
                vertex.label = vertex.label, vertex.size = mlg.number*3, 
                vertex.color = palette(1),  ...)
    legend(-1.55,1,bty = "n", cex = 0.75, 
           legend = populations, title = "Populations", fill = palette(1), 
           border = NULL)
  }
  
  # Save variables and return plot.
  # E(mst)$width <- E(mst)$width
  V(mst)$size  <- mlg.number
  V(mst)$color <- palette(1)
  V(mst)$label <- vertex.label
  return(list(graph = mst, populations = populations, colors = palette(1)))
}

#==============================================================================#
# Calculate Bruvo's distance from a bruvomat object.
#
# Public functions utilizing this function:
# # bruvo.msn, bruvo.dist, bruvo.boot
#
# Internal functions utilizing this function:
# # singlepop_msn
#==============================================================================#

bruvos_distance <- function(bruvomat, funk_call = match.call(), add = TRUE, loss = TRUE){
  x      <- bruvomat@mat
  ploid  <- bruvomat@ploidy
  replen <- bruvomat@replen
  x[is.na(x)] <- 0
  # Dividing the data by the repeat length of each locus.
  x <- x / rep(replen, each=ploid*nrow(x))
  x <- matrix(as.integer(round(x)), ncol=ncol(x))
  # Getting the permutation vector.
  perms <- .Call("permuto", ploid, PACKAGE = "poppr")
  # Calculating bruvo's distance over each locus. 
  distmat <- .Call("bruvo_distance", x, perms, ploid, add, loss, PACKAGE = "poppr")
  # If there are missing values, the distance returns 100, which means that the
  # comparison is not made. These are changed to NA.
  distmat[distmat == 100] <- NA
  # Obtaining the average distance over all loci.
  avg.dist.vec <- apply(distmat, 1, mean, na.rm=TRUE)
  # presenting the information in a lower triangle distance matrix.
  dist.mat <- matrix(ncol=nrow(x), nrow=nrow(x))
  dist.mat[which(lower.tri(dist.mat)==TRUE)] <- avg.dist.vec
  dist.mat <- as.dist(dist.mat)
  attr(dist.mat, "Labels") <- bruvomat@ind.names
  attr(dist.mat, "method") <- "Bruvo"
  attr(dist.mat, "call")   <- funk_call
  return(dist.mat)
}

#==============================================================================#
# A function for creating a population hierarchy using a formula and data frame
# 
# hier = a nested formula such as ~ A/B/C where C is nested within B, which is
# nested within A.
#
# df = a data frame containing columns corresponding to the variables in hier.
#
# example:
# df <- data.frame(list(a = letters, b = LETTERS, c = 1:26))
# newdf <- make_hierarchy(~ a/b/c, df)
# df[names(newdf)] <- newdf # Add new columns.
#
# Public functions utilizing this function:
#
# # poppr.amova, setpop, gethierarchy
#
# Internal functions utilizing this function:
# # none
#==============================================================================#
make_hierarchy <- function(hier, df, expand_label = FALSE){
  newlevs <- attr(terms(hier), "term.labels")
  levs <- all.vars(hier)
  if (length(levs) > 1){
    newlevs <- gsub(":", "_", newlevs)
  }
  if (!all(levs %in% names(df))){
    stop(hier_incompatible_warning(levs, df))
  }
  newdf <- df[levs[1]]
  if (!expand_label){
    newlevs <- levs
  }
  lapply(1:length(levs), function(x) newdf[[newlevs[x]]] <<- as.factor(pop_combiner(df, levs[1:x])))
  return(newdf)
}

#==============================================================================#
# Function for creating the structure data frame needed for ade4's AMOVA
# implementation.
#
# Public functions utilizing this function:
#
# # poppr.amova
#
# Internal functions utilizing this function:
# # none
#==============================================================================#

make_ade_df <- function(hier, df, expanded = FALSE){
  if (expanded){
    levs <- attr(terms(hier), "term.labels")
  } else {
    levs <- all.vars(hier)
  }
  if(length(levs) <= 1){
    # stop("Only one level present")
    return(NULL)
  }
  levs <- gsub(":", "_", levs)
  if(!all(levs %in% names(df))){
    stop(hier_incompatible_warning(levs, df))
  }
  smallest  <- df[[levs[length(levs)]]]
  smallinds <- !duplicated(smallest)
  newdf     <- df[smallinds, , drop = FALSE]
  newdf     <- newdf[-length(levs)]
  if (length(newdf) > 1){
    factlist <- lapply(newdf, function(x) factor(x, unique(x)))
  } else {
    factlist        <- list(factor(newdf[[1]], unique(newdf[[1]])))
    names(factlist) <- names(newdf)
  }  
  return(rev(data.frame(factlist)))
}

# Function for determining if a genind object has any heterozygous sites.
# Public functions utilizing this function:
#
# # none...yet
#
# Internal functions utilizing this function:
# # none

check_Hs <- function(x){
  res <- any(x@tab > 0 & x@tab < 1, na.rm = TRUE)
  return(res)
}

## Obtains specific haplotypes from a genind object.
## 
## Arguments:
##   inds      - indexes of the first and last characters defining the haplotype
##   x         - the loci object
##   loci_cols - the columns defining the loci
##   ind_names - the sample names
##   hap_fac   - a vector defining the haplotype. 
##
## Since the inds and hap_fac args are not intuitive, I should demonstrate.
## Let's say we have a two individuals with two loci:
##      1      2
##   A: 10/20  30/40
##   B: 15/15  25/10
##
## They have the haplotypes of 
##   A1: 10 30 
##   A2: 20 40 
##   B1: 15 25
##   B2: 15 10
##
## Each genotype at a locus is 5 characters long. To get the haplotypes, we need
## to avoid the separator.
## inds in this case will be c(1, 2) for the first haplotype and c(4, 5) for the
## second haplotype and hap_fac will be c(1,1,1,2,2,2) for both.
##
# Public functions utilizing this function:
# # none
#
# Internal functions utilizing this function:
# # separate_haplotypes

hap2genind <- function(inds, x, loci_cols, ind_names, hap_fac){
  new_ind_names <- paste(ind_names, hap_fac[inds[2]], sep = ".")
  new_pop       <- rep(hap_fac[inds[2]], nrow(x))
  x2 <- lapply(x[loci_cols], substr, inds[1], inds[2])
  x2 <- data.frame(x2, stringsAsFactors = F)
  x2 <- df2genind(x2, ploidy = 1, pop=new_pop, ind.names=new_ind_names)
  return(x2)
}

## Secondary function. Transforms the genind object into a loci object and
## collects information necessary to collect individual haplotypes.
# Public functions utilizing this function:
# # none
#
# Internal functions utilizing this function:
# # pool_haplotypes
separate_haplotypes <- function(x){
  ploidy        <- ploidy(x)
  allele_list   <- unlist(lapply(x@all.names, nchar))
  allele_length <- sum(allele_list)/length(allele_list)
  if (!all(allele_list == allele_length)){
    stop("not all alleles are of equal length.")
  }
  inds        <- indNames(x)
  x.loc       <- as.loci(x)
  loci_cols   <- attr(x.loc, "locicol")
  pop_col     <- which(names(x.loc) == "population")
  geno_length <- allele_length*ploidy + ploidy - 1
  sep         <- which(1:geno_length %% (allele_length + 1) == 0)
  sep         <- c(0, sep, geno_length + 1)
  hap_fac     <- rep(1:ploidy, each = allele_length + 1)
  allele_inds <- lapply(1:ploidy, function(i) c(sep[i] + 1, sep[i + 1] - 1))
  haplist     <- lapply(allele_inds, hap2genind, x.loc, loci_cols, inds, hap_fac)
  return(haplist)
}

#==============================================================================#
# Haplotype pooling. 
# The following functions are necessary to account for within sample variation. 
# They will separate the haplotypes of a genind object and repool them so that
# there are n*k individuals in the new data set where n is the number of
# individuals and k is the ploidy. 
# Public functions utilizing this function:
# # poppr.amova
#
# Internal functions utilizing this function:
# # none
#==============================================================================#
## Main Function. Lengthens the population hierarchy as well.
pool_haplotypes <- function(x, dfname = "population_hierarchy"){
  ploidy        <- ploidy(x)
  df            <- other(x)[[dfname]]
  df$Individual <- indNames(x)
  df            <- df[rep(1:nrow(df), ploidy), , drop = FALSE]
  newx          <- repool(separate_haplotypes(x))
  the_names     <- indNames(newx)
  the_names     <- substr(the_names, 1, nchar(the_names) - 2)
  df <- df[df$Individual %in% the_names, ]
  pop(newx)     <- df$Individual
  other(newx)[[dfname]] <- df
  return(newx)
}

#==============================================================================#
# The function locus_table_pegas is the internal workhorse. It will process a 
# summary.loci object into a nice table utilizing the various diversity indices 
# provided by vegan. Note that it has a catch which will remove all allele names
# that are any amount of zeroes and nothing else. The reason for this being that
# these alleles actually represent missing data. 
# The wrapper for this is locus_table, which will take in a genind object and 
# send it into locus_table_pegas. 
# Note: lev argument has only the options of "allele" or "genotype"
# Public functions utilizing this function:
# # locus_table
#
# Internal functions utilizing this function:
# # none
#==============================================================================#
locus_table_pegas <- function(x, index = "simpson", lev = "allele", type = "codom"){
  unique_types <- x[[lev]]
  # Removing any zero-typed alleles that would be present with polyploids.
  zero_names   <- grep("^0+?$", names(unique_types))
  if (length(zero_names) > 0 & type == "codom"){
    unique_types <- unique_types[-zero_names]
  }
  
  N       <- length(unique_types)
  H       <- vegan::diversity(unique_types)
  G       <- vegan::diversity(unique_types, "inv")
  Simp    <- vegan::diversity(unique_types, "simp")
  nei     <- (N/(N-1)) * Simp
  
  if (index == "simpson"){
    idx        <- Simp
    names(idx) <- "1-D"
  } else if (index == "shannon"){
    idx        <- H
    names(idx) <- "H"
  } else {
    idx        <- G
    names(idx) <- "G"
  }
  
  E.5        <- (G - 1)/(exp(H) - 1)
  names(N)   <- lev
  return(c(N, idx, Hexp = nei, Evenness = E.5))
}

#==============================================================================#
# Function to plot phylo objects the way I want to.
#
# Public functions utilizing this function:
# bruvo.boot
#
# Private functions utilizing this function:
# # nei.boot any.boot
#==============================================================================#
poppr.plot.phylo <- function(tree, type = "nj"){
  ARGS <- c("nj", "upgma")
  type <- match.arg(type, ARGS)
  barlen <- min(median(tree$edge.length), 0.1)
  if (barlen < 0.1) barlen <- 0.01
  if (type == "nj"){
    tree <- ladderize(tree)
  }
  plot.phylo(tree, cex = 0.8, font = 2, adj = 0, xpd = TRUE, 
             label.offset = 0.0125)
  nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8, 
             font = 3, xpd = TRUE)
  if (type == "nj"){
    add.scale.bar(lwd = 5, length = barlen)
  } else {
    axisPhylo(3)
  }
}


#==============================================================================#
# Function to do something with Rodger's distance
#
# Public functions utilizing this function:
# rodger.dist
#
# Private functions utilizing this function:
# # none
#==============================================================================#
# From adegenet dist.genpop
dcano <- function(mat) {
  daux       <- mat%*%t(mat)
  vec        <- diag(daux)
  daux       <- -2*daux + vec[col(daux)] + vec[row(daux)]
  diag(daux) <- 0
  # if (any(daux == Inf)){
  #   daux <- infinite_vals_replacement(daux, warning)
  # }
  daux       <- sqrt(daux*0.5)
  return(daux)
}

#==============================================================================#
# tabulate the amount of missing data per locus. 
#
# Public functions utilizing this function:
# none
#
# Private functions utilizing this function:
# # percent_missing
#==============================================================================#

number_missing_locus <- function(x, divisor){
  missing_result <- colSums(1 - propTyped(x, by = "both"))
  return(missing_result/divisor)
}

#==============================================================================#
# tabulate the amount of missing data per genotype. 
#
# Public functions utilizing this function:
# none
#
# Private functions utilizing this function:
# # percent_missing
#==============================================================================#

number_missing_geno <- function(x, divisor){
  missing_result <- rowSums(1 - propTyped(x, by = "both"))
  return(missing_result/divisor)
}

#==============================================================================#
# Replace infinite values with the maximum finite value of a distance matrix.
#
# Public functions utilizing this function:
# nei.dist
#
# Private functions utilizing this function:
# # none
#==============================================================================#

infinite_vals_replacement <- function(D, warning){
  if (warning){
    warning("Infinite values detected.")
  }
  maxval      <- max(D[!D == Inf])
  D[D == Inf] <- maxval*10
  return(D)
}

#==============================================================================#
# Given a tree function and a distance function, this will generate an
# automatic tree generating function. This is useful for functions such as
# boot.phylo.
#
# Public functions utilizing this function:
# anyboot
#
# Private functions utilizing this function:
# # none
#==============================================================================#
tree_generator <- function(tree, distance, quiet = TRUE, ...){
  TREEFUNK <- match.fun(tree)
  DISTFUNK <- match.fun(distance)
  distargs <- formals(distance)
  otherargs <- list(...)
  #print(otherargs)
  matchargs <- names(distargs)[names(distargs) %in% names(otherargs)]
  distargs[matchargs] <- unlist(otherargs[matchargs])
  #print(distargs)
  if (!quiet) cat("\nTREE....... ", tree,"\nDISTANCE... ", distance)
  treedist <- function(x){
    distargs[[1]] <- x
    #print(distargs)
    TREEFUNK(do.call(DISTFUNK, distargs))
  }
  return(treedist)
}

#==============================================================================#
# This will retrieve a genetic matrix based on genpop status or not.
#
# Public functions utilizing this function:
# *.dist
#
# Private functions utilizing this function:
# # none
#==============================================================================#
get_gen_mat <- function(x){
  if (is.genpop(x) && x@type == "codom"){
    MAT  <- makefreq(x, missing = "mean", quiet = TRUE)$tab
  } else {
    MAT  <- x@tab
  }
  return(MAT)
}

#==============================================================================#
# This will retrieve the labels for the distance matrix from "gen" objects
#
# Public functions utilizing this function:
# *.dist
#
# Private functions utilizing this function:
# # none
#==============================================================================#
get_gen_dist_labs <- function(x){
  if (is.genind(x)){
    labs <- indNames(x)
  } else if (is.genpop(x)){
    labs <- x@pop.names
  } else if (is(x, "bootgen")){
    labs <- names(x)
  } else {
    labs <- rownames(x)
  }
  return(labs)
}

#==============================================================================#
# This will give attributes to genetic distance matrices. 
#
# Public functions utilizing this function:
# *.dist
#
# Private functions utilizing this function:
# # none
#==============================================================================#
make_attributes <- function(d, nlig, labs, method, matched_call){
  attr(d, "Size")   <- nlig
  attr(d, "Labels") <- labs
  attr(d, "Diag")   <- FALSE
  attr(d, "Upper")  <- FALSE
  attr(d, "method") <- method
  attr(d, "call")   <- matched_call
  class(d) <- "dist"
  return(d)
}

#==============================================================================#
# Function used to update colors in poppr msn 
#
# Public functions utilizing this function:
# plot_poppr_msn
#
# Private functions utilizing this function:
# # none
#==============================================================================#
update_poppr_graph <- function(graphlist, PALETTE){
  PALETTE <- match.fun(PALETTE)
  lookup  <- data.frame(old    = graphlist$populations, 
                        update = PALETTE(length(graphlist$populations)), 
                        stringsAsFactors = FALSE)
  if (nrow(lookup) > 1){
    colorlist                    <- V(graphlist$graph)$pie.color
    V(graphlist$graph)$pie.color <- lapply(colorlist, update_colors, lookup)
  } else {
    colorlist <- V(graphlist$graph)$color
    V(graphlist$graph)$color <- rep(PALETTE(1), length(colorlist))
  }
  graphlist$colors <- lookup[[2]]
  names(graphlist$colors) <- graphlist$populations
  return(graphlist)
}

#==============================================================================#
# Function used to update colors in poppr msn 
#
# Public functions utilizing this function:
# none
#
# Private functions utilizing this function:
# # update_poppr_graph
#==============================================================================#
update_colors <- function(colorvec, lookup){
  #   x <- vapply(1:length(colorvec), update_single_color, colorvec, lookup, colorvec)
  #   return(x)
  pops     <- lookup[[1]]
  update   <- lookup[[2]]
  names(update) <- pops
  matching <- match(names(colorvec), pops)
  return(update[matching[!is.na(matching)]])
}
#==============================================================================#
# Function used to update colors in poppr msn 
#
# Public functions utilizing this function:
# none
#
# Private functions utilizing this function:
# # update_colors
#==============================================================================#
update_single_color <- function(x, lookup, colorvec){
  update   <- lookup[[2]]
  original <- lookup[[1]]
  return(update[original %in% colorvec[x]])
}

#==============================================================================#
# Function used to obtain information about local ploidy in a genind data set
# for polyploids. When polyploids are imported, the entire data set is coded in
# such a way that the whole data set is considered to be the ploidy of the higest
# observed datum. The missing data are coded as "0". This function will simply
# take the difference between the maximum ploidy and the number of zeroes present. 
#
# Public functions utilizing this function:
# info_table
#
# Private functions utilizing this function:
# # none
#==============================================================================#
get_local_ploidy <- function(x){
  ploidy <- x@ploidy
  stopifnot(ploidy > 2)
  stopifnot(test_zeroes(x))
  zerocol <- which(as.numeric(x@all.names[[1]]) == 0)
  locs <- names(x@loc.names)
  locmat <- vapply(1:nLoc(x), function(z) as.integer(round(ploidy - x[loc = locs[z]]@tab[, zerocol]*ploidy)), 
                   integer(nInd(x)))
  return(locmat)
}
#==============================================================================#
# Test if a data set is comprised of microsatellites.
#
# Public functions utilizing this function:
# none
#
# Private functions utilizing this function:
# # test_zeroes
#==============================================================================#
test_microsat <- function(x){
  allnames <- unlist(lapply(x@all.names, as.numeric))
  if (any(is.na(allnames))){
    return(FALSE)
  } else {
    return(TRUE)
  }
}
#==============================================================================#
# Test if a polyploid microsatellite data set represent missing data as "0"
#
# Public functions utilizing this function:
# none
#
# Private functions utilizing this function:
# # get_local_ploidy
#==============================================================================#
test_zeroes <- function(x){
  if (test_microsat(x)){
    allnames <- unlist(lapply(x@all.names, as.numeric))
    if (any(allnames == 0) & x@ploidy > 2){
      return(TRUE)
    }
  }
  return(FALSE)
}

#==============================================================================#
# Internal plotting function for mlg.table
#
# Public functions utilizing this function:
# none
#
# Private functions utilizing this function:
# # print_mlg_barplot
#==============================================================================#
mlg_barplot <- function(mlgt){

  # create a data frame that ggplot2 can read.
  mlgt.df <- as.data.frame(list(MLG   = colnames(mlgt), 
                                count = as.vector(mlgt)), 
                           stringsAsFactors = FALSE)

  # Organize the data frame by count in descending order.
  rearranged <- order(mlgt.df$count, decreasing = TRUE)
  mlgt.df <- mlgt.df[rearranged, , drop = FALSE]
  mlgt.df[["MLG"]] <- factor(mlgt.df[["MLG"]], unique(mlgt.df[["MLG"]]))

  # plot it
  return(ggplot(mlgt.df, aes_string(x = "MLG", y = "count")) + 
         geom_bar(aes_string(fill = "count"), position="identity", stat = "identity"))
}

#==============================================================================#
# Internal plotting function for mlg.table
#
# Public functions utilizing this function:
# mlg.table
#
# Private functions utilizing this function:
# # none
#==============================================================================#

print_mlg_barplot <- function(n, mlgtab, quiet=quiet) {
  if(!quiet) cat("|",n,"\n")

  # Gather all nonzero values
  mlgt <- mlgtab[n, mlgtab[n, ] > 0, drop=FALSE]

  # controlling for the situation where the population size is 1.
  if (sum(mlgtab[n, ]) > 1){ 
    print(mlg_barplot(mlgt) +
            theme_classic() %+replace%
            theme(axis.text.x=element_text(size=10, angle=-45, hjust=0, vjust=1)) + 
            labs(title=paste("Population:", n, "\nN =", sum(mlgtab[n, ]),
                             "MLG =", length(mlgt))))
  }
}

#==============================================================================#
# Internal function for resampling loci for genotype accumulation curve.
#
# Public functions utilizing this function:
# genotype_curve
#
# Private functions utilizing this function:
# # none
#==============================================================================#
get_sample_mlg <- function(size, samp, nloci, gen, progbar){
  if (!is.null(progbar)){
    setTxtProgressBar(progbar, size/(nloci-1))
  }
  out <- vapply(1:samp, function(x) nrow(unique(sample(gen, size))), integer(1))
  return(out)
}

