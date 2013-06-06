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
#' rot pathogen, \emph{Aphanomyces euteiches} collected from two different
#' fields in NW Oregon and W Washington, USA. 
#' @format a \code{\link{genind}} object with two popualations containing a 
#' data frame in the \code{other} slot called \code{population_hierarchy}.
#' This data frame gives indices of the populations and subpopulations for the
#' data set.
#' @references Grunwald, NJ and Hoheisel, G.-A. 2006. Hierarchical Analysis 
#' of Diversity, Selfing, and Genetic Differentiation in Populations of the 
#' Oomycete \emph{Aphanomyces euteiches}. Phytopathology 96:1134-1141
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
#' 99.9\% clonal reproduction over 10,000 generations. Populations were assigned
#' post-hoc and are simply present for the purposes of demonstrating a minimum
#' spanning network with Bruvo's distance.  
#' @format a \code{\link{genind}} object with 50 individuals, 10 loci, and four 
#' popualations. 
#' @references Bo Peng and Christopher Amos (2008) Forward-time simulations of 
#' nonrandom mating populations using simuPOP. \emph{bioinformatics}, 24 (11): 
#' 1408-1409. 
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
.file.type <- function(pop, quiet=TRUE, missing="ignore", cutoff=0.05, keep=1,
                            clonecorrect=FALSE, hier=c(1), dfname="hier"){
  if (!is.genind(pop)){
    x <- pop
    if(toupper(.readExt(x)) == "CSV"){
      try(pop <- read.genalex(x), silent=quiet)
      try(pop <- read.genalex(x, region=TRUE), silent=quiet)
      try(pop <- read.genalex(x, geo=TRUE), silent=quiet)
      try(pop <- read.genalex(x, geo=TRUE, region=TRUE), silent=quiet)
      # try(pop <- read.aflp(x), silent=quiet)
    }
    else{
      try(pop <- import2genind(x, quiet=quiet), silent=quiet)
    }
    stopifnot(is.genind(pop))
    pop@call[2] <- x
    popcall <- pop@call
    pop <- missingno(pop, type=missing, cutoff=cutoff, quiet=quiet)
    pop@call <- popcall
    if (clonecorrect == TRUE){
      poplist <- clonecorrect(pop, hier=hier, dfname=dfname, keep=keep)
      pop <- poplist
      pop@call <- popcall
      #poplist <- .pop.divide(pop)
    }
  }
  else if (is.genind(pop)) {
    x <- as.character(pop@call)[2]
    popcall <- pop@call
    pop <- missingno(pop, type=missing, cutoff=cutoff, quiet=quiet)
    if (clonecorrect == TRUE){
      poplist <- clonecorrect(pop, hier=hier, dfname=dfname, keep=keep)
      pop <- poplist
      pop@call <- popcall
      #poplist <- .pop.divide(pop)
    }
  }
  return(list(X=x, GENIND=pop))
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
  res <- -which(duplicated(x@tab[, 1:ncol(x@tab)]))
  # conditional for the case that all individuals are unique.
  if(is.na(res[1])){
    res <- which(!duplicated(x@tab[, 1:ncol(x@tab)]))
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
#==============================================================================#

geno.na <- function(pop){
  pop2 <- -unique(which(is.na(pop@tab), arr.ind=T)[,1])  
  if(is.na(pop2[1])){
    return(unique(which(!is.na(pop@tab), arr.ind=T)[,1]))
  }
  else return(pop2)
}

#==============================================================================#
# loci.na will find the loci in the population that contain na's and remove
# them.
# Public functions utilizing this function:
# # none
#
# Internal functions utilizing this function:
# # percent_missing
#==============================================================================#

loci.na <- function(pop) {
  pop2 <- -unique(which(is.na(pop@tab), arr.ind=T)[,2])  
  if(is.na(pop2[1])){
    return(unique(which(!is.na(pop@tab), arr.ind=T)[,2]))
  }
  else return(pop2)
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
  if(toupper(type) == "LOCI"){
    misslist <- loci.na(pop)
    if(all(misslist > 0)){
      return(misslist)
    }
    poplen <- nInd(pop)
    filter <- vapply(-misslist, function(x) 
      length(which(is.na(pop@tab[, x])))/poplen, 1) > cutoff
    if(is.na(filter[1])){
      filter <- 1:length(misslist)
    }
  }
  else{
    misslist <- geno.na(pop)
    if(all(misslist > 0)){
      return(misslist)
    }
    poplen <- nLoc(pop)
    filter <- vapply(-misslist, function(x)
      length(unique( pop@loc.fac[which(is.na(pop@tab[x, ]))] )) / poplen, 1) > cutoff
  }
  if(all(filter %in% FALSE)){
    filter <- 1:length(misslist)
    misslist <- 1:length(misslist)
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
#==============================================================================#

ia.pval <- function(index="index", sampled, observed){
  if(all(is.nan(sampled[[index]]))){
    return(NA)
  }
  pval <- mean(ifelse(!is.na(sampled[[index]]) & sampled[[index]] >= observed,1,0))
  return(pval)
}

#==============================================================================#
# this is simply a function to print out information to the screen depending on
# what the user decides.
# Public functions utilizing this function:
# # none
#
# Internal functions utilizing this function:
# # .ia
#==============================================================================#

.quiet <- function(quiet="minimal", IarD=NULL, pop=pop, N=NULL){
  if (quiet != TRUE){
    if(quiet == FALSE){
      if (!is.null(N)){
        cat("Now Analyzing Population: ", paste(pop,",", sep=""),"N:",N,"\n")
      }
      else{
        cat("|", pop,"\n")
      }
    }
    else if(quiet == "noisy"){
      cat("Population: ", pop,"\n")
      if(!is.null(IarD)){
        cat("Index of Association: ", IarD[1],"\n")
        cat("Standardized Index of Association (rbarD): ", IarD[2],"\n")
      }  
    }
    else{
      cat("|", pop ,"\n")
    }
  }
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
    df[[paste("h",x, sep="")]] <- "NA"
    df[[paste("h",x, sep="")]] <- vapply(strsplit(df[[1]],sep), 
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
    if(is.numeric(sublist) & is.numeric(blacklist) | class(sublist) == class(blacklist)){
      sublist <- sublist[!sublist %in% blacklist]
    }
    
    # if the sublist is numeric and blacklist is a character. eg s=1:10, b="USA"
    else if(is.numeric(sublist) & class(blacklist) == "character"){
      sublist <- sublist[sublist %in% which(!pop@pop.names %in% blacklist)]
    }
    else{

      # no sublist specified. Ideal situation
      if(all(pop@pop.names %in% sublist)){
        sublist <- sublist[-blacklist]
      }

      # weird situation where the user will specify a certain sublist, yet index
      # the blacklist numerically. Interpreted as an index of populations in the
      # whole data set as opposed to the sublist.
      else{
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

mlg.matrix <- function(pop){
  
  # getting the genotype counts
  countvec2 <- mlg.vector(pop)
  
  if(!is.null(pop@pop)){
  
    # creating a new population matrix. Rows are the population indicator and 
    # columns are the genotype indicator.
    mlg.mat <- matrix(ncol=length(unique(countvec2)),nrow=length(levels(pop@pop)),
                    data=0)
    # populating (no, pun intended.) the matrix with genotype counts.
  
    lapply(levels(pop@pop),function(z){
                           # This first part gets the index for the row names. 
                           count <- as.numeric(paste(unlist
                                    (strsplit(z,""))[2:nchar(z)],
                                                       collapse=""))
                           sapply(countvec2[which(pop@pop==z)], 
                                  function(a) mlg.mat[count, a] <<-
                                              mlg.mat[count, a] + 1)
                         })
    rownames(mlg.mat) <- pop@pop.names
  }
  else{

    # if there are no populations to speak of.

    mlg.mat <- t(as.matrix(
                  vector(length=length(unique(countvec2)), mode="numeric")))
    sapply(countvec2, function(a) mlg.mat[a] <<- mlg.mat[a] + 1)
    rownames(mlg.mat) <- "Total"
  }

  colnames(mlg.mat) <- paste("MLG",seq(ncol(mlg.mat)), sep=".")
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
	numLoci <- ncol(pop@tab)
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
	vardpair.vector <- .Call("pairwise_covar", vard.vector)
	# The sum of the variances necessary for the calculation of Ia is calculated
	sigVarj <- sum(vard.vector)
	rm(vard.vector)
	# Finally, the Index of Association and the standardized Index of associati-
	# on are calculated.
	Ia <- (varD/sigVarj)-1
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
  }
  else{    
    temp.d.vector <- vapply(seq(numLoci), 
                          function(x) as.vector(dist(pop@tab[,x])), 
                          temp.d.vector[,1])
    # checking for missing data and imputing the comparison to zero.
    if(any(is.na(temp.d.vector))){
      temp.d.vector[which(is.na(temp.d.vector))] <- 0
    }
  }
  if(ploidy(pop) > 1){
    # multiplying by two is the proper way to evaluate P/A diploid data because
    # one cannot detect heterozygous loci (eg, a difference of 1).
    temp.d.vector <- temp.d.vector*ploidy(pop)
    d.vector <- as.vector(colSums(temp.d.vector))
    d2.vector <- as.vector(colSums(temp.d.vector^2))
    D.vector <- as.vector(rowSums(temp.d.vector))
  }
  else{
    d.vector <- as.vector(colSums(temp.d.vector))
    d2.vector <- d.vector
    D.vector <- as.vector(rowSums(temp.d.vector))
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
  }
  else{
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

.ia <- function(pop,sample=0,method=1,quiet="minimal",namelist=NULL,missing="ignore",
                    hist=TRUE){
  METHODS = c("multilocus", "permute alleles", "parametric bootstrap",
      "non-parametric bootstrap")
  if(pop@type!="PA"){
    type <- pop@type
    popx <- seploc(pop)
  }
  else {
    type <- pop@type
    popx <- pop
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
  .quiet(quiet=quiet, IarD=IarD, pop=namelist$population)
  names(IarD) <- c("Ia", "rbarD")
  # no sampling, it will simply return two named numbers.
  if (sample==0){
    Iout <- IarD
    result <- NULL
  }
  # sampling will perform the iterations and then return a data frame indicating
  # the population, index, observed value, and p-value. It will also produce a 
  # histogram.
  else{
    Iout <- NULL 
    idx <- as.data.frame(list(Index=names(IarD)))
    samp <- .sampling(popx, sample, missing, quiet=quiet, type=type, method=method)
    samp2 <- rbind(samp, IarD)
    p.val <- ia.pval(index="Ia", samp2, IarD[1])
    p.val[2] <- ia.pval(index="rbarD", samp2, IarD[2])
    if(hist == TRUE){
      if(require(ggplot2)){
        poppr.plot(samp, observed=IarD, pop=namelist$population,
                          file=namelist$File, pval=p.val, N=nrow(pop@tab))
      }
      else{      
        permut.histogram(samp, IarD, p.val[1], pop=namelist$population, 
                        file=namelist$File)
      }
    }
    result <- 1:4
    result[c(1,3)] <- IarD
    result[c(2,4)] <- p.val
    names(result) <- c("Ia","p.Ia","rbarD","p.rD")
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
# DEPRECIATED
#==============================================================================#
.pairwise.differences <- function(pop,numLoci,np, missing){  
  temp.d.vector <- matrix(nrow=np, ncol=numLoci, data=as.numeric(NA))
  if( missing == "MEAN" )
    temp.d.vector <- matrix(nrow=np, ncol=numLoci,
                            data=vapply(vapply(pop, pairwisematrix, 
                                        temp.d.vector[,1], np),round.poppr,1))  
  else    
    temp.d.vector <- vapply(pop, pairwisematrix, temp.d.vector[,1], np)
  d.vector <- as.vector(colSums(temp.d.vector))
  d2.vector <- as.vector(colSums(temp.d.vector^2))
  D.vector <- as.vector(rowSums(temp.d.vector))
  vectors <- list(d.vector=d.vector, d2.vector=d2.vector, D.vector=D.vector)
  return(vectors)
}
#==============================================================================#
# pairwisematrix performs a pairwise comparison over all individuals per locus
# and returns a vector that will make its way into the final matrix for d.
# the conditional for this is that each locus must not be completely fixed for
# one allele. In that case, the resulting pairwise differences will all be zero.
#
#
# DEPRECIATED 
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
# DEPRECIATED
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
# DEPRECIATED
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
  numLoci <- length(pop)
  numIsolates <- length(pop[[1]]@ind.names)
  np <- choose(numIsolates, 2)
  if (np < 2) {
    return(as.numeric(c(NaN, NaN)))
  }
  V <- pair_diffs(pop, numLoci, np)
  varD <- ((sum(V$D.vector^2) - ((sum(V$D.vector))^2)/np))/np
  vard.vector <- ((V$d2.vector - ((V$d.vector^2)/np))/np)
  vardpair.vector <- .Call("pairwise_covar", vard.vector)
  sigVarj <- sum(vard.vector)
  rm(vard.vector)
  Ia <- (varD/sigVarj) - 1
  rbarD <- (varD - sigVarj)/(2 * sum(vardpair.vector))
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
  temp.d.vector <- vapply(pop, function(x) .Call("pairdiffs", x@tab)*(ploid/2), 
                          temp.d.vector[, 1])
  d.vector <- colSums(temp.d.vector)
  d2.vector <- colSums(temp.d.vector^2)
  D.vector <- rowSums(temp.d.vector)
  return(list(d.vector = d.vector, d2.vector = d2.vector, D.vector = D.vector))
}

#==============================================================================#
# Internal counter...probably depreciated.
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
#==============================================================================#
phylo.bruvo.dist <- function(ssr.matrix, replen=c(2), ploid=2, add = TRUE, loss = TRUE){
  # Preceeding functions should take care of this:
  # ssr.matrix <- genind2df(pop, sep="/", usepop=FALSE)
  # ssr.matrix[is.na(ssr.matrix)] <- paste(rep(0, ploid), collapse="/")
  # Bruvo's distance needs a matrix with the number of columns equal to the
  # number of loci multiplied by the polidy. 

  ssr.matrix <- apply(ssr.matrix, 1, strsplit, "/")
  # Getting the values into numeric form.
  ssr.matrix <- apply(as.matrix(t(sapply(ssr.matrix, unlist))), 2, as.numeric)
  # Dividing each column by the repeat length and changing the values to integers.
  ssr.matrix <- apply(ssr.matrix / rep(replen, each=ploid*nrow(ssr.matrix)), 2, as.integer)
  perms <- .Call("permuto", ploid)
  distmat <- .Call("bruvo_distance", ssr.matrix, perms, ploid, add, loss)
  distmat[distmat == 100] <- NA
  avg.dist.vec <- apply(distmat, 1, mean, na.rm=TRUE)
  # presenting the information in a lower triangle distance matrix.
  dist.mat <- matrix(ncol=nrow(ssr.matrix), nrow=nrow(ssr.matrix))
  dist.mat[which(lower.tri(dist.mat)==TRUE)] <- avg.dist.vec
  dist.mat <- as.dist(dist.mat)
  return(dist.mat)
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

adjustcurve <- function(weights, glim = c(0,0.8), correction = 3, show=FALSE){
  w <- weights
  maxg <- max(glim)
  ming <- 1-(min(glim)/maxg)
  if (correction < 0){
    adj <-  (w^abs(correction))/(1/ming) 
    adj <- (adj + 1-ming) / ((1 / maxg))
  }
  else{
    adj <-  (1 - (((1-w)^abs(correction))/(1/ming)) )
    adj <- adj / (1/maxg)
  }
  if (show == FALSE){
    return(adj)
  }
  else{
    cols <- grey(sort(adj))
    hist(w, col=cols, border=NA, breaks=w, ylim=0:1, xlab="Observed Value", 
         ylab="Grey Adjusted", 
         main=paste("Grey adjustment\n min:", min(glim), "max:", max(glim), 
                    "adjust:",abs(correction)))
    points(x=w, y=adj, col=grey(rev(adj)), pch=20)
    if (correction < 0){
      text(bquote(frac(bgroup("(",frac(scriptstyle(x)^.(abs(correction)),
                                       .(ming)^-1),")") + .(1-ming), 
                       .(maxg)^-1)) , 
           x=0.25,y=0.75, col="red")
    }
    else{
      text(bquote(frac(1-bgroup("(",frac((1-scriptstyle(x))^.(abs(correction)),
                                         .(ming)^-1),")"), 
                       .(maxg)^-1)) , 
           x=0.15,y=0.75, col="red")
    }
    lines(x=0:1, y=c(min(glim),min(glim)), col="yellow")
    lines(x=0:1, y=c(max(glim),max(glim)), col="yellow")    
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
  if(length(vec) > 1){
    lens <- vapply(2:length(vec), function(x) abs(vec[x] - vec[x - 1]), 1)
    if(all(lens == 1)){
      return(1)
    }
    else{
      return(min(lens[lens > 1]))
    }
  }
  else{
    return(1)
  }
}
