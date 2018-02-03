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
process_file <- function(input, quiet=TRUE, missing="ignore", cutoff=0.05, 
                         keep=1, clonecorrect=FALSE, strata=1){
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
      poplist    <- clonecorrect(input, strata = strata, keep = keep)
      input      <- poplist
      input@call <- popcall
    }
  } else if (is.genind(input)) {
    x         <- as.character(match.call()[2])
    popcall   <- input@call
    input     <- missingno(input, type=missing, cutoff=cutoff, quiet=quiet)
    if (clonecorrect == TRUE){
      poplist    <- clonecorrect(input, strata = strata, keep = keep)
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
  if (is.genclone(x) | is(x, "snpclone")){
    is_duplicated <- duplicated(x@mlg[])
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
    locfac              <- locFac(pop)
    names(missing_loci) <- levels(locfac)
    missing_loci        <- missing_loci[missing_loci <= cutoff]
    misslist            <- 1:length(locfac)
    filter              <- locfac %in% names(missing_loci)
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

round.poppr <- Vectorize(function(x){
  ix <- as.integer(x)
  is_even <- ix %% 2 == 0
  if (is_even) {
    if (x - ix == 0.5)
      x <- round(x) + 1
    else if (-x + ix == 0.5)  
      x <- round(x) - 1  
  } else {
    x <- round(x)
  }
  return(x)
})
#==============================================================================#
# Subsetting the population and returning the indices.
# 
# Public functions utilizing this function:
# ## mlg.crosspop poppr.msn
#
# Internal functions utilizing this function:
# ## none
#==============================================================================#
sub_index <- function(pop, sublist="ALL", blacklist=NULL){
  numList <- seq(nInd(pop))
  if (is.null(pop(pop))){
    return(numList)
  }
  if(toupper(sublist[1]) == "ALL"){
    if (is.null(blacklist)){
      return(numList)
    } else {
      # filling the sublist with all of the population names.
      sublist <- popNames(pop) 
    }
  }
  # Treating anything present in blacklist.
  if (!is.null(blacklist)){
    # If both the sublist and blacklist are numeric or character.
    if (is.numeric(sublist) & is.numeric(blacklist) | class(sublist) == class(blacklist)){
      sublist <- sublist[!sublist %in% blacklist]
    } else if (is.numeric(sublist) & class(blacklist) == "character"){
    # if the sublist is numeric and blacklist is a character. eg s=1:10, b="USA"
      sublist <- sublist[sublist %in% which(!popNames(pop) %in% blacklist)]
    } else {
      # no sublist specified. Ideal situation
      if(all(popNames(pop) %in% sublist)){
        sublist <- sublist[-blacklist]
      } else {
      # weird situation where the user will specify a certain sublist, yet index
      # the blacklist numerically. Interpreted as an index of populations in the
      # whole data set as opposed to the sublist.
        warning("Blacklist is numeric. Interpreting blacklist as the index of the population in the total data set.")
        sublist <- sublist[!sublist %in% popNames(pop)[blacklist]]
      }
    }
  }

  # subsetting the population. 
  if (is.numeric(sublist)){
    sublist <- popNames(pop)[sublist]
  } else{
    sublist <- popNames(pop)[popNames(pop) %in% sublist]
  }
  # sublist <- (1:length(pop@pop))[pop@pop %in% sublist]
  sublist <- pop(pop) %in% sublist
  if (sum(sublist) == 0){
    warning("All items present in Sublist are also present in the Blacklist.\nSubsetting not taking place.")
    return(seq(nInd(pop)))
  } 
  #cat("Sublist:\n", sublist,"\n")
  return(numList[sublist])
}


#==============================================================================#
# Internal function to create mlg.table.
# 
# Public functions utilizing this function:
# # mlg.table mlg.crosspop
#
# Internal functions utilizing this function:
# # none
#==============================================================================#
mlg.matrix <- function(x){
  visible <- "original"
  if (is.genclone(x) | is(x, "snpclone")){
    mlgvec <- x@mlg[]
    if (is(x@mlg, "MLG")){
      visible <- visible(x@mlg)
    }
  } else {
    mlgvec <- mlg.vector(x)
  }
  
  if (!is.null(pop(x))){
    mlg.mat <- table(pop(x), mlgvec)
  } else {
    mlg.mat <- t(as.matrix(table(mlgvec)))
    rownames(mlg.mat) <- "Total"
  }
  names(attr(mlg.mat, "dimnames")) <- NULL
  if (visible == "custom"){
    return(mlg.mat)
  }
  if (is.null(colnames(mlg.mat))){
    mlgs <- length(unique(mlgvec))
    colnames(mlg.mat) <- 1:mlgs
  }
  colnames(mlg.mat) <- paste("MLG", colnames(mlg.mat), sep=".")
  return(unclass(mlg.mat))
}
# #==============================================================================#
# # DEPRECATED
# #==============================================================================#
# old.mlg.matrix <- function(x){
#   mlgvec <- mlg.vector(x)
#   mlgs   <- length(unique(mlgvec))
#   
#   if (!is.null(x@pop)){
#     # creating a new population matrix. Rows are the population indicator and 
#     # columns are the genotype indicator.
#     mlg.mat <- matrix(ncol=mlgs, nrow=length(levels(x@pop)), data=0L)
#     # populating (no, pun intended.) the matrix with genotype counts.
#     lapply(levels(x@pop), function(z){
#                            # This first part gets the index for the row names. 
#                            count <- as.numeric(substr(z, 2, nchar(z)))
#                            sapply(mlgvec[which(x@pop==z)], 
#                                   function(a) mlg.mat[count, a] <<-
#                                               mlg.mat[count, a] + 1L)
#                          })
#     rownames(mlg.mat) <-popNames(x)
#   } else {
#     # if there are no populations to speak of.
#     mlg.mat <- t(as.matrix(vector(length=mlgs, mode="numeric")))
#     sapply(mlgvec, function(a) mlg.mat[a] <<- mlg.mat[a] + 1)
#     rownames(mlg.mat) <- "Total"
#   }
#   colnames(mlg.mat) <- paste("MLG", 1:mlgs, sep=".")
#   return(mlg.mat)
# }
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
	V <- .PA.pairwise.differences(pop, numLoci, np, missing=missing)
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

.PA.pairwise.differences <- function(pop, numLoci, np, missing){  
  temp.d.vector <- matrix(nrow = np, ncol = numLoci, data = NA_real_)
  if( missing == "MEAN" ){
    # this will round all of the values if the missing indicator is "mean"  
    temp.d.vector <- vapply(seq(numLoci), 
                            function(x) as.vector(dist(pop@tab[, x])), 
                            temp.d.vector[, 1])
    # since the replacement was with "mean", the missing data will not produce
    # a binary distance. The way we will handle this is to replace numbers that
    # are not one or zero with a rounded value. 
    tempz <- !temp.d.vector %in% 0:1
    temp.d.vector[tempz] <- vapply(temp.d.vector[tempz], round.poppr, 1)
  } else {    
    temp.d.vector <- vapply(seq(numLoci), 
                          function(x) as.vector(dist(pop@tab[, x])), 
                          temp.d.vector[, 1])
    # checking for missing data and imputing the comparison to zero.
    if(any(is.na(temp.d.vector))){
      temp.d.vector[which(is.na(temp.d.vector))] <- 0
    }
  }
  if (max(ploidy(pop)) > 1){
    # multiplying by two is the proper way to evaluate P/A diploid data because
    # one cannot detect heterozygous loci (eg, a difference of 1).
    temp.d.vector <- temp.d.vector * max(ploidy(pop))
    d.vector  <- as.vector(colSums(temp.d.vector))
    d2.vector <- as.vector(colSums(temp.d.vector^2))
    D.vector  <- as.vector(rowSums(temp.d.vector))
  } else {
    d.vector  <- as.vector(colSums(temp.d.vector))
    d2.vector <- d.vector
    D.vector  <- as.vector(rowSums(temp.d.vector))
  }
  vectors <- list(d.vector = d.vector, d2.vector = d2.vector, D.vector = D.vector)
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
                missing="ignore", hist=TRUE, index = "rbarD"){
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
    IarD <- as.numeric(c(NA, NA))
    names(IarD) <- c("Ia", "rbarD")
    if(sample==0){
      return(IarD)
    }
    else{
      IarD <- as.numeric(rep(NA, 4))
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
  if (sample == 0){
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
      print(poppr.plot(samp, observed=IarD, pop=namelist$population, index = index,
                       file=namelist$File, pval=p.val, N=nrow(pop@tab)))
    }
    result <- 1:4
    result[c(1, 3)] <- IarD
    result[c(2, 4)] <- p.val
    names(result)  <- c("Ia","p.Ia","rbarD","p.rD")
    iaobj <- list(index = final(Iout, result), samples = samp)
    class(iaobj) <- "ialist"
    return(iaobj)
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
# # DEPRECATED
# #==============================================================================#
# .pairwise.differences <- function(pop, numLoci, np, missing){  
#   temp.d.vector <- matrix(nrow=np, ncol=numLoci, data=as.numeric(NA))
#   if( missing == "MEAN" )
#     temp.d.vector <- matrix(nrow=np, ncol=numLoci,
#                             data=vapply(vapply(pop, pairwisematrix, 
#                                         temp.d.vector[, 1], np), round.poppr, 1))  
#   else    
#     temp.d.vector <- vapply(pop, pairwisematrix, temp.d.vector[, 1], np)
#   d.vector  <- as.vector(colSums(temp.d.vector))
#   d2.vector <- as.vector(colSums(temp.d.vector^2))
#   D.vector  <- as.vector(rowSums(temp.d.vector))
#   vectors   <- list(d.vector=d.vector, d2.vector=d2.vector, D.vector=D.vector)
#   return(vectors)
# }
#==============================================================================#
# pairwisematrix performs a pairwise comparison over all individuals per locus
# and returns a vector that will make its way into the final matrix for d.
# the conditional for this is that each locus must not be completely fixed for
# one allele. In that case, the resulting pairwise differences will all be zero.
#
#
# DEPRECATED 
#==============================================================================#
# pairwisematrix <- function(pop, np){
#   temp.d.vector <- vector(mode="numeric", length=np)
#   if ( ncol(pop@tab) != 1 )
#     temp.d.vector <- as.numeric(colSums(.pairwise.diffs(t(pop@tab)), na.rm=TRUE))
#   return(temp.d.vector)
# }
#==============================================================================#
# The original function pairwise.diffs can be found here
# https://stat.ethz.ch/pipermail/r-help/2004-August/055324.html
#
#
# DEPRECATED
#==============================================================================#
# .pairwise.diffs <- function(x){
#   stopifnot(is.matrix(x))
# 
#   # create column combination pairs
#   prs <- cbind(rep(1:ncol(x), each = ncol(x)), 1:ncol(x))
#   col.diffs <- prs[prs[, 1] < prs[, 2], , drop = FALSE]
# 
#   # do pairwise differences 
#   result <- abs(x[, col.diffs[, 1]] - x[, col.diffs[, 2], drop = FALSE])
# 
#   return(result)
# }
#==============================================================================#
# To calculate rbarD, the pairwise variances for each locus needs to be
# calculated. 
#
#
# DEPRECATED
#==============================================================================#
# .pairwise.variances <- function(vard.vector, pair.alleles){  
#   # Here the roots of the products of the variances are being produced and
#   # the sum of those values is taken. 
#   vardpair.vector <- vector(length=pair.alleles)
#   vardpair.vector <- sqrt(combn(vard.vector, 2, prod))
#   return(vardpair.vector)
# }
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
  numIsolates <- nInd(pop[[1]])
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
# This creates the results from the pairwise difference matrix used by .Ia.Rd
# to calculate the index of association.
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
  temp.d.vector <- pair_matrix(pop, numLoci, np)
  d.vector  <- colSums(temp.d.vector)
  d2.vector <- colSums(temp.d.vector^2)
  D.vector  <- rowSums(temp.d.vector)
  return(list(d.vector = d.vector, d2.vector = d2.vector, D.vector = D.vector))
}

#==============================================================================#
# This creates a pairwise difference matrix via the C function pairdiffs in
# src/poppr_distance.c
# 
# Public functions utilizing this function:
# # none
#
# Internal functions utilizing this function:
# # pair_diffs, pair.ia
#
#==============================================================================#
pair_matrix <- function(pop, numLoci, np)
{
  temp.d.vector <- matrix(nrow = np, ncol = numLoci, data = as.numeric(NA))
  temp.d.vector <- vapply(pop, function(x) .Call("pairdiffs", tab(x), 
                                                 PACKAGE = "poppr")/2, 
                          FUN.VALUE = temp.d.vector[, 1])
  temp.d.vector <- ceiling(temp.d.vector)
  return(temp.d.vector)
}
#==============================================================================#
# Internal counter...probably DEPRECATED.
#==============================================================================#
# .new_counter <- function() {
#   i <- 0
#   function() {
#     i <<- i + 1
#     i
#   }
# }

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
# phylo.bruvo.dist <- function(ssr.matrix, replen=c(2), ploid=2, add = TRUE, loss = TRUE){
#   # Preceeding functions should take care of this:
#   # ssr.matrix <- genind2df(pop, sep="/", usepop=FALSE)
#   # ssr.matrix[is.na(ssr.matrix)] <- paste(rep(0, ploid), collapse="/")
#   # Bruvo's distance needs a matrix with the number of columns equal to the
#   # number of loci multiplied by the polidy. 
#   indnames <- rownames(ssr.matrix)
#   ssr.matrix <- apply(ssr.matrix, 1, strsplit, "/")
#   # Getting the values into numeric form.
#   ssr.matrix <- apply(as.matrix(t(sapply(ssr.matrix, unlist))), 2, as.numeric)
#   # Dividing each column by the repeat length and changing the values to integers.
#   ssr.matrix <- apply(ssr.matrix / rep(replen, each=ploid*nrow(ssr.matrix)), 2, round)
#   ssr.matrix <- apply(ssr.matrix, 2, as.integer)
#   perms <- .Call("permuto", ploid, PACKAGE = "poppr")
#   distmat <- .Call("bruvo_distance", ssr.matrix, perms, ploid, add, loss, PACKAGE = "poppr")
#   distmat[distmat == 100] <- NA
#   avg.dist.vec <- apply(distmat, 1, mean, na.rm=TRUE)
#   # presenting the information in a lower triangle distance matrix.
#   dist.mat <- matrix(ncol=nrow(ssr.matrix), nrow=nrow(ssr.matrix))
#   dist.mat[which(lower.tri(dist.mat)==TRUE)] <- avg.dist.vec
#   dist.mat <- as.dist(dist.mat)
#   attr(dist.mat, "labels") <- indnames
#   return(dist.mat)
# }

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
  if (!is.finite(minx) || !is.finite(maxx)){
    warning("non-finite values found for distances, returning 0.5")
    return(rep(0.5, length(x)))
  }
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

adjustcurve <- function(weights, glim = c(0, 0.8), correction = 3, show=FALSE, 
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
    graphics::plot(xlims, 0:1, type = "n", ylim = 0:1, xlim = xlims, xlab = "", ylab = "")
    graphics::rasterImage(wq_raster, xlims[1], 0, xlims[2], 1)
    graphics::points(x = sort(weights), y = sort(adj), col = grDevices::grey(rev(sort(adj))), pch=20)
    title(xlab="Observed Value", ylab="Grey Adjusted", 
          main=paste("Grey adjustment\n min:", 
                     min(glim), 
                     "max:", max(glim), 
                     "adjust:", abs(correction)))
    if (correction < 0){
      graphics::text(bquote(frac(bgroup("(", frac(scriptstyle(x)^.(abs(correction)),
                                       .(ming)^-1),")") + .(1-ming), 
                       .(maxg)^-1)) , 
           x = min(weights) + (0.25*max(weights)), y=0.75, col="red")
    } else {
      graphics::text(bquote(frac(1-bgroup("(", frac((1-scriptstyle(x))^.(abs(correction)),
                                         .(ming)^-1),")"), 
                       .(maxg)^-1)) , 
           x= min(weights) + (0.15*max(weights)), y=0.75, col="red")
    }
    graphics::lines(x=xlims, y=c(min(glim), min(glim)), col="yellow")
    graphics::lines(x=xlims, y=c(max(glim), max(glim)), col="yellow")    
  } else {
    with_quantiles <- sort(weights)
    wq_raster      <- t(grDevices::as.raster(as.matrix(grDevices::gray(sort(adj)), nrow = 1)))
    no_quantiles   <- seq(min(weights), max(weights), length = 1000)
    nq_raster      <- adjustcurve(no_quantiles, glim, correction, show = FALSE)
    nq_raster      <- t(grDevices::as.raster(as.matrix(grDevices::gray(nq_raster), nrow = 1)))
    graphics::layout(matrix(1:2, nrow = 2))
    graphics::plot.new()
    graphics::rasterImage(wq_raster, 0, 0.5, 1, 1)
    graphics::polygon(c(0, 1, 1), c(0.5, 0.5, 0.8), col = "white", border = "white", lwd = 2)
    graphics::axis(3, at = c(0, 0.25, 0.5, 0.75, 1), labels = round(quantile(with_quantiles), 3))
    graphics::text(0.5, 0, labels = "Quantiles From Data", font = 2, cex = 1.5, adj = c(0.5, 0))
    graphics::plot.new()
    graphics::rasterImage(nq_raster, 0, 0.5, 1, 1)
    graphics::polygon(c(0, 1, 1), c(0.5, 0.5, 0.8), col = "white", border = "white", lwd = 2)
    graphics::axis(3, at = c(0, 0.25, 0.5, 0.75, 1), labels = round(quantile(no_quantiles), 3))
    graphics::text(0.5, 0, labels = "Quantiles From Smoothing", font = 2, cex = 1.5, adj = c(0.5, 0))
    # Return top level plot to defau lts.
    graphics::layout(matrix(c(1), ncol=1, byrow=T))
    graphics::par(mar=c(5, 4, 4, 2) + 0.1) # number of lines of margin specified.
    graphics::par(oma=c(0, 0, 0, 0)) # Figure margins
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
# This now has a few caveats: 
#  1. If the parent branch contains a polytomy, then only the branch that is
#  equal to the negative branch will be fixed.
#  2. If there are branches that cannot be fixed (i.e., the absolute value of 
#  the negative branch is greater than the value of the sibling), then it will 
#  not be fixed.
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
    # Check for polytomies
    if (length(fork) > 2){
      # Fix only siblings of equal magnitude.
      forkinds <- abs(fork) == -min(fork)
      fork     <- fork[forkinds]
      fixed_lengths <- abs(fork) + min(fork)
      # Check that branches are actually fixed.
      if (all(fixed_lengths >= 0)){
        index.table[, "length"][the_parents][forkinds] <- abs(fork) + min(fork)
      }
    } else { # No polytomies
      fixed_lengths <- abs(fork) + min(fork)
      # Check that branches are actually fixed.
      if (all(fixed_lengths >= 0)){
        index.table[, "length"][the_parents] <- abs(fork) + min(fork)
      }
    }
  }
  # replacing the branch length for each negative value in the total table
  name_match <- match(rownames(index.table), rownames(all.lengths))
  all.lengths[, "length"][name_match] <- index.table[, "length"]
  # replacing the branch lengths to the original tree
  tre$edge.length <- all.lengths[, "length"]
  return(tre)
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

bruvos_distance <- function(bruvomat, funk_call = match.call(), add = TRUE, 
                            loss = TRUE, by_locus = FALSE){
  
  
  x      <- bruvomat@mat
  ploid  <- bruvomat@ploidy
  if (getOption("old.bruvo.model") && ploid > 2 && (add | loss)){
    msg <- paste("The option old.bruvo.model has been set to TRUE, which does",
                 "not represent every ordered combinations of alleles in the",
                 "genome addition or loss models. This could result in",
                 "potentially incorrect results.",
                 "\n\n To use every ordered combination of alleles for",
                 "estimating short genotypes, enter the following command in",
                 "your R console:",
                 "\n\n\toptions(old.bruvo.model = FALSE)\n")
    warning(msg, call. = FALSE, immediate. = TRUE)
  }
  replen <- bruvomat@replen
  x[is.na(x)] <- 0

  # Dividing the data by the repeat length of each locus.
  x <- x / rep(replen, each = ploid * nrow(x))
  x <- matrix(as.integer(round(x)), ncol = ncol(x))

  # Getting the permutation vector.
  perms <- .Call("permuto", ploid, PACKAGE = "poppr")

  # Calculating bruvo's distance over each locus. 
  distmat <- .Call("bruvo_distance", 
                   x,     # data matrix
                   perms, # permutation vector (0-indexed)
                   ploid, # maximum ploidy
                   add,   # Genome addition model switch
                   loss,  # Genome loss model switch
                   getOption("old.bruvo.model"), # switch to use unordered genotypes
                   PACKAGE = "poppr")

  # If there are missing values, the distance returns 100, which means that the
  # comparison is not made. These are changed to NA.
  distmat[distmat == 100] <- NA

  if (!by_locus){
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
  } else {
    n    <- nrow(x)
    cols <- seq(ncol(distmat))
    labs <- bruvomat@ind.names
    meth <- "Bruvo"
    return(lapply(cols, function(i) make_attributes(distmat[, i], n, labs, meth, funk_call)))
  }

}


#==============================================================================#
# match repeat lengths to loci present in data
#
# Public functions utilizing this function:
# ## initialize,bruvomat-method
#
# Internal functions utilizing this function:
# ## none
#==============================================================================#
match_replen_to_loci <- function(gid_loci, replen){
  # unnamed loci with the same length
  if (is.null(names(replen))){
    if (is.character(gid_loci)){
      names(replen) <- gid_loci # give them names
    }
    return(replen)
  } else if (!all(gid_loci %in% names(replen))){ # names don't match up
    unmatched_repeats <- names(replen[!names(replen) %in% gid_loci])
    stop(unmatched_loci_warning(unmatched_repeats, gid_loci), call. = FALSE)
  } else {
    return(replen[gid_loci])
  }
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
# ## none
#
# Internal functions utilizing this function:
# ## none
# 
# DEPRECATED
#==============================================================================#
# make_hierarchy <- function(hier, df, expand_label = FALSE){
#   newlevs <- attr(terms(hier), "term.labels")
#   levs <- all.vars(hier)
#   if (length(levs) > 1){
#     newlevs <- gsub(":", "_", newlevs)
#   }
#   if (!all(levs %in% names(df))){
#     stop(hier_incompatible_warning(levs, df))
#   }
#   newdf <- df[levs[1]]
#   if (!expand_label){
#     newlevs <- levs
#   }
#   lapply(1:length(levs), function(x) newdf[[newlevs[x]]] <<- as.factor(pop_combiner(df, levs[1:x])))
#   return(newdf)
# }

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
# # poppr.amova
#
# Internal functions utilizing this function:
# # none

check_Hs <- function(x){
  if (all(ploidy(x) == 1)) return(FALSE)
  res <- sweep(tab(x), 1, ploidy(x), function(x, y) any(x < y))
  return(res)
}

#==============================================================================#
## Obtains specific haplotypes from a genind object.
#  
# Arguments:
#    index - name or index of locus
#    locus - a data frame of genotypes
#    sep_location - a matrix with the location of each separator
#    geno_lengths - a matrix with the char lengths of each genotype
#    first - when TRUE, the first haplotype is extracted, when FALSE, the second
#
# Returns:
#   a character vector with a single locus haplotype for each sample. 
#   
# Public functions utilizing this function:
# # none
#
# Internal functions utilizing this function:
# # separate_haplotypes
#==============================================================================#
# get_haplotype <- function(index, locus, sep_location, geno_lengths, first = TRUE){
#   if (first){
#     res <- substr(locus[[index]], 
#                   start = 1, 
#                   stop = sep_location[, index] - 1)
#   } else {
#     res <- substr(locus[[index]], 
#                   start = sep_location[, index] + 1,
#                   stop = geno_lengths[, index])
#   }
#   return(res)
# }

#==============================================================================#
# Arguments:
#   x a diploid genind object. 
# 
# Returns:
#   list:
#     loci_data    - a data frame of genotypes
#     sep_location - a matrix with the location of each separator
#     geno_lengths - a matrix with the char lengths of each genotype
#
# External functions utilizing this function:
# # none
#
# Internal functions utilizing this function:
# # separate_haplotypes
#==============================================================================#
# haplotype_detector <- function(x){
#   x.loc        <- genind2df(x, sep = "/", usepop = FALSE)
#   facstr       <- function(i) nchar(as.character(i))
#   sep_location <- lapply(x.loc, function(i) regexpr("/", i))
#   sep_location <- matrix(unlist(sep_location, use.names = FALSE), 
#                          nrow = nrow(x.loc))
#   dimnames(sep_location) <- dimnames(x.loc)
#   geno_lengths <- vapply(x.loc, facstr, integer(nrow(x.loc)))
#   dimnames(geno_lengths) <- dimnames(sep_location)
#   haplist <- list(loci_data = x.loc, 
#                   sep_location = sep_location, 
#                   geno_lengths = geno_lengths)
#   return(haplist)
# }

#==============================================================================#
# Arguments:
#   x a diploid genind object. 
# 
# Returns:
#   a matrix with n*2 rows where haplotypes are interleaved. This is to be 
#   converted into a genind object. To explain what it looks like:
#   
#   if you have your genotypes arranged like so:
#      L1   L2
#   A  x/y  a/b
#   B  z/z  c/a
#   
#   The result of this function would be to return a matrix that looks like:
#        L1 L2
#   A.1  x  a
#   A.2  y  b
#   B.1  z  c
#   B.2  z  a
#     
# External functions utilizing this function:
# # none
#
# Internal functions utilizing this function:
# # separate_haplotypes
#==============================================================================#
separate_haplotypes <- function(x){
  # haps      <- haplotype_detector(x)
  # the_loci  <- colnames(haps$sep_location)
  # ind_names <- rownames(haps$sep_location)
  # h1 <- lapply(the_loci, get_haplotype, haps$loci_data, haps$sep_location, 
  #              haps$geno_lengths, first = TRUE)
  # h2 <- lapply(the_loci, get_haplotype, haps$loci_data, haps$sep_location, 
  #              haps$geno_lengths, first = FALSE)
  # newnames <- paste(rep(ind_names, each = 2), 1:2, sep = ".")
  # outlength <- length(newnames)
  # outvec <- seq(outlength)
  # index1 <- outvec %% 2 == 1
  # index2 <- !index1
  # matdimnames <- list(newnames, the_loci)
  # hapmat <- matrix(NA_character_, nrow = outlength,
  #                  ncol = length(the_loci), dimnames = matdimnames)
  # hapmat[index1, ] <- unlist(h1)
  # hapmat[index2, ] <- unlist(h2)
  # return(hapmat)
  if (max(ploidy(x)) > 2){
    x <- recode_polyploids(x, addzero = TRUE)
  }
  df <- genind2df(x, sep = "/", usepop = FALSE)
  df <- apply(df, 1, strsplit, "/")
  df <- lapply(df, lapply, function(i){ i[i == "0"] <- NA; i }) # replace zeroes as missing
  df <- lapply(df, data.frame, stringsAsFactors = FALSE, check.names = FALSE)
  dplyr::bind_rows(df)
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
pool_haplotypes <- function(x){
  ploidy <- max(ploidy(x))
  if (is.null(strata(x))){
    strata(x) <- data.frame(pop = pop(x))
  }
  addStrata(x)  <- data.frame(Individual = indNames(x))
  df            <- strata(x)
  df            <- df[rep(seq(nrow(df)), each = ploidy), , drop = FALSE]
  newdf         <- separate_haplotypes(x)
  if (ploidy > 2){
    is_typed <- apply(newdf, 1, function(i) sum(is.na(i))) < nLoc(x)
    newdf <- newdf[is_typed, , drop = FALSE]
    df    <- df[is_typed, , drop = FALSE]
  }
  newx          <- df2genind(newdf, ploidy = 1, strata = df)
  setPop(newx)  <- ~Individual
  return(newx)
}

#==============================================================================#
# The function locus_table_pegas is the internal workhorse. It will process a
# summary.loci object into a nice table utilizing the various diversity indices
# provided by vegan. Note that it has a catch which will remove all allele names
# that are any amount of zeroes and nothing else. The reason for this being that
# these alleles actually represent missing data. The wrapper for this is
# locus_table, which will take in a genind object and send it into
# locus_table_pegas. 
# Note: lev argument has only the options of "allele" or "genotype"
# 
# UPDATE 2015-07-16: I had noticed in issue #47 that the definition of Hexp was
# completely wrong. It was correcting by the number of classes as opposed to the
# number of samples, artificially increasing the diversity. It has been
# corrected with one exception: polyploids will not be able to utilize this
# correction because the allelic dosage is ambiguous. This is corrected by
# utilizing Müller's index, which is equivalent to the unbiased Simpson's index.
# 
# Public functions utilizing this function:
# # locus_table
#
# Internal functions utilizing this function:
# # none
#==============================================================================#
locus_table_pegas <- function(x, index = "simpson", lev = "allele", ploidy = NULL, 
                              type = "codom", hexp_only = FALSE){
  unique_types <- x[[lev]]
  # Removing any zero-typed alleles that would be present with polyploids.
  unique_types <- remove_zeroes(unique_types, type)

  n    <- sum(unique_types)
  Simp <- vegan::diversity(unique_types, "simp")
  nei  <- (n/(n-1)) * Simp
  
  if (hexp_only){
    return(nei)
  }
  
  N <- length(unique_types) # Resetting N to be the number of unique types.
  H <- vegan::diversity(unique_types)
  G <- vegan::diversity(unique_types, "inv")
  
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
  E.5      <- (G - 1)/(exp(H) - 1)
  names(N) <- lev
  return(c(N, idx, Hexp = nei, Evenness = E.5)) 

}
#==============================================================================#
# This function takes a vector of allele counts and searches for the names
# that appear to be zero-length alleles. These alleles are then removed from
# analysis.
# 
# Public functions utilizing this function:
# # none
#
# Internal functions utilizing this function:
# # locus_table_pegas
# # get_hexp_from_loci
#==============================================================================#
remove_zeroes <- function(unique_types, type){
  zero_names   <- grep("^0+?$", names(unique_types))
  if (length(zero_names) > 0 & type == "codom"){
    unique_types <- unique_types[-zero_names]
  }
  return(unique_types)
}

#==============================================================================#
# This function replaces the old definition of Hexp in poppr and will actually
# return a measure of expected heterozygosity. The only catch is that if an 
# incoming population has varying ploidy or ploidy greater than two, the
# returned value is the same as the result of Hs (except missing values do not
# corrupt the whole calculation). This is then treated in the poppr function by
# multiplying by the sample size correction to give a corrected simpson's index.
# 
# Public functions utilizing this function:
# # poppr
#
# Internal functions utilizing this function:
# # none
#==============================================================================#
get_hexp_from_loci <- function(loci, type = "codom", ploidy = NULL){
  loci <- vapply(summary(loci), FUN = locus_table_pegas, FUN.VALUE = numeric(1),
                 ploidy = ploidy, type = type, hexp_only = TRUE)    
  return(mean(loci, na.rm = TRUE))
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
poppr.plot.phylo <- function(tree, type = "nj", root = FALSE){
  barlen <- min(median(tree$edge.length), 0.1)
  if (barlen < 0.1) barlen <- 0.01
  if (!root && type != "upgma"){
    tree <- ladderize(tree)
  }
  nodelabs <- round(tree$node.label, 2)
  plot.phylo(tree, cex = 0.8, font = 2, adj = 0, xpd = TRUE, 
             label.offset = 0.0125)
  nodelabels(nodelabs, adj = c(1.3, -0.5), frame = "n", cex = 0.8, 
             font = 3, xpd = TRUE)
  if (type != "upgma"){
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
  distargs <- as.list(formals(distance))
  otherargs <- list(...)
  #print(otherargs)
  matchargs <- names(distargs)[names(distargs) %in% names(otherargs)]
  distargs[matchargs] <- otherargs[matchargs]
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
  if (suppressWarnings(is.genpop(x)) && x@type == "codom"){
    MAT  <- makefreq(x, missing = "mean", quiet = TRUE)
  } else {
    MAT  <- tab(x, freq = TRUE)
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
    labs <- popNames(x)
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
# Input: 
#  - d a vector to be coerced into a distance matrix.
#  - nlig number of samples
#  - labs the names of the samples
#  - method the name of the method used to create the distances
#  - matched_call the call that created the matrix.
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
# Parse incoming palettes
# 
# Color palettes could come in the form of functions, vectors, or named vectors.
# This internal helper will parse them and return a named vector of colors. 
#
# Public functions utilizing this function:
# ## bruvo.msn poppr.msn
#
# Internal functions utilizing this function:
# ## none
#==============================================================================#
palette_parser <- function(inPAL, npop, pnames){
  PAL <- try(match.fun(inPAL, descend = FALSE), silent = TRUE)
  if ("try-error" %in% class(PAL)){
    if (all(pnames %in% names(inPAL))){
      color <- inPAL[pnames]
    } else if (npop == length(inPAL)){
      color <- stats::setNames(inPAL, pnames)
    } else if (npop < length(inPAL)){
      warning("Number of populations fewer than number of colors supplied. Discarding extra colors.")
      color <- stats::setNames(inPAL[1:npop], pnames)
    } else {
      warning("insufficient color palette supplied. Using topo.colors().")
      color <- stats::setNames(topo.colors(npop), pnames)
    }
  } else {
    color   <- stats::setNames(PAL(npop), pnames)
  }
  return(color)
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
  updates <- setNames(PALETTE(length(graphlist$populations)), NULL)
  lookup  <- data.frame(old    = graphlist$populations, 
                        update = updates, 
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
# ## none
#
# Private functions utilizing this function:
# ## none
# 
# DEPRECATED
#==============================================================================#
# update_single_color <- function(x, lookup, colorvec){
#   update   <- lookup[[2]]
#   original <- lookup[[1]]
#   return(update[original %in% colorvec[x]])
# }

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
  if (any(as.numeric(unlist(x@all.names, use.names = FALSE)) == 0)){
    x <- recode_polyploids(x, newploidy = TRUE)    
  }
  tabx   <- tab(x)
  cols   <- split(colnames(tabx), locFac(x))
  locmat <- vapply(cols, function(i) rowSums(tabx[, i, drop = FALSE]), numeric(nInd(x)))
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
    
    allnames  <- as.numeric(unlist(alleles(x), use.names = FALSE))
    
    if (any(allnames == 0)){
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
#' @param mlgt a table with populations in rows and MLGs in columns
#'
#' @param color logical. Should the output be colored by population?
#' @param background logical. should the data be plotted against the background
#'  data?
#' @noRd
#' @keywords internal
#'
#' @importFrom dplyr %>% 
#' @importFrom dplyr arrange_ group_by_ ungroup
#' @importFrom dplyr n
#' @importFrom rlang "!!"
mlg_barplot <- function(mlgt, color = FALSE, background = FALSE){
  names(dimnames(mlgt)) <- c("Population", "MLG")
  mlgt.df <- as.data.frame.table(mlgt, responseName = "count", stringsAsFactors = FALSE)
  use_old_dplyr <- utils::packageVersion("dplyr") <= package_version("0.5.0") | getOption("poppr.old.dplyr")
  # Ensure that the population is a factor
  if (use_old_dplyr) {
    mlgt.df <- mlgt.df %>%
      dplyr::mutate_(.dots = list(Population = ~factor(Population, unique(Population)))) %>%
      dplyr::filter_("count > 0")
  } else {
    qPop <- quote(Population)
    mlgt.df <- mlgt.df %>%
      dplyr::mutate(Population = factor(!!qPop, unique(!!qPop))) %>%
      dplyr::filter(!!quote(count > 0))
  }
    
  if (color | background) {
    # summarize with the total counts, and merge with original data
    if (use_old_dplyr) {
      mlgt.df <- mlgt.df %>% 
        dplyr::group_by_("MLG") %>%
        dplyr::summarize_(.dots = list(n = ~sum(count))) %>%
        dplyr::full_join(mlgt.df, by = "MLG") %>%
        dplyr::arrange_(~dplyr::desc(n)) %>%
        dplyr::mutate_(.dots = list(MLG = ~factor(MLG, levels = unique(MLG))))
    } else {
      qMLG <- quote(MLG)
      mlgt.df <- mlgt.df %>% 
        dplyr::group_by(!!qMLG) %>%
        dplyr::summarize(n = sum(!!quote(count))) %>%
        dplyr::full_join(mlgt.df, by = "MLG") %>%
        dplyr::arrange(dplyr::desc(!!quote(n))) %>%
        dplyr::mutate(MLG = factor(!!qMLG, levels = unique(!!qMLG)))
    }

    the_breaks <- pretty(mlgt.df$n)
  } else {
    if (use_old_dplyr) {
      mlgt.df <- mlgt.df %>%
        dplyr::group_by_("Population") %>%
        dplyr::do_(~dplyr::arrange(., dplyr::desc(count)))
    } else {
      mlgt.df <- mlgt.df %>%
        dplyr::group_by(!!quote(Population)) %>%
        dplyr::arrange(dplyr::desc(!!quote(count)), .by_group = TRUE)
    }

    the_breaks <- pretty(mlgt.df$count)
  }
  if (use_old_dplyr) {
    mlgt.df <- mlgt.df %>%
      dplyr::ungroup() %>%
      unique() %>%
      dplyr::filter_("count > 0") %>%
      dplyr::mutate_(.dots = list(order = ~seq(nrow(.)))) %>%
      dplyr::mutate_(.dots = list(order = ~factor(order, unique(order))))
  } else {
    mlgt.df <- mlgt.df %>%
      dplyr::ungroup() %>%
      unique() %>%
      dplyr::filter(!!quote(count > 0)) %>%
      dplyr::mutate(order = seq(n())) %>%
      dplyr::mutate(order = factor(!!quote(order), unique(!!quote(order))))
  }

  if (color | background) {
    if (use_old_dplyr) {
      mlgt.df <- mlgt.df %>%
        dplyr::select_(~Population, ~MLG, ~count, ~order, ~n)
    } else {
      mlgt.df <- mlgt.df %>%
        dplyr::select(!!quote(Population), !!quote(MLG), !!quote(count), !!quote(order), !!quote(n))
    }

  }
  the_breaks <- the_breaks[the_breaks %% 1 == 0]
  # Conditionals for plotting
  x    <- if (color) "MLG" else "order"
  baes <- if (color) aes_string(fill = "Population") else aes_string()
  
  
  the_plot <- ggplot(mlgt.df, aes_string(x = x, y = "count")) 
  if (background){
    the_plot <- the_plot + 
      geom_bar(aes_string(y = "n"), color = "grey25", alpha = 0, stat = "identity")
  }
  the_plot <- the_plot + 
    geom_bar(baes, stat = "identity") + 
    theme(panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(size = 10, angle = 90, 
                                     hjust = 1, vjust = 1)) +
    xlab("MLG") + 
    scale_y_continuous(expand = c(0, 0), breaks = the_breaks)
  if (!color || (color && background)){
    xbreak <- if (color) mlgt.df$MLG else mlgt.df$order
    the_plot <- the_plot +
      scale_x_discrete(labels = mlgt.df$MLG, breaks = xbreak) +
      facet_wrap(~Population, scales = "free_x", shrink = TRUE, drop = TRUE)
  }
  return(the_plot)
}

#==============================================================================#
# Internal function for creating title for index of association histogram.
#
# Public functions utilizing this function:
# # none
#
# Private functions utilizing this function:
# # poppr.plot
#==============================================================================#
make_poppr_plot_title <- function(samp, file = NULL, N = NULL, pop = NULL){
  plot_title <- ifelse(is.null(pop), "", paste0("Population:", pop))
  plot_title <- ifelse(is.null(N), paste0(plot_title, ""), 
                       paste(plot_title, paste("N:", N), sep = "\n"))
  plot_title <- ifelse(is.null(file), paste0(plot_title, ""), 
                       paste(plot_title, paste("Data:", file), sep = "\n"))
  perms      <- paste("Permutations:", length(samp))
  plot_title <- paste(plot_title, perms, sep = "\n")
  return(plot_title)
}

#==============================================================================#
# fill a single genotype with zeroes if the number of alleles is maxploid.
#
# Public functions utilizing this function:
# # none
#
# Private functions utilizing this function:
# # fill_zero_locus
#==============================================================================#
fill_zero <- function(x, maxploid, mat = FALSE){
  if (length(x) < maxploid){
    if (!mat){
      zeroes <- paste(rep(0, maxploid - length(x)), collapse = "/")
      res    <- paste(x, collapse = "/")
      res    <- paste(zeroes, res, sep = "/")     
    } else {
      res <- c(rep(0.0, maxploid - length(x)), as.numeric(x))
    }
 
  } else {
    if (!mat){
      res <- paste(x, collapse = "/")
    } else {
      res <- as.numeric(x)
    }
  }
  return(res)
}

#==============================================================================#
# Fill short genotypes in a character vector with zeroes.
#
# Public functions utilizing this function:
# # none
#
# Private functions utilizing this function:
# # generate_bruvo_mat
#==============================================================================#
fill_zero_locus <- function(x, sep = "/", maxploid, mat = FALSE){
  x <- strsplit(x, sep)
  if (mat){
    result <- numeric(maxploid)
  } else {
    result <- character(1)
  }
  return(t(vapply(x, fill_zero, result, maxploid, mat)))
}

#==============================================================================#
# This will fill in zeroes from a data frame with uneven ploidy to force it to
# have uniform ploidy.
#
# Example: 
#               locus_1     locus_2
# sample_1     25/91/20 15/33/10/55
# sample_2  71/29/34/69          45
# sample_3           24    23/83/25
# sample_4     12/18/94           4
# sample_5           68    15/54/57
# sample_6     36/45/91    32/74/22
# sample_7           72 22/27/70/75
# sample_8       100/54          92
# sample_9        62/50 17/11/62/50
# sample_10       41/31    17/30/57
#
# Will become
#
#               locus_1     locus_2
# sample_1   0/25/91/20 15/33/10/55
# sample_2  71/29/34/69    0/0/0/45
# sample_3     0/0/0/24  0/23/83/25
# sample_4   0/12/18/94     0/0/0/4
# sample_5     0/0/0/68  0/15/54/57
# sample_6   0/36/45/91  0/32/74/22
# sample_7     0/0/0/72 22/27/70/75
# sample_8   0/0/100/54    0/0/0/92
# sample_9    0/0/62/50 17/11/62/50
# sample_10   0/0/41/31  0/17/30/57
#
# or, if mat = TRUE
#
#           locus_1.1 locus_1.2 locus_1.3 locus_1.4 locus_2.1 locus_2.2 locus_2.3 locus_2.4
# sample_1          0        25        91        20        15        33        10        55
# sample_2         71        29        34        69         0         0         0        45
# sample_3          0         0         0        24         0        23        83        25
# sample_4          0        12        18        94         0         0         0         4
# sample_5          0         0         0        68         0        15        54        57
# sample_6          0        36        45        91         0        32        74        22
# sample_7          0         0         0        72        22        27        70        75
# sample_8          0         0       100        54         0         0         0        92
# sample_9          0         0        62        50        17        11        62        50
# sample_10         0         0        41        31         0        17        30        57
#
#
# Public functions utilizing this function:
# # none
#
# Private functions utilizing this function:
# # none
#==============================================================================#
generate_bruvo_mat <- function(x, maxploid, sep = "/", mat = FALSE){
  if (mat){
    result <- matrix(numeric(nrow(x)*maxploid), ncol = maxploid, nrow = nrow(x))
  } else {
    result <- character(nrow(x))
  }
  res <- vapply(x, fill_zero_locus, result, sep, maxploid, mat)
  if (length(dim(res)) > 2){
    redim    <- dim(res)
    dim(res) <- c(redim[1], redim[2]*redim[3])
    xcols  <- colnames(x)
    maxseq <- seq(maxploid)
    colnames(res) <- vapply(xcols, FUN = paste, 
                            FUN.VALUE = character(maxploid), maxseq, sep = ".")
  } else {
    colnames(res) <- colnames(x)
  }
  if (!mat){
    res[grep("NA", res)] <- NA_character_
  }
  rownames(res) <- rownames(x)
  return(res)
}
#==============================================================================#
# Function to subset the custom MLGs by the computationally derived MLGs in the
# data set. This is necessary due to the fact that minimum spanning networks
# will clone correct before calculations, but this is performed on the visible
# multilocus genotypes. for custom MLGs that may not be monophyletic, this 
# results in observed networks that may be incorrect. 
# 
# A solution to this would simply be to label the multilocus genotypes with 
# their custom labels, but collapse them with the computationally derived 
# labels.
# 
# In order to parse these out, we have three possible situations we can think
# of:
# 
#  1. computational MLGs match the custom MLGs: pretty easy, simply return the
#     non-duplicated mlgs
#  2. There are more computational MLG classes than custom MLGs. This is also
#     fairly simple: return the custom MLGs censored by the computational MLGs
#  3. More custom MLGs than computational MLGs. For labelling purposes, 
#     the custom MLGs that occupy a single MLG should be concatenated in a 
#     string.
#  4. A mix of 2 and 3. Same strategy as 3.
#  
#  Input: a genclone or snpclone object with an MLG object in the @mlg slot
#  Output: A string of clone-censored custom MLGs
#
# Public functions utilizing this function:
# # bruvo.msn poppr.msn plot_poppr_msn
#
# Private functions utilizing this function
# # singlepop_msn
#==============================================================================#
correlate_custom_mlgs <- function(x, by = "original", subset = TRUE){
  if (!is.genclone(x) & !is(x, "snpclone")){
    stop("needs a genclone or snpclone object")
  }
  if (!is(x@mlg, "MLG")){
    stop("the @mlg slot needs to be of class MLG. This procedure is meaningless otherwise.")
  }
  the_mlgs <- x@mlg@mlg
  customs  <- as.character(the_mlgs[["custom"]])
  ncustom  <- nlevels(the_mlgs[["custom"]])
  mlg_dup  <- !duplicated(the_mlgs[[by]])
  ndup     <- sum(mlg_dup)
  
  # Create contingency table with custom genotypes in rows and computed in
  # columns
  cont_table <- table(customs, the_mlgs[[by]])
  
  # Create true/false table of MLG identity
  i_table <- cont_table > 0
  
  # Count up the number of custom MLGs contained within each computed MLG.
  check_less_custom <- colSums(i_table)
  
  if ((ndup == ncustom | ndup > ncustom) & all(check_less_custom == 1)){
    if (!subset){
      return(customs)
    }
    res        <- customs[mlg_dup]
    names(res) <- the_mlgs[[by]][mlg_dup]
    return(res)
  }
  
  cust_names <- rownames(cont_table)
  res <- apply(cont_table, 2, function(i) paste(cust_names[i > 0], collapse = "\n"))
  if (!subset){
    res <- res[as.character(as.character(the_mlgs[[by]]))]
  }
  return(res)
}


#==============================================================================#
# Function to boostrap a single population from a mlg.matrix
# 
# Public functions utilizing this function:
# ## diversity_boot
# 
# Internal functions utilizing this function:
# ## none
#==============================================================================#
boot_per_pop <- function(x, rg = multinom_boot, n, mle = NULL, H = TRUE, 
                         G = TRUE, lambda = TRUE, E5 = TRUE, ...){
  xi  <- extract_samples(x)
  res <- boot::boot(xi, boot_stats, R = n, sim = "parametric", ran.gen = rg, 
                    mle = mle, H = H, G = G, lambda = lambda, E5 = E5, ...)
  return(res)
}

#==============================================================================#
# multinomial sampler for bootstrapping
# 
# Public functions utilizing this function:
# ## diversity_boot
# 
# Internal functions utilizing this function:
# ## boot_per_pop
#==============================================================================#
multinom_boot <- function(x, mle = NULL){
  if (is.null(mle) || mle < 2){
    res <- stats::rmultinom(1, length(x), prob = tabulate(x))
  } else {
    res <- stats::rmultinom(1, mle, prob = tabulate(x))
  }
  extract_samples(res)
}

#==============================================================================#
# subsampler for rarefaction bootstrapping
# 
# Public functions utilizing this function:
# ## diversity_boot
# 
# Internal functions utilizing this function:
# ## boot_per_pop
#==============================================================================#
rare_sim_boot <- function(x, mle = 10){
  if (length(x) < mle){
    sample(x)
  } else {
    sample(x, mle)    
  }
}

#==============================================================================#
# wrapper to calculate statistics from bootstrapped samples
# 
# Public functions utilizing this function:
# ## none
# 
# Internal functions utilizing this function:
# ## boot_per_pop
#==============================================================================#
boot_stats <- function(x, i, H = TRUE, G = TRUE, lambda = TRUE, E5 = TRUE, ...){
  xi  <- tabulate(x[i])
  res <- diversity_stats(xi, H, G, lambda, E5, ...)
  return(res)
}

#==============================================================================#
# Sets up a vector of mlg counts for bootstrapping. Input is a vector where each
# element represents the count of a specific mlg, output is a vector where each
# element represents the mlg assignment of a particular sample. 
# 
# Public functions utilizing this function:
# ## none
# 
# Internal functions utilizing this function:
# ## boot_per_pop, multinom_boot
#==============================================================================#
extract_samples <- function(x) rep(1:length(x), x)

#==============================================================================#
# calculates confidence interval for a single population and a single statistic
# given the result of a bootstrap. Note, around_estimate is always true, and 
# will center the CI around the observed estimate. The structure for changing it
# is set up, but currently not used.
# 
# Public functions utilizing this function:
# ## none
# 
# Internal functions utilizing this function:
# ## get_ci
#==============================================================================#
get_boot_ci <- function(index, x, type = "normal", conf = 0.95, 
                        around_estimate = TRUE, ...){
  if (length(unique(x$t[, index])) == 1){
    return(c(NA_real_, NA_real_))
  } else if (around_estimate){
    res <- boot::norm.ci(conf = conf, t0 = x$t0[index], var.t0 = var(x$t[, index]))[-1]
  } else {
    res <- boot::norm.ci(x, conf = conf, index = index)[1, ]
    # res <- boot::boot.ci(x, conf, type, index, ...)[[type]][1, ]
  }
  return(utils::tail(res, 2))
}

#==============================================================================#
# calculates confidence intervals for all statistics per population. If 
# bci = TRUE, then the confidence interval is calculated using norm.ci or 
# boot.ci. Otherwise, the CI is calculated using quantiles from the bootstrapped
# data.
# 
# Public functions utilizing this function:
# ## none
# 
# Internal functions utilizing this function:
# ## get_all_ci
#==============================================================================#
get_ci <- function(x, lb, ub, bci = TRUE, btype = "normal", center = TRUE){
  if (bci){
    lenout <- length(x$t0)
    conf   <- diff(c(lb, ub))
    res    <- vapply(1:lenout, get_boot_ci, numeric(2), x, btype, conf, center)
    if (!is.null(dim(res))) rownames(res) <- paste(c(lb, ub)*100, "%")
  } else {
    res <- apply(x$t, 2, quantile, c(lb, ub), na.rm = TRUE)
    if (all(all.equal(res[1, ], res[2, ]) == TRUE)){
      res[] <- NA_real_
    }
  }
  return(res)
}

#==============================================================================#
# compiles an array of confidence intervals per statistic, per population.
# 
# Public functions utilizing this function:
# ## diversity_ci
# 
# Internal functions utilizing this function:
# ## none
#==============================================================================#
get_all_ci <- function(res, ci = 95, index_names = c("H", "G", "Hexp", "E.5"),
                       center = TRUE, btype = "normal", bci = TRUE){
  lower_bound  <- (100 - ci)/200
  upper_bound  <- 1 - lower_bound
  n_indices    <- length(index_names)
  funval       <- matrix(numeric(n_indices*2), nrow = 2)
  CI           <- vapply(res, FUN = get_ci, FUN.VALUE = funval, 
                         lower_bound, upper_bound, bci = bci, btype = btype,
                         center = center)
  dCI          <- dimnames(CI)
  dimnames(CI) <- list(CI    = dCI[[1]], 
                       Index = index_names,
                       Pop   = dCI[[3]])
  
  return(CI)
}

#==============================================================================#
# Takes the output of diversity_ci and converts it to a data frame.
# 
# Public functions utilizing this function:
# ## diversity_ci
# 
# Internal functions utilizing this function:
# ## none
#==============================================================================#
pretty_info <- function(obs, est, CI, boots = NULL){
  pretty_ci <- t(apply(round(CI, 3), 2:3, 
                       function(x){
                         if (all(is.na(x))) return(NA_character_)
                         paste0("(", paste(x, collapse = ", "), ")")
                       }))
  colnames(est) <- paste(colnames(est), "est", sep = ".")
  out <- vector(mode = "list", length = ncol(est)*3)
  colnames(pretty_ci) <- paste(colnames(pretty_ci), "ci", sep = ".")
  names(out)[(1:length(out)) %% 3 != 0] <- intersp(colnames(obs), colnames(est))
  names(out)[(1:length(out)) %% 3 == 0] <- colnames(pretty_ci)
  for (i in names(out)){
    out[[i]] <- obs[, 1]
  }
  out <- data.frame(out)
  out[colnames(obs)] <- obs
  out[colnames(est)] <- est
  out[colnames(pretty_ci)] <- pretty_ci
  class(out) <- c("popprtable", "data.frame")
  return(out)
}

#==============================================================================#
# Takes two vectors and intersperses their elements
# 
# Public functions utilizing this function:
# ## none
# 
# Internal functions utilizing this function:
# ## pretty_info
#==============================================================================#
intersp <- function(v1, v2){
  v1l <- length(v1)
  v2l <- length(v2)
  the_evens <- evens(v1l + v2l)
  stopifnot(v1l == v2l)
  ov <- vector(length = v1l + v2l, mode = class(c(v1[1], v2[1])))
  ov[the_evens]  <- v2
  ov[!the_evens] <- v1
  return(ov)
}

#==============================================================================#
# Returns a vector of even elements
# 
# Public functions utilizing this function:
# ## none
# 
# Internal functions utilizing this function:
# ## intersp
#==============================================================================#
evens <- function(x){
  if (length(x) > 1){
    if (is.numeric(x)){
      res <- x %% 2 == 0
    } else {
      res <- 1:length(x) %% 2 == 0
    }
  } else {
    res <- seq(x) %% 2 == 0
  }
  return(res)
}

#==============================================================================#
# Plots the results of bootstrapping and confidence interval calculation in
# facetted barplots
# 
# Public functions utilizing this function:
# ## diversity_ci
# 
# Internal functions utilizing this function:
# ## none
#==============================================================================#
#' Plot the results of boot_ci
#'
#' @param res a list of "boot" class objects
#' @param orig a table of the observed statistics
#' @param CI 
#'
#' @return prints a ggplot2 object
#' @noRd
#' @keywords internal
boot_plot <- function(res, orig, CI){
  statnames <- colnames(orig)
  popnames  <- rownames(orig)
  orig <- as.data.frame.table(orig, responseName = "value")
  if (!is.null(CI)) {
    cidf <- as.data.frame.table(CI, responseName = "v", stringsAsFactors = FALSE) %>%
      stats::reshape(idvar = c("Pop", "Index"), timevar = "CI", direction = "wide")
    names(cidf) <- gsub("^v\\.", "", names(cidf))
    orig <- merge(orig, cidf)
    colnames(orig)[4:5] <- c("lb", "ub")
  }
  samp <- vapply(res, "[[", FUN.VALUE = res[[1]]$t, "t")
  dimnames(samp) <- list(NULL, 
                         Index = statnames,
                         Pop = popnames)
  sampmelt <- as.data.frame.table(samp, responseName = "value")
  pl <- ggplot(sampmelt, aes_string(x = "Pop", y = "value", group = "Pop")) + 
    geom_boxplot() + 
    geom_point(aes_string(color = "Pop", x = "Pop", y = "value"), 
               size = rel(4), pch = 16, data = orig) +
    xlab("Population") + labs(color = "Observed") +
    facet_wrap(~Index, scales = "free_y") + myTheme
  if (!is.null(CI)){
    pl <- pl +
      geom_errorbar(aes_string(color = "Pop", x = "Pop", ymin = "lb", ymax = "ub"),
                    data = orig)
  }
  print(pl)
}

#==============================================================================#
# Extract the observed bootstrap statistics from a boot object. 
# 
# Public functions utilizing this function:
# ## diversity_ci
# 
# Internal functions utilizing this function:
# ## none
#==============================================================================#
get_boot_stats <- function(bootlist){
  npop   <- length(bootlist)
  bstats <- bootlist[[1]]$t0
  nstat  <- length(bstats)
  resmat <- matrix(nrow = npop, ncol = nstat,
                   dimnames = list(Pop = names(bootlist), Index = names(bstats)))
  resmat[] <- t(vapply(bootlist, FUN = "[[", FUN.VALUE = bstats, "t0"))
  return(resmat)
}

#==============================================================================#
# mean and sd methods for boot 
# 
# Public functions utilizing this function:
# ## none
# 
# Internal functions utilizing this function:
# ## get_boot_se
#==============================================================================#
mean.boot <- function(x, ...) apply(x$t, 2, mean, na.rm = TRUE)
sd.boot <- function(x, na.rm = TRUE) apply(x$t, 2, stats::sd, na.rm)

#==============================================================================#
# retrieves standard error or mean from list of boot objects.
# 
# Public functions utilizing this function:
# ## diversity_ci
# 
# Internal functions utilizing this function:
# ## none
#==============================================================================#
get_boot_se <- function(bootlist, res = "sd"){
  npop   <- length(bootlist)
  bstats <- bootlist[[1]]$t0
  nstat  <- length(bstats)
  THE_FUN <- match.fun(paste0(res, ".boot")) # mean.boot sd.boot
  resmat <- matrix(nrow = npop, ncol = nstat,
                   dimnames = list(Pop = names(bootlist), Index = names(bstats)))
  resmat[] <- t(vapply(bootlist, FUN = THE_FUN, FUN.VALUE = bstats))
  if (res == "sd"){
    colnames(resmat) <- paste(colnames(resmat), res, sep = ".")  
  }
  return(resmat)
}

#==============================================================================#
# Given a maximum and minimum value, this makes an n x 2 matrix of windows that
# encompass these values.
# 
# Public functions utilizing this function:
# ## win.ia
# 
# Internal functions utilizing this function:
# ## none
#==============================================================================#
make_windows <- function(maxp, minp = 1L, window = 100L){
  nwin   <- ceiling(maxp/window)
  minwin <- ceiling(minp/window)
  winmat <- matrix(window * seq(nwin), nrow = nwin, ncol = 2)[minwin:nwin, , drop = FALSE]
  winmat[, 1] <- winmat[, 1] - window + 1
  return(winmat)
}
#==============================================================================#
# add tied edges to minimum spanning tree
#
# Public functions utilizing this function:
# ## bruvo.msn poppr.msn
#
# Internal functions utilizing this function:
# ## singlepop_msn
#==============================================================================#
add_tied_edges <- function(mst, distmat, tolerance = .Machine$double.eps ^ 0.5){
  tied_edges <- .Call("msn_tied_edges", as.matrix(mst[]), as.matrix(distmat), tolerance)
  if (length(tied_edges) > 0){
    the_edges   <- tied_edges[c(TRUE, TRUE, FALSE)]
    the_weights <- tied_edges[c(FALSE, FALSE, TRUE)]
    mst <- add.edges(mst, dimnames(mst[])[[1]][the_edges], weight = the_weights)
  }
  return(mst)
}
#==============================================================================#
#' Correct minor allele frequencies derived from rraf (INTERNAL)
#' 
#' \strong{This is an internal function. The documentation is for use with 
#' \code{\link{rraf}}, \code{\link{pgen}}, and \code{\link{psex}}. Do not 
#' attempt to use this function directly.} Minor alleles are often lost when
#' calculating allele frequencies from a round-robin approach, resulting in
#' zero-valued allele frequencies (Arnaud-Haond et al. 2007, Parks and Werth
#' 1993). This can be problematic when calculating values for \code{\link{pgen}}
#' and \code{\link{psex}}. This function gives options for giving a value to
#' these zero-valued frequencies.
#' 
#' 
#' @param rraf \emph{internal} a list or matrix produced from \code{\link{rraf}}
#'   (with uncorrected MAF)
#' @param rrmlg \emph{internal} a matrix containing multilocus genotypes per 
#'   locus derived from \code{\link{rrmlg}}
#' @param e a numeric epsilon value to use for all missing allele frequencies.
#' @param sum_to_one when \code{TRUE}, the original frequencies will be reduced 
#'   so that all allele frequencies will sum to one. \strong{Default: 
#'   \code{FALSE}}
#' @param d the unit by which to take the reciprocal. \code{div = "sample"} will
#'   be 1/(n samples), \code{div = "mlg"} will be 1/(n mlg), and \code{div = 
#'   "rrmlg"} will be 1/(n mlg at that locus). This is overridden by \code{e}.
#' @param mul a multiplier for div. Default is \code{mult = 1}. This parameter
#'   is overridden by \code{e}
#' @param mlg \emph{internal} the number of MLGs in the sample. Only required if
#'   \code{d = "mlg"}.
#' @param pop \emph{internal} a vector of factors that define the population 
#'   definition for each observation in \code{rrmlg}. This must be supplied if 
#'   \code{rraf} is a matrix.
#' @param locfac \emph{internal} a vector of factors that define the columns 
#'   belonging to the loci.
#'   
#' @details Arguments of interest to the user are: 
#' \itemize{
#'  \item \strong{e}
#'  \item \strong{sum_to_one}
#'  \item \strong{d}
#'  \item \strong{m}
#' }
#' By default (\code{d = "sample", e = NULL, sum_to_one = FALSE, mul = 1}), this
#' will add 1/(n samples) to all zero-value alleles. The basic formula is
#' \strong{1/(d * m)} unless \strong{e} is specified. If \code{sum_to_one =
#' TRUE}, then the frequencies will be scaled as x/sum(x) AFTER correction,
#' indicating that the allele frequencies will be reduced. See the examples for
#' details. The general pattern of correction is that the value of the MAF will
#' be \emph{rrmlg > mlg > sample} 
#'   
#' @return a matrix or vector the same type as rraf
#' @author Zhian N. Kamvar
#' @noRd
#' @references
#' 
#' Arnaud-Haond, S., Duarte, C. M., Alberto, F., & Serrão, E. A. 2007.
#' Standardizing methods to address clonality in population studies.
#' \emph{Molecular Ecology}, 16(24), 5115-5139.
#' 
#' Parks, J. C., & Werth, C. R. 1993. A study of spatial features of clones in a
#' population of bracken fern, \emph{Pteridium aquilinum} (Dennstaedtiaceae).
#' \emph{American Journal of Botany}, 537-544.
#' 
#' @seealso \code{\link{rraf}}, 
#'   \code{\link{pgen}}, 
#'   \code{\link{psex}}, 
#'   \code{\link{rrmlg}}
#'   
#' @examples
#' \dontrun{
#' 
#' data(Pram)
#' #-------------------------------------
#' 
#' # If you set correction = FALSE, you'll notice the zero-valued alleles
#' 
#' rraf(Pram, correction = FALSE)
#' 
#' # By default, however, the data will be corrected by 1/n
#' 
#' rraf(Pram)
#' 
#' # Of course, this is a diploid organism, we might want to set 1/2n
#' 
#' rraf(Pram, mul = 1/2)
#' 
#' # To set MAF = 1/2mlg
#' 
#' rraf(Pram, d = "mlg", mul = 1/2)
#' 
#' # Another way to think about this is, since these allele frequencies were
#' # derived at each locus with different sample sizes, it's only appropriate to
#' # correct based on those sample sizes.
#' 
#' rraf(Pram, d = "rrmlg", mul = 1/2)
#' 
#' # If we were going to use these frequencies for simulations, we might want to
#' # ensure that they all sum to one. 
#' 
#' rraf(Pram, d = "mlg", mul = 1/2, sum_to_one = TRUE) 
#' 
#' #-------------------------------------
#' # When we calculate these frequencies based on population, they are heavily
#' # influenced by the number of observed mlgs. 
#' 
#' rraf(Pram, by_pop = TRUE, d = "rrmlg", mul = 1/2)
#' 
#' # This can be fixed by specifying a specific value
#' 
#' rraf(Pram, by_pop = TRUE, e = 0.01)
#' 
#' }
#' 
#==============================================================================#
rare_allele_correction <- function(rraf, rrmlg, e = NULL, sum_to_one = FALSE, 
                                    d = c("sample", "mlg", "rrmlg"), mul = 1, 
                                    mlg = NULL, pop = NULL, locfac = NULL){
  
  if (identical(parent.frame(), globalenv())){
    msg <- paste0("\n\n\n",
                  "    !!! rare_allele_correction() is for internal use only !!!",
                  "\n\n",
                  "    Input types are not checked within this function and may\n",
                  "    result in an error. USE AT YOUR OWN RISK.\n\n\n")
    warning(msg, immediate. = TRUE)
  }
  d <- match.arg(d, c("sample", "mlg", "rrmlg"))

  if (is.list(rraf)){
    if (is.null(e)){
      e <- get_minor_allele_replacement(rrmlg, d, mul, mlg)
    }
    if (length(e) == 1){
      e <- setNames(rep(e, ncol(rrmlg)), colnames(rrmlg))
    }
    res        <- lapply(names(rraf), replace_zeroes, rraf, e, sum_to_one)
    names(res) <- names(rraf)
  } else if (is.matrix(rraf)){
    
    # split matrix by population and locus
    poplist <- apply(rraf, 1, split, locfac)
    
    # loop over populations and call this function again on the list of loci
    res <- lapply(names(poplist), function(i){
      prraf  <- poplist[[i]]
      prrmlg <- rrmlg[pop == i, ]
      rare_allele_correction(prraf, prrmlg, mlg = mlg, e = e, d = d, mul = mul, 
                              sum_to_one = sum_to_one)
    })
    
    # make this list of loci a matrix again
    res           <- t(vapply(res, unlist, rraf[1, , drop = TRUE]))
    dimnames(res) <- dimnames(rraf)
  } 
  return(res)
}
#==============================================================================#
# round robin clone-correct allele frequency estimates
# 
# This will take in
# - i the index of the locus.
# - loclist a list of genind objects each one locus.
# - mlgs a matrix from rrmlg
#
# Public functions utilizing this function:
# ## rraf
#
# Internal functions utilizing this function:
# ## none
#==============================================================================#
rrcc <- function(i, loclist, mlgs){
  cc  <- !duplicated(mlgs[, i])
  mat <- tab(loclist[[i]], freq = TRUE)
  res <- colMeans(mat[cc, , drop = FALSE], na.rm = TRUE)
  return(res)
}
#==============================================================================#
#' round robin clone-correct allele frequency estimates by population
#' 
#' @param i the index of the locus.
#' @param loclist a list of genind objects each one locus.
#' @param mlgs a matrix from rrmlg
#' @param pnames the population names.
#'
#' @note
#' Public functions utilizing this function:
#' ## rraf
#'
#' Internal functions utilizing this function:
#' ## none
#' @noRd
#==============================================================================#
rrccbp <- function(i, loclist, mlgs, pnames){
  
  mat  <- tab(loclist[[i]], freq = TRUE)
  npop <- length(pnames)
  res  <- matrix(numeric(ncol(mat)*npop), nrow = npop)
  rownames(res) <- pnames
  colnames(res) <- colnames(mat)
  pops <- pop(loclist[[i]])
  for (p in pnames){
    psub <- which(pops %in% p)
    cc   <- psub[!duplicated(mlgs[psub, i])]
    out  <- colMeans(mat[cc, , drop = FALSE], na.rm = TRUE)
    res[p, ] <- out
  }
  
  return(res)
}

#' Treat the optional "G" argument for psex
#'
#' @param G either NULL or an integer vector that can be named or not
#' @param N an integer or integer vector
#' @param dat a data set
#' @param population a population name
#' @param method passed from psex
#' 
#' @note 
#' Public functions: psex
#' Private functions: none
#'
#' @return a value for G
#' @noRd
treat_G <- function(G, N, dat, population, method){
  if (is.null(G)){
    return(N)
  } else if (length(G) == 1){
    return(G)
  } else if (!is.null(names(G)) && all(names(G) %in% popNames(dat))){
    if (method == "multiple"){
      G <- G[population]
    } else {
      G <- G[pop(dat)]
    }
  } else if (length(G) > 1){
    stop("G must be a named vector of integers.")
  } else {
    stop("G must be NULL or an integer vector of length 1 or nPop(gid)", 
         call. = FALSE)
  }
  G
}

#==============================================================================#
# replace missing allele frequencies with a specified value
# 
# @param i the name or index of a locus
# @param loci a list of allele frequencies by locus derived from the round-robin
#   method
# @param e a value to replace zero-value frequencies with
# @param sum_to_one a boolean indicating whether or not allele frequencies at a
#   single locus should sum to one.
#
# Public functions utilizing this function:
# ## rare_allele_correction
#
# Internal functions utilizing this function:
# ## none
#==============================================================================#

replace_zeroes <- function(i, loci, e, sum_to_one = FALSE){
  locus           <- loci[[i]]
  missing_alleles <- locus <= .Machine$double.eps^0.5
  locus[missing_alleles] <- e[i]
  if (sum_to_one){
    locus <- locus/sum(locus)
  }
  return(locus)
}

#==============================================================================#
# prepare the value to set the minor allele
#
# @param rrmlg a matrix derived from rrmlg
# @param d one of "sample", "mlg", or "rrmlg"
# @param m a multiplier
# @param mlg an integer or NULL specifying the number of unique multilocus
#   genotypes in the sample.
#
# Public functions utilizing this function:
# ## rare_allele_correction
#
# Internal functions utilizing this function:
# ## none
#==============================================================================#

get_minor_allele_replacement <- function(rrmlg, d, m, mlg = NULL){
  if (d == "sample"){
    e <- (1/nrow(rrmlg)) * m
  } else if (d == "mlg"){
    e <- (1/mlg) * m
  } else {
    clones <- !apply(rrmlg, 2, duplicated)
    e      <- setNames((1/colSums(clones)) * m, colnames(rrmlg))
  }
  return(e)
}

#==============================================================================#
# Pairwise index of association calculation
# 
# Arguments
# - pair a character vector of length 3 containing the pair of loci and the 
#   iteration.
# - V a choose(n, 2) x m matrix containing the distances at all loci. 
# - np choose(n, 2)
# - progbar either NULL or a progressbar object.
# - iterations the number of loci combinations.
#
# Public functions utilizing this function:
# ## pair.ia
#
# Internal functions utilizing this function:
# ## none
#==============================================================================#
ia_pair_loc <- function(pair, V, np, progbar, iterations){
  if (!is.null(progbar)){
    setTxtProgressBar(progbar, as.numeric(pair[3])/iterations)
  }
  newV <- V[, pair[-3]]
  V    <- list(d.vector  = colSums(newV), 
               d2.vector = colSums(newV * newV), 
               D.vector  = rowSums(newV)
  )
  return(ia_from_d_and_D(V, np))
}
#==============================================================================#
# Get the index of association from sums of distances over loci and samples
# 
# This will take in a list called that contains 3 vectors:
# 1. d.vector : a vector containing sums of distances per locus
# 2. d2.vector : like d.vector, but containing sums of squares
# 3. D.vector : a vector containing the distance over all samples.
#
# Public functions utilizing this function:
# ## none
#
# Internal functions utilizing this function:
# ## ia_pair_loc
#==============================================================================#
ia_from_d_and_D <- function(V, np){
  varD <- ((sum(V$D.vector^2) - ((sum(V$D.vector))^2)/np))/np
  vard.vector <- ((V$d2.vector - ((V$d.vector^2)/np))/np)
  vardpair.vector <- .Call("pairwise_covar", vard.vector, PACKAGE = "poppr")
  sigVarj <- sum(vard.vector)
  rm(vard.vector)
  Ia <- (varD/sigVarj) - 1
  rbarD <- (varD - sigVarj)/(2 * sum(vardpair.vector))
  return(c(Ia, rbarD))
}


#' scans genlight objects for genotypes that have no alternate homozygous sites
#' and manually adds a vector of zeroes in for bitwise functions to work.
#'
#' @param x a genlight or snpclone object
#'
#' @return a genlight object
#' @noRd
#'
fix_uneven_diploid <- function(x){
  snplens <- x$gen %>% lapply(slot, "snp") %>% vapply(length, integer(1))
  if (any(snplens == 1)){
    replacement <- as.raw(rep(0, length(x$gen[[1]]$snp[[1]])))
    for (i in which(snplens == 1)){
      x$gen[[i]]$snp[[2]] <- replacement
    }
  }
  return(x)
}

#' MLG index handler for the "[[" methods
#' 
#' This handles the indices for MLGs within the subset methods since they don't 
#' know what the sample names are.
#' @param i a vector used to subset data
#' @param gid a genclone or snpclone object
#'   
#' @return a vector of integers or logicals
#' @noRd
#' 
#' @seealso \code{\link{initialize,genclone-method}} 
#'   \code{\link{initialize,snpclone-method}}
handle_mlg_index <- function(i, gid){
  if (is.logical(i) || is.numeric(i)){
    return(i)
  }
  i <- if (is.factor(i)) as.character(i) else i
  i <- match(i, indNames(gid))         # match to the index names
  i <- i[!is.na(i)]                    # remove missing indices
  i <- if (length(i) == 0) TRUE else i # ignore if no result
  return(i)
}

#' Population index handler for the "[[" methods
#'
#' If the user specifies populations, it takes priority
#'
#' @param pops vector of population indices or names
#' @param gid a genclone or snpclone object
#'
#' @return a vector of logicals for each sample
#' @noRd
#' @seealso \code{\link{initialize,genclone-method}} 
#'   \code{\link{initialize,snpclone-method}}
handle_pops_index <- function(pops, gid){
  if (!is.null(pops) && !is.null(pop(gid))){
    pops <- if (is.factor(pops)) as.character(pops) else pops
    pops <- if (!is.character(pops)) popNames(gid)[pops] else pops
    newi <- pop(gid) %in% pops
  } 
  return(newi)
}

#' Should poppr be quiet?
#' 
#' If it's not an interactive session or it's in a knitr document, messages in
#' poppr should be suppressed UNLESS poppr.debug is set to TRUE
#'
#' @return TRUE or FALSE
#' @noRd
#'
should_poppr_be_quiet <- function(quiet){
  # Suppress the noise if it's not interactive or in knitr
  # This is thanks to @jimhester
  # https://github.com/hadley/dplyr/commit/c8beb59217620614b36cd82df0a7e89c556fb374
  in_knitr    <- !is.null(getOption("knitr.in.progress"))
  in_script   <- !interactive()
  poppr_debug <- getOption("poppr.debug")
  if ((in_script || in_knitr) && !poppr_debug){
    quiet <- TRUE
  }
  return(quiet)
}
#' Set the population from a formula or vector
#'
#' @param gid a genind/genclone/genlight/snpclone object
#' @param pop a formula or factor specifying population
#'
#' @return the object with the population set in the pop slot
#' @noRd
set_pop_from_strata_or_vector <- function(gid, pop) {
  if (is.language(pop)){ # incoming is a formula, e.g. ~Country/Year
    setPop(gid) <- pop
  } else {               # incoming is a vector
       pop(gid) <- pop
  }
  gid
}
# poppr's theme for ggplot2 (mainly rotating x axis labels)
myTheme <- theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


#' Calculate Psex for a given genotype
#'
#' @param n_encounters integer the number of observations of an MLG
#' @param p_genotype numeric the probability of the given MLG
#' @param sample_ids list the names of the given MLG
#' @param n_samples integer the number of samples observed
#'
#' @return a vector of probabilities
#' @noRd
#'
#' @examples
#' make_psex(10, 0.005, n_samples = 100)
make_psex <- function(n_encounters, p_genotype, sample_ids = NULL, n_samples){
  encounters <- seq(n_encounters) - 1
  out        <- dbinom(encounters, n_samples, p_genotype)
  names(out) <- sample_ids
  return(out)
}

cromulent_replen <- function(gid, replen){
  the_loci <- locNames(gid)
  if (length(replen) != nLoc(gid) && is.null(names(replen))){
    msg <- mismatched_repeat_length_warning(replen, nLoc(gid))
    stop(msg, call. = FALSE)
  } else if (length(replen) < nLoc(gid)){
    msg <- mismatched_repeat_length_warning(replen, nLoc(gid))
    stop(msg, call. = FALSE)
  } else if (length(replen) > nLoc(gid)){
    msg <- trimmed_repeats_warning(replen, the_loci)
    warning(msg, call. = FALSE, immediate. = TRUE)
  } 
  new_replen <- match_replen_to_loci(the_loci, replen)
  return(new_replen)
}


#' calculate the distance between two circles
#'
#' @param c1 radius of circle 1
#' @param c2 radius of circle 2
#'
#' @return the distance between two centers on one axis
#' @noRd
#' @details 
#' https://stackoverflow.com/a/14830596/2752888
cdist <- function(c1, c2){
  a <- (c1 + c2) ^ 2
  b <- (c1 - c2) ^ 2
  sqrt(a - b)
}

# spread the circles out 
#' Arrange circles adjacent to each other
#'
#' @param radii a vector indicating the radii of the circles to be arranged.
#'
#' @return a vector of positions on which to plot the centers of circles
#' @noRd
make_adjacent_circles <- function(radii){
  res <- vapply(seq(radii), function(i) {
    if (i == 1)
      0.0
    else
      cdist(radii[i], radii[i - 1])
  }, numeric(1))
  cumsum(res)
}

#' Make labels for circles
#' 
#' This function takes in a range of numbers and figures out three circles that
#' exist in these numbers to create a legend using the range and approximate 
#' mean of the numbers
#'
#' @param mlg_number a vector of integers
#'
#' @return 1, 2, or 3 integers.
#' @noRd
#'
#' @examples
#' set.seed(40)
#' make_circle_labs(sample(100, 50, replace = TRUE))
make_circle_labs <- function(mlg_number){
  labs   <- range(mlg_number)
  if (diff(labs) > 1) {
    m    <- mean(labs)
    m    <- mlg_number[which.min(abs(mlg_number - m))]
    labs <- c(labs[1], m, labs[2])
  } else if (diff(labs) == 0) {
    labs <- labs[1]
  }
  labs
}

#' Get the correct side of the legend box
#'
#' @param a the result of a "legend" call
#' @param side either "top", "bottom", "right", or "left"
#'
#' @return a vector of at max length 4
#' @noRd
#'
#' @examples
#' plot(seq(-1, 1), seq(-1, 1), asp = 1)
#' tr <- legend("topright", fill = "blue", legend = "blue")
#' tl <- legend("topleft", fill = "blue", legend = "blue")
#' bl <- legend("bottomleft", fill = "blue", legend = "blue")
#' br <- legend("bottomright", fill = "blue", legend = "blue")
#' ce <- legend("center", fill = "blue", legend = "blue")
#' points(x = get_legend_side(ce, 3:4), 
#'        y = get_legend_side(ce, 1:2), 
#'        col = c("red", "blue"))
#' points(x = get_legend_side(tl, 3:4), 
#'        y = get_legend_side(tl, 1:2), 
#'        col = c("red", "blue"))
#' points(x = get_legend_side(bl, 3:4), 
#'        y = get_legend_side(bl, 1:2), 
#'        col = c("red", "blue"))
#' points(x = get_legend_side(br, 3:4), 
#'        y = get_legend_side(br, 1:2), 
#'        col = c("red", "blue"))
#' points(x = get_legend_side(tr, 3:4), 
#'        y = get_legend_side(tr, 1:2), 
#'        col = c("red", "blue"))
get_legend_side <- function(a, side = NULL) {
  res <- c(
    top    = a$rect$top,
    bottom = a$rect$top  - a$rect$h,
    right  = a$rect$left + a$rect$w,
    left   = a$rect$left
  )
  return( if (is.null(side)) res else res[side] )
}

#' Create a circle legend
#'
#' @param a the output of a "legend" call. Defaults to `NULL`.
#' @param mlg_number a vector of integers
#' @param scale a number by which to multiply the node sizes
#' @param cex character expansion for the text
#' @param radmult multiplier for the radius (specifically for [plot_poppr_msn])
#' @param xspace the defined xspacer (currently is wonky :/)
#' @param font font for the title
#' @param pos a numeric position for the title or NULL
#' @param x the position for the leftmost part of the legend
#' @param y the position of the topmost part of the legend
#' @param txt The position of the text
#' 
#' @details
#'   If a legend is provided in `a`, then the position of the circle legend will
#'   be below it. If not, the arguments left, top, and txt will control the
#'   position of the text within the box. 
#'
#' @return a legend with circles
#'
#' @noRd
make_circle_legend <-
  function(a = NULL,
           mlg_number,
           scale = 5,
           cex = 0.75,
           radmult = 1,
           xspace = NULL,
           font = 1,
           pos = NULL,
           x = -1.55,
           y = 1,
           txt = -1.3) {
    
  
  # Create example circles for comparison
  labs   <- make_circle_labs(mlg_number)
  rads   <- (sqrt(labs) * scale)/200
  
  if (!is.null(a)){
    txt    <- a$text$x[1]    
    left   <- get_legend_side(a, "left")   
    bottom <- get_legend_side(a, "bottom")   
  } else {
    left   <- x
    bottom <- y
  }
  
  # Get the space between legend elements
  yspace <- if (!is.null(a)) diff(a$text$y) else 0.05
  yspace <- if (length(yspace) > 0L) min(abs(yspace)) else 0.05
  xspace <- if (is.null(xspace)) 0.25 * yspace else 2 * xspace * yspace
  
  # Create positions of circles horizontally
  diam <- max(rads) * 2
  big_circles <- diam > abs(txt - left)/2
  circlx <- if (big_circles) left + diam + xspace else txt
  circlx <- rep(circlx, length(rads))
  
  # shift the y position of the circles
  cpos   <- make_adjacent_circles(rads)
  circly <- bottom - ((2.5 * yspace) + max(cpos) + max(rads)/2)
  circly <- cpos + circly
  
  # Create the circle legend
  xpos <- if (is.null(pos)) left else pos
  tadj <- if (is.null(pos)) c(0, 0.5) else c(0.5, 0.5)
  graphics::text(x     = xpos, 
                 y     = bottom - yspace, 
                 label = "Samples/Node",
                 adj   = tadj,
                 font  = font,
                 cex   = cex)
  
  graphics::symbols(x       = circlx - (rads * radmult + xspace), 
                    y       = circly, 
                    circles = rads * radmult, 
                    add     = TRUE, 
                    inches  = FALSE, 
                    asp     = 1)
  
  graphics::text(x      = circlx, 
                 y      = circly, 
                 labels = labs, 
                 adj    = c(0, 0.5),
                 cex    = cex)
}