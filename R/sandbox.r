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


new.read.genalex <- function(genalex, ploidy=2, geo=FALSE, region=FALSE){
  # The first two lines from a genalex file contain all of the information about
  # the structure of the file (except for ploidy and geographic info)
  gencall <- match.call()
  all.info <- strsplit(readLines(genalex, n=2), ",")
  if (any(all.info[[1]]=="")){
    num.info <- as.numeric(all.info[[1]][-which(all.info[[1]]=="")])
  }
  else {
    num.info <- as.numeric(all.info[[1]])
  }
  if (any(all.info[[2]]=="")){
    pop.info <- all.info[[2]][c(-1,-2,-3,-which(all.info[[2]]==""))]  
  }
  else {
    pop.info <- all.info[[2]][c(-1,-2,-3)]
  }
  
  # glob.info is for the number of loci, individuals, and populations.
  glob.info <- num.info[1:3]
  if(any(is.na(glob.info))){
    stop("Something is wrong with your csv file. Perhaps it is not comma delimited?\n")
  }
  #   cat("Global Information:",glob.info,"\n")
  #   cat("Populations:",pop.info,"\n")
  #   cat("num.info:",num.info,"\n")
  
  gena <- read.csv(genalex, header=TRUE, skip=2)
  
  # Removing all null columns 
  if(!is.na(which(is.na(gena[1, ]))[1])){
    gena <- gena[, -which(is.na(gena[1, ]))]
  }
  
  #----------------------------------------------------------------------------#
  # Checking for extra information such as Regions or XY coordinates
  #----------------------------------------------------------------------------#
  
  # Creating vectors that correspond to the different information fields.
  # If the regions are true, then the length of the pop.info should be equal to
  # the number of populations "npop"(glob.info[3]) plus the number of regions
  # which is the npop+4th entry in the vector. 
  # Note that this strategy will only work if the name of the first region does
  # not match any of the populations. 
  
  if (region==TRUE & length(pop.info) == glob.info[3] + num.info[glob.info[3]+4]){
    # Info for the number of columns the loci can take on.
    loci.adj <- c(glob.info[1], glob.info[1]*ploidy)
    
    # First question, do you have two or four extra columns? Two extra would
    # indicate no geographic data. Four extra would indicate geographic data.
    # Both of these indicate that, while a regional specification exists, a 
    # column indicating the regions was not specified, so it needs to be created
    if(((ncol(gena) %in% (loci.adj + 4)) & (geo == TRUE)) | (ncol(gena) %in% (loci.adj + 2))){
      pop.vec <- gena[, 2]
      ind.vec <- gena[, 1]
      xy <- gena[, c((ncol(gena)-1), ncol(gena))]
      
      # Get the indices for the regions
      region.inds <- ((glob.info[3]+5):length(num.info))
      # Get the number of individuals per region
      reg.inds <- num.info[region.inds]
      # Get the names of the regions
      reg.names <- all.info[[2]][region.inds]
      # Paste them all into a single vector.
      reg.vec <- rep(reg.names, reg.inds)
      if(geo == TRUE){
        geoinds <- c((ncol(gena)-1), ncol(gena))
        xy <- gena[, geoinds]
        gena <- gena[, -geoinds]
      }
      else{
        xy <- NULL
      }
      gena <- gena[, c(-1,-2)]
    }
    #
    # The Regions are specified in one of the first two columns.
    #
    else{
      pop.vec <- ifelse(any(gena[, 1] == pop.info[1]), 1, 2)
      reg.vec <- ifelse(pop.vec == 2, 1, 2)
      orig.ind.vec <- NULL
      # Regional Vector    
      reg.vec <- gena[, reg.vec]
      # Population Vector
      pop.vec <- gena[, pop.vec]    
      if(geo == TRUE){
        geoinds <- c((ncol(gena)-1), ncol(gena))
        xy <- gena[, geoinds]
        gena <- gena[, -geoinds]
      }
      else{
        xy <- NULL
      }
      # Individual Vector
      ind.vec <- gena[, ncol(gena)]
      # removing the non-genotypic columns from the data frame
      gena <- gena[, c(-1,-2,-ncol(gena))]
    }
  }
  
  #
  # There are no Regions specified, but there are geographic coordinates
  #
  else if (geo == TRUE & length(pop.info) == glob.info[3]){
    reg.vec <- NULL
    pop.vec <- gena[, 2]
    ind.vec <- gena[, 1]
    xy <- gena[, c((ncol(gena)-1), ncol(gena))]
    gena <- gena[, c(-1,-2,-(ncol(gena)-1),-ncol(gena))]
  }
  #
  # There are no Regions or geographic coordinates
  #
  else{
    reg.vec <- NULL
    pop.vec <- gena[, 2]
    ind.vec <- gena[, 1]
    xy <- NULL
    gena <- gena[, c(-1,-2)]
  }
  
  #----------------------------------------------------------------------------#
  # The genotype matrix has been isolated at this point. Now this will
  # reconstruct the matrix in a way that adegenet likes it.
  #----------------------------------------------------------------------------#
  
  clm <- ncol(gena)
  gena.mat <- as.matrix(gena)
  # Checking for greater than haploid data.
  if (glob.info[1] == clm/ploidy){
    # Missing data in genalex is coded as "0" for non-presence/absence data.
    # this converts it to "NA" for adegenet.
#     if(any(gena.mat =="0")){
#       #gena.mat <- as.matrix(gena)
#       gena[gena.mat == "0"] <- NA
#       #gena <- as.data.frame(gena.mat)
#     }
    type <- 'codom'
    loci <- which((1:clm) %% ploidy==1)
    gena2 <- gena[, loci]
    lapply(loci, function(x) gena2[, ((x-1)/ploidy)+1] <<-
             pop_combiner(gena, hier = x:(x+ploidy-1), sep = "/"))
    #res <- list(Gena=gena2, Glob.info=glob.info, ploidy=ploidyy)
    res.gid <- df2genind(gena2, sep="/", ind.names=ind.vec, pop=pop.vec,
                         ploidy=ploidy, type=type)
  }
  # Checking for AFLP data.
  else if (glob.info[1] == clm & all(gena.mat %in% as.integer(-1:1))) {
    # Missing data in genalex is coded as "-1" for presence/absence data.
    # this converts it to "NA" for adegenet.
    if(any(gena.mat == -1L)){
      #gena.mat <- as.matrix(gena)
      gena[gena.mat == -1L] <- NA
      #gena <- as.data.frame(gena.mat)
    }
    type <- 'PA'
    #res <- list(Gena=gena, Glob.info=glob.info, Ploid=ploidy)
    res.gid <- df2genind(gena, ind.names=ind.vec, pop=pop.vec,
                         ploidy=ploidy, type=type)
  }
  # Checking for haploid microsattellite data or SNP data
  else if (glob.info[1] == clm & !all(gena.mat %in% as.integer(-1:1))) {
    if(any(gena.mat == "0")){
      #gena.mat <- as.matrix(gena)
      gena[gena.mat == "0"] <- NA
      #gena <- as.data.frame(gena.mat)
    }
    type <- 'codom'
    #res <- list(Gena=gena, Glob.info=glob.info, Ploid=1)
    res.gid <- df2genind(gena, ind.names=ind.vec, pop=pop.vec,
                         ploidy= 1, type=type)
  }
  else {
    stop("Something went wrong. Check your geo and region flags to make sure they are set correctly. Otherwise, the problem may lie within the data structure itself.")
  }
  if (any(duplicated(ind.vec))){
    # ensuring that all names are unique
    res.gid@ind.names <- paste("ind",1:length(ind.vec))
    res.gid@other[["original_names"]] <- ind.vec
  }
  
  res.gid@other[["population_hierarchy"]] <- as.data.frame(list(Pop=pop.vec))
  res.gid@call <- gencall
  res.gid@call[2] <- basename(genalex)
  if(region==TRUE){
    res.gid@other[["population_hierarchy"]]$Region <- reg.vec
    #return(res.gid)
  }
  if(geo==TRUE){
    res.gid@other[["xy"]] <- xy
    #return(res.gid)
  }
  return(res.gid)
}



.new.permut.shuff <- function(mat, ploidy = 2){
  # gathering the counts of allelic states (represented as column numbers)
  bucket <- as.vector(colSums(mat, na.rm=T) * ploidy)
  bucketlist <- rep(1:length(bucket), bucket)
  # reshuffling the allelic states and making a matrix with rows for each
  # allele and columns for each individual.
  samplist <- sample(bucketlist)
  newmat <- mat  
  newmat[!is.na(mat)] <- 0
  
  # Treating Missing Data
  #
  # go on to the next locus if all are NA
  if(all(is.na(mat))){
    next
  }
  # placing the individuals with missing data in the sample list so that it is
  # the same length as the number of rows * 2 to avoid errors.
  if(any(is.na(mat))){
    temp <- vector(mode="numeric", length=nrow(mat) * ploidy)
    miss <- which(is.na(mat[, 1]))
    temp[c(-miss * ploidy, -((miss * ploidy)-1))] <- samplist
    samplist <- temp
    #return(samplist)
  }
  samplist <- matrix(samplist, nrow = ploidy)
  # A function that will repopulate the new matrix with the shuffled alleles
  # x is a list of even numbers, newmat is an empty matrix, and samplist is 
  # a list of shuffled alleleic states.
  populate.mat <- function(x, newmat, samplist, ploidy){
    samps <- table(samplist[, x])/ploidy
    newmat[x, as.numeric(names(samps))] <<- samps
  }
  lapply(1:nrow(mat), populate.mat, newmat, samplist, ploidy)
  return(newmat)
}

testing_funk <- function(){
  cat("This test worked...maybe.\n")
}





new.poppr <- function(pop,total=TRUE,sublist=c("ALL"),blacklist=c(NULL), sample=0,
  method=1,missing="ignore", quiet="minimal",clonecorrect=FALSE,hier=c(1),dfname="hier",
  hist=TRUE, minsamp=10){
  METHODS = c("multilocus", "permute alleles", "parametric bootstrap",
      "non-parametric bootstrap")
	x <- .file.type(pop, missing=missing, clonecorrect=clonecorrect, hier=hier, 
                  dfname=dfname, quiet=TRUE)	
  # The namelist will contain information such as the filename and population
  # names so that they can easily be ported around.
  namelist <- NULL
  callpop <- match.call()
  if(!is.na(grep("system.file", callpop)[1])){
    popsplt <- unlist(strsplit(pop, "/"))
    namelist$File <- popsplt[length(popsplt)]
  }
  else if(is.genind(pop)){
    namelist$File <- x$X
  }
  else{
    namelist$File <- basename(x$X)
  }
  #poplist <- x$POPLIST
  pop <- popsub(x$GENIND, sublist=sublist, blacklist=blacklist)
  poplist <- .pop.divide(pop)
  # Creating the genotype matrix for vegan's diversity analysis.
  pop.mat <- mlg.matrix(pop)
  if (total==TRUE & !is.null(poplist)){
    poplist$Total <- pop
    pop.mat <- rbind(pop.mat, colSums(pop.mat))
  }
  sublist <- names(poplist)
  Iout <- NULL
  result <- NULL
  origpop <- x$GENIND
  rm(x)
  total <- toupper(total)
  missing <- toupper(missing)
  type <- pop@type
  # For presence/absences markers, a different algorithm is applied. 
  if(type=="PA"){
    .Ia.Rd <- .PA.Ia.Rd
  }
  #''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''#
  # Creating an indicator for multiple subpopulations.
  # MPI = Multiple Population Indicator
  #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,#
  if (is.null(poplist)){
    MPI <- NULL
  }
  else{
    MPI <- 1
  }

  #''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''#
  # ANALYSIS OF MULTIPLE POPULATIONS.
  #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,# 
 
	if (!is.null(MPI)){
    
    #''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''#
    # Calculations start here.
    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,#

    MLG.vec <- vapply(sublist, function(x) mlg(poplist[[x]], quiet=TRUE), 1)
    N.vec <- vapply(sublist, function(x) length(poplist[[x]]@ind.names), 1)
    # Shannon-Weiner vegan:::diversity index.
    H <- vegan:::diversity(pop.mat)
    # E_1, Pielou's evenness.
    # J <- H / log(rowSums(pop.mat > 0))
    # inverse Simpson's index aka Stoddard and Taylor: 1/lambda
    G <- vegan:::diversity(pop.mat, "inv")
    Hexp <- (N.vec/(N.vec-1))*vegan:::diversity(pop.mat, "simp")
    # E_5
    E.5 <- (G-1)/(exp(H)-1)
    # rarefaction giving the standard errors. This will use the minimum pop size
    # above a user-defined threshold.
    raremax <- ifelse(is.null(nrow(pop.mat)), sum(pop.mat), 
                      ifelse(min(rowSums(pop.mat)) > minsamp, 
                             min(rowSums(pop.mat)), minsamp))
    
    N.rare <- rarefy(pop.mat, raremax, se=TRUE)
    IaList <- NULL
    invisible(lapply(sublist, function(x) 
                              IaList <<- rbind(IaList, 
                                               .new.ia(poplist[[x]], 
                                                       sample=sample, 
                                                       method=method, 
                                                       quiet=quiet, 
                                                       missing=missing, 
                                                       namelist=list(File=namelist$File, population = x),
                                                       hist=hist
              ))))

    #''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''#
    # Making the data look pretty.
    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,#
    Iout <- as.data.frame(list(Pop=sublist, N=N.vec, MLG=MLG.vec, 
                                eMLG=round(N.rare[1, ], 3), 
                                SE=round(N.rare[2, ], 3), 
                                H=round(H, 3), 
                                G=round(G,3),
                                Hexp=round(Hexp, 3),
                                E.5=round(E.5,3),
                                IaList,
                                File=namelist$File))
    rownames(Iout) <- NULL
    return(final(Iout, result))
	}
  #''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''#
  # ANALYSIS OF SINGLE POPULATION. This is for if there are no subpopulations to
  # be analyzed. For details of the functions utilized, see the notes above.
  #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,#
  else { 
    MLG.vec <- mlg(pop, quiet=TRUE)
    N.vec <- length(pop@ind.names)
    # Shannon-Weiner vegan:::diversity index.
    H <- vegan:::diversity(pop.mat)
    # E_1, Pielou's evenness.
    # J <- H / log(rowSums(pop.mat > 0))
    # inverse Simpson's index aka Stoddard and Taylor: 1/lambda
    G <- vegan:::diversity(pop.mat, "inv")
    Hexp <- (N.vec/(N.vec-1))*vegan:::diversity(pop.mat, "simp")
    # E_5
    E.5 <- (G-1)/(exp(H)-1)
    # rarefaction giving the standard errors. No population structure means that
    # the sample is equal to the number of individuals.
    N.rare <- rarefy(pop.mat, sum(pop.mat), se=TRUE)
    IaList <- .new.ia(pop, sample=sample, method=method, quiet=quiet, missing=missing,
                      namelist=(list(File=namelist$File, population="Total")),
                      hist=hist)
    Iout <- as.data.frame(list(Pop="Total", N=N.vec, MLG=MLG.vec, 
                          eMLG=round(N.rare[1, ], 3), 
                          SE=round(N.rare[2, ], 3),
                          H=round(H, 3), 
                          G=round(G,3), 
                          Hexp=round(Hexp, 3), 
                          E.5=round(E.5,3), 
                          as.data.frame(t(IaList)),
                          File=namelist$File))
    rownames(Iout) <- NULL
    return(final(Iout, result))
  }
}


new.poppr.all <- function(filelist, ...) {
	result <- NULL
	for(a in filelist){
    cat("| File: ",basename(a),"\n")
		result <- rbind(result, new.poppr(a, ...))
	}
	result <- extract.info(result)
	return(result)
}

.new.ia <- function(pop,sample=0,method=1,quiet="minimal",namelist=NULL,missing="ignore",
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
    sampstatIa <- c(range(samp[[1]]), median(samp[[1]]), mean(samp[[1]]), var(samp[[1]]), sd(samp[[1]]))
    names(sampstatIa) <- c("Ia.min", "Ia.max", "Ia.med", "Ia.mean", "Ia.var", "Ia.sd")
    sampstatRd <- c(range(samp[[2]]), median(samp[[2]]), mean(samp[[2]]), var(samp[[2]]), sd(samp[[2]]))
    names(sampstatRd) <- c("Rd.min", "Rd.max", "Rd.med", "Rd.mean", "Rd.var", "Rd.sd")
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
    result <- c(result[1:2], sampstatIa, result[3:4], sampstatRd)
  } 
  return(final(Iout, result))
}

.new.Ia.Rd <- function (pop, missing = NULL) 
{
  vard.vector <- NULL
  numLoci <- length(pop)
  numIsolates <- length(pop[[1]]@ind.names)
  np <- choose(numIsolates, 2)
  if (np < 2) {
    return(as.numeric(c(NaN, NaN)))
  }
  V <- new.pair_diffs(pop, numLoci, np)
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
new.pair_diffs <- function(pop, numLoci, np)
{
  ploid <- ploidy(pop[[1]])
  temp.d.vector <- matrix(nrow = np, ncol = numLoci, data = as.numeric(NA))
  temp.d.vector <- mclapply(pop, function(x) .Call("pairdiffs", x@tab)*(ploid/2))
  return(temp.d.vector)
  d.vector <- colSums(temp.d.vector)
  d2.vector <- colSums(temp.d.vector^2)
  D.vector <- rowSums(temp.d.vector)
  return(list(d.vector = d.vector, d2.vector = d2.vector, D.vector = D.vector))
}

getinds <- function(x){
  # x is a vector indicating the size of loci in a data set (must be continuous).
  # The sum of x is equal to the number of observations.
  to <- cumsum(x)
  from <- c(1, to[-length(to)]+1)
  indices <- rbind(from, to)
  rownames(indices) <- NULL
  return(indices)
}


brian.ia <- function(pop, sample=0, valuereturn = FALSE, method=1, 
                     quiet="minimal", missing="ignore",hist=TRUE){
  METHODS = c("multilocus", "permute alleles", "parametric bootstrap",
      "non-parametric bootstrap")
  namelist <- NULL
  namelist$population <- ifelse(length(levels(pop@pop)) > 1 | 
                                is.null(pop@pop), "Total", pop@pop.names)
  namelist$File <- as.character(pop@call[2])
  popx <- pop
  missing <- toupper(missing)
  type <- pop@type
  if(type=="PA"){
    .Ia.Rd <- .PA.Ia.Rd
  }
  else {
    popx <- seploc(popx)
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
    if (valuereturn == TRUE){
      return(list(index = final(Iout, result), samples = samp))
    }
  }  
  return(final(Iout, result))
}



new.jack.ia <- function(pop, iterations, type = type, divisor = 2){ 
  if(!is.list(pop)){
    if(type=="PA"){
      .Ia.Rd <- .PA.Ia.Rd
    }
  }
  inds <- nInd(pop[[1]])
  half <- round(inds/divisor)
  np <- choose(half, 2)
	funrun <- function(pop, inds = 10, half = 5, np = 10){
	
	  samp <- sample(inds, half)
    comb <- combn(inds, 2) # pairwise vector
    jackvec <- colSums(matrix(comb %in% samp, nrow = 2)) > 1
    V <- jack.Ia.Rd(pop)
    V <- list(d.vector = colSums(V[jackvec, ]), 
              d2.vector = colSums(V[jackvec, ]^2), 
              D.vector = rowSums(V[jackvec, ])
              )
    IarD <- jack.calc(V, np)
	}
	sample.data <- vapply(1:iterations, function(x) funrun(pop, inds = inds, half = half, np = np), c(pi, pi))
	rownames(sample.data) <- c("Ia", "rbarD")
	return(data.frame(t(sample.data)))
}

jack.Ia.Rd <- function (pop){
  vard.vector <- NULL
  numLoci <- length(pop)
  numIsolates <- length(pop[[1]]@ind.names)
  np <- choose(numIsolates, 2)
  if (np < 2) {
    return(as.numeric(c(NaN, NaN)))
  }
  V <- jack.pair_diffs(pop, numLoci, np)
  return(V)
}

jack.calc <- function(V, np){
  varD <- ((sum(V$D.vector^2) - ((sum(V$D.vector))^2)/np))/np
  vard.vector <- ((V$d2.vector - ((V$d.vector^2)/np))/np)
  vardpair.vector <- .Call("pairwise_covar", vard.vector)
  sigVarj <- sum(vard.vector)
  rm(vard.vector)
  Ia <- (varD/sigVarj) - 1
  rbarD <- (varD - sigVarj)/(2 * sum(vardpair.vector))
  return(c(Ia, rbarD))
}

jack.pair_diffs <- function(pop, numLoci, np)
{
  ploid <- ploidy(pop[[1]])
  temp.d.vector <- matrix(nrow = np, ncol = numLoci, data = as.numeric(NA))
  temp.d.vector <- vapply(pop, function(x) .Call("pairdiffs", x@tab)*(ploid/2), 
                          temp.d.vector[, 1])
  return(temp.d.vector)
}

# Function jackcomp is not fully thought out yet. Right now, it is unclear what 
# the application for this could be. The idea behind it is that if you randomly
# sampled half of the individuals from your population, how does that affect the
# index of association? This
#
#

jackcomp <- function(pop, sample = 999, quiet = TRUE, method = 1, divisor = 2){
  inds <- nInd(pop)
  half <- round(inds/divisor)
  cat("Creating Null Distribution...\t")
  null <- brian.ia(pop, sample = sample, valuereturn = TRUE, quiet = quiet, 
    hist = FALSE, method = method)
  cat("Alternative Distribution...\t")
  alt <- new.jack.ia(seploc(pop), sample, type = pop@type, divisor)
  library(reshape)
  mnull <- melt(null$sample)
  malt <- melt(alt)
  mnull$Distribution <- "null"
  malt$Distribution <- "alternative"
  dat <- rbind(mnull, malt)
  cat("Creating Plots\n")
  obs.df <- data.frame(list(variable = names(null$samples), value = null$index[c(1,3)]))
  distplot <- ggplot(dat, aes_string(x = "value", fill = "Distribution")) + geom_histogram(alpha = 0.5, position = "identity") + geom_rug(alpha = 0.5, aes_string(color = "Distribution")) + facet_grid(" ~ variable", scales = "free_x") + theme_classic() + ggtitle(paste("Data:", deparse(substitute(pop)), "\n", inds, "Individuals,", half, "Sampled for Alt.")) + geom_vline(data = obs.df, aes_string(xintercept = "value", group = "variable"), linetype = "dashed")
  #print(distplot)
  return(list(observed = null$index, null_samples = null$samples, alt_samples = alt, plot = distplot))
}


new.bootia <- function(pop, inds){
  print("heyy!")
}

bootia <- function(pop, inds, quiet = TRUE, method = 1){
  pop <- pop[sample(inds, replace = TRUE), ]
  pop@ind.names <- paste("ind", 1:inds)
  return(.ia(pop, quiet = quiet))
}

bootcomp <- function(pop, sample = 999, quiet = TRUE, method = 1){
  inds <- nInd(pop)
  cat("Creating Null Distribution...\t")
  null <- brian.ia(pop, sample = sample, valuereturn = TRUE, quiet = quiet, 
    hist = FALSE, method = method)
  cat("Alternative Distribution...\t")
  alt <- vapply(1:sample, function(x) bootia(pop, inds, quiet, method), c(pi, pi))
  library(reshape)
  mnull <- melt(null$sample)
  malt <- melt(data.frame(t(alt)))
  mnull$Distribution <- "null"
  malt$Distribution <- "alternative"
  dat <- rbind(mnull, malt)
  cat("Creating Plots\n")
  obs.df <- data.frame(list(variable = names(null$samples), value = null$index[c(1,3)]))
  distplot <- ggplot(dat, aes_string(x = "value", fill = "Distribution")) + geom_histogram(alpha = 0.5, position = "identity") + geom_rug(alpha = 0.5, aes_string(color = "Distribution")) + facet_grid(" ~ variable", scales = "free_x") + theme_classic() + ggtitle(paste("Data:", deparse(substitute(pop)), "\n", inds, "Individuals,", inds, "Sampled for Alt.")) + geom_vline(data = obs.df, aes_string(xintercept = "value", group = "variable"), linetype = "dashed")
  #print(distplot)
  return(list(observed = null$index, null_samples = null$samples, alt_samples = data.frame(t(alt)), plot = distplot))
}


total.shuffler <- function(pop, method){
  return(lapply(pop, pop.sampler, method))
}
  


pop.sampler <- function(pop, method=1){
  pop_names <- names(table(pop@pop))
  total_inds <- 1:length(pop@pop)
  repop <- function(x, pop, pop_names, total_inds, method){
    pop_vec <- total_inds[pop@pop %in% pop_names[x]]
    if(all(is.na(pop@tab[pop_vec, ]))){
      pop@tab[pop_vec, ] <<- pop@tab[pop_vec, ]
      return(1)
    }
    cols <- which(colSums(pop[pop_vec, ]@tab, na.rm=TRUE) > 0)
    newpop <- pop
    if(method != 2){
      newpop@tab[pop_vec, -cols] <- 0
    }
    newpop@tab[pop_vec, cols] <- .locus.shuffler(popsub(pop, x), method)@tab
    if(method == 1){
      replacements <- is.na(rowSums(newpop@tab[pop_vec, cols]))
      if(any(replacements == TRUE)){
        #cat(pop_vec,"\t",replacements,"\n")
        newpop@tab[pop_vec[replacements], -cols] <- NA
      }
    }
    pop@tab[pop_vec, ] <<- newpop@tab[pop_vec, ]
    return(0)
  }
  invisible(lapply(1:length(pop_names), repop, pop, pop_names, total_inds, method))
#   nas <- function(x, pop, pop_names, total_inds, method){
#     pop_vec <- total_inds[pop@pop %in% pop_names[x]]
#     if(any(is.na(pop@tab[pop_vec, ]))){
#       newpop <- pop
#       replacena <- function(locus, newpop, pop_vec){
#         locus <- newpop@loc.fac %in% locus
#         truth <- is.na(rowSums(newpop@tab[pop_vec, locus]))
#         if(any(truth == TRUE)){
#           #cat(pop_vec,"\n",as.numeric(locus),"\n", truth, "\n")
#           newpop@tab[pop_vec[truth], locus] <<- NA
#         }
#       }
#       invisible(lapply(names(pop@loc.names), replacena, newpop, pop_vec))
#       pop@tab[pop_vec, ] <<- newpop@tab[pop_vec, ]
#     }
#   }
#   if(method == 1){
#     invisible(lapply(1:length(pop_names), nas, pop, pop_names, total_inds, method))
#   }
  return(pop)
}





################################################################################
#################### Zhian's Functions above ###################################
################################################################################
################################################################################
################################################################################
#################### DO NOT CROSS THIS LINE!!!! ################################
################################################################################
################################################################################
################################################################################
#################### Javier's Functions below ##################################
################################################################################

genoid.bruvo.boot <- function(pop, replen = 1, add = TRUE, loss = TRUE, 
                              sample = 100,  tree = "upgma", showtree = TRUE, 
                              cutoff = NULL,  quiet = FALSE, ...){
  # This attempts to make sure the data is true microsatellite data. It will
  # reject snp and aflp data. 
  if(pop@type != "codom" | all(is.na(unlist(lapply(pop@all.names, as.numeric))))){
    stop("\nThis dataset does not appear to be microsatellite data. Bruvo's Distance can only be applied for true microsatellites.")
  }
  ploid <- ploidy(pop)
  # Bruvo's distance depends on the knowledge of the repeat length. If the user
  # does not provide the repeat length, it can be estimated by the smallest
  # repeat difference greater than 1. This is not a preferred method. 
  if (length(replen) != length(pop@loc.names)){
    replen <- vapply(pop@all.names, function(x) guesslengths(as.numeric(x)), 1)
    #    replen <- rep(replen[1], numLoci)
    warning("\n\nRepeat length vector for loci is not equal to the number of loci represented.\nEstimating repeat lengths from data:\n", immediate.=TRUE)
    cat(replen,"\n\n")
  }
  # This controlls for the user correcting missing data using "mean". 
  if(any(!pop@tab,10 %in% c(0,((1:ploid)/ploid), 1, NA))){
    pop@tab[!pop@tab %in% c(0,((1:ploid)/ploid), 1, NA)] <- NA
  }
  # Converting the genind object into a matrix with each allele separated by "/"
  bar <- as.matrix(genind2df(pop, sep="/", usepop = FALSE))
  # The bruvo algorithm will ignore missing data, coded as 0.
  bar[bar %in% c("", NA)] <- paste(rep(0, ploid), collapse="/")
  stopifnot(require(phangorn))
  # Steps: Create initial tree and then use boot.phylo to perform bootstrap
  # analysis, and then place the support labels on the tree.
  if(tree == "upgma"){
    newfunk <- match.fun(upgma)
  }
  else if(tree == "nj"){
    newfunk <- match.fun(nj)
  }
  tre <- newfunk(phylo.bruvo.dist(bar, replen=replen, ploid=ploid, add = add, loss = loss))
  if (any (tre$edge.length < 0)){
    tre$edge.length[tre$edge.length < 0] <- 0
  }
  if(quiet == FALSE){
    cat("\nBootstrapping... (note: calculation of node labels can take a while even after the progress bar is full)\n\n")
  }
  bp <- boot.phylo(tre, bar, FUN = function (x) newfunk(phylo.bruvo.dist(x, replen = replen, ploid = ploid, add = add, loss = loss)), B = sample, quiet = quiet, ...)
  tre$node.labels <- round(((bp / sample)*100))
#  if (!is.null(cutoff)){
#    if (cutoff < 1 | cutoff > 100){
#      cat("Cutoff value must be between 0 and 100.\n")
#      cutoff<- as.numeric(readline(prompt = "Choose a new cutoff value between 0 and 100:\n"))
#    }
    tre$node.labels[tre$node.labels < cutoff] <- NA
#  }
  tre$tip.label <- pop@ind.names
#  if(showtree == TRUE){
#    plot(tre, show.node.label=TRUE)
#  }
#  if(tree=="upgma"){
#    axisPhylo(3)
#  }
  return(tre)
}









javier<-function(x){
  cat ("http://www.youtube.com/watch?v=1-ctsxVXvO0")
}
