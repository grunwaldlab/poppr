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

# Keeping this for compatability
new.read.genalex <- function(genalex, ploidy=2, geo=FALSE, region=FALSE){
  return(read.genalex(genalex, ploidy, geo, region))
}

#==============================================================================#
# Create a population hierarchy in a genind object if one is not there.
#==============================================================================#
genind_hierarchy <- function(x, df = NULL, dfname = "population_hierarchy"){
  if (is.null(df)){
    df <- data.frame(Pop = pop(x))
  } 
  other(x)[[dfname]] <- df
  return(x)
}


#==============================================================================#
# Idea, create a table representing missing data in the data set. This will 
# show missing data across loci and populations. 
# Args:
# 
#   x - a genind object
#   percent - should the data be presented as a percent or raw value? Note that 
#             this represents percentage with respect to the population. 
#   plot - should this be plotted
#   df - return a data frame in long format
#   returnplot - return the plot as the second element of a list
#   low - color representing no missing data
#   high - color representing missing data
#   plotlab - print the missing values on the plot.
#   scaled - if percent = TRUE, then this will make the scale from 0 to 100, if 
#   this is not turned on (default), then the highest missing value will be the max.
#==============================================================================#

# tabulate the amount of missing data per locus. 
number_missing_locus <- function(x, divisor){
  missing_result <- colSums(1 - propTyped(x, by = "both"))
  return(missing_result/divisor)
}

# tabulate the amount of missing data per genotype. 
number_missing_geno <- function(x, divisor){
  missing_result <- rowSums(1 - propTyped(x, by = "both"))
  return(missing_result/divisor)
}

# main function



#==============================================================================#
# This function will simply return a table denoting private alleles within
# populations. The rows will indicate populations with private alleles and the
# columns will be the specific alleles. 
#==============================================================================#
private_alleles <- function(gid){
  if (!is.genind(gid) & !is.genpop(gid)) stop(paste(gid, "is not a genind or genpop object."))
    if (is.genind(gid) & !is.null(pop(gid)) | is.genpop(gid) & nrow(gid@tab) > 1){
      if (is.genind(gid)){
          gid.pop <- truenames(genind2genpop(gid, quiet = TRUE))
        } else {
          gid.pop <- truenames(gid)
        }
        privates <- gid.pop[, colSums(ifelse(gid.pop > 0, 1, 0), na.rm = TRUE) < 2]
        privates <- privates[rowSums(privates) > 0, ]
        return(privates)
      } else {
    stop("There are no populations detected")
  }
}

#==============================================================================#
# The following two functions serve as summary functions across loci. It was
# noted that there are no functions that really give explicit summaries
# accross loci, so these were written to do so. 
#
# The function locus_table_pegas is the internal workhorse. It will process a 
# summary.loci object into a nice table utilizing the various diversity indices 
# provided by vegan. Note that it has a catch which will remove all allele names
# that are any amount of zeroes and nothing else. The reason for this being that
# these alleles actually represent missing data. 
# The wrapper for this is locus_table, which will take in a genind object and 
# send it into locus_table_pegas. 
# Note: lev argument has only the options of "allele" or "genotype"
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

locus_table <- function(x, index ="simpson", lev = "allele", population = "ALL", information = TRUE){
  INDICES <- c("shannon", "simpson", "invsimpson")
  index   <- match.arg(index, INDICES)
  x       <- popsub(x, population, drop = FALSE)
  x.loc   <- summary(as.loci(x))
  outmat  <- vapply(x.loc, locus_table_pegas, numeric(4), index, lev, x@type)
  loci    <- colnames(outmat)
  divs    <- rownames(outmat)
  outmat  <- t(outmat)
  dimlist <- list(`locus` = loci, `summary` = divs)
  attr(outmat, "dimnames") <- dimlist
  if (information){
    if (index == "simpson"){
      msg <- "Simpson index"
    } else if (index == "shannon"){
      msg <- "Shannon-Wiener index"
    } else {
      msg <- "Stoddard and Taylor index"
    }
    cat("\n", divs[1], "= Number of observed", paste0(divs[1], "s"))
    cat("\n", divs[2], "=", msg)
    cat("\n", divs[3], "= Nei's 1978 expected heterozygosity\n")
    cat("------------------------------------------\n")
  }
  return(outmat)
}


pair_ia <- function(pop){

  if(pop@type == "codom"){
    pop_loci <- seploc(pop)
    loci_pairs <- combn(pop@loc.names, 2)
    pair_ia_vector <- apply(loci_pairs, 2, function(x) .Ia.Rd(pop_loci[x]))
    colnames(pair_ia_vector) <- apply(loci_pairs, 2, paste, collapse = ":")
  } else {
    loci_pairs <- combn(1:nLoc(pop), 2)
    pair_ia_vector <- apply(loci_pairs, 2, function(x) .PA.Ia.Rd(pop[, x], missing = "ignore"))
    colnames(pair_ia_vector) <- apply(combn(pop@loc.names, 2), 2, paste, collapse = ":")
  }
  rownames(pair_ia_vector) <- c("Ia", "rbarD")    
  return(pair_ia_vector)
}


poppr_pair_ia <- function(pop){
  if(is.null(pop(pop))){
    return(pair_ia(pop))
  }
  pops <- seppop(pop, drop = FALSE)
  loci_pairs <- choose(nLoc(pop), 2)
  res_mat <- matrix(0.5, 2, loci_pairs)
  pops_array <- vapply(pops, pair_ia, res_mat)
  return(pops_array)
}


testing_funk <- function(){
  cat("This test worked...maybe.\n")
}





new.poppr <- function(pop,total=TRUE,sublist=c("ALL"),blacklist=c(NULL), sample=0,
  method=1,missing="ignore", quiet=FALSE,clonecorrect=FALSE,hier=c(1),dfname="hier",
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
    # Shannon-Weiner vegan::diversity index.
    H <- vegan::diversity(pop.mat)
    # E_1, Pielou's evenness.
    # J <- H / log(rowSums(pop.mat > 0))
    # inverse Simpson's index aka Stoddard and Taylor: 1/lambda
    G <- vegan::diversity(pop.mat, "inv")
    Hexp <- (N.vec/(N.vec-1))*vegan::diversity(pop.mat, "simp")
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
    # Shannon-Weiner vegan::diversity index.
    H <- vegan::diversity(pop.mat)
    # E_1, Pielou's evenness.
    # J <- H / log(rowSums(pop.mat > 0))
    # inverse Simpson's index aka Stoddard and Taylor: 1/lambda
    G <- vegan::diversity(pop.mat, "inv")
    Hexp <- (N.vec/(N.vec-1))*vegan::diversity(pop.mat, "simp")
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

.new.ia <- function(pop,sample=0,method=1,quiet=FALSE,namelist=NULL,missing="ignore",
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
  if (!quiet){
    cat("|", namelist$population ,"\n")
  }
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
      poppr.plot(samp, observed=IarD, pop=namelist$population,
                        file=namelist$File, pval=p.val, N=nrow(pop@tab))
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
  temp.d.vector <- vapply(pop, function(x) .Call("pairdiffs", x@tab)*(ploid/2), numeric(np))
  return(temp.d.vector)
  d.vector <- colSums(temp.d.vector)
  d2.vector <- colSums(temp.d.vector^2)
  D.vector <- rowSums(temp.d.vector)
  return(list(d.vector = d.vector, d2.vector = d2.vector, D.vector = D.vector))
}

getinds <- function(x){
  # x is a vector indicating the size of loci in a data set (must be continuous).
  # The sum of x is equal to the number of observations.
  indices <- matrix(ncol = 2, nrow = length(x))
  to   <- cumsum(x)
  from <- c(1, to[-length(to)]+1)
  indices[] <- as.integer(c(from, to))
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
      poppr.plot(samp, observed=IarD, pop=namelist$population,
                        file=namelist$File, pval=p.val, N=nrow(pop@tab))
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



jack.ia <- function(pop, iterations, type = type, divisor = 2){ 
  if(!is.list(pop)){
    inds <- nInd(pop)
  }
  else{
    inds <- nInd(pop[[1]])
  }
  half <- round(inds/divisor)
  np <- choose(half, 2)    
  V <- jack.Ia.Rd(pop)
	funrun <- function(V, inds = 10, half = 5, np = 10){
	
	  samp <- sample(inds, half)
    comb <- combn(inds, 2) # pairwise vector
    jackvec <- colSums(matrix(comb %in% samp, nrow = 2)) > 1

    V <- list(d.vector = colSums(V[jackvec, ]), 
              d2.vector = colSums(V[jackvec, ]^2), 
              D.vector = rowSums(V[jackvec, ])
              )
    IarD <- jack.calc(V, np)
	}
	sample.data <- vapply(1:iterations, function(x) funrun(V, inds = inds, half = half, np = np), c(pi, pi))
	rownames(sample.data) <- c("Ia", "rbarD")
	return(data.frame(t(sample.data)))
}

jack.Ia.Rd <- function (pop){

  vard.vector <- NULL
  if(!is.list(pop)){
    numLoci <- nLoc(pop)
    numIsolates <- nInd(pop)
  } else {
    numLoci <- length(pop)
    numIsolates <- length(pop[[1]]@ind.names)
  }
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
  temp.d.vector <- matrix(nrow = np, ncol = numLoci, data = as.numeric(NA))
  if(!is.list(pop)){
    ploid <- 2
    temp.d.vector <- vapply(seq(numLoci), 
                      function(x) as.vector(dist(pop@tab[,x])), 
                      temp.d.vector[,1])
    # checking for missing data and imputing the comparison to zero.
    if(any(is.na(temp.d.vector))){
      temp.d.vector[which(is.na(temp.d.vector))] <- 0
    }
  } 
  else {
    ploid <- ploidy(pop[[1]])
    temp.d.vector <- vapply(pop, function(x) .Call("pairdiffs", x@tab)*(ploid/2), 
                            temp.d.vector[, 1])
  }
  return(temp.d.vector)
}

# Function jackcomp is not fully thought out yet. Right now, it is unclear what 
# the application for this could be. The idea behind it is that if you randomly
# sampled half of the individuals from your population, how does that affect the
# index of association? This
#
#

jackcomp <- function(pop, sample = 999, quiet = TRUE, method = 1, divisor = 1.587302){
  inds <- nInd(pop)
  half <- round(inds/divisor)
  cat("Creating Null Distribution...\t")
  null <- brian.ia(pop, sample = sample, valuereturn = TRUE, quiet = quiet, 
    hist = FALSE, method = method)
  cat("Alternative Distribution...\t")
  if(pop@type == "codom"){
    popx <- seploc(pop)
  }
  else{
    popx <- pop
  }
  alt <- jack.ia(popx, sample, type = pop@type, divisor)
  # library(reshape)
  mnull <- melt(null$sample)
  malt <- melt(alt)
  mnull$Distribution <- "null"
  malt$Distribution <- "alternative"
  dat <- rbind(mnull, malt)
  cat("Creating Plots\n")
  obs.df <- data.frame(list(variable = names(null$samples), value = null$index[c(1,3)]))
  distplot <- jackbootplot(dat, obs.df) + 
              ggtitle(paste("Data:", deparse(substitute(pop)), "\n", inds, 
                            "Individuals,", half, "Sampled for Alt.")
                     ) +
              theme(strip.text = element_text(size = rel(2)))
  #print(distplot)
  return(list(observed = null$index, null_samples = null$samples, alt_samples = alt, plot = distplot))
}

jackbootplot <- function(df, obs.df){
  Indexfac <- c("I[A]","bar(r)[d]")
  levels(df$variable) <- Indexfac
  obs.df$variable <- Indexfac
  dfIanull <- df[df$variable == Indexfac[1] & df$Distribution == "null", ]
  dfrbarDnull <- df[df$variable == Indexfac[2] & df$Distribution == "null", ]
  dfIaalt <- df[df$variable == Indexfac[1] & df$Distribution != "null", ]
  dfrbarDalt <- df[df$variable == Indexfac[2] & df$Distribution != "null", ]
  distplot <- ggplot(df, aes_string(x = "value", fill = "Distribution")) +
              geom_histogram(alpha = 0.5, position = "identity", 
                             data = dfIanull)+#,
                             # binwidth=diff(range(dfIanull$value))/30) +
              geom_histogram(alpha = 0.5, position = "identity", 
                             data = dfrbarDnull)+#,
                             # binwidth=diff(range(dfrbarDnull$value))/30) +
              geom_histogram(alpha = 0.5, position = "identity", 
                             data = dfIaalt)+#,
                             # binwidth=diff(range(dfIaalt$value))/30) +
              geom_histogram(alpha = 0.5, position = "identity", 
                             data = dfrbarDalt)+#,
                             # binwidth=diff(range(dfrbarDalt$value))/30) +

              geom_rug(alpha = 0.5, aes_string(color = "Distribution")) +
              facet_grid(" ~ variable", scales = "free_x", labeller = label_parsed) + theme_classic() +
              geom_vline(data = obs.df, aes_string(xintercept = "value", 
                                                   group = "variable"), 
                         linetype = "dashed") 
  return(distplot)      
}

#==============================================================================#
# This will perform a bootstrap analysis using the resulting distance matrix
# produced from pair_diffs(). Since the matrix is values of pairwise comparisons,
# the method of bootstrapping involves a few steps in order to resort the matrix.
#==============================================================================#
old.bootia <- function(pop, inds, iterations, comb){
  V <- jack.Ia.Rd(pop) # calculate matrix. rows = inds, cols = loci
  
  bootit <- function(V, inds, comb){
    #set.seed(9000)
    sam <- sample(inds, replace = TRUE) 
    #print(head(sam))
    
    
    # Indicating which columns contain ONLY values included in `sam`.
    uneak <- colSums(matrix(comb %in% sam, nrow = 2)) > 1
    samcomb <- comb[, uneak, drop = FALSE]
    
    comsam <- combn(sam, 2)
    compast <- vapply(1:ncol(comsam), function(x) paste(sort(comsam[, x]), collapse = ""), "a")
    sampast <- vapply(1:ncol(samcomb), function(x) paste(sort(samcomb[, x]), collapse = ""), "a")
    inds <- vapply(sampast, function(x) sum(compast %in% x), 1)
    samdump <- rep(as.numeric(names(inds)), inds)  
    
    # Now we have to subset the vector. Adding the duplicated comparisons will
    # not complete the entire data set, so the rest has to be filled in with 
    # zeroes, since that's the value one would obtain from comparisons to self.
    V[1:length(samdump), ] <- V[samdump, ]
    V[-(1:length(samdump)), ] <- 0
    
    # Calculate necessary vectors and send them to be calculated into Ia and rd.
    V <- list(d.vector = colSums(V), 
          d2.vector = colSums(V^2), 
          D.vector = rowSums(V)
          )
    return(jack.calc(V, ncol(comb)))
  }
  comb <- combn(inds, 2) # matrix of pairwise indexes.
  
  # Labelling the columns of this matrix provides tractability. When we subset
  # it, we will know where those subsetted columns came from. 
  colnames(comb) <- 1:ncol(comb)
	sample.data <- vapply(1:iterations, function(x) bootit(V, inds, comb), c(pi, pi))
	rownames(sample.data) <- c("Ia", "rbarD")
	return(data.frame(t(sample.data)))
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
  #if(pop@type == "codom"){
  #  popx <- seploc(pop)
  #}
  #else{
  #  popx <- pop
  #}
  #alt <- bootia(popx, inds, sample)
  alt <- data.frame(t(alt))
  # library(reshape)
  mnull <- melt(null$sample)
  malt <- melt(alt)
  mnull$Distribution <- "null"
  malt$Distribution <- "alternative"
  dat <- rbind(mnull, malt)
  cat("Creating Plots\n")
  obs.df <- data.frame(list(variable = names(null$samples), value = null$index[c(1,3)]))
  distplot <-distplot <- jackbootplot(dat, obs.df) + 
              ggtitle(paste("Data:", deparse(substitute(pop)), "\n", inds, 
                            "Individuals, Method: Bootstrap")
                     ) +
              theme(strip.text = element_text(size = rel(2)))
  #print(distplot)
  return(list(observed = null$index, null_samples = null$samples, alt_samples = alt, plot = distplot))
}

#' @importFrom reshape2 melt

jackbootcomp <- function(pop, sample = 999, quiet = TRUE, method = 1, divisor = 1.587302){
  inds <- nInd(pop)
  cat("Creating Null Distribution...\t")
  null <- brian.ia(pop, sample = sample, valuereturn = TRUE, quiet = quiet, 
    hist = FALSE, method = method)
  if(pop@type == "codom"){
    popx <- seploc(pop)
  }
  else{
    popx <- pop
  }
  cat("Alternative Distribution (Bootstrap)...\t")
  altboot <- vapply(1:sample, function(x) bootia(pop, inds, quiet, method), c(pi, pi))
  altboot <- data.frame(t(altboot))
  #altboot <- bootia(popx, inds, sample)
  cat("Alternative Distribution (Jack Knife)...\t")
  half <- round(inds/divisor)
  altjack <- jack.ia(popx, sample, type = pop@type, divisor)  
  # library(reshape)
  mnull <- melt(null$sample)
  bmalt <- melt(altboot)
  jmalt <- melt(altjack)
  mnull$Distribution <- "null"
  bmalt$Distribution <- "bootstrap"
  jmalt$Distribution <- "jack knife"
  dat <- rbind(mnull, bmalt, jmalt)
  cat("Creating Plots\n")
  obs.df <- data.frame(list(variable = names(null$samples), value = null$index[c(1,3)]))
  distplot <-distplot <- jackbootplot(dat, obs.df) + 
              ggtitle(paste("Data:", deparse(substitute(pop)), "\n", inds, 
                            "Individuals,", half, "Sampled for Jack knife")
                     ) +
              theme(strip.text = element_text(size = rel(2)))
  #print(distplot)
  return(list(observed = null$index, null_samples = null$samples, 
              alt_samples_boot = altboot, alt_samples_jack = altjack, 
              plot = distplot)
        )
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
  # stopifnot(require(phangorn))
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
