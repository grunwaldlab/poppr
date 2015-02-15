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


read_allele_columns <- function(x, ploidy = 2){
  clm <- ncol(x)
  if (length(ploidy) == 1){
    if (clm %% ploidy == 0){
      ploidy <- rep(ploidy, clm/ploidy)
    } else {
      stop("ploidy specified does not match number of columns.")
    }
  }
  loci  <- first_allele_col(ploidy)
  x2    <- x[, loci]
  lapply(1:length(loci), function(a) x2[, a] <<-
           pop_combiner(x, hier = loci[a]:(loci[a]+ploidy[a]-1), sep = "/"))
  
  x2 <- data.frame(lapply(x2, add_zeroes,  ploidy = max(ploidy)))
  x2[apply(x2, 2, function(na) grepl("NA", na))] <- NA
  return(x2)
}

first_allele_col <- function(ploidies) cumsum(ploidies) - ploidies + 1

# Function that will add extra zeroes onto genotypes that are deficient in the 
# right number of alleles.
add_zeroes <- function(x, ploidy = 3){
  extras <- ploidy - vapply(strsplit(x, "/"), length, 69)
  vapply(1:length(extras), zero_adder, character(1), extras, x)
}

zero_adder <- function(index, extras, df){
  if(extras[index] > 0){
    return(paste(c(rep(0, extras[index]), df[index]), collapse = "/"))
  } else {
    return(df[index])
  }
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

# gen2polysat <- function(gen, newploidy = gen@ploidy){
#   if (!require(polysat)){
#     stop("User needs polysat installed")
#   }
#   gen   <- recode_polyploids(gen, newploidy)
#   gendf <- genind2df(gen, sep = "/", usepop = FALSE)
#   gendf <- lapply(gendf, strsplit, "/")
#   gendf <- lapply(gendf, lapply, as.numeric)
#   ambig <- new("genambig", samples = indNames(gen), loci = locNames(gen))
#   lapply(names(gendf), function(x) Genotypes(ambig, loci = x) <<- gendf[[x]])
#   return(ambig)
# }


new_graph_pops <- function(graph, dat, color){
  cmlg       <- unique(mlg.vector(dat))
  mlg.cp     <- mlg.crosspop(dat, mlgsub = 1:length(cmlg), quiet = TRUE)
  mlg.cp     <- mlg.cp[rank(cmlg)]
  mlg.color  <- lapply(mlg.cp, function(x) color[dat@pop.names %in% names(x)])
  V(graph)$pie <- mlg.cp
  V(graph)$pie.color <- mlg.color
  return(graph)
}


new.poppr.plot <- function(sample, pval = c("0.05", "0.05"), pop="pop", 
                           observed = c(Ia = 0, rbarD = 0), file="file",
                           N=NA, index = "rbarD"){
  srange <- range(sample[[index]])
  labs <- c(Ia = "I[A]", rbarD = "bar(r)[d]")
  binw <- diff(srange/30)
  if (all(observed[index] > srange)){
    xval <- observed[index]
    just <- 1.1
  } else {
    maxobsmin <- c(srange[1], observed[index], srange[2])
    xjust <- which.max(abs(diff(maxobsmin)))
    xval <- srange[xjust]
    just <- ifelse(xjust < 2, 0, 1)     
  }
  obslab <- paste0("observed:\n", labs[index], ":", signif(observed[index], 2))
  thePlot <- ggplot(sample, aes_string(x = index))
  thePlot <- thePlot + 
    geom_histogram(binwidth = binw, position = "identity") +
    geom_rug() + 
    geom_vline(xintercept = observed[index], color = "blue", linetype = 2) +
    annotate(geom = "text", x = xval, y = Inf, color = "blue",
             label = obslab, vjust = 2, hjust = just, parse = TRUE)
  
  if (index == "rbarD"){
    thePlot <- thePlot + xlab(expression(paste(bar(r)[d])))
  } else if (index == "Ia"){
    thePlot <- thePlot + xlab(expression(paste(I[A])))
  }
  return(thePlot)
}

#==============================================================================#
# An attempt at making the shuffling schemes handle polyploid data better.
#==============================================================================#

new.permut.shuff <- function(mat, ploidy = 2L, ploidvec = NULL, zerocol = 1L){
  bucket     <- round(colSums(mat, na.rm = TRUE)*ploidy)[-zerocol]
  bucketlen <- (1:ncol(mat))[-zerocol]
  bucketlist <- as.integer(sample(rep(bucketlen, bucket)))
  nas <- !is.na(ploidvec)
  ploidvec[nas] <- sample(ploidvec[nas])
  mat <- .Call("new_permute_shuff", mat, bucketlist - 1L, 1/ploidy, ploidy, 
                     ploidvec, zerocol - 1L, PACKAGE = "poppr")
  # mat <- list(mat = mat, bucket = bucket, bucketlist = bucketlist - 1L, 
  #             allele = 1/ploidy, ploidy = ploidy, ploidvec = ploidvec, 
  #             zerocol = zerocol - 1L)
  return(mat)
}


new.all.shuff <- function(pop, type=type, method=1, ploidvec= NULL, zerocol = 1L){
  METHODS = c("permute alleles", "parametric bootstrap",
              "non-parametric bootstrap", "multilocus")
  if(type=="PA"){
    if(method == 1 | method == 4){
      pop@tab <- vapply(1:ncol(pop@tab),
                        function(x) sample(pop@tab[, x]), pop@tab[, 1])
    }
    else if(method == 2){
      paramboot <- function(x){
        one <- mean(pop@tab[, x], na.rm=TRUE)
        zero <- 1-one
        return(sample(c(1,0), length(pop@tab[, x]), prob=c(one, zero), replace=TRUE))
      }
      pop@tab <- vapply(1:ncol(pop@tab), paramboot, pop@tab[, 1])
    }
    else if(method == 3){
      pop@tab <- vapply(1:ncol(pop@tab),
                        function(x) sample(pop@tab[, x], replace=TRUE), pop@tab[, 1])
    }
  } 
  else {
    pop <- lapply(1:length(pop), function(x){
      new.locus.shuffler(pop[[x]], method = method, ploidvec[, x], zerocol)
    })
  }
  return(pop)
}

new.locus.shuffler <- function(pop, method=1, ploidvec = NULL, zerocol = 1L){
  ploid <- ploidy(pop)
  if (method == 1){
    pop@tab <- new.permut.shuff(pop@tab, ploidy = ploid, ploidvec, zerocol)    
  } else if (method == 2){

    uploids <- unique(ploidvec)
    weights <- colMeans(pop@tab[, -zerocol, drop = FALSE], na.rm = TRUE)
    
    if (length(uploids) > 1){
      ploidprobs <- table(ploidvec)
      newploids <- sample(unique(ploidvec), length(ploidvec), 
                          prob = ploidprobs/sum(ploidprobs), replace = TRUE)
      pop@tab[, zerocol] <- newploids/ploid
      lapply(unique(newploids), function(x){
        pop@tab[newploids == x, -zerocol] <<- t(rmultinom(sum(newploids == x), 
                                                size = x, prob = weights))/ploid
        
      })
    
    } else {
      pop@tab[, -zerocol] <- t(rmultinom(length(ploidvec), size = uploids, 
                               prob = weights))/ploid
    }
  } else if (method == 3){
    uploids <- unique(ploidvec)
    if (length(uploids) > 1){
      newploids <- sample(unique(ploidvec), length(ploidvec), replace = TRUE)
      pop@tab[, zerocol] <- newploids/ploid
      lapply(unique(newploids), function(x){
        pop@tab[newploids == x, -zerocol] <<- t(rmultinom(sum(newploids == x), 
                                                size = x, prob = rep(1, ploid)))/ploid
        
      })
    } else {
      pop@tab[, -zerocol] <- t(rmultinom(length(ploidvec), size = uploids, 
                               prob = rep(1, ploid)))/ploid
    }
  } else if (method == 4){
    pop@tab <- pop@tab[sample(nInd(pop)), ]
  }
  return(pop)
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
  pops       <- seppop(pop, drop = FALSE)
  loci_pairs <- choose(nLoc(pop), 2)
  res_mat    <- matrix(0.5, 2, loci_pairs)
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
	x <- process_file(pop, missing=missing, clonecorrect=clonecorrect, hier=hier, 
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
  vardpair.vector <- .Call("pairwise_covar", vard.vector, PACKAGE = "poppr")
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
  temp.d.vector <- vapply(pop, function(x) .Call("pairdiffs", x@tab, PACKAGE = "poppr")*(ploid/2), numeric(np))
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

#==============================================================================#
# bootjack is a function that will calculate values of the index of association
# after sampling with or without replacement. The purpose of this function is
# to help determine distributions of I_A under non-random mating by creating 
# exact copies of samples by sampling with replacement and also to derive a
# distribution of the data by sampling a subset of the data (63% by default)
# without replacement. 
#
# Since the data itself is not being changed, we can use the distances observed.
# These distances are calculated per-locus and presented in a matrix (V) with the
# number of columns equal to the number of loci and the number of rows equal to
# choose(N, 2) where N is the number of samples in the data set. Sampling has
# been optimized with this strategy by first creating a lower triangle distance
# matrix containing indices from 1 to choose(N, 2). This is then converted to a
# square matrix (mat) and subset after sampling with or without replacement. The
# remaining indices in the distance matrix is used to subset the distance-per-
# locus matrix (V). Calculations are then performed on this matrix. It must be
# noted that for sampling with replacement, duplicated indices must be
# supplemented with rows of zeroes to indicate no distance.  
#==============================================================================#
bootjack <- function(gid, iterations = 999, half = NULL, progbar = NULL){
  if (!is.list(gid)){
    N <- nInd(gid)
    numLoci <- nLoc(gid)
  } else {
    N <- nInd(gid[[1]])
    numLoci <- length(gid)
  }
  # Step 1: make a distance matrix defining the indices for the pairwise 
  # distance matrix of loci.
  np  <- choose(N, 2)
  dis <- 1:np
  dis <- make_attributes(dis, N, 1:N, "dist", call("dist"))
  mat <- as.matrix(dis)
  # Step 2: calculate the pairwise distances for each locus. 
  V   <- jack.pair_diffs(gid, numLoci, np)
  if (!is.null(half)){
    np <- choose(half, 2)
    sample.data <- vapply(1:iterations, jackrun, numeric(2), V, mat, N, half, 
                          np, progbar, iterations)
  } else {
    sample.data <- vapply(1:iterations, bootrun, numeric(2), V, mat, N, np, 
                          progbar, iterations)
  }
  rownames(sample.data) <- c("Ia", "rbarD")
  return(data.frame(t(sample.data)))
}

#------------------------------------------------------------------------------#
# Variables for bootrun and jackrun
#
# count      = the current iteration
# V          = the matrix of pairwise distances per locus
# mat        = square matrix containing indices for V
# N          = number of observations in data set
# half       = the number of observations to be sampled without replacement
# np         = choose(N, 2)
# progbar    = a progress bar object
# iterations = the number of total samples
#------------------------------------------------------------------------------#
bootrun <- function(count, V, mat, N, np, progbar, iterations){
  if (!is.null(progbar)){
    setTxtProgressBar(progbar, count/iterations)
  }
  # Step 3: sample individuals with replacement and subset distance matrix
  # from step 1.
  inds    <- sample(N, replace = TRUE)
  newInds <- as.vector(as.dist(mat[inds, inds])) 
  
  # Step 4: Find the number of zeroes in from the distance matrix. This will be 
  # the number of duplicated genotypes.
  zeroes <- length(inds == 0)

  # Step 5: subset the pairwise locus matrix.
  newV <- V[newInds, ]
  V    <- list(d.vector  = colSums(newV), 
               d2.vector = colSums(newV * newV), 
               D.vector  = c(rowSums(newV), rep(0, zeroes))
              )
  # Step 6: Calculate the index of association.
  return(jack.calc(V, np))
}

jackrun <- function(count, V, mat, N, half, np, progbar, iterations){
  if (!is.null(progbar)){
    setTxtProgressBar(progbar, count/iterations)
  }

  inds    <- sample(N, half)
  newInds <- as.vector(as.dist(mat[inds, inds])) 

  newV <- V[newInds, ]
  V    <- list(d.vector  = colSums(newV), 
               d2.vector = colSums(newV * newV), 
               D.vector  = rowSums(newV)
              )
  return(jack.calc(V, np))
}

jack.calc <- function(V, np){
  varD <- ((sum(V$D.vector^2) - ((sum(V$D.vector))^2)/np))/np
  vard.vector <- ((V$d2.vector - ((V$d.vector^2)/np))/np)
  vardpair.vector <- .Call("pairwise_covar", vard.vector, PACKAGE = "poppr")
  sigVarj <- sum(vard.vector)
  rm(vard.vector)
  Ia <- (varD/sigVarj) - 1
  rbarD <- (varD - sigVarj)/(2 * sum(vardpair.vector))
  return(c(Ia, rbarD))
}

jack.pair_diffs <- function(pop, numLoci, np){
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
    temp.d.vector <- vapply(pop, function(x) .Call("pairdiffs", x@tab, PACKAGE = "poppr")*(ploid/2), 
                            temp.d.vector[, 1])
  }
  return(temp.d.vector)
}


jackbootplot <- function(df, obs.df, index = c("rd", "ia", "both"), obsline = TRUE){
  ARGS                <- c("rd", "ia", "both")
  index               <- match.arg(index, ARGS)
  Indexfac            <- c("I[A]","bar(r)[d]")
  levels(df$variable) <- Indexfac
  obs.df$variable     <- Indexfac

  if (index == "rd"){
    df     <- df[df$variable == "bar(r)[d]", ]
    obs.df <- obs.df[obs.df$variable == "bar(r)[d]", ]
  } else if (index == "ia"){
    df     <- df[df$variable == "I[A]", ]
    obs.df <- obs.df[obs.df$variable == "I[A]", ]
  }

  distplot <- ggplot(df, aes_string(y = "value", x = "Distribution")) +
                geom_violin(aes_string(fill = "Distribution")) + 
                geom_boxplot(alpha = 0.25, width = 0.125)
  if (index == "both"){
    distplot <- distplot + facet_grid("variable ~ .", scales = "free", 
                                      labeller = label_parsed) 
  } 
  if (obsline){
    aesth    <- aes_string(yintercept = "value", group = "variable") 
    distplot <- distplot + geom_hline(aesth, data = obs.df, linetype = 2)
  }
  return(distplot)      
}

#' @importFrom reshape2 melt

jackbootcomp <- function(pop, sample = 999, quiet = FALSE, method = 1, 
                         jack = 0.63, plotindex = "rd"){
  ARGS <- c("rd", "ia", "both")
  plotindex <- match.arg(plotindex, ARGS)
  
  inds <- nInd(pop)
  message("Null Distribution...")

  null <- ia(pop, sample = sample, valuereturn = TRUE, quiet = quiet, 
             hist = FALSE, method = method)
  if(pop@type == "codom"){
    popx <- seploc(pop)
  } else {
    popx <- pop
  }

  message("Alternative Distribution (Bootstrap)...")
  if(!quiet){
    bootbar <- txtProgressBar(style = 3)
  } else {
    bootbar <- NULL
  }
  jackbar <- bootbar

  altboot <- bootjack(popx, sample, half = NULL, bootbar)

  message("Alternative Distribution (Jack Knife)...")
  half    <- round(inds*jack)
  altjack <- bootjack(popx, sample, half, jackbar)  
  datlist <- list(null = null$sample, bootstrap = altboot, `jack knife` = altjack)
  dat           <- melt(datlist, measure.vars = c("Ia", "rbarD"))
  names(dat)[3] <- "Distribution"
  message("Creating Plots")
  obs.df   <- data.frame(list(variable = names(null$samples), 
                              value = null$index[c(1,3)]))
  distplot <- jackbootplot(dat, obs.df, plotindex) + 
              ggtitle(paste("Data:", deparse(substitute(pop)), "\n", inds, 
                            "Individuals,", half, "Sampled for Jack knife")
                     )
  if (plotindex == "both"){
    distplot <- distplot + theme(strip.text = element_text(size = rel(2)))    
  } else if (plotindex == "rd"){
    distplot <- distplot + ylab(expression(bar(r)[d]))
  } else {
    distplot <- distplot + ylab(expression(I[A]))
  }
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



javier<-function(x){
  cat ("http://www.youtube.com/watch?v=1-ctsxVXvO0")
}
