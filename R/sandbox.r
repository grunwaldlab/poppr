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








testing_funk <- function(){
  cat("This test worked...maybe.\n")
}


# Randomly inserts NA's into alleles. So far this is one allele per locus. I could
# possibly write a function that would simulate systematic missing data later.
# This function would take a list of loci from a single population and return it
# with missing data.
missing_randomizer <- function(pop){
  loci_pos <- sample(length(pop), 1)
  allele_pos <- sample(nrow(pop[[1]]@tab), 1)
  pop[[loci_pos]]@tab[allele_pos,] <- NA
  return(pop)
}
# This repeats missing randomizer to create a dataset with any number of missing
# datapoints.
sprinkle_missing <- function(pop, reps=1){
  vapply(1:reps, function(x) pop <<- missing_randomizer(pop), pop)
  print(lapply(pop, function(x) which(is.na(x@tab))))
  return(pop)
}


# This is their pairwisefst function. It demonstrates what to do with separate
# population structures.
new_pairwise.fst <- function (x, pop = NULL, res.type = c("dist", "matrix"), truenames = TRUE) 
{
    if (!is.genind(x)) 
        stop("x is not a valid genind object")
    #``````````````````````````````````````````````````````````````````````````#
    # Replacing the population with another selected population vector, allowing
    # for popsubulations. 
    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,#    
    if (!is.null(pop)) {
        pop(x) <- pop
    }
    temp <- pop(x)
    if (is.null(temp)) 
        stop("no grouping factor (pop) provided")
    if (length(levels(temp)) < 2) {
        warning("There is only one pop - returning NULL")
        return(NULL)
    }
    res.type <- match.arg(res.type)
    f1 <- function(pop1, pop2) {
        #``````````````````````````````````````````````````````````````````````#
        # Lining up the individuals of the populations and repooling them for
        # Pairwise comparison
        #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,#
        n1 <- nrow(pop1@tab)
        n2 <- nrow(pop2@tab)
        temp <- repool(pop1, pop2)
        #``````````````````````````````````````````````````````````````````````#
        # This returns a single mean value of the heterozygosities of these
        # populations. This represents the "expected" value. 
        #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,#
        b <- weighted.mean(Hs(temp), c(n1, n2))
        pop(temp) <- NULL
        #``````````````````````````````````````````````````````````````````````#
        # returning the means of both populations. and then the f-statistic
        #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,#
        a <- Hs(temp)
        return((a - b)/a)
    }
    #``````````````````````````````````````````````````````````````````````````#
    # Separating out the popsubulations. temp <- pop(x) is not necessary
    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,#    
    lx <- seppop(x, treatOther = FALSE)
    temp <- pop(x)
    levPop <- levels(temp)
    #``````````````````````````````````````````````````````````````````````````#
    # getting all possible combinations using combn. This might be where I can
    # Modify it.
    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,#    
    allPairs <- combn(1:length(levPop), 2)
    if (!is.matrix(allPairs)) {
        allPairs <- matrix(allPairs, nrow = 2)
    }
    vecRes <- numeric()
    for (i in 1:ncol(allPairs)) {
        vecRes[i] <- f1(lx[[allPairs[1, i]]], lx[[allPairs[2, 
            i]]])
    }
    squelres <- dist(1:length(levPop))
    res <- vecRes
    #``````````````````````````````````````````````````````````````````````````#
    # Converting this into a distance triangle. dis
    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,#
    attributes(res) <- attributes(squelres)
    if (res.type == "matrix") {
        res <- as.matrix(res)
        if (truenames) {
            lab <- x@pop.names
        }
        else {
            lab <- names(x@pop.names)
        }
        colnames(res) <- rownames(res) <- lab
    }
    return(res)
}


.new.sampling <- function(pop,iterations,quiet="noisy",missing="ignore",type=type, method=1){ 
  METHODS = c("multilocus", "permute alleles", "parametric bootstrap",
      "non-parametric bootstrap")
  if(!is.list(pop)){
    if(type=="PA"){
      .Ia.Rd <- .PA.Ia.Rd
    }
  }
	sample.data <- NULL
	for (c in 1:iterations){
    IarD <- .Ia.Rd(total.shuffler(pop, method=method), missing=missing)
		sample.data <- rbind(sample.data, as.data.frame(list( 
                                      Ia=IarD[1],
                                      rbarD=IarD[2]
                                      )))
    if (quiet != TRUE){
      if(quiet == "noisy"){
        cat("Sample: ",c,"\n")
        cat("Index of Association: ", IarD[1],"\n")
        cat("Standardized Index of Association (rbarD): ", IarD[2],"\n")
      }
      else{
        if(c%%50 != 0){
          cat(".")
          if (c == iterations){
            cat("\n")
          }
#          else if( c == 666 ){
#            cat("\\m/")
#          }
        }
        else{
          cat(".")
          cat("\n")
        }      
      }
    }
	}
	return(sample.data)
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



getinds <- function(x){
  # x is a vector indicating the size of loci in a data set (must be continuous).
  # The sum of x is equal to the number of observations.
  to <- cumsum(x)
  from <- c(1, to[-length(to)]+1)
  indices <- rbind(from, to)
  rownames(indices) <- NULL
  return(indices)
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


new.diss.dist <- function(pop){
  ploid <- ploidy(pop)
  ind.names <- pop@ind.names
  inds <- nInd(pop)
  np <- choose(inds, 2)
  dist.vec <- matrix(data = 0, nrow=inds, ncol=inds)
  if(pop@type == "PA"){
    dist.vec[lower.tri(dist.vec)] <- .Call("pairdiffs",pop@tab*ploid)/ploid
  }
  else{
    pop <- seploc(pop)
    numLoci <- length(pop)
    temp.d.vector <- matrix(nrow = np, ncol = numLoci, data = as.numeric(NA))
    temp.d.vector <- vapply(pop, function(x) .Call("pairdiffs",x@tab*ploid)/ploid, 
                            temp.d.vector[, 1])
    dist.vec[lower.tri(dist.vec)] <- rowSums(temp.d.vector)
  }
  colnames(dist.vec) <- ind.names
  rownames(dist.vec) <- ind.names
  loci <- ifelse(is.list(pop), length(pop), nLoc(pop))
  return(as.dist(dist.vec/(loci*ploid)))
}



new.poppr.msn <- function (pop, distmat, palette = topo.colors, 
                           sublist = "All", blacklist = NULL, vertex.label = "MLG", 
                           gscale=TRUE, glim = c(0,0.8), gadj = 3, gweight = 1, 
                           wscale=TRUE, ...){
  if(!require(igraph)){
    stop("You must have the igraph library installed to use this function.\n")
  }
  if(class(distmat) != "dist"){
    if(is.matrix(distmat)){
      if(any(nInd(pop) != dim(distmat))){
        stop("The size of the distance matrix does not match the size of the data.\n")
      }
      distmat <- as.dist(distmat)
    }
    else{
      stop("The distance matrix is neither a dist object nor a matrix.\n")
    }
  }
  if(nInd(pop) != attr(distmat, "Size")){
    stop("The size of the distance matrix does not match the size of the data.\n")
  }
  gadj <- ifelse(gweight == 1, gadj, -gadj)
  # Storing the MLG vector into the genind object
  pop$other$mlg.vec <- mlg.vector(pop)
  #bclone <- as.matrix(distmat)[!duplicated(pop$other$mlg.vec), !duplicated(pop$other$mlg.vec)]
  
  singlepop <- function(pop, vertex.label){
    cpop <- pop[.clonecorrector(pop), ]
    mlg.number <- table(pop$other$mlg.vec)[rank(cpop$other$mlg.vec)]
    rownames(bclone) <- cpop$pop
    colnames(bclone) <- cpop$pop
    #bclone <- discreet.dist(cpop)
    mclone <- as.dist(bclone)
    #attr(bclone, "Labels") <- paste("MLG.", cpop$other$mlg.vec, sep="")
    g <- graph.adjacency(as.matrix(bclone),weighted=TRUE,mode="undirected")
    mst <- (minimum.spanning.tree(g,algorithm="prim",weights=E(g)$weight))
    if(!is.na(vertex.label[1]) & length(vertex.label) == 1){
      if(toupper(vertex.label) == "MLG"){
        vertex.label <- paste("MLG.", cpop$other$mlg.vec, sep="")
      }
      else if(toupper(vertex.label) == "INDS"){
        vertex.label <- cpop$ind.names
      }
    }
    
    if(gscale == TRUE){
      E(mst)$color <- gray(adjustcurve(E(mst)$weight, glim=glim, correct=gadj, 
                                       show=FALSE))
    }
    else{
      E(mst)$color <- rep("black", length(E(mst)$weight))
    }
    
    edgewidth <- 2
    if(wscale==TRUE){
      edgewidth <- 1/(E(mst)$weight)
      if(any(E(mst)$weight < 0.08)){
        edgewidth <- 1/(E(mst)$weight + 0.08)
      }
    }
    populations <- ifelse(is.null(pop(pop)), NA, pop$pop.names)
    plot(mst, edge.width = edgewidth, edge.color = E(mst)$color,  
         vertex.label = vertex.label, vertex.size = mlg.number*3, 
         vertex.color = palette(1),  ...)
    legend(-1.55,1,bty = "n", cex = 0.75, 
           legend = populations, title = "Populations", fill = palette(1), 
           border = NULL)
    E(mst)$width <- edgewidth
    V(mst)$size <- mlg.number
    V(mst)$color <- palette(1)
    V(mst)$label <- vertex.label
    return(list(graph = mst, populations = populations, colors = palette(1)))
  }
  bclone <- as.matrix(distmat)
  # The clone correction of the matrix needs to be done at this step if there
  # is only one or no populations. 
  if(is.null(pop(pop)) | length(pop@pop.names) == 1){
    bclone <- bclone[!duplicated(pop$other$mlg.vec), !duplicated(pop$other$mlg.vec)]
    return(singlepop(pop, vertex.label))
  }
  # This will subset both the population and the matrix. 
  if(sublist[1] != "ALL" | !is.null(blacklist)){
    sublist_blacklist <- sub_index(pop, sublist, blacklist)
    bclone <- bclone[sublist_blacklist, sublist_blacklist]
    pop <- popsub(pop, sublist, blacklist)
  }
  
  # This will clone correct the incoming matrix. 
  bclone <- bclone[!duplicated(pop$other$mlg.vec), !duplicated(pop$other$mlg.vec)]
  
  if(is.null(pop(pop)) | length(pop@pop.names) == 1){
    return(singlepop(pop, vertex.label))
  }
  # Obtaining population information for all MLGs
  mlg.cp <- mlg.crosspop(pop, mlgsub=1:mlg(pop, quiet=TRUE), quiet=TRUE)
  names(mlg.cp) <- paste("MLG.",sort(unique(pop$other$mlg.vec)),sep="")
  cpop <- pop[.clonecorrector(pop), ]
  
  # This will determine the size of the nodes based on the number of individuals
  # in the MLG. Subsetting by the MLG vector of the clone corrected set will
  # give us the numbers and the population information in the correct order.
  # Note: rank is used to correctly subset the data
  mlg.number <- table(pop$other$mlg.vec)[rank(cpop$other$mlg.vec)]
  mlg.cp <- mlg.cp[rank(cpop$other$mlg.vec)]
  #bclone <- discreet.dist(cpop)
  rownames(bclone) <- cpop$pop
  colnames(bclone) <- cpop$pop

  g <- graph.adjacency(as.matrix(bclone), weighted=TRUE, mode="undirected")
  mst <- (minimum.spanning.tree(g,algorithm="prim",weights=E(g)$weight))
  
  if(!is.na(vertex.label[1]) & length(vertex.label) == 1){
    if(toupper(vertex.label) == "MLG"){
      vertex.label <- paste("MLG.", cpop$other$mlg.vec, sep="")
    }
    else if(toupper(vertex.label) == "INDS"){
      vertex.label <- cpop$ind.names
    }
  }
  ###### Color schemes #######  
  # The pallete is determined by what the user types in the argument. It can be 
  # rainbow, topo.colors, heat.colors ...etc.
  palette <- match.fun(palette)
  color <- palette(length(pop@pop.names))
  
  ###### Edge adjustments ######
  # Grey Scale Adjustment weighting towards more diverse or similar populations.
  if(gscale == TRUE){
    E(mst)$color <- gray(adjustcurve(E(mst)$weight, glim=glim, correct=gadj, 
                                     show=FALSE))
  }
  else{
    E(mst)$color <- rep("black", length(E(mst)$weight))
  }
  
  # Width scale adjustment to avoid extremely large widths.
  # by adding 0.08 to entries, the max width is 12.5 and the min is 0.9259259
  edgewidth <- 2
  if(wscale==TRUE){
    edgewidth <- 1/(E(mst)$weight)
    if(any(E(mst)$weight < 0.08)){
      edgewidth <- 1/(E(mst)$weight + 0.08)
    }
  }
  
  # This creates a list of colors corresponding to populations.
  mlg.color <- lapply(mlg.cp, function(x) color[pop@pop.names %in% names(x)])

  plot(mst, edge.width = edgewidth, edge.color = E(mst)$color, 
       vertex.size = mlg.number*3, vertex.shape = "pie", vertex.pie = mlg.cp, 
       vertex.pie.color = mlg.color, vertex.label = vertex.label, ...)
  legend(-1.55 ,1 ,bty = "n", cex = 0.75, legend = pop$pop.names, 
         title = "Populations", fill=color, border=NULL)
  E(mst)$width <- edgewidth
  V(mst)$size <- mlg.number
  V(mst)$shape <- "pie"
  V(mst)$pie <- mlg.cp
  V(mst)$pie.color <- mlg.color
  V(mst)$label <- vertex.label
  return(list(graph = mst, populations = pop$pop.names, colors = color))
}



new.bruvo.msn <- function (pop, replen=c(1), palette = topo.colors,
                           sublist = "All", blacklist = NULL, vertex.label = "MLG", 
                           gscale=TRUE, glim = c(0,0.8), gadj = 3, gweight = 1, wscale=TRUE, ...){
  stopifnot(require(igraph))
  gadj <- ifelse(gweight == 1, gadj, -gadj)
  # Storing the MLG vector into the genind object
  pop$other$mlg.vec <- mlg.vector(pop)
  
  singlepop <- function(pop, vertex.label){
    cpop <- pop[.clonecorrector(pop), ]
    mlg.number <- table(pop$other$mlg.vec)[rank(cpop$other$mlg.vec)]
    bclone <- bruvo.dist(cpop, replen=replen)
    mclone<-as.dist(bclone)
    #attr(bclone, "Labels") <- paste("MLG.", cpop$other$mlg.vec, sep="")
    g <- graph.adjacency(as.matrix(bclone),weighted=TRUE,mode="undirected")
    mst <- (minimum.spanning.tree(g,algorithm="prim",weights=E(g)$weight))
    if(!is.na(vertex.label[1]) & length(vertex.label) == 1){
      if(toupper(vertex.label) == "MLG"){
        vertex.label <- paste("MLG.", cpop$other$mlg.vec, sep="")
      }
      else if(toupper(vertex.label) == "INDS"){
        vertex.label <- cpop$ind.names
      }
    }
    
    if(gscale == TRUE){
      E(mst)$color <- gray(adjustcurve(E(mst)$weight, glim=glim, correct=gadj, 
                                       show=FALSE))
    }
    else{
      E(mst)$color <- rep("black", length(E(mst)$weight))
    }
    
    edgewidth <- 2
    if(wscale == TRUE){
      edgewidth <- 1/(E(mst)$weight)
      if(any(E(mst)$weight < 0.08)){
        edgewidth <- 1/(E(mst)$weight + 0.08)
      }
    }
    populations <- ifelse(is.null(pop(pop)), NA, pop$pop.names)
    plot(mst, edge.width = edgewidth, edge.color = E(mst)$color,  
         vertex.label = vertex.label, vertex.size = mlg.number*3, 
         vertex.color = palette(1),  ...)
    legend(-1.55,1,bty = "n", cex = 0.75, 
           legend = populations, title = "Populations", fill = palette(1), 
           border = NULL)
    E(mst)$width <- edgewidth
    V(mst)$size <- mlg.number
    V(mst)$color <- palette(1)
    V(mst)$label <- vertex.label
    return(list(graph = mst, populations = populations, colors = palette(1)))
  }
  if(is.null(pop(pop)) | length(pop@pop.names) == 1){
    return(singlepop(pop, vertex.label))
  }
  if(sublist[1] != "ALL" | !is.null(blacklist)){
    pop <- popsub(pop, sublist, blacklist)
  }
  if(is.null(pop(pop)) | length(pop@pop.names) == 1){
    return(singlepop(pop, vertex.label))
  }
  # Obtaining population information for all MLGs
  mlg.cp <- mlg.crosspop(pop, mlgsub=1:mlg(pop, quiet=TRUE), quiet=TRUE)
  names(mlg.cp) <- paste("MLG.",sort(unique(pop$other$mlg.vec)),sep="")
  cpop <- pop[.clonecorrector(pop), ]
  # This will determine the size of the nodes based on the number of individuals
  # in the MLG. Subsetting by the MLG vector of the clone corrected set will
  # give us the numbers and the population information in the correct order.
  # Note: rank is used to correctly subset the data
  mlg.number <- table(pop$other$mlg.vec)[rank(cpop$other$mlg.vec)]
  mlg.cp <- mlg.cp[rank(cpop$other$mlg.vec)]
  bclone <- bruvo.dist(cpop, replen=replen)

  ###### Create a graph #######
  g <- graph.adjacency(as.matrix(bclone), weighted = TRUE, mode = "undirected")
  mst <- (minimum.spanning.tree(g, algorithm = "prim", weights = E(g)$weight))
  
  if(!is.na(vertex.label[1]) & length(vertex.label) == 1){
    if(toupper(vertex.label) == "MLG"){
      vertex.label <- paste("MLG.", cpop$other$mlg.vec, sep="")
    }
    else if(toupper(vertex.label) == "INDS"){
      vertex.label <- cpop$ind.names
    }
  }
  ###### Color schemes #######  
  # The pallete is determined by what the user types in the argument. It can be 
  # rainbow, topo.colors, heat.colors ...etc.
  palette <- match.fun(palette)
  color <- palette(length(pop@pop.names))
  if(gscale == TRUE){
    E(mst)$color <- gray(adjustcurve(E(mst)$weight, glim=glim, correct=gadj, 
                                     show=FALSE))
  }
  else{
    E(mst)$color <- rep("black", length(E(mst)$weight))
  }
  
  edgewidth <- 2
  if(wscale == TRUE){
    edgewidth <- 1/(E(mst)$weight)
    if(any(E(mst)$weight < 0.08)){
      edgewidth <- 1/(E(mst)$weight + 0.08)
    }
  }
  # This creates a list of colors corresponding to populations.
  mlg.color <- lapply(mlg.cp, function(x) color[pop@pop.names %in% names(x)])
  
  plot(mst, edge.width = edgewidth, edge.color = E(mst)$color, 
       vertex.size = mlg.number*3, vertex.shape = "pie", vertex.pie = mlg.cp, 
       vertex.pie.color = mlg.color, vertex.label = vertex.label, ...)
  legend(-1.55 ,1 ,bty = "n", cex = 0.75, legend = pop$pop.names, 
         title = "Populations", fill = color, border = NULL)
  E(mst)$width <- edgewidth
  V(mst)$size <- mlg.number
  V(mst)$shape <- "pie"
  V(mst)$pie <- mlg.cp
  V(mst)$pie.color <- mlg.color
  V(mst)$label <- vertex.label
  return(list(graph = mst, populations = pop$pop.names, colors = color))
}


new.genind2genalex <- function(pop, filename="genalex.csv", quiet=FALSE, geo=FALSE, geodf="xy"){
  if(!is.genind(pop)) stop("A genind object is needed.")
  if(is.null(pop@pop)){
    pop(pop) <- rep("Pop", nInd(pop))
  }
  popcall <- match.call()
  #topline is for the number of loci, individuals, and populations.
  topline <- c(nLoc(pop), nInd(pop), length(pop@pop.names))
  popsizes <- table(pop@pop)
  # The sizes of the populations correspond to the second line, which is the pop
  # names. 
  topline <- c(topline, popsizes)
  secondline <- c("", "", "", pop@pop.names)
  ploid <- ploidy(pop)
  # Constructing the locus names. GenAlEx separates the alleles of the loci, so
  # There is one locus name for every p ploidy columns you have.
  if(ploid > 1 & pop@type == "codom"){
    locnames <- unlist(strsplit(paste(pop@loc.names, 
                                      paste(rep(" ", ploidy(pop)-1), 
                                            collapse="/"), sep="/"),"/"))
  }
  else{
    locnames <- pop@loc.names
  }
  thirdline <- c("Ind","Pop", locnames)
  
  # This makes sure that you don't get into a stacking error when stacking the
  # first three rows.
  if(length(thirdline) > length(topline)){
    lenfac <- length(thirdline) - length(topline)
    topline <- c(topline, rep("", lenfac))
    secondline <- c(secondline, rep("", lenfac))
  }
  else if(length(thirdline) < length(topline)){
    lenfac <- length(topline) - length(thirdline)
    thirdline <- c(thirdline, rep("", lenfac))
  }
  infolines <- rbind(topline, secondline, thirdline)
  
  # converting to a data frame
  if(any(!pop@tab %in% c(0, ((1:ploid)/ploid), 1, NA))){
    pop@tab[!pop@tab %in% c(0, ((1:ploid)/ploid), 1, NA)] <- NA
  }
  if(!quiet) cat("Extracting the table ... ")
  df <- genind2df(pop, oneColPerAll=TRUE)
  
  # making sure that the individual names are included.
  if(all(pop@ind.names == "") | is.null(pop@ind.names)){
    pop@ind.names <- paste("ind", 1:nInd(pop), sep="")
  }
  df <- cbind(pop@ind.names, df)
  # setting the NA replacement. This doesn't work too well. 
  replacement <- ifelse(pop@type =="PA","-1","0")
  if(!quiet) cat("Writing the table to",filename,"... ")
  
  if(geo == TRUE & !is.null(pop$other[[geodf]])){
    replacemat <- matrix("", 3, 3)
    replacemat[3, 2:3] <- c("X", "Y")
    infolines <- cbind(infolines, replacemat)
    df2 <- data.frame(list("Space" = rep("", nInd(pop))))
    gdf <- as.matrix(pop@other[[geodf]])
    if(nrow(gdf) < nInd(pop)){
      gdf <- rbind(gdf, matrix("", nInd(pop) - nrow(gdf), 2))
    }
    df <- cbind(df, df2, gdf)
  }
  else{
    popcall <- popcall[2]
    warning(paste("There is no data frame or matrix in ",
                  paste(substitute(popcall), collapse=""), 
                  "@other called ",geodf,
                  ".\nThe xy coordinates will not be represented in the resulting file.", sep=""))
  }
  
  write.table(infolines, file=filename, quote=TRUE, row.names=FALSE, 
              col.names=FALSE, sep=",")
  write.table(df, file=filename, quote=TRUE, na=replacement, append=TRUE, 
              row.names=FALSE, col.names=FALSE, sep=",")
  if(!quiet) cat("Done.\n")
}

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
    else{
      pop.vec <- ifelse(any(gena[, 1]==pop.info[1]), 1, 2)
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
  else if (geo == TRUE & length(pop.info) == glob.info[3]){
    reg.vec <- NULL
    pop.vec <- gena[, 2]
    ind.vec <- gena[, 1]
    xy <- gena[, c((ncol(gena)-1), ncol(gena))]
    gena <- gena[, c(-1,-2,-(ncol(gena)-1),-ncol(gena))]
  }
  else{
    reg.vec <- NULL
    pop.vec <- gena[, 2]
    ind.vec <- gena[, 1]
    xy <- NULL
    gena <- gena[, c(-1,-2)]
  }
  clm <- ncol(gena)
  # Checking for diploid data.
  if (glob.info[1] == clm/2){
    # Missing data in genalex is coded as "0" for non-presence/absence data.
    # this converts it to "NA" for adegenet.
    #gena <- as.data.frame(gsub("  0", "NA", as.matrix(gena)))
    if(any(gena =="0")){
      gena.mat <- as.matrix(gena)
      gena.mat[which(gena=="0")] <- NA
      gena <- as.data.frame(gena.mat)
    }
    type='codom'
    gena2 <- gena[, which((1:clm)%%2==1)]
    loci <- which((1:clm)%%2==1)
    lapply(loci, function(x) gena2[, ((x-1)/2)+1] <<-
             paste(gena[, x],"/",gena[, x+1], sep=""))
    res <- list(Gena=gena2, Glob.info=glob.info, Ploid=ploidy)
    res.gid <- df2genind(res$Gena, sep="/", ind.names=ind.vec, pop=pop.vec,
                         ploidy=res$Ploid, type=type)
  }
  # Checking for AFLP data.
  else if (glob.info[1] == clm & any(gena == 1)) {
    # Missing data in genalex is coded as "-1" for presence/absence data.
    # this converts it to "NA" for adegenet.
    if(any(gena==-1)){
      gena.mat <- as.matrix(gena)
      gena.mat[which(gena==-1)] <- NA
      gena <- as.data.frame(gena.mat)
    }
    type='PA'
    res <- list(Gena=gena, Glob.info=glob.info, Ploid=ploidy)
    res.gid <- df2genind(res$Gena, ind.names=ind.vec, pop=pop.vec,
                         ploidy=res$Ploid, type=type)
  }
  # Checking for haploid microsattellite data.
  else if (glob.info[1] == clm & all(!gena %in% c(1,-1))) {
    if(any(gena =="0")){
      gena.mat <- as.matrix(gena)
      gena.mat[which(gena=="0")] <- NA
      gena <- as.data.frame(gena.mat)
    }
    type = 'codom'
    res <- list(Gena=gena, Glob.info=glob.info, Ploid=1)
    res.gid <- df2genind(res$Gena, ind.names=ind.vec, pop=pop.vec,
                         ploidy=res$Ploid, type=type)
  }
  else {
    stop("Did you forget to set the geo or region flags?")
  }
  if (any(duplicated(ind.vec))){
    # preserving the names
    orig.ind.vec <- ind.vec
    # ensuring that all names are unique
    ind.vec <- 1:length(ind.vec)
    res.gid@other[["original_names"]] <- orig.ind.vec
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
  cat ("hello dawg")
}