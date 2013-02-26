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

# Creates Minimum spanning trees using Bruvo's distance an igraph minimun number 
# spanning trees. In test, made by Javier "Oh, great master" Tabima. Hey Zhian,
# do you like my new nickname??
bruvo.mstree <- function (pop, replen=1, vertex.size=8, interactive=FALSE){
  library(igraph)
  if (!is.genind(pop)) {
    stop("This function requires a genind object.")
  }
  bd <- bruvo.dist(pop,replen=replen)  
  g <- graph.adjacency(as.matrix(bd),weighted=TRUE,mode="undirected")
  V(g)$pop <- as.character(pop$pop)
  color.list <- colors()[1:500%%5 == 0]
  fc <- factor(V(g)$pop)
  levels(fc) <- color.list[1:length(levels(fc))]
  V(g)$color <- as.character(fc)
  if (interactive==TRUE){
    derp <- tkplot(minimum.spanning.tree(g),edge.color="black",edge.width=E(g)$weight,vertex.size=vertex.size)
    return(derp)
  }
  else{
    plot.igraph(minimum.spanning.tree(g),edge.color="black",edge.width=E(g)$weight,vertex.size=vertex.size)
  }  
}

################################################################################
# 
# Inserted methods for subsetting and dealing with null or single populations.
# There was a little trickiness dealing with naming the nodes, but I made a hack
# to fix it. It's not pretty, but it gets the job done. Currently the user has
# a few choices, to label each node based on the MLGs, they select vertex.label
# = "mlg", if they want individuals, they select vertex.label = "inds", if they
# want something else, they add a vector into vertex.label and get whatever
# value is in that vector. NA is a popular option. 
#
# Function for manual:
# set.seed(2013)
# bruvo.msn(nancycats, replen=rep(1, 9), vertex.label="inds", sublist=c(8:9), 
#   vertex.label.cex=0.7, vertex.label.dist=0.4, palette=heat.colors)
################################################################################

test.bruvo.msn <- function (pop, replen=c(1), palette = topo.colors,
                            sublist = "All", blacklist = NULL, vertex.label = "MLG", ...){
  stopifnot(require(igraph))
  
  # Storing the MLG vector into the genind object
  pop$other$mlg.vec <- mlg.vector(pop)
  
  singlepop <- function(pop, vertex.label){
    cpop <- pop[.clonecorrector(pop), ]
    mlg.number <- table(pop$other$mlg.vec)[rank(cpop$other$mlg.vec)]
    bclone <- bruvo.dist(cpop, replen=replen)
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
    l <- layout.drl
    plot(mst, edge.color="black", edge.width=2, vertex.label = vertex.label,
         vertex.size=mlg.number*3, vertex.color = palette(1), layout=l, edge.label=E(mst)$weight, ...)
    legend(-1.55,1,bty = "n", cex=0.75, legend=pop$pop.names, title="Populations",
           fill=palette(1), border=NULL)
    return(invisible(1))
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
  ###### Change names to MLGs #######
  #attr(bclone, "Labels") <- paste("MLG.", cpop$other$mlg.vec, sep="")
  ###### Create a graph #######
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
  
  # This creates a list of colors corresponding to populations.
  mlg.color <- lapply(mlg.cp, function(x) color[pop@pop.names %in% names(x)])
  #print(mlg.color)
  plot(mst, edge.color="black", edge.width=2, vertex.size=mlg.number*3,
       vertex.shape="pie", vertex.pie=mlg.cp, vertex.pie.color=mlg.color,
       vertex.label = vertex.label, edge.label=E(mst)$weight, ...)
  legend(-1.55,1,bty = "n", cex=0.75, legend=pop$pop.names, title="Populations",
         fill=color, border=NULL)
}



discreet.dist <- function(pop){
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
  return(as.dist(dist.vec))
}





poppr.msn <- function (pop, distmat, palette = topo.colors,
                            sublist = "All", blacklist = NULL, vertex.label = "MLG", ...){
  stopifnot(require(igraph))
  
  # Storing the MLG vector into the genind object
  pop$other$mlg.vec <- mlg.vector(pop)
  mat <- as.matrix(distmat)[!duplicated(pop$other$mlg), !duplicated(pop$other$mlg)]
  cpop <- pop[.clonecorrector(pop),]
  rownames(mat) <- cpop$pop
  colnames(mat) <- cpop$pop
  mlg.cp <- mlg.crosspop(pop, mlgsub=1:mlg(pop, quiet=TRUE), quiet=TRUE)
  names(mlg.cp) <- paste("MLG.",sort(unique(pop$other$mlg.vec)),sep="")
  mlg.number <- table(pop$other$mlg.vec)[rank(cpop$other$mlg.vec)]
  mlg.cp <- mlg.cp[rank(cpop$other$mlg.vec)]
  g <- graph.adjacency(as.matrix(mat),weighted=TRUE,mode="undirected")
  mst <- (minimum.spanning.tree(g,algorithm="prim",weights=E(g)$weight))
  if(!is.na(vertex.label[1]) & length(vertex.label) == 1){
    if(toupper(vertex.label) == "MLG"){
      vertex.label <- paste("MLG.", cpop$other$mlg.vec, sep="")
    }
    else if(toupper(vertex.label) == "INDS"){
      vertex.label <- cpop$ind.names
    }
  }
  
  palette <- match.fun(palette)
  color <- palette(length(pop@pop.names))
  
  # This creates a list of colors corresponding to populations.
  mlg.color <- lapply(mlg.cp, function(x) color[pop@pop.names %in% names(x)])
  #print(mlg.color)
  plot(mst, edge.color="black", edge.width=2, vertex.size=mlg.number*3,
       vertex.shape="pie", vertex.pie=mlg.cp, vertex.pie.color=mlg.color,
       vertex.label = vertex.label, ...)
  legend(-1.55,1,bty = "n", cex=0.75, legend=pop$pop.names, title="Populations",
         fill=color, border=NULL)
}