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
# This will make a histogram of permutations of a dataset and show a red line
# indicating p = 0.05 and a green line indicating the observed value. 
# 
#==============================================================================#
permut.histogram.maker <- function(index="index",sampled,observed,pop){
  if(all(is.nan(sampled[[index]]))){
    warning("Cannot create histogram due to NaN's")
    next
  }
	xmin <- min(observed, sampled[[index]])
	xmax <- max(observed, sampled[[index]])
  # R does not like it when the max value for xlim is infinity, so this will fix
  # the culprit where it occurs, usually in the rbarD caluclation. I will be
  # surprised if I see it in the Ia calculation.
  if(!is.nan(xmax) & xmax==Inf & index=="rbarD"){
    xmax <- 1
  }
  if(is.nan(xmin)){
    xmin <- -0.5
  }
  if(is.nan(xmax)){
    xmax <- ifelse(index=="Ia", 18, 1)
  }
	hist(sampled[[index]], xlim=c(xmin, xmax),
        main=c("Population:",as.character(pop)),
        xlab=sprintf("%s (%d permutations)",names(sampled[index]),
                      length(sampled[[index]])), col="grey")
	abline(v=observed, col="green")
  perc95 <- quantile(sampled[[index]], 0.95, na.rm=TRUE)[[1]]
  abline(v=perc95, col="red")
}
#==============================================================================#
# histogram will create two histograms showing the sampled distributions of the
# index of association and the standardized index of association along with an
# informative legend between them.
#==============================================================================#
permut.histogram <- function(sampled, observed, pval, pop="pop", file="file"){
  par(xpd=TRUE,mfrow=c(3,1))
  permut.histogram.maker(index="Ia",sampled,observed[1],pop)
  if (!is.na(pval)){
    pval <- ifelse(pval==0, sprintf("< %g", 1/length(sampled[["Ia"]])), pval)
    frame()
    legend('center', col=c("red", "green"), lty=c(1,1), c("p = 0.05",
              paste("Observed\n","(p-value  ",pval,")", sep="")), 
                        horiz=TRUE, title=paste("File : ",file,sep=""))
  }
  permut.histogram.maker(index="rbarD",sampled,observed[2],pop)
}

poppr.plot <- function(sample, pval = c("0.05", "0.05"), pop="pop", 
                        observed = observed, file="file", N=NA){
  #   if (!is.na(pval[1])){
  #     pval[1] <- ifelse(pval[1]==0, sprintf("< %g", 1/length(sample[["Ia"]])), pval[1])
  #   }
  #   if (!is.na(pval[2])){
  #     pval[2] <- ifelse(pval[2]==0, sprintf("< %g", 1/length(sample[["rbarD"]])), pval[2])
  #   }
  #````````````````````````````````````````````````````````````````````````````#
  # In the case that the sample contains all NaNs, a graph cannot be displayed.
  # In place of the graph will be an orange box with a warning message.
  #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,#
  
  if (all(is.nan(sample$rbarD))){
    warning(paste("The data from ",file,", population: ",pop," contains only missing values and cannot be displayed graphically", sep=""))
    background_stuff <- theme(panel.grid.major.y = element_line(size=0)) +
      theme(panel.grid.minor.y = element_line(size=0)) +
      theme(panel.grid.major.x = element_line(size=0)) +
      theme(panel.background = element_rect(fill="grey95")) +
      #theme(plot.background = element_rect(size=1, color="black")) +
      theme(axis.ticks.y = element_line(size=0)) +
      theme(axis.text.y = element_text(size=0)) +
      theme(axis.ticks.x = element_line(size=0)) +
      theme(axis.text.x = element_text(size=0)) +
      theme(axis.title.y = element_text(size=rel(0)))+
      theme(axis.title.x = element_text(size=rel(0)))
    oops <- ggplot(as.data.frame(list(x=-10:9)), aes_string(x = "x")) + 
      geom_histogram(binwidth=1, fill="orange") + 
      geom_text(aes(label="Warning:", x=0, y=0.8), color="black", size=rel(15)) + 
      geom_text(aes(label="Data contains only NaNs and\ncannot be displayed graphically", 
                    x=0, y=0.5, hjust=0.5), color="black", size=rel(10)) +
      labs(title=paste("Population: ", pop, "; N: ", N, "\nPermutations: ", 
                       length(sample$Ia), "\nFile: ", file, sep="")) + 
      theme(plot.title = element_text(vjust=1, size=rel(2), face="bold")) +
      background_stuff
    print(oops)
  }
  else {
    #``````````````````````````````````````````````````````````````````````````#
    # Normal Cases
    #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,#
    if(any(is.nan(sample[["rbarD"]])))
      sample[["rbarD"]][which(is.nan(sample[["rbarD"]]))] <- mean(sample[["rbarD"]], na.rm=TRUE)
    if(any(is.nan(sample[["Ia"]])))
      sample[["Ia"]][which(is.nan(sample[["Ia"]]))] <- mean(sample[["Ia"]], na.rm=TRUE)
    
    # Transforming the data into a data frame that ggplot can understand and put
    # into facets. It will have the following rows: Value, Index, and Quant.
    Indexfac <- factor(1:2, levels=1:2, labels=c("I[A]","bar(r)[d]"))
    infodata <- as.data.frame(list(Value=c(sample[["Ia"]], sample[["rbarD"]]), 
                                   Index=rep(Indexfac, each=length(sample[["Ia"]]))))
    
    # Setting up the observed values for formatting the overlaying lines.
    # It will have the following rows: Observed, Index, Quant, Sam, P, min, max
    
    obsdata <- data.frame(list(Observed=observed[1:2], Index=Indexfac))
    obsdata[["P"]] <- paste("Observed \n(p-value: ", pval,")",sep="")
    obsdata[["median"]] <- c(median(infodata[infodata[["Index"]] == Indexfac[1], "Value"]),
                        median(infodata[infodata[["Index"]] == Indexfac[2], "Value"]))


    obsdata[["label"]] <- paste("Observed: ",obsdata[["Observed"]], sep="")
    if(any(is.na(observed))){
      warning(paste("The Index of Association values from ",file,", population: ",pop," contain missing values and cannot be displayed graphically", sep=""))
      iard_plot <- ggplot(infodata, aes_string(x = "Value")) + 
        # Giving the data over to the histogram creating function and removing
        # all of the lines from each bar, so it's displayed as a solid area.
        geom_histogram(linetype="blank", #alpha=0.8, 
                       data = infodata[infodata[["Index"]] == Indexfac[1], ], 
                       position="identity",
                       binwidth=diff(range(infodata[infodata[["Index"]] == Indexfac[1], "Value"]))/30) + 
        geom_histogram(linetype="blank", #alpha=0.8, 
                       data = infodata[infodata[["Index"]] == Indexfac[2], ], 
                       position="identity",
                       binwidth=diff(range(infodata[infodata[["Index"]] == Indexfac[2], "Value"]))/30) + 
        geom_rug() + 
        # The label for the observed line is a bit more difficult to code as
        # it has the ability to appear anywhere on the chart. Here, I'm
        # forcing it to flip to one side or the other based on which side of
        # the mean that the observed value falls on.
        geom_text(aes_string(label="label", x=0, y=Inf, vjust=1.5),
                  data=obsdata, angle=0, color="red") + 
        
        # Splitting the data into separate plots with free x-axes.
        facet_grid(". ~ Index", scales="free_x", labeller=label_parsed) +
        
        # Title of the plot.
        labs(title=paste("Population: ", pop, "; N: ", N, "\nPermutations: ", 
                         length(sample$Ia), "\nFile: ", file, sep="")) + 
        theme_classic() %+replace%
        theme(plot.title = element_text(vjust=1, size=rel(2), face="bold")) +
        theme(panel.background = element_rect(fill="grey98")) +
        
        # Making the Index titles bigger. 
        theme(strip.text.x = element_text(size=rel(3), face="bold"))
    }
    else{
      iard_plot <- ggplot(infodata, aes_string(x = "Value")) + 
        # Giving the data over to the histogram creating function and removing
        # all of the lines from each bar, so it's displayed as a solid area.
        geom_histogram(linetype="blank", #alpha=0.8, 
                       data = infodata[infodata[["Index"]] == Indexfac[1], ], 
                       position="identity",
                       binwidth=diff(range(infodata[infodata[["Index"]] == Indexfac[1], "Value"]))/30) + 
        geom_histogram(linetype="blank", #alpha=0.8, 
                       data = infodata[infodata[["Index"]] == Indexfac[2], ], 
                       position="identity",
                       binwidth=diff(range(infodata[infodata[["Index"]] == Indexfac[2], "Value"]))/30) + 
        geom_rug() + 
        # Positioning the observed line and labeling it.
        geom_vline(aes_string(xintercept = "Observed"), data = obsdata, color="blue", 
                   show_guide=TRUE, linetype="dashed") +
        
        # The label for the observed line is a bit more difficult to code as
        # it has the ability to appear anywhere on the chart. Here, I'm
        # forcing it to flip to one side or the other based on which side of
        # the mean that the observed value falls on.
        #geom_text(aes_string(label=paste("Observed \n(p-value: ", P,")",sep=""),
        geom_text(aes_string(label="P", 
                      x="Observed",y=Inf,vjust=2,hjust="ifelse(Observed > median, 1.01, -0.01)"),
                  data=obsdata, angle=0, color="blue") + 
        
        # Splitting the data into separate plots with free x-axes.
        facet_grid(". ~ Index", scales="free_x", labeller=label_parsed) +
        
        # Title of the plot.
        labs(title=paste("Population: ", pop, "; N: ", N, "\nPermutations: ", 
                         length(sample$Ia), "\nFile: ", file, sep="")) + 
        theme_classic() %+replace%
        theme(plot.title = element_text(vjust=1, size=rel(2), face="bold")) +
        theme(panel.background = element_rect(fill="grey98")) +
        
        # Making the Index titles bigger. 
        theme(strip.text.x = element_text(size=rel(3), face="bold"))
    }
    print(iard_plot)
  }
} 

#==============================================================================#
#
#' Create a minimum spanning network of selected populations using a distance
#' matrix.
#'
#' @param pop a \code{\link{genind}} object
#'
#' @param distmat a distance matrix that has been derived from your data set.
#' 
#' @param palette a \code{function} defining the color palette to be used to
#' color the populations on the graph. It defaults to \code{\link{topo.colors}},
#' but you can easily create new schemes by using \code{\link{colorRampPalette}}
#' (see examples for details)
#' 
#' @param sublist a \code{vector} of population names or indexes that the user
#' wishes to keep. Default to "ALL".
#'
#' @param blacklist a \code{vector} of population names or indexes that the user
#' wishes to discard. Default to \code{NULL}
#' 
#' @param vertex.label a \code{vector} of characters to label each vertex. There
#' are two defaults: \code{"MLG"} will label the nodes with the multilocus genotype
#' from the original data set and \code{"inds"} will label the nodes with the 
#' representative individual names.
#' 
#' @param gscale "grey scale". If this is \code{TRUE}, this will scale the color
#' of the edges proportional to the observed distance, with the lines becoming 
#' darker for more related nodes. See \code{\link{greycurve}} for details.
#' 
#' @param glim "grey limit". Two numbers between zero and one. They determine 
#' the upper and lower limits for the \code{\link{gray}} function. Default is 0
#' (black) and 0.8 (20\% black). See \code{\link{greycurve}} for details.
#' 
#' @param gadj "grey adjust". a positive \code{integer} greater than zero that
#' will serve as the exponent to the edge weight to scale the grey value to
#' represent that weight. See \code{\link{greycurve}} for details.
#' 
#' @param gweight "grey weight". an \code{integer}. If it's 1, the grey scale
#' will be weighted to emphasize the differences between closely related nodes.
#' If it is 2, the grey scale will be weighted to emphasize the differences
#' between more distantly related nodes. See \code{\link{greycurve}} for details.
#' 
#' @param wscale "width scale". If this is \code{TRUE}, the edge widths will be
#' scaled proportional to the inverse of the observed distance , with the lines 
#' becoming thicker for more related nodes.
#' 
#' @param ... any other arguments that could go into plot.igraph
#'
#' @return 
#' \item{graph}{a minimum spanning network with nodes corresponding to MLGs within
#' the data set. Colors of the nodes represent population membership. Width and
#' color of the edges represent distance.}
#' \item{populations}{a vector of the population names corresponding to the 
#' vertex colors}
#' \item{colors}{a vector of the hexadecimal representations of the colors used
#' in the vertex colors}
#' 
#' @note The edges of these graphs may cross each other if the graph becomes too
#' large. 
#'
#' @seealso \code{\link{nancycats}}, \code{\link{upgma}}, \code{\link{nj}},
#' \code{\link{nodelabels}}, \code{\link{na.replace}}, \code{\link{missingno}}, 
#' \code{\link{bruvo.msn}}, \code{\link{greycurve}}.
#' 
#' @export
#' @author Javier F. Tabima, Zhian N. Kamvar
#' @examples
#' 
#' # Load the data set and calculate the distance matrix for all individuals.
#' data(Aeut)
#' A.dist <- diss.dist(Aeut)
#' 
#' # Graph it.
#' A.msn <- poppr.msn(Aeut, A.dist, gadj=15, vertex.label=NA)
#' 
#' \dontrun{
#' # Set subpopulation structure.
#' Aeut.sub <- splitcombine(Aeut, method=2, hier=c("Pop", "Subpop"))
#' 
#' # Plot respective to the subpopulation structure
#' As.msn <- poppr.msn(Aeut.sub, A.dist, gadj=15, vertex.label=NA)
#' 
#' # Show only the structure of the Athena population.
#' As.msn <- poppr.msn(Aeut.sub, A.dist, gadj=15, vertex.label=NA, sublist=1:10)
#' 
#' # Let's look at the structure of the microbov data set
#' data(microbov)
#' micro.dist <- diss.dist(microbov)
#' micro.msn <- poppr.msn(microbov, diss.dist(microbov), vertex.label=NA)
#' 
#' Let's plot it and show where individuals have < 15% of their genotypes 
#' different.
#' 
#' plot(micro.msn$graph, edge.label = ifelse(E(micro.msn$graph)$weight < 0.15, 
#' round(E(micro.msn$graph)$weight, 3), NA), vertex.size=2, edge.label.color="red")
#' 
#' }
#' 
#==============================================================================#
poppr.msn <- function (pop, distmat, palette = topo.colors, 
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
      E(mst)$color <- gray(adjustcurve(E(mst)$weight, glim=glim, correction=gadj, 
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
    E(mst)$color <- gray(adjustcurve(E(mst)$weight, glim=glim, correction=gadj, 
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
