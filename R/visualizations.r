
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
        geom_rug(alpha = 0.5) + #theme(legend.position = "none") + 
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
        geom_rug(alpha = 0.5) + #theme(legend.position = "none") + 
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
#'   color the populations on the graph. It defaults to
#'   \code{\link{topo.colors}}, but you can easily create new schemes by using
#'   \code{\link{colorRampPalette}} (see examples for details)
#'   
#' @param sublist a \code{vector} of population names or indexes that the user 
#'   wishes to keep. Default to "ALL".
#'   
#' @param blacklist a \code{vector} of population names or indexes that the user
#'   wishes to discard. Default to \code{NULL}
#'   
#' @param vertex.label a \code{vector} of characters to label each vertex. There
#'   are two defaults: \code{"MLG"} will label the nodes with the multilocus
#'   genotype from the original data set and \code{"inds"} will label the nodes
#'   with the representative individual names.
#'   
#' @param gscale "grey scale". If this is \code{TRUE}, this will scale the color
#'   of the edges proportional to the observed distance, with the lines becoming
#'   darker for more related nodes. See \code{\link{greycurve}} for details.
#'   
#' @param glim "grey limit". Two numbers between zero and one. They determine 
#'   the upper and lower limits for the \code{\link{gray}} function. Default is
#'   0 (black) and 0.8 (20\% black). See \code{\link{greycurve}} for details.
#'   
#' @param gadj "grey adjust". a positive \code{integer} greater than zero that 
#'   will serve as the exponent to the edge weight to scale the grey value to 
#'   represent that weight. See \code{\link{greycurve}} for details.
#'   
#' @param gweight "grey weight". an \code{integer}. If it's 1, the grey scale 
#'   will be weighted to emphasize the differences between closely related
#'   nodes. If it is 2, the grey scale will be weighted to emphasize the
#'   differences between more distantly related nodes. See
#'   \code{\link{greycurve}} for details.
#'   
#' @param wscale "width scale". If this is \code{TRUE}, the edge widths will be 
#'   scaled proportional to the inverse of the observed distance , with the
#'   lines becoming thicker for more related nodes.
#'   
#' @param showplot logical. If \code{TRUE}, the graph will be plotted. If
#'   \code{FALSE}, it will simply be returned.
#'   
#' @param include.ties logical. If \code{TRUE}, the graph will include all
#'   edges that were arbitrarily passed over in favor of another edge of equal
#'   weight. If \code{FALSE}, which is the default, one edge will be arbitrarily 
#'   selected when two or more edges are tied, resulting in a pure minimum spanning
#'   network. 
#'
#' @param threshold numeric. If greater than the default value of 0.0, this will
#'   be passed to \code{\link{mlg.filter}} prior to creating the msn.
#'
#' @param clustering.algorithm string. If \code{threshold} is greater than 0, this
#'   this will also be passed to \code{\link{mlg.filter}} prior to creating the msn.
#'   For both of these arguments, see \code{\link{mlg.filter}} for more details.
#'
#' @param ... any other arguments that could go into plot.igraph
#'   
#' @return \item{graph}{a minimum spanning network with nodes corresponding to 
#'   MLGs within the data set. Colors of the nodes represent population 
#'   membership. Width and color of the edges represent distance.} 
#'   \item{populations}{a vector of the population names corresponding to the 
#'   vertex colors} \item{colors}{a vector of the hexadecimal representations of
#'   the colors used in the vertex colors}
#'   
#' @note The edges of these graphs may cross each other if the graph becomes too
#'   large.
#'   
#' @seealso \code{\link{nancycats}}, \code{\link{upgma}}, \code{\link{nj}}, 
#'   \code{\link{nodelabels}}, \code{\link{na.replace}},
#'   \code{\link{missingno}}, \code{\link{bruvo.msn}}, \code{\link{greycurve}}.
#' 
#' @export
#' @aliases msn.poppr
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
#' plot.igraph(micro.msn$graph, edge.label = ifelse(E(micro.msn$graph)$weight < 0.15, 
#' round(E(micro.msn$graph)$weight, 3), NA), vertex.size=2, edge.label.color="red")
#' 
#' }
#' 
#==============================================================================#
poppr.msn <- function (pop, distmat, palette = topo.colors, 
                       sublist = "All", blacklist = NULL, vertex.label = "MLG", 
                       gscale=TRUE, glim = c(0,0.8), gadj = 3, gweight = 1, 
                       wscale=TRUE, showplot = TRUE, 
                       include.ties = FALSE, threshold = 0.0, clustering.algorithm="farthest_neighbor", ...){
  if (class(distmat) != "dist"){
    if (is.matrix(distmat)){
      if (any(nInd(pop) != dim(distmat))){
        stop("The size of the distance matrix does not match the size of the data.\n")
      }
      distmat <- as.dist(distmat)
    } else {
      stop("The distance matrix is neither a dist object nor a matrix.\n")
    }
  }
  if (nInd(pop) != attr(distmat, "Size")){
    stop("The size of the distance matrix does not match the size of the data.\n")
  }
  if(!is.genclone(pop)){
    pop <- as.genclone(pop)
  }
  gadj <- ifelse(gweight == 1, gadj, -gadj)
  # Updating MLG with filtered data
  if(threshold > 0){
    filter.stats <- mlg.filter(pop,threshold,distance=distmat,algorithm=clustering.algorithm,stats="ALL")
    pop$mlg <- filter.stats[[1]]
    bclone <- filter.stats[[3]]
  } else {
    bclone <- as.matrix(distmat)
  }

  # The clone correction of the matrix needs to be done at this step if there
  # is only one or no populations. 
  if (is.null(pop(pop)) | length(pop@pop.names) == 1){
    mlgs <- pop@mlg
    if(threshold > 0){
      cpop <- pop[if(length(-which(duplicated(pop$mlg))==0)) which(!duplicated(pop$mlg)) else -which(duplicated(pop$mlg)) ,]
    }
    else {
      cpop <- pop[.clonecorrector(pop), ]
      # This will clone correct the incoming matrix. 
      bclone <- bclone[!duplicated(mlgs), !duplicated(mlgs)]
    }
    rownames(bclone) <- cpop$ind.names
    colnames(bclone) <- cpop$ind.names
    return(singlepop_msn(pop, vertex.label, distmat = bclone, gscale = gscale, 
                         glim = glim, gadj = gadj, wscale = wscale, 
                         palette = palette, include.ties = include.ties,
                         threshold=threshold, clustering.algorithm=clustering.algorithm, ...))
  }
  # This will subset both the population and the matrix. 
  if(sublist[1] != "ALL" | !is.null(blacklist)){
    sublist_blacklist <- sub_index(pop, sublist, blacklist)
    bclone <- bclone[sublist_blacklist, sublist_blacklist]
    pop    <- popsub(pop, sublist, blacklist)
  }
  mlgs <- pop@mlg
  if(threshold > 0){
    cpop <- pop[if(length(-which(duplicated(pop$mlg))==0)) which(!duplicated(pop$mlg)) else -which(duplicated(pop$mlg)) ,]
  }
  else {
    cpop <- pop[.clonecorrector(pop), ]
    # This will clone correct the incoming matrix. 
    bclone <- bclone[!duplicated(mlgs), !duplicated(mlgs)]
  }
  rownames(bclone) <- cpop$ind.names
  colnames(bclone) <- cpop$ind.names
  cmlg <- cpop@mlg
  
  if (is.null(pop(pop)) | length(pop@pop.names) == 1){
    return(singlepop_msn(pop, vertex.label, distmat = bclone, gscale = gscale, 
                         glim = glim, gadj = gadj, wscale = wscale, 
                         palette = palette, showplot = showplot, include.ties = include.ties,
                         threshold=threshold, clustering.algorithm=clustering.algorithm, ...))
  }
  # Obtaining population information for all MLGs
  subs <- sort(unique(mlgs))
  mlg.cp <- mlg.crosspop(pop, mlgsub = subs, quiet=TRUE)

  names(mlg.cp) <- paste0("MLG.", sort(unique(mlgs)))
  # This will determine the size of the nodes based on the number of individuals
  # in the MLG. Sub-setting by the MLG vector of the clone corrected set will
  # give us the numbers and the population information in the correct order.
  # Note: rank is used to correctly subset the data
  mlg.number <- table(mlgs)[rank(cmlg)]
  mlg.cp     <- mlg.cp[rank(cmlg)]
  bclone <- as.matrix(bclone)
  
  g   <- graph.adjacency(bclone, weighted=TRUE, mode="undirected")
  if(length(cpop@mlg) > 1){
    mst <- minimum.spanning.tree(g, algorithm="prim", weights=E(g)$weight)
    # Add any relevant edges that were cut from the mst while still being tied for the title of optimal edge
    if(include.ties){
      tied_edges <- .Call("msn_tied_edges",as.matrix(mst[]),as.matrix(bclone),(.Machine$double.eps ^ 0.5))
      if(length(tied_edges) > 0){
        mst <- add.edges(mst, dimnames(mst[])[[1]][tied_edges[c(TRUE,TRUE,FALSE)]], weight=tied_edges[c(FALSE,FALSE,TRUE)])
      }
    }
  } else {
    mst <- minimum.spanning.tree(g)
  }
  
  if (!is.na(vertex.label[1]) & length(vertex.label) == 1){
    if(toupper(vertex.label) == "MLG"){
      vertex.label <- paste("MLG.", cmlg, sep="")
    }
    else if(toupper(vertex.label) == "INDS"){
      vertex.label <- cpop$ind.names
    }
  }
  ###### Color schemes #######  
  # The pallete is determined by what the user types in the argument. It can be 
  # rainbow, topo.colors, heat.colors ...etc.
  palette <- match.fun(palette)
  color   <- palette(length(pop@pop.names))
  
  if(length(cpop@mlg) > 1){
    ###### Edge adjustments ######
    mst <- update_edge_scales(mst, wscale, gscale, glim, gadj)
  }
 
  # This creates a list of colors corresponding to populations.
  mlg.color <- lapply(mlg.cp, function(x) color[pop@pop.names %in% names(x)])
  if (showplot){
    plot.igraph(mst, edge.width = E(mst)$width, edge.color = E(mst)$color, 
         vertex.size = mlg.number*3, vertex.shape = "pie", vertex.pie = mlg.cp, 
         vertex.pie.color = mlg.color, vertex.label = vertex.label, ...)
    legend(-1.55 ,1 ,bty = "n", cex = 0.75, legend = pop$pop.names, 
           title = "Populations", fill=color, border=NULL)
  }
  V(mst)$size      <- mlg.number
  V(mst)$shape     <- "pie"
  V(mst)$pie       <- mlg.cp
  V(mst)$pie.color <- mlg.color
  V(mst)$label     <- vertex.label
  return(list(graph = mst, populations = pop$pop.names, colors = color))
}


#==============================================================================#
#' Create a table summarizing missing data or ploidy information of a genind or
#' genclone object
#' 
#' @param gen a \linkS4class{genind} or \linkS4class{genclone} object.
#'   
#' @param type \code{character}. What information should be returned. Choices
#'   are "missing" (Default) and "ploidy". See Description.
#'   
#' @param percent \code{logical}. (ONLY FOR \code{type = 'missing'}) If
#'   \code{TRUE} (default), table and plot will represent missing data as a
#'   percentage of each cell. If \code{FALSE}, the table and plot will represent
#'   missing data as raw counts. (See details)
#'   
#' @param plot \code{logical}. If \code{TRUE}, a simple heatmap will be 
#'   produced. If \code{FALSE} (default), no heatmap will be produced.
#'   
#' @param df \code{logical}. If \code{TRUE}, the data will be returned as a long
#'   form data frame. If \code{FALSE} (default), a matrix with samples in rows
#'   and loci in columns will be returned.
#'   
#' @param returnplot \code{logical}. If \code{TRUE}, a list is returned with two
#'   elements: \code{table} - the normal output and \code{plot} - the ggplot 
#'   object. If \code{FALSE}, the table is returned.
#'   
#' @param low \code{character}. What color should represent no missing data or 
#'   lowest observed ploidy? (default: "blue")
#'   
#' @param high \code{character}. What color should represent the highest amount 
#'   of missing data or observed ploidy? (default: "red")
#'   
#' @param plotlab \code{logical}. (ONLY FOR \code{type = 'missing'}) If
#'   \code{TRUE} (default), values of missing data greater than 0\% will be
#'   plotted. If \code{FALSE}, the plot will appear unappended.
#'   
#' @param scaled \code{logical}. (ONLY FOR \code{type = 'missing'}) This is for
#'   when \code{percent = TRUE}. If \code{TRUE} (default), the color specified
#'   in \code{high} will represent the highest observed value of missing data.
#'   If \code{FALSE}, the color specified in \code{high} will represent 100\%.
#'   
#' @return a matrix, data frame (\code{df = TRUE}), or a list (\code{returnplot 
#'   = TRUE}) representing missing data per population (\code{type = 'missing'})
#'   or ploidy per individual (\code{type = 'ploidy'}) in a \linkS4class{genind}
#'   or \linkS4class{genclone} object.
#' 
#' @details 
#'   Missing data is accounted for on a per-population level.\cr
#'   Ploidy is accounted for on a per-individual level.
#'   
#'   \subsection{For type = 'missing'}{
#'   This data is potentially useful for identifying areas of systematic missing
#'   data. There are a few caveats to be aware of. \itemize{ \item
#'   \strong{Regarding counts of missing data}: Each count represents the number
#'   of individuals with missing data at each locus. The last column, "mean" can
#'   be thought of as the average number of individuals with missing data per
#'   locus. \item \strong{Regarding percentage missing data}: This percentage is
#'   \strong{relative to the population and locus}, not ot the enitre data set.
#'   The last colum, "mean" represents the average percent of the population
#'   with missing data per locus. }} 
#'   \subsection{For type = 'ploidy'}{
#'   This option is useful for data that has been imported with mixed ploidies.
#'   It will summarize the relative levels of ploidy per individual per locus.
#'   This is simply based off of observed alleles and does not provide any
#'   further estimates.}
#'   
#' @export
#' @keywords missing ploidy
#' @author Zhian N. Kamvar
#' @examples
#' data(nancycats)
#' nancy.miss <- info_table(nancycats, plot = TRUE, type = "missing")
#' data(Pinf)
#' Pinf.ploid <- info_table(Pinf, plot = TRUE, type = "ploidy")
#' 
#==============================================================================#
info_table <- function(gen, type = c("missing", "ploidy"), percent = TRUE, plot = FALSE, 
                       df = FALSE, returnplot = FALSE, low = "blue", 
                       high = "red", plotlab = TRUE, scaled = TRUE){
  datalabel <- as.character(match.call()[2])
  ARGS      <- c("missing", "ploidy")
  type      <- match.arg(type, ARGS)

  if (type == "missing"){

    valname    <- "Missing"
    pops       <- seppop(gen, drop = FALSE)
    pops$Total <- gen
    inds       <- 1
    if (percent){
      inds <- c(table(pop(gen)), nInd(gen))
    }
    data_table <- matrix(0, nrow = nLoc(gen) + 1, ncol = length(pops))
    data_table[1:nLoc(gen), ] <- vapply(pops, number_missing_locus, numeric(nLoc(gen)), 1)
    data_table[-nrow(data_table), ] <- t(apply(data_table[-nrow(data_table), ], 1, "/", inds))
    data_table[nrow(data_table), ]  <- colMeans(data_table[-nrow(data_table), ])
    rownames(data_table)         <- c(gen@loc.names, "Mean")
    colnames(data_table)         <- names(pops)
    dimnames(data_table) <- list(Locus = c(gen@loc.names, "Mean"), Population = names(pops))
    if (all(data_table == 0)){
      cat("No Missing Data Found!")
      return(NULL)
    }
    if (plot){
      data_df      <- melt(data_table, value.name = valname)
      leg_title    <- valname
      data_df[1:2] <- data.frame(lapply(data_df[1:2], 
                                       function(x) factor(x, levels = unique(x))))
      plotdf <- textdf <- data_df
      if (percent) {
        plotdf$Missing <- round(plotdf$Missing*100, 2)
        textdf$Missing <- paste(plotdf$Missing, "%")
        miss           <- "0 %"
        title          <- paste("Percent missing data per locus and population of", 
                                datalabel)
        leg_title      <- paste("Percent", leg_title)
        if(!scaled){
          lims <- c(0, 100)
        }
      } else {
        textdf$Missing <- round(textdf$Missing, 2)
        miss           <- 0
        title          <- paste("Missing data per locus and population of", 
                                datalabel)
      }
      if (scaled | !percent){
        lims <- c(0, max(plotdf$Missing))
      }
      linedata <- data.frame(list(yint = 1.5, xint = nrow(data_table) - 0.5))
      textdf$Missing <- ifelse(textdf$Missing == miss, "", textdf$Missing)
      plotdf$Missing[plotdf$Locus == "Mean" & plotdf$Population == "Total"] <- NA

      outplot <- ggplot(plotdf, aes_string(x = "Locus", y = "Population")) + 
        geom_tile(aes_string(fill = valname)) +
        labs(list(title = title, x = "Locus", y = "Population")) +
        labs(fill = leg_title) + 
        scale_fill_gradient(low = low, high = high, na.value = "white", 
                            limits = lims) +
        geom_hline(aes_string(yintercept = "yint"), data = linedata) + 
        geom_vline(aes_string(xintercept = "xint"), data = linedata) 
      if (plotlab){
        outplot <- outplot + geom_text(aes_string(label = valname), 
                                       data = textdf)
      }
      outplot <- outplot +
        theme_classic() + 
        scale_x_discrete(expand = c(0, -1)) + 
        scale_y_discrete(expand = c(0, -1), 
                         limits = rev(unique(plotdf$Population))) + 
        theme(axis.text.x = element_text(size = 10, angle = -45, 
                                         hjust = 0, vjust = 1)) 
      print(outplot)
    }

  } else if (type == "ploidy"){

    valname <- "Observed_Ploidy"
    if (gen@ploidy <= 2){
      warning("This function is meant for polyploid data.")
      data_table <- matrix(gen@ploidy, nrow = nInd(gen), ncol = nLoc(gen))
      missing <- propTyped(gen, "both") == 0
      data_table[missing] <- NA
    } else {
      data_table <- get_local_ploidy(gen)
    }
    dimnames(data_table) <- list(Samples = indNames(gen), Loci = locNames(gen))
    if (plot){
      data_df <- melt(data_table, value.name = valname)
      data_df[1:2] <- data.frame(lapply(data_df[1:2], 
                                       function(x) factor(x, levels = unique(x))))
      vars <- aes_string(x = "Loci", y = "Samples", fill = valname)

      mytheme <- theme_classic() +  
                 theme(axis.text.x = element_text(size = 10, angle = -45, 
                                                  hjust = 0, vjust = 1)) 

      title <- paste("Observed ploidy of", datalabel)
      outplot <- ggplot(data_df) + geom_tile(vars) + 
                   scale_fill_gradient(low = low, high = high) + 
                   scale_x_discrete(expand = c(0, -1)) + 
                   scale_y_discrete(expand = c(0, -1), 
                                    limits = rev(unique(data_df$Samples))) + 
                   labs(list(title = title, x = "Locus", y = "Sample", 
                             fill = "Observed\nPloidy")) +
                   mytheme 

      print(outplot)
    }
  } 
  if (df){
    if(!exists("data_df")){
      data_df <- melt(data_table, value.name = valname)
    }
    data_table <- data_df
  } else {
    if (type == "missing"){
      data_table <- t(data_table)
    }
    class(data_table) <- c("locustable", "matrix")
  }
  if (returnplot & exists("outplot")){
    data_table <- list(table = data_table, plot = outplot)
  }
  return(data_table)
}

#==============================================================================#
#' Display a greyscale gradient adjusted to specific parameters
#' 
#' This function has one purpose. It is for deciding the appropriate scaling for
#' a grey palette to be used for edge weights of a minimum spanning network.
#'
#'
#' @param data a sequence of numbers to be converted to greyscale.
#' 
#' @param glim "grey limit". Two numbers between zero and one. They determine 
#' the upper and lower limits for the \code{\link{gray}} function. Default is 0
#' (black) and 0.8 (20\% black). 
#' 
#' @param gadj "grey adjust". a positive \code{integer} greater than zero that
#' will serve as the exponent to the edge weight to scale the grey value to
#' represent that weight.
#' 
#' @param gweight "grey weight". an \code{integer}. If it's 1, the grey scale
#' will be weighted to emphasize the differences between closely related nodes.
#' If it is 2, the grey scale will be weighted to emphasize the differences
#' between more distantly related nodes. 
#'
#' @param scalebar When this is set to \code{TRUE}, two scalebars will be
#' plotted. The purpose of this is for adding a scale bar to minimum spanning
#' networks produced in earlier versions of poppr. 
#' 
#' @return A plot displaying a grey gradient from 0.001 to 1 with minimum and 
#' maximum values displayed as yellow lines, and an equation for the correction 
#' displayed in red. 
#' 
#' @author Zhian N. Kamvar
#'
#' @examples
#' # Normal grey curve with an adjustment of 3, an upper limit of 0.8, and
#' # weighted towards smaller values.
#' greycurve()
#' \dontrun{
#' # 1:1 relationship grey curve.
#' greycurve(gadj=1, glim=1:0)
#' 
#' # Grey curve weighted towards larger values.
#' greycurve(gweight=2)
#' 
#' # Same as the first, but the limit is 1.
#' greycurve(glim=1:0)
#' 
#' # Setting the lower limit to 0.1 and weighting towards larger values.
#' greycurve(glim=c(0.1,0.8), gweight=2)
#' }
#' @export
#==============================================================================#
greycurve <- function(data = seq(0, 1, length = 1000), glim = c(0,0.8), 
                      gadj = 3, gweight = 1, scalebar = FALSE){
  gadj <- ifelse(gweight == 1, gadj, -gadj)
  adjustcurve(data, glim, correction=gadj, show = TRUE, scalebar = scalebar)
}


#==============================================================================#
#' Plot minimum spanning networks produced in poppr.
#' 
#' This function allows you to take the output of poppr.msn and bruvo.msn and 
#' customize the plot by labeling groups of individuals, size of nodes, and 
#' adjusting the palette and scale bar.
#' 
#' @param x a \code{\linkS4class{genind}} or \code{\linkS4class{genclone}}
#'   object from which \code{poppr_msn} was derived.
#'   
#' @param poppr_msn a \code{list} produced from either \code{\link{poppr.msn}}
#'   or \code{\link{bruvo.msn}}. This list should contain a graph, a vector of
#'   population names and a vector of hexadecimal color definitions for each
#'   popualtion.
#'   
#' @inheritParams greycurve
#' 
#' @inheritParams poppr.msn
#'   
#' @param inds a character vector indicating which individual names to label
#'   nodes with. See details.
#'   
#' @param quantiles \code{logical}. When set to \code{TRUE} (default), the scale
#'   bar will be composed of the quantiles from the observed edge weights. When
#'   set to \code{FALSE}, the scale bar will be composed of a smooth gradient
#'   from the minimum edge weight to the maximum edge weight.
#'   
#' @param nodelab an \code{integer} specifying the smallest size of node to
#'   label. See details.
#'   
#' @param cutoff a number indicating the longest distance to display in your
#'   graph. This is performed by removing edges with weights greater than this
#'   number.
#'   
#' @param palette a function or character corresponding to a specific palette 
#'   you want to use to delimit your populations. The default is whatever
#'   palette was used to produce the original graph.
#'   
#' @param layfun a function specifying the layout of nodes in your graph. It 
#'   defaults to \code{\link[igraph]{layout.auto}}.
#'   
#' @param beforecut if \code{TRUE}, the layout of the graph will be computed
#'   before any edges are removed with \code{cutoff}. If \code{FALSE} (Default),
#'   the layout will be computed after any edges are removed.
#'   
#' @param ... any other parameters to be passed on to
#'   \code{\link[igraph]{plot.igraph}}.
#'   
#' @details The previous incarnation of msn plotting in poppr simply plotted the
#'   minimum spanning network with the legend of populations, but did not 
#'   provide a scale bar and it did not provide the user a simple way of 
#'   manipulating the layout or labels. This function allows the user to 
#'   manipulate many facets of graph creation, making the creation of minimum 
#'   spanning networks ever so slightly more user friendly. Note that this 
#'   function will only plot individual names, not MLG names since the naming 
#'   convention for those are arbitrary. 
#'   
#'   This function must have both the source data and the output msn to work. 
#'   The source data must contain the same population structure as the graph. 
#'   Every other parameter has a default setting.
#'   
#'   \subsection{Parameter details}{ \itemize{ \item \code{inds} This will take
#'   in the name of a query individual in your data set and will use that to
#'   query any other individuals that share multilocus genotypes and label their
#'   node on the graph. The default is to label all the nodes, but you can set
#'   it to a name that doesn't exist to label none of the nodes. \item
#'   \code{nodelab} If a node is not labeled by individual, this will label the
#'   size of the nodes greater than or equal to this value. If you don't want to
#'   label the size of the nodes, simply set this to a very high number. \item
#'   \code{cutoff} This is useful for when you want to investigate groups of
#'   multilocus genotypes separated by a specific distance or if you have two
#'   distinct populations and you want to physically separate them in your
#'   network. \item \code{beforecut} This is an indicator useful if you want to
#'   maintain the same position of the nodes before and after removing edges
#'   with the \code{cutoff} argument. This works best if you set a seed before
#'   you run the function.}}
#' 
#' @seealso \code{\link[igraph]{layout.auto}} \code{\link[igraph]{plot.igraph}}
#' \code{\link{poppr.msn}} \code{\link{bruvo.msn}} \code{\link{greycurve}}
#' \code{\link[igraph]{delete.edges}} \code{\link{palette}}
#' 
#' @author Zhian N. Kamvar
#' @export
#' 
#' @examples
#' # Using a data set of the Aphanomyces eutieches root rot pathogen.
#' data(Aeut)
#' adist <- diss.dist(Aeut, percent = TRUE)
#' amsn <- poppr.msn(Aeut, adist, showplot = FALSE)
#' 
#' # Default
#' library(igraph) # To get all the layouts.
#' set.seed(500)
#' plot_poppr_msn(Aeut, amsn, gadj = 15, beforecut = TRUE)
#' 
#' # Removing link between populations and labelling no individuals
#' set.seed(500)
#' plot_poppr_msn(Aeut, amsn, inds = "none", gadj = 15, beforecut = TRUE, cutoff = 0.2)
#' 
#' # Labelling individual #57 because it is an MLG that crosses popualtions
#' # Showing clusters of MLGS with at most 5% variation
#' # Notice that the Mt. Vernon population appears to be more clonal
#' set.seed(50) 
#' plot_poppr_msn(Aeut, amsn, gadj = 15, cutoff = 0.05, inds = "57")
#' 
#' 
#' \dontrun{
#' data(partial_clone)
#' pcmsn <- bruvo.msn(partial_clone, replen = rep(1, 10))
#' plot_poppr_msn(partial_clone, pcmsn, palette = rainbow, inds = "sim 20")
#' 
#' # Something pretty
#' data(microbov)
#' mdist <- diss.dist(microbov, percent = TRUE)
#' micmsn <- poppr.msn(microbov, mdist, showplot = FALSE)
#' 
#' plot_poppr_msn(microbov, micmsn, palette = "terrain.colors", inds = "n", 
#'   quantiles = FALSE)
#' plot_poppr_msn(microbov, micmsn, palette = "terrain.colors", inds = "n", 
#'   cutoff = 0.3, quantiles = FALSE)
#' }
#==============================================================================#
#' @importFrom igraph layout.auto delete.edges
plot_poppr_msn <- function(x, poppr_msn, gscale = TRUE, gadj = 3, glim = c(0, 0.8),
                           gweight = 1, wscale = TRUE, inds = "ALL", quantiles = TRUE, 
                           nodelab = 2, cutoff = NULL, palette = NULL,
                           layfun = layout.auto, beforecut = FALSE, ...){
  if (!is.genind(x)){
    stop(paste(substitute(x), "is not a genind or genclone object."))
  }
  if (!identical(names(poppr_msn), c("graph", "populations", "colors"))){
    stop("graph not compatible")
  }
  if (!is.null(palette)){
    poppr_msn <- update_poppr_graph(poppr_msn, palette)
  }
  # Making sure incoming data matches so that the individual names match.
  x <- popsub(x, sublist = poppr_msn$populations)
  
  if (beforecut){
    LAYFUN <- match.fun(layfun)
    lay <- LAYFUN(poppr_msn$graph)
  } else {
    lay <- match.fun(layfun)
  }
  # delete.edges      <- match.fun(igraph::delete.edges)
  if (!is.null(cutoff) && !is.na(cutoff)){
    if (all(cutoff < E(poppr_msn$graph)$weight)){
      msg <- paste0("Cutoff value (", cutoff, ") is below the minimum observed",
                    " distance. Edges will not be removed.")
      warning(msg)
    } else {
      E_above_cutoff  <- E(poppr_msn$graph)[E(poppr_msn$graph)$weight >= cutoff]
      poppr_msn$graph <- delete.edges(poppr_msn$graph, E_above_cutoff)
    }
  }
  # Adjusting color scales. This will replace any previous scaling contained in
  # poppr_msn.
  weights <- E(poppr_msn$graph)$weight
  wmin    <- min(weights)
  wmax    <- max(weights)
  gadj    <- ifelse(gweight == 1, gadj, -gadj)
  poppr_msn$graph <- update_edge_scales(poppr_msn$graph, wscale, gscale, glim, gadj)
  
  # Highlighting only the names of the submitted genotypes and the matching
  # isolates.
  x.mlg <- mlg.vector(x)
  labs  <- unique(x.mlg)
  # The labels in the graph are organized by MLG, so we will use that to extract
  # the names we need.
  if (length(inds) == 1 & toupper(inds[1]) == "ALL"){
    x.input <- unique(x.mlg)
  } else {
    x.input <- unique(x.mlg[x@ind.names %in% inds])
  }
  # Combine all the names that match with each particular MLG in x.input.
  combined_names <- vapply(x.input, function(mlgname)
                           paste(rev(x@ind.names[x.mlg == mlgname]),
                                 collapse = "\n"),
                           character(1))
  # Remove labels that are not specified.
  labs[which(!labs %in% x.input)] <- NA
  labs[!is.na(labs)] <- combined_names
  if (any(is.na(labs))){
    sizelabs <- V(poppr_msn$graph)$size
    sizelabs <- ifelse(sizelabs >= nodelab, sizelabs, NA)
    labs     <- ifelse(is.na(labs), sizelabs, labs)
  }
  # Change the size of the vertices to a log scale.
  vsize <- log(V(poppr_msn$graph)$size, base = 1.15) + 3
  
  # Plotting parameters.
  def.par <- par(no.readonly = TRUE)
  # Setting up the matrix for plotting. One Vertical panel of width 1 and height
  # 5 for the legend, one rectangular panel of width 4 and height 4.5 for the
  # graph, and one horizontal panel of width 4 and height 0.5 for the greyscale.
  layout(matrix(c(1,2,1,3), ncol = 2, byrow = TRUE),
         widths = c(1, 4), heights= c(4.5, 0.5))
  # mar = bottom left top right
  
  ## LEGEND
  par(mar = c(0, 0, 1, 0) + 0.5)
  too_many_pops   <- as.integer(ceiling(length(x$pop.names)/30))
  pops_correction <- ifelse(too_many_pops > 1, -1, 1)
  yintersperse    <- ifelse(too_many_pops > 1, 0.51, 0.62)
  plot(c(0, 2), c(0, 1), type = 'n', axes = F, xlab = '', ylab = '',
       main = 'POPULATION')
  legend("topleft", bty = "n", cex = 1.2^pops_correction,
         legend = poppr_msn$populations, fill = poppr_msn$color, border = NULL,
         ncol = too_many_pops, x.intersp = 0.45, y.intersp = yintersperse)
  
  ## PLOT
  par(mar = c(0,0,0,0))
  plot.igraph(poppr_msn$graph, vertex.label = labs, vertex.size = vsize, 
              layout = lay, ...)
  
  ## SCALE BAR
  if (quantiles){
    scales <- sort(weights)
  } else {
    scales <- seq(wmin, wmax, l = 1000)
  }
  greyscales <- gray(adjustcurve(scales, show=FALSE, glim=glim, correction=gadj))
  legend_image <- as.raster(matrix(greyscales, nrow=1))
  par(mar = c(0, 1, 0, 1) + 0.5)
  plot.new()
  rasterImage(legend_image, 0, 0.5, 1, 1)
  polygon(c(0, 1, 1), c(0.5, 0.5, 0.8), col = "white", border = "white", lwd = 2)
  axis(3, at = c(0, 0.25, 0.5, 0.75, 1), labels = round(quantile(scales), 3))
  text(0.5, 0, labels = "DISTANCE", font = 2, cex = 1.5, adj = c(0.5, 0))
  
  # Return top level plot to defaults.
  layout(matrix(c(1), ncol=1, byrow=T))
  par(mar=c(5,4,4,2) + 0.1) # number of lines of margin specified.
  par(oma=c(0,0,0,0)) # Figure margins
}


#==============================================================================#
#' Produce a genotype accumulation curve
#' 
#' GA curves are useful for determinining the minimum number of loci necessary 
#' to discriminate between individuals in a population. This function will
#' randomly sample loci without replacement and count the number of multilocus
#' genotypes observed.
#' 
#' @param gen a \code{\linkS4class{genclone}} or \code{\linkS4class{genind}}
#'   object.
#'   
#' @param sample an \code{integer} defining the number of times loci will be
#'   resampled.
#'   
#' @param quiet if \code{FALSE}, a progress bar will be displayed. If
#'   \code{TRUE}, nothing is printed to screen as the function runs.
#'   
#' @param thresh a number from 0 to 1. This will draw a line at this fraction of
#'   multilocus genotypes.
#'   
#' @return a matrix of integers showing the results of each randomization.
#'   Columns represent the number of loci sampled and rows represent an
#'   independent sample.
#'   
#' @author Zhian N. Kamvar
#' @export
#' @examples
#' data(nancycats)
#' nan_geno <- genotype_curve(nancycats)
#' \dontrun{
#' # With AFLP data, it is often necessary to include more markers for resolution
#' data(Aeut)
#' Ageno <- genotype_curve(Aeut)
#' 
#' # Many microsatellite data sets have hypervariable markers
#' data(microbov)
#' mgeno <- geotype_curve(microbov)
#' 
#' # This data set has been pre filtered
#' data(monpop)
#' mongeno <- genotype_curve(monpop)}
#==============================================================================#

genotype_curve <- function(gen, sample = 100, quiet = FALSE, thresh = 0.9){
  datacall <- match.call()
  if (!is.genind(gen)){
    stop(paste(datacall[2], "must be a genind object"))
  }
  genloc <- as.loci(gen)
  if (!is.null(pop(gen))){
    genloc <- genloc[-1]
  }
  nloci  <- nLoc(gen)
  if (!quiet){
    cat("Calculating genotype accumulation for", nloci - 1, "loci...\n")
    progbar <- txtProgressBar(style = 3)
  } else {
    progbar <- NULL
  }
  out <- vapply(1:(nloci-1), get_sample_mlg, integer(sample), sample, nloci, genloc, progbar)
  colnames(out) <- 1:(nloci-1)
  threshdf <- data.frame(x = mlg(gen, quiet = TRUE)*thresh)
  outmelt <- melt(out, value.name = "MLG", varnames = c("sample", "NumLoci"))
  aesthetics <- aes_string(x = "factor(NumLoci)", y = "MLG", group = "NumLoci")
  outplot <- ggplot(outmelt) + geom_boxplot(aesthetics) + 
             labs(list(title = paste("Genotype accumulation curve for", datacall[2]), 
                       y = "Number of multilocus genotypes",
                       x = "Number of loci sampled")) 
  if (!is.null(thresh)){
    outbreaks <- sort(c(seq(min(out), max(out), length.out = 5), threshdf$x))
    outplot <- outplot + geom_hline(aes_string(yintercept = "x"), 
                                    data = threshdf, color = "red", linetype = 2) + 
                         annotate("text", x = 1, y = threshdf$x, vjust = -1, 
                                  label = paste0(thresh*100, "%"), 
                                  color = "red", hjust = 0) +
                         scale_y_continuous(breaks = outbreaks)
  }
  print(outplot)
  return(out)
}
