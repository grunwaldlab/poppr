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
    oops <- ggplot(as.data.frame(list(x=-10:9)), aes(x)) + 
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
    if(any(is.nan(sample$rbarD)))
      sample$rbarD[which(is.nan(sample$rbarD))] <- mean(sample$rbarD, na.rm=TRUE)
    if(any(is.nan(sample$Ia)))
      sample$Ia[which(is.nan(sample$Ia))] <- mean(sample$Ia, na.rm=TRUE)

	  # Transforming the data into a data frame that ggplot can understand and put
	  # into facets. It will have the following rows: Value, Index, and Quant.

	  infodata <- as.data.frame(list(Value=c(sample$Ia, sample$rbarD), 
		  Index=rep(c("Ia","rbarD"), each=length(sample$Ia))))

	  # Setting up the observed values for formatting the overlaying lines.
	  # It will have the following rows: Observed, Index, Quant, Sam, P, min, max

	  obsdata <- data.frame(list(Observed=observed[1:2], Index=c("Ia","rbarD")))
	  obsdata$P <- pval
	  obsdata$median <- c(median(unlist(subset(infodata,Index="Ia",select=Value))), 
      median(unlist(subset(infodata,Index="rbarD",select=Value))))
    if(any(is.na(observed))){
      warning(paste("The Index of Association values from ",file,", population: ",pop," contain missing values and cannot be displayed graphically", sep=""))
	    derp <- ggplot(infodata, aes(Value)) + 
		    # Giving the data over to the histogram creating function and removing
		    # all of the lines from each bar, so it's displayed as a solid area.
		    geom_histogram(linetype="blank", alpha=0.8, 
          data=subset(infodata, Index=="Ia"), 
          position="identity",
          binwidth=diff(range(subset(infodata, Index=="Ia", select=Value)))/30) + 
        geom_histogram(linetype="blank", alpha=0.8, 
          data=subset(infodata, Index=="rbarD"), 
          position="identity",
          binwidth=diff(range(subset(infodata,Index=="rbarD",select=Value)))/30) + 

		    # The label for the observed line is a bit more difficult to code as
		    # it has the ability to appear anywhere on the chart. Here, I'm
		    # forcing it to flip to one side or the other based on which side of
		    # the mean that the observed value falls on.
		    geom_text(aes(label=paste("Observed: ",Observed,sep=""), 
			    x=0,y=Inf,vjust=1.5),
          data=obsdata, angle=0, color="red") + 

		    # Splitting the data into separate plots with free x-axes.
		    facet_grid(.~Index, scales="free_x") +

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
      derp <- ggplot(infodata, aes(Value)) + 
		  # Giving the data over to the histogram creating function and removing
		  # all of the lines from each bar, so it's displayed as a solid area.
		  geom_histogram(linetype="blank", alpha=0.8, 
        data=subset(infodata, Index=="Ia"), 
        position="identity",
        binwidth=diff(range(subset(infodata, Index=="Ia", select=Value)))/30) + 
      geom_histogram(linetype="blank", alpha=0.8, 
        data=subset(infodata, Index=="rbarD"), 
        position="identity",
        binwidth=diff(range(subset(infodata,Index=="rbarD",select=Value)))/30) + 
		  # Positioning the observed line and labeling it.
		  geom_vline(aes(xintercept=Observed), data=obsdata, color="blue", 
			  show_guide=TRUE, linetype="dashed") +

		  # The label for the observed line is a bit more difficult to code as
		  # it has the ability to appear anywhere on the chart. Here, I'm
		  # forcing it to flip to one side or the other based on which side of
		  # the mean that the observed value falls on.
		  geom_text(aes(label=paste("Observed \n(p-value: ", P,")",sep=""), 
			  x=Observed,y=Inf,vjust=2,hjust=ifelse(Observed > median, 1.01, -0.01)),
        data=obsdata, angle=0, color="blue") + 

		  # Splitting the data into separate plots with free x-axes.
		  facet_grid(.~Index, scales="free_x") +

		  # Title of the plot.
		  labs(title=paste("Population: ", pop, "; N: ", N, "\nPermutations: ", 
			  length(sample$Ia), "\nFile: ", file, sep="")) + 
      theme_classic() %+replace%
    	theme(plot.title = element_text(vjust=1, size=rel(2), face="bold")) +
		  theme(panel.background = element_rect(fill="grey98")) +

		  # Making the Index titles bigger. 
		  theme(strip.text.x = element_text(size=rel(3), face="bold"))
    }
	  print(derp)
  }
} 

