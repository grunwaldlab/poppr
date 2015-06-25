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

#' @method print ialist
#' @export
print.ialist <- function(x, ...){
  cat("Index\n")
  print(x$index)
  cat("Samples\n")
  if (nrow(x$samples) > 13){
    print(head(x$samples))
    cat("...\n")
    print(tail(x$samples))    
  } else {
    print(x$samples)
  }
}

#' @method plot ialist
#' @export
plot.ialist <- function(x, y = NULL, ..., index = "rbarD"){
  poppr.plot(x$samples, pval = x$index[c(2, 4)], file = substitute(x),
             observed = x$index[c(1, 3)], index = index, ...)
}

#' @method print amova
#' @export
print.amova <- function(x, full = FALSE, ...) 
{
  if (full == TRUE) 
    print(x)
  else print(x[-((length(x) - 2):length(x))])
}

#' @method print popprtable
#' @export
print.popprtable <- function(x, ...){
  call <- list(...)
  if (length(call > 0) && names(call) %in% "digits"){
    print.data.frame(x, ...)
  } else {
    print.data.frame(x, digits = 3, ...)
  }
}

#' @method print locustable
#' @export
print.locustable <- function(x, ...){
  call <- list(...)
  if (length(call > 0) && names(call) %in% "digits"){
    print.table(x, ...)
  } else {
    print.table(x, digits = 2, zero.print = ".", ...)
  }
}

#' @method print pairia
#' @export
print.pairia <- function(x, ...){
  print.locustable(x, ...)
}

#' @method plot pairia
#' @export
plot.pairia <- function(x, ..., index = "rbarD", low = "blue", high = "red"){
  df.index <- x[, index]
  theLoci  <- strsplit(rownames(x), ":")
  lnames   <- unique(unlist(theLoci))
  theTitle <- ifelse(index == "rbarD", 
                     expression(paste(bar(r)[d])), 
                     expression(paste(I[A])))
  L1 <- factor(vapply(theLoci, "[[", character(1), 1), lnames)
  L2 <- factor(vapply(theLoci, "[[", character(1), 2), rev(lnames))
  df <- data.frame(value = df.index, L1 = L1, L2 = L2)
  basic_plot <- ggplot(df, aes_string(x = "L1", y = "L2", fill = "value")) +
    geom_tile()
  basic_plot <- basic_plot + 
    scale_fill_gradient(low = low, high = high, na.value = "white") +
    scale_x_discrete(expand = c(0, -1)) +
    scale_y_discrete(expand = c(0, -1)) +
    theme(axis.title = element_blank(), title = element_text(size = rel(2))) + 
    # theme(legend.title.align = 0.5) + 
    myTheme +
    labs(fill = theTitle)
  print(basic_plot)
}