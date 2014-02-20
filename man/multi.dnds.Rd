\name{multi.dnds}
\alias{multi.dnds}
\title{multi.dnds}
\usage{
multi.dnds(x, quiet = FALSE)
}
\arguments{
  \item{x}{a \code{list} of multiple alignment \code{fasta}
  files (preferably full path, fasta|fas extension).}

  \item{quiet}{a \code{boolean} for verbosity.}
}
\value{
a table with the dN/dS calculation for all the files in the
list.
}
\description{
dN/dS (as Ka/Ks) wrapper for multiple sets of FASTA files
and the \code{\link{dn.ds}} function
}
\note{
To create the list of files use the
\code{\link{list.files}} function. (See examples for help)
}
\examples{
# Create the list of FASTA files
#list.pairs <- list.files(path="fasta_file_folder/", pattern="*.fasta$", full.names=TRUE)
#
# Run the function
# multi.dnds(list.pairs)
}
\author{
Javier F. Tabima
}
\seealso{
\code{\link{dn.ds}}
}

