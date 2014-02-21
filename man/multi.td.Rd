\name{multi.td}
\alias{multi.td}
\title{multi.td}
\usage{
multi.td(x, cutoff = 0.05, quiet = FALSE)
}
\arguments{
  \item{x}{a \code{list} of multiple alignment \code{fasta}
  files (preferably full path, fasta|fas extension).}

  \item{cutoff}{a \code{integer} of the cutoff value for
  the p-statistic.}

  \item{quiet}{a \code{boolean} for verbosity.}
}
\value{
a table with the Tajima's D calculation for all the files
in the list.
}
\description{
Tajima-D wrapper for multiple sets of FASTA files and the
\code{\link{td}} function
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
# multi.td(list.pairs)
}
\author{
Javier F. Tabima
}
\seealso{
\code{\link{td}}
}

