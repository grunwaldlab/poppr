\name{multi.nd}
\alias{multi.nd}
\title{multi.nd}
\usage{
multi.nd(x)
}
\arguments{
  \item{x}{a \code{list} of multiple alignment \code{fasta}
  files (preferably full path, fasta|fas extension).}
}
\value{
a table with the nucleotide diversity calculation for all
the files in the list.
}
\description{
Nucleotide diversity wrapper for multiple sets of FASTA
files and the \code{\link{n.diversity}} function
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
# multi.nd(list.pairs)
}
\author{
Javier F. Tabima
}
\seealso{
\code{\link{td}}
}

