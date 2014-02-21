\name{n.diversity}
\alias{n.diversity}
\title{n.diversity}
\usage{
n.diversity(x)
}
\arguments{
  \item{x}{an alignment \code{fasta} file (preferably full
  path, fasta|fas extension).}
}
\value{
a table with the Nucleotide Diversity calculation results.
}
\description{
Nucleotide diversity calculation for sequence data. Wrapper
of \code{pegas} function \code{nuc.div} to get a data frame
out of the data with the nucleotide diversity estimate.
}
\note{
The input FASTA file MUST BE AN ALIGNMENT! Be sure to use
the full path. The FASTA files are read by
\code{\link{ape}} function \code{read.dna} and then
processed and ran trough \code{pegas} function
\code{\link{nuc.div}}.
}
\examples{
# Analyze one of the included FASTA alignment files
# n.diversity("data_1.fasta",format="fasta")
}
\author{
Javier F. Tabima
}
\references{
Nei, M. (1987) Molecular evolutionary genetics. New York:
Columbia University Press. Nei, M and Li, W. 1979.
"Mathematical Model for Studying Genetic Variation in Terms
of Restriction Endonucleases". PNAS 76
}
\seealso{
\code{\link{nuc.div}}
}

