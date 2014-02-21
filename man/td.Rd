\name{td}
\alias{td}
\title{td}
\usage{
td(x)
}
\arguments{
  \item{x}{an alignment \code{fasta} file (preferably full
  path, fasta|fas extension).}
}
\value{
a table with the Tajima's D calculation results: Tajimas D
value, the p-values for normality and poisson, the files
and if the p-value (for normality) is significant or not.
}
\description{
Tajima-D wrapper of \code{pegas} function
\code{\link{tajima.test}} to get a data frame out of the
data with the Tajima D value, the p-values for normality
and poisson, the files and if the p-value (for normality)
is significant or not.
}
\note{
The input FASTA file MUST BE AN ALIGNMENT! Be sure to use
the full path. The FASTA files are read by
\code{\link{ape}} function \code{read.dna} and then
processed and ran trough \code{pegas} function
\code{\link{tajima.test}}.
}
\examples{
# Analyze one of the included FASTA alignment files
# td("data_1.fasta",format="fasta")
}
\author{
Javier F. Tabima
}
\references{
Tajima, F. (Nov 1989). "Statistical method for testing the
neutral mutation hypothesis by DNA polymorphism.". Genetics
123 (3): 585-95. PMC 1203831. PMID 2513255.
}
\seealso{
\code{\link{tajima.test}}
}

