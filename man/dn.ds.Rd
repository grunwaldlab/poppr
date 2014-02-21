\name{dn.ds}
\alias{dn.ds}
\title{dn.ds}
\usage{
dn.ds(x)
}
\arguments{
  \item{x}{an alignment \code{fasta} file (preferably full
  path, fasta|fas extension).}
}
\value{
a table with the dN/dS (Ka/Ks) value and the estimate for
selection (Positive, Negative or Neutral)
}
\description{
dN/dS (as Ka/Ks estimation) estimation of selection.
Wrapper of \code{seqinr} function \code{\link{kaks}}
}
\note{
The input FASTA file MUST BE AN ALIGNMENT! Be sure to use
the full path. The FASTA files are read by
\code{\link{seqinr}} function \code{read.alignment}.
}
\examples{
# Analyze one of the included FASTA alignment files
# dn.ds("data_1.fasta",format="fasta")
}
\author{
Javier F. Tabima
}
\references{
Li WH, Wu CI, Luo CC. 1985. A new method for estimating
synonymous and nonsynonymous rates of nucleotide
substitution considering the relative likelihood of
nucleotide and codon changes. Mol Biol Evol 2(2):150-174.
Zhang, Z., Li, J., Zhao, X., Wang, J., Wong, G.K. and Yu,
J. (2006) KaKs_Calculator: calculating Ka and Ks through
model selection and model averaging, Genomics Proteomics
Bioinformatics, 4(4): 259-263.
}
\seealso{
\code{\link{kaks}}
}

