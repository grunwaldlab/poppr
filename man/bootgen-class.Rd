\docType{class}
\name{bootgen-class}
\alias{bootgen-class}
\title{Bootgen object}
\description{
An internal object used for bootstrapping. Not intended for
user interaction.
}
\section{Slots}{

\describe{
\item{\code{alllist}}{a list with numeric vectors, each representing a different
  locus where each element in the vector represents the index for a specific
  allele.}

\item{\code{names}}{a vector containing names of the observed samples}
}}
\section{Extends}{
  Virtual Class \code{"\linkS4class{gen}"}.
}
\author{
Zhian N. Kamvar
}

