% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.r
\name{initialize,genclone-method}
\alias{initialize,genclone-method}
\alias{[,genclone,ANY,ANY,ANY-method}
\alias{show,genclone-method}
\alias{print,genclone-method}
\title{Methods used for the genclone object}
\usage{
\S4method{initialize}{genclone}(.Object, ..., mlg, mlgclass = TRUE)

\S4method{[}{genclone,ANY,ANY,ANY}(x, i, j, ..., mlg.reset = FALSE, drop = FALSE)

\S4method{show}{genclone}(object)

\S4method{print}{genclone}(x, ..., fullnames = TRUE)
}
\arguments{
\item{.Object}{a character, "genclone"}

\item{...}{passed on to the \code{\link[adegenet:genind-class]{genind}} object.}

\item{mlg}{a vector where each element assigns the multilocus genotype of 
that individual in the data set.}

\item{mlgclass}{a logical value specifying whether or not to translate the 
mlg object into an MLG class object.}

\item{x}{a genclone object}

\item{i}{vector of numerics indicating number of individuals desired}

\item{j}{a vector of numerics corresponding to the loci desired.}

\item{mlg.reset}{logical. Defaults to \code{FALSE}. If \code{TRUE}, the mlg
vector will be reset}

\item{drop}{set to \code{FALSE}}

\item{object}{a genclone object}

\item{fullnames}{\code{logical}. If \code{TRUE}, then the full names of the
populations will be printed. If \code{FALSE}, then only the first and last
three population names are displayed.}
}
\description{
Default methods for subsetting genclone objects.
}
\author{
Zhian N. Kamvar
}
\keyword{internal}
