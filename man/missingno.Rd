\name{missingno}
\alias{missingno}
\title{How to deal with missing data in a genind object.}
\usage{
missingno(pop, type = "loci", cutoff = 0.05, quiet = FALSE)
}
\arguments{
  \item{pop}{a \code{\linkS4class{genclone}} or
  \code{\linkS4class{genind}} object.}

  \item{type}{a \code{character} string: can be "ignore",
  "zero", "mean", "loci", or "geno" (see \code{Details} for
  definitions).]}

  \item{cutoff}{\code{numeric}. A number from 0 to 1
  indicating the allowable rate of missing data in either
  genotypes or loci. This will be ignored for \code{type}
  values of \code{"mean"} or \code{"zero"}.}

  \item{quiet}{if \code{TRUE}, it will print to the screen
  the action performed.}
}
\value{
a \code{\linkS4class{genclone}} or
\code{\linkS4class{genind}} object.
}
\description{
missingno gives the user four options to deal with missing
data.
}
\details{
These methods provide a way to deal with systematic missing
data and to give a wrapper for \code{adegenet}'s \code{
\link{na.replace}} function.  ALL OF THESE ARE TO BE USED
WITH CAUTION.

\subsection{Treatment types}{ \itemize{
\item{\code{"ignore"} - does not remove or replace missing
data.} \item{\code{"loci"} - removes all loci containing
missing data in the entire data set. }
\item{\code{"genotype"} - removes any
genotypes/isolates/individuals with missing data.}
\item{\code{"mean"} - replaces all NA's with the mean of
the alleles for the entire data set.} \item{\code{"zero"}
or \code{"0"} - replaces all NA's with "0". Introduces more
diversity.}}}
}
\note{
\emph{"wild missingno appeared!"}
}
\examples{
data(nancycats)

nancy.locina <- missingno(nancycats, type = "loci")

## Found 617 missing values.
## 2 loci contained missing values greater than 5\%.
## Removing 2 loci : fca8 fca45

nancy.genona <- missingno(nancycats, type = "geno")

## Found 617 missing values.
## 38 genotypes contained missing values greater than 5\%.
## Removing 38 genotypes : N215 N216 N188 N189 N190 N191 N192 N302 N304 N310
## N195 N197 N198 N199 N200 N201 N206 N182 N184 N186 N298 N299 N300 N301 N303
## N282 N283 N288 N291 N292 N293 N294 N295 N296 N297 N281 N289 N290

# Replacing all NA with "0" (see na.replace in the adegenet package).
nancy.0 <- missingno(nancycats, type = "0")

## Replaced 617 missing values

# Replacing all NA with the mean of each column (see na.replace in the
# adegenet package).
nancy.mean <- missingno(nancycats, type = "mean")

## Replaced 617 missing values
}
\author{
Zhian N. Kamvar
}
\seealso{
\code{\link{na.replace}}, \code{\link{poppr}},
\code{\link{poppr.amova}}, \code{\link{nei.dist}},
\code{\link{aboot}}
}

