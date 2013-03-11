\name{poppr}
\alias{poppr}
\title{Produce a basic summary table for population genetic analyses.}
\usage{
  poppr(pop, total = TRUE, sublist = c("ALL"),
    blacklist = c(NULL), sample = 0, method = 1,
    missing = "ignore", cutoff = 0.05, quiet = "minimal",
    clonecorrect = FALSE, hier = c(1),
    dfname = "population_hierarchy", hist = TRUE,
    minsamp = 10)
}
\arguments{
  \item{pop}{a \code{\link{genind}} object OR any fstat,
  structure, gentix, genpop, or genalex formatted file.}

  \item{total}{default \code{TRUE}. Should indecies be
  calculated for the combined populations represented in
  the entire file?}

  \item{sublist}{a list of character strings or integers to
  indicate specific population names (located in
  \code{$pop.names} within the \code{\link{genind}} object)
  Defaults to "ALL".}

  \item{blacklist}{a list of character strings or integers
  to indicate specific populations to be removed from
  analysis. Defaults to NULL.}

  \item{sample}{an integer indicating the number of
  permutations desired to obtain p-values. Sampling will
  shuffle genotypes at each locus to simulate a panmictic
  population using the observed genotypes. Calculating the
  p-value includes the observed statistics, so set your
  sample number to one off for a round p-value (eg.
  \code{sample = 999} will give you p = 0.001 and
  \code{sample = 1000} will give you p = 0.000999001).}

  \item{method}{an integer from 1 to 4 indicating the
  method of sampling desired. see \code{popsample} for
  details.}

  \item{missing}{how should missing data be treated?
  \code{"zero"} and \code{"mean"} will set the missing
  values to those documented in \code{\link{na.replace}}.
  \code{"loci"} and \code{"geno"} will remove any loci or
  genotypes with missing data, respectively (see
  \code{\link{missingno}} for more information.}

  \item{quiet}{Should the function print anything to the
  screen while it is performing calculations? \code{TRUE}
  prints nothing, \code{"minimal"} (defualt) will print the
  population name and dots indicating permutation progress,
  \code{FALSE} will print the number of individuals in the
  population. \code{"noisy"} will print out the individual
  indecies as they are produced.}

  \item{clonecorrect}{default \code{FALSE}. must be used
  with the \code{hier} and \code{dfname} parameters, or the
  user will potentially get undesiered results. see
  \code{\link{clonecorrect}} for details.}

  \item{hier}{a \code{numeric or character list}. This is
  the list of vectors within a data frame (specified in
  \code{dfname}) in the 'other' slot of the
  \code{\link{genind}} object. The list should indicate the
  population hierarchy to be used for clone correction.}

  \item{dfname}{a \code{character string}. This is the name
  of the data frame or list containing the vectors of the
  population hierarchy within the \code{other} slot of the
  \code{\link{genind}} object.}

  \item{hist}{\code{logical} if \code{TRUE} a histogram
  will be produced for each population.}

  \item{minsamp}{an \code{integer} indicating the minimum
  number of individuals to resample for rarefaction
  analysis.}
}
\value{
  \item{Pop}{A vector indicating the pouplation factor}
  \item{N}{An integer vector indicating the number of
  individuals/isolates in the specified population.}
  \item{MLG}{An integer vector indicating the number of
  multilocus genotypes found in the specified poupulation,
  (see: \code{\link{mlg}})} \item{eMLG}{The expected number
  of MLG at the lowest common sample size (set by the
  parameter \code{minsamp}.} \item{SE}{The standard error
  for the rarefaction analysis} \item{H}{Shannon-Weiner
  Diversity index} \item{G}{Stoddard and Taylor's Index}
  \item{Hexp}{Expected heterozygosity or Nei's 1987
  genotypic diversity corrected for sample size.}
  \item{E.5}{Evenness} \item{Ia}{A numeric vector giving
  the value of the Index of Association for each population
  factor, (see \code{\link{ia}}).} \item{p.Ia}{A numeric
  vector indicating the p-value for Ia from the number of
  reshufflings indicated in \code{sample}. Lowest value is
  1/n where n is the number of observed values.}
  \item{rbarD}{A numeric vector giving the value of the
  Standardized Index of Association for each population
  factor, (see \code{\link{ia}}).} \item{p.rD}{A numeric
  vector indicating the p-value for rbarD from the number
  of reshufflings indicated in \code{sample}. Lowest value
  is 1/n where n is the number of observed values.}
  \item{File}{A vector indicating the name of the original
  data file.}
}
\description{
  This function allows the user to quickly view indecies of
  distance, heterozygosity, and inbreeding to aid in the
  decision of a path to further analyze a specified
  dataset. It natively takes \code{\link{genind}} formatted
  files, but can convert any raw data formats that adegenet
  can take (fstat, structure, gentix, and genpop) as well
  as genalex files exported into a csv format (see
  \code{\link{read.genalex}} for details).
}
\note{
  All values are rounded to three significant digits for
  the final table.
}
\examples{
data(nancycats)
poppr(nancycats)

\dontrun{
poppr(nancycats, sample=99, total=FALSE, quiet=FALSE)

# Note: this is a larger data set that could take a couple of minutes to run
# on slower computers.
data(H3N2)
poppr(H3N2, total=FALSE, sublist=c("Austria", "China", "USA"),
				clonecorrect=TRUE, hier="country", dfname="x")
}
}
\seealso{
  \code{\link{clonecorrect}}, \code{\link{poppr.all}},
  \code{\link{ia}}, \code{\link{missingno}},
  \code{\link{mlg}}
}

