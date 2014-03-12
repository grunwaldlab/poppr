\name{read.genalex}
\alias{read.genalex}
\title{Importing data from genalex formatted *.csv files.}
\usage{
read.genalex(genalex, ploidy = 2, geo = FALSE, region = FALSE,
  genclone = TRUE, sep = ",")
}
\arguments{
  \item{genalex}{a *.csv file exported from genalex}

  \item{ploidy}{indicate the ploidy of the dataset}

  \item{geo}{indicates the presence of geographic data in
  the file. This data will be included in a data frame
  labeled \code{xy} in the\code{\link{other}} slot.}

  \item{region}{indicates the presence of regional data in
  the file.}

  \item{genclone}{should the output be a genclone object?
  Defaults to \code{TRUE}}

  \item{sep}{A character specifying the column separator of
  the data. Defaults to ",".}
}
\value{
A \code{\link{genind}} object.
}
\description{
read.genalex will read in a genalex-formatted file that has
been exported in a comma separated format and will parse
most types of genalex data. The output is a
\code{\link{genind}} object.
}
\note{
This function cannot handle raw allele frequency data.

The resulting genind object will have a data frame in the
\code{other} slot called \code{population_hierarchy}. This
will contain a column for your population data and a column
for your Regional data if you have set the flag.

In the case that there are duplicated names within the
file, this function will assume separate individuals and
rename each one to a sequence of integers from 1 to the
number of individuals. A vector of the original names will
be saved in the \code{other} slot under
\code{original_names}.

FOR POLYPLOID (> 2n) DATA SETS: Adegenet's genind object
has an all-or-none approach to missing data. If a sample
has missing data at a particular locus, then the entire
locus is considered missing. This works for diploids and
haploids where allelic dosage is unambiguous. For
polyploids this poses a problem as much of the data set
would be transformed into missing data.  With this
function, I have created a workaround.

When importing polyploid data sets, missing data is scored
as "0" and kept within the genind object as an extra
allele. This will break most analyses relying on allele
frequencies. All of the functions in poppr will work
properly with these data sets as multilocus genotype
analysis is agnostic of ploidy and we have written both
Bruvo's distance and the index of association in such a way
as to be able to handle polyploids presented in this
manner.
}
\examples{
\dontrun{
Aeut <- read.genalex(system.file("files/rootrot.csv", package="poppr"))

genalex2 <- read.genalex("genalex2.csv", geo=TRUE)
# A genalex file with geographic coordinate data.

genalex3 <- read.genalex("genalex3.csv", region=TRUE)
# A genalex file with regional information.

genalex4 <- read.genalex("genalex4.csv", region=TRUE, geo=TRUE)
# A genalex file with both regional and geographic information.
}
}
\author{
Zhian N. Kamvar
}
\seealso{
\code{\link{clonecorrect}}, \code{\link{genind}}
}

