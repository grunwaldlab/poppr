\name{missing_table}
\alias{missing_table}
\title{Create a table summarizing missing data a genind or genclone object}
\usage{
missing_table(x, percent = TRUE, plot = FALSE, df = FALSE,
  returnplot = FALSE, low = "blue", high = "red", plotlab = TRUE,
  scaled = TRUE)
}
\arguments{
  \item{x}{a \linkS4class{genind} or \linkS4class{genclone}
  object.}

  \item{percent}{\code{logical}. If \code{TRUE} (default),
  table and plot will represent missing data as a
  percentage of each cell. If \code{FALSE}, the table and
  plot will represent missing data as raw counts. (See
  details)}

  \item{plot}{\code{logical}. If \code{TRUE}, a simple
  heatmap will be produced. If \code{FALSE} (default), no
  heatmap will be produced.}

  \item{df}{\code{logical}. If \code{TRUE}, the data will
  be returned as a long form data frame. If \code{FALSE}
  (default), a matrix with populations in rows and loci in
  columns will be returned.}

  \item{returnplot}{\code{logical}. If \code{TRUE}, a list
  is returned with two elements: \code{table} - the normal
  output and \code{plot} - the ggplot object. If
  \code{FALSE}, the missing table is returned.}

  \item{low}{\code{character}. What color should represent
  no missing data?  (default: "blue")}

  \item{high}{\code{character}. What color should represent
  the highest amount of missing data? (default: "red")}

  \item{plotlab}{\code{logical}. If \code{TRUE} (default),
  values of missing data greater than 0\% will be plotted.
  If \code{FALSE}, the plot will appear unappended.}

  \item{scaled}{\code{logical}. This is for when
  \code{percent = TRUE}. If \code{TRUE} (default), the
  color specified in \code{high} will represent the highest
  observed value of missing data. If \code{FALSE}, the
  color specified in \code{high} will represent 100\%.}
}
\value{
a matrix, data frame (\code{df = TRUE}), or a list
(\code{returnplot = TRUE}) representing missing data in a
\linkS4class{genind} or \linkS4class{genclone} object.
}
\description{
Create a table summarizing missing data a genind or
genclone object
}
\details{
This data is potentially useful for identifying areas of
systematic missing data. There are a few caveats to be
aware of.  \itemize{ \item \strong{For counts of missing
data}: Each count represents the number of individuals with
missing data at each locus. The last column, "mean" can be
thought of as the average number of individuals with
missing data per locus. \item \strong{For percentage
missing data}: This percentage is \strong{relative to the
population and locus}, not ot the enitre data set. The last
colum, "mean" represents the average percent of the
population with missing data per locus. }
}
\examples{
data(nancycats)
nancy.miss <- missing_table(nancycats, plot = TRUE)
}
\author{
Zhian N. Kamvar
}

