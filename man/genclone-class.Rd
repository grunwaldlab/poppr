\docType{class}
\name{genclone-class}
\alias{genclone-class}
\title{Genclone class}
\description{
Genclone is an S4 class that extends the
\code{\linkS4class{genind-class}} from the
\emph{\link{adegenet}} package. It will have all of the
same attributes as the \code{\linkS4class{genind-class}},
but it will contain two extra slots that will help retain
information about population hierarchies and multilocus
genotypes.
}
\details{
The genclone class will allow for more optimized methods of
clone correcting and analyzing data over multiple levels of
population hierarchy.

Previously, for hierarchical analysis to work in a
\code{\link{genind}} object, the user had to place a data
frame in the \code{\link{other}} slot of the object. The
suggested name of the data frame was
\code{population_hierarchy}, and this was used to be able
to store the hierarchical information inside the object so
that the user did not have to keep track of that
information. This method worked, but it became apparent
that this method of doing things became a bit confusing to
the user as the method for changing the population of an
object became:

\code{pop(object) <-
other(object)$population_hierarchy$population_name}

That is a lot to keep track of. The new hierarchy slot will
allow the user to be able to insert a population hierarchy
into the slot and then change it with one function and a
formula:

\code{sethierarchy(object, ~Population/Subpopulation)}

making this become slightly more intuitive and tractable.

Calculations of multilocus genotypes is rapid, but
unfortunately, the assignments can change if the number of
individuals in the population is different. This means that
if you analyze the multilocus genotypes for two
subpopulations, they will have different multilocus
genotype assignments even though they may share multilocus
genotypes. The new slot allows us to be able to assign the
multilocus genotypes and retain that information no matter
how we subset the data set.
}
\section{Slots}{

\describe{
\item{\code{mlg}}{a vector representing multilocus genotypes for the data set.}

\item{\code{hierarchy}}{a data frame containing hierarchical levels.}
}}
\section{Extends}{
  Class \code{"\linkS4class{genind-class}"}, directly.
}
\author{
Zhian N. Kamvar
}

