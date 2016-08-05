\name{runsNAltraj}
\alias{runsNAltraj}
\alias{summaryNAltraj}
\alias{plotNAltraj}
\title{Highlighting the Patterns in Missing Values in Trajects}
\description{
  \code{runsNAltraj} performs a runs test to detect any autocorrelation
  in the location of missing relocations, for each burst of an object of
  class \code{ltraj}.
  
  \code{summaryNAltraj} returns a summary of the number and proportion
  of missing values for each burst of an object of class
  \code{ltraj}.

  \code{plotNAltraj} plots the missing values in an object of class
  \code{ltraj} against the time.
}
\usage{
runsNAltraj(x, nrep = 500, plotit = TRUE, \dots) 

summaryNAltraj(x)

plotNAltraj(x, \dots)
}
\arguments{
  \item{x}{ An object of class \code{ltraj}}
  \item{nrep}{Number of randomisations}
  \item{plotit}{logical.  Whether the results should be plotted on a
    graph}
  \item{\dots}{Further arguments to be passed to the generic function
    \code{plot}}
}
\details{
  The statistics used here for the test is the number
  of runs in the sequence of relocations. For example, the sequence
  reloc-NA-NA-reloc-reloc-reloc-NA-NA-NA-reloc contains 5 runs, 3 runs
  of successful relocations and 2 runs of missing values. Under the
  hypothesis of random distribution of the missing values in the
  sequence, the theoretical expectation and standard deviation of the
  number of runs is known. The runs test is a randomization test that
  compares the standardized value of the number of runs
  (i.e. (value-expectation)/(standard deviation)) to the distribution of
  values obtained after randomizing the distribution of the NA in the
  sequence. Thus, a negative value of this standardized number of runs
  indicates that the missing values tend to be clustered together in the
  sequence.
}
\value{
  \code{runsNAltraj} returns a list of objects of class
  \code{randtest} (if a burst does not contain any missing value, the
  corresponding component is NULL).
}
\note{
  In the versions of adehabitatLT prior to 0.3.21, a bug occurred in the
  calculation of the P-value (the test actually presented the value
  1-P). This bug is now corrected.
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}
\seealso{\code{\link{ltraj}} for additional information about objects of
  class \code{ltraj}, \code{\link{setNA}} for additional information
  about missing values in such objects}
\examples{

## Two relocations are theoretically separated by
## 10 minutes (600 seconds)
data(puechcirc)
puechcirc

## plot the missing values
plotNAltraj(puechcirc)


## Test for an autocorrelation pattern in the missing values
(runsNAltraj(puechcirc))

}
\keyword{hplot}
\keyword{NA}
