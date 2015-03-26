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
\value{
  \code{runsNAltraj} returns a list of objects of class
  \code{randtest} (if a burst does not contain any missing value, the
  corresponding component is NULL).
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
