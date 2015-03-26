\name{c.ltraj}
\alias{c.ltraj}
\title{Combine Bursts of Relocations in Objects of Class "ltraj"}
\description{
  This function combines several objects of class \code{ltraj}.
}
\usage{
\method{c}{ltraj}(...)
}
\arguments{
  \item{\dots}{objects of class \code{ltraj} to be combined}
}
\value{
  An object of class \code{ltraj}.
}
\author{
  Clement Calenge \email{clement.calenge@oncfs.gouv.fr}
}
\seealso{\code{\link{ltraj}} for further information on the class
  \code{ltraj}, \code{\link{Extract.ltraj}} to extract or replace 
  bursts of relocations, \code{\link{plot.ltraj}} and
  \code{\link{trajdyn}} for graphical
  displays, \code{\link{gdltraj}} to specify a time period}
\examples{

data(puechcirc)

(i <- puechcirc[1])
(j <- puechcirc[3])

(toto <- c(i,j))

}
\keyword{spatial}
