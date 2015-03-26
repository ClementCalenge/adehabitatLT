\name{na.omit.ltraj}
\alias{na.omit.ltraj}
\title{Removes the missing values in a trajectory}
\description{
  \code{na.omit.ltraj} can be used to remove missing relocations from a
  trajectory.
}
\usage{
\method{na.omit}{ltraj}(object, \dots)
}
\arguments{
  \item{object}{an object of class \code{ltraj}}
  \item{\dots}{additionnal arguments to be passed to or from other methods}
}
\value{
  An object of class \code{ltraj}
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}
\seealso{
  \code{\link{setNA}} to place the missing values in the trajectory
}
\examples{
data(puechcirc)
puechcirc

na.omit(puechcirc)

}
\keyword{NA}
