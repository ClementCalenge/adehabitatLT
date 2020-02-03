\name{is.regular}
\alias{is.regular}
\title{Regular Trajectories}
\description{
  \code{is.regular} tests whether a trajectory is regular (i.e. constant
  time lag between successive relocations).
}
\usage{
is.regular(ltraj)
}
\arguments{
  \item{ltraj}{an object of class \code{ltraj}}
}
\value{
  \code{is.regular} returns a logical value
}
\author{ Clement Calenge \email{clement.calenge@ofb.gouv.fr} }
\seealso{\code{\link{ltraj}}}
\examples{
data(capreotf)
is.regular(capreotf)
plotltr(capreotf, "dt")

data(albatross)
is.regular(albatross)
plotltr(albatross, "dt")
}
\keyword{programming}
