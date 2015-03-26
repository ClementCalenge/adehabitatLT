\name{hbrown}
\alias{hbrown}
\title{ Estimates the value of h for a Brownian motion }
\description{
  \code{hbrown} estimates the scaling factor \code{h} (used in the
  Brownian motion, see \code{help(simm.brown)}) from a trajectory.
}
\usage{
hbrown(x)
}
\arguments{
  \item{x}{an object of class \code{ltraj}}
}
\value{
  a vector with one estimate per burst of the object of class
  \code{ltraj}.
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}
\seealso{\code{\link{simm.brown}}}
\examples{
toto <- simm.brown(1:200, h=4)
hbrown(toto)

toto <- simm.brown(1:200, h=20)
hbrown(toto)


}
\keyword{spatial}
