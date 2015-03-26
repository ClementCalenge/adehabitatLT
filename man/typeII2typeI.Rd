\name{typeII2typeI}
\alias{typeII2typeI}
\title{Change the Type of a Trajectory}
\description{
  This function transforms a trajectory of type II (time recorded) into a
  trajectory of type I (time not recorded).
}
\usage{
typeII2typeI(x)
}
\arguments{
  \item{x}{a object of class "ltraj" of type II}
}
\value{
  An object of class "ltraj"
}
\author{  Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }
\seealso{ \code{\link{as.ltraj}} for additional information on the
  objects of class "ltraj" }
\examples{
data(puechcirc)
puechcirc
typeII2typeI(puechcirc)
}
\keyword{programming}
