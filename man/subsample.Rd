\name{subsample}
\alias{subsample}
\title{Subsample a Trajectory}
\description{
  This function subsamples a regular trajectory (i.e. changes the time
  lag between successive relocations).
}
\usage{
subsample(ltraj, dt, nlo = 1,
          units = c("sec", "min", "hour", "day"), ...)
}
\arguments{
  \item{ltraj}{an object of class \code{ltraj}}
  \item{dt}{numeric value.  The new time lag (should be a multiple of
    the time lag in \code{ltraj})}
  \item{nlo}{an integer, or a vector of integers (with length equal to
    the number of bursts in \code{ltraj}), indicating the position of
    the first location of the new bursts in the old bursts.  For
    example, if the previous time lag is equal to 300 seconds and the
    new time lag is 900 seconds, the new bursts may begin at the first,
    second or third relocations of the old bursts in \code{ltraj}.}
  \item{units}{character string. The time units of \code{dt}}
  \item{\dots}{additional arguments to be passed to other functions}
}
\value{
  An object of class \code{ltraj}
}
\author{ Clement Calenge \email{clement.calenge@ofb.gouv.fr} }
\seealso{\code{\link{ltraj}} for additional information on objects of
  class \code{ltraj}, \code{\link{is.regular}} for regular trajectories. }
\examples{
data(capreotf)
plot(capreotf)

toto <- subsample(capreotf, dt = 900)
plot(toto)
}
\keyword{programming}
