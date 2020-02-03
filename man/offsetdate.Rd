\name{offsetdate}
\alias{offsetdate}
\title{Date Handling in an Object of Class 'ltraj'}
\description{
  This functions allows to set an offset value from the date in an
  object of class \code{ltraj} of type II (time recorded).
}
\usage{
offsetdate(ltraj, offset, units = c("sec", "min", "hour", "day"))
}
\arguments{
  \item{ltraj}{an object of class \code{ltraj} of type II (time
    recorded)}
  \item{offset}{a numeric value indicating the offset to be deducted
    from the date}
  \item{units}{a character string indicating the time units for
    \code{offset}}
}
\details{
  The use of offset is a convenient way to define reference dates in an
  object of class \code{ltraj}.  For example, if the animal is monitored
  every night, from 18H00 to 06H00, the fact that the beginning and the
  end of the monitoring do not correspond to the same day may cause
  difficulties to handle the trajectory.  Though these difficulties are
  not unsurmountable, it is often convenient to deduct an offset to the
  trajectory, so that the first relocation is collected at 0H and the
  last one at 12H00 the same day (i.e., in this example, an offset of
  18 hours).

}
\value{
  an object of class \code{ltraj}
}
\author{  Clement Calenge \email{clement.calenge@ofb.gouv.fr} }
\seealso{ \code{\link{ltraj}} for additional information on objects of
  class \code{ltraj}}
\examples{
data(puechcirc)

plotltr(puechcirc, "dt")

toto <- offsetdate(puechcirc, 17, "hour")

plotltr(puechcirc, "dt")

}
\keyword{programming}
