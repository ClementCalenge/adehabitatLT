\name{set.limits}
\alias{set.limits}
\title{Define the Same Time Limits for several Bursts in a
  Regular Trajectory}
\description{
  This function sets the same time limits for several
  bursts in a regular trajectory.
}
\usage{
set.limits(ltraj, begin, dur, pattern,
           units = c("sec", "min", "hour", "day"),
           tz = "", ...)
}
\arguments{
  \item{ltraj}{an object of class \code{ltraj}}
  \item{begin}{a character string which is used to determine the
    time of beginning of the study period (see below) }
  \item{dur}{the duration of the study period}
  \item{pattern}{a character string indicating the conversion
    specifications for \code{begin} (see below)}
  \item{units}{a character string indicating the time units of
    \code{dur}}
  \item{tz}{A timezone specification to be used for the conversion of
    \code{begin}.  System-specific, but \code{""} is the current time
    zone, and \code{"GMT"} is UTC (see \code{help(strptime)})}
  \item{\dots}{additional arguments to be passed to other functions}
}
\details{
  Some studies are intended to compare regular trajectories of the same
  duration collected at different period.  For example, the aim may be to
  identify the differences/similarities between different days (each one
  corresponding to a burst of relocation) in the pattern of movements of
  an animal between 05H00 and 08H00, with a time lag of 5 minutes.  In
  such cases, it is often convenient that the relocations of the bursts
  are paired (e.g. the fifth relocation correspond to the position of
  the animal at 5H30 for all bursts).
  
  The function \code{set.limits} is intended to ensure that the time of
  beginning, the end, and the duration of the trajectory is the same for
  all bursts of the object \code{ltraj}.  If relocations are collected
  outside the limits, they are removed (and so is the corresponding
  metadata in the attribute \code{infolocs}).  If the actual time limits
  of the burst cover a shorter period than those specified, missing values
  are added to the trajectory (and in the corresponding
  metadata in the attribute \code{infolocs}).

  Note that "time of beginning" is not a synonym for "date".  That is,
  two trajectories of the same animal, both beginning at 05H00 and
  ending at 08H00, have the same time of beginning, but are necessarily
  not sampled on the same day, which implies that they correspond to
  different dates.  For this reason, the time of beginning is indicated
  to the function \code{set.limits} by a character string, and the
  parameter \code{pattern} should indicate the conversion specifications.
  These conversions specifications are widely documented on the help
  page of the function \code{strptime}.  For example, to indicate that the
  trajectory begins at 5H00, the value for \code{begin} should be
  \code{"05:00"} and the value for \code{pattern} should be
  \code{"\%H:\%M"}. If the trajectory should begin on january 10th, the
  value for \code{begin} should be \code{"01:10"} and \code{pattern}
  should be \code{"\%m:\%d"}.  Note that the only conversion
  specifications allowed in this function are \code{\%S} (seconds),
  \code{\%M} (minutes), \code{\%H} (hours), \code{\%d} (day), \code{\%m}
  (month), \code{\%Y} (year with century), \code{\%w} (weekday), and
  \code{\%j} (yearday). See \code{help(strptime)} for additional
  information on these convention specifications.

}
\value{
  an object of class \code{ltraj}
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}
\seealso{\code{\link{ltraj}} for additional information on objects of
  class \code{ltraj}, \code{\link{sett0}} for additional
  information on regular trajectories, and \code{\link{sd2df}} for
  additionnal information about regular trajectories of the same
  duration.  See also \code{\link{strptime}} for further information
  about conversion specifications for dates.}
\examples{

## load data on the ibex
data(ibex)
ibex

## The monitoring of the 4 ibex should start and end at the same time
## define the time limits

ib2 <- set.limits(ibex, begin="2003-06-01 00:00", dur=14,
                  units="day", pattern="\%Y-\%m-\%d \%H:\%M", tz="Europe/Paris")
ib2
is.sd(ib2)

## All the trajectories cover the same study period
## Relocations are collected at the same time. This dataset can now be
## used for studies of interactions between animals

}
\keyword{programming}
