\name{sett0}
\alias{sett0}
\title{Round the Timing of Collection of Relocations to Obtain Regular
  Trajectory}
\description{
  This function rounds the timing of collection of relocations in an
  object of class \code{ltraj} to obtain a regular trajectory, based on
  a reference date.
}
\usage{
sett0(ltraj, date.ref, dt, correction.xy = c("none", "cs"),
      tol = dt/10, units = c("sec", "min", "hour", "day"), ...)
}
\arguments{
  \item{ltraj}{an object of class \code{ltraj}}
  \item{date.ref}{an object of class \code{POSIXt} containing either one
    reference date (the same for all animals) or n reference dates, where
    n is the number of bursts in \code{ltraj} (see below)}
  \item{dt}{the time lag between relocations}
  \item{correction.xy}{the correction for the coordinates.
    \code{"none"} (default), does not performs any correction.
    \code{"cs"} performs a correction based on the hypothesis that the
    animal moves at constant speed (see below).}
  \item{tol}{the tolerance, which measures the imprecision in the timing
    of data collection (see below)}
  \item{units}{the time units for \code{dt} and \code{tol}}
  \item{\dots}{additional arguments to be passed to the function
    \code{rec}}
}
\details{
  Trajectories are stored in adehabitatLT as lists of "bursts" of
  successive relocations with the timing of relocation.  Regular
  trajectories are characterized by a constant time lag \code{dt}
  between successive relocations (don't mix animals located every 10
  minutes and animals located every day in a regular trajectory).
  
  However, in many cases, the actual time lag in the data may not be
  equal to the theoretical time lag \code{dt}: there may be some
  negligible imprecision in the time of collection of the data (e.g. an
  error of a few seconds on a time lag of one hour).
  
  But many functions of adehabitatLT require exact regular
  trajectories. \code{sett0} allows to round the date so that all the
  successive relocations are separated exactly by \code{dt}.  The
  function \code{sett0} requires that the imprecision is at most equal
  to \code{tol}. To proceed, it is necessary to pass a reference date as
  argument.

  The reference date is chosen so that the rest of the division of
  (date.relocations - reference.date) by \code{dt} is equal to zero.
  For example, if it is known that one of the relocations of the
  trajectory should have been collected on January 16th 1996 at 18H00,
  and if the theoretical time lag between two relocations is of one
  hour, the date of reference could be (for example) the August 1st 2017
  at 05H00, because these two dates are separated by an exact number of
  hours.  Alternatively, the August 1st 2007 at 05H30 is an uncorrect
  reference date, because the number of hours separating these two dates
  is not an integer.
  
  Note that this rounding adds an error on the relocation.  For example,
  the position of a moving animal at 17H57 is not the same as its
  position at 18H00.  If the time imprecision in the data collection is
  negligible (e.g. a few seconds, while \code{dt} is equal to an hour),
  this "noise" in the relocations can be ignored, but if it is more
  important, a correction on the relocation is needed.  The function
  \code{sett0} may correct the relocations based on the hypothesis of
  constant speed (which is not necessarily biologically relevant, see
  examples).
  
  Note finally that missing values can be present in the trajectory.
  Indeed, there are modes of data collection that fail to locate the
  animal at some dates.  These failures should appear as missing values
  in the regular trajectory.  It is often convenient to use the function
  \code{setNA} before the function \code{sett0} to set the missing
  values in a (nearly) regular trajectory.
}
\value{
  an object of class \code{ltraj} containing a regular trajectory.
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }
\seealso{\code{\link{ltraj}} for additional information on objects of
  class \code{ltraj}, \code{\link{is.regular}} for regular trajectories,
  \code{\link{setNA}} to place missing values in the trajectory and
  \code{\link{cutltraj}} to cut a trajectory into several bursts based
  on a criteria. 
}
\examples{
\dontrun{
#########################################################################
##
##
## Transform a GPS monitoring on 4 ibex into a regular trajectory
##

data(ibexraw)
is.regular(ibexraw)

## the data are not regular: see the distribution of dt (in hours)
## according to the date

plotltr(ibexraw, "dt/3600")

## The relocations have been collected every 4 hours, and there are some
## missing data

## The reference date: the hour should be exact (i.e. minutes=0):
refda <- strptime("00:00", "\%H:\%M", tz="Europe/Paris")
refda

## Set the missing values
ib2 <- setNA(ibexraw, refda, 4, units = "hour") 

## now, look at dt for the bursts:
plotltr(ib2, "dt")

## dt is nearly regular: round the date:

ib3 <- sett0(ib2, refda, 4, units = "hour") 

plotltr(ib3, "dt")
is.regular(ib3)

## ib3 is now regular
}
}
\keyword{programming}
