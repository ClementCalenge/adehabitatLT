\name{setNA}
\alias{setNA}
\title{ Place Missing Values in Objects of Class 'ltraj' }
\description{
  This function places missing values in an (approximately) regular
  trajectory, when a relocation should have been collected, but is
  actually missing.
}
\usage{
setNA(ltraj, date.ref, dt, tol = dt/10,
      units = c("sec", "min", "hour", "day"), \dots)
}
\arguments{
  \item{ltraj}{an object of class \code{ltraj}}
  \item{date.ref}{an object of class \code{POSIXt} (see below)}
  \item{dt}{the time lag between relocations}
  \item{tol}{the tolerance, which measures the imprecision in the timing
    of data collection (see below)}
  \item{units}{a character string indicating the time units for
    \code{dt} and \code{tol}}
  \item{\dots}{additional arguments to be passed to the function
    \code{rec}}
}
\details{
  During the field study, the collection of the relocations of a
  trajectory may sometimes fail, which results into missing values.  The
  class \code{ltraj} deal with these missing values, so that it is
  recommended to store the missing values in the data *before* the
  creation of the object of class \code{ltraj}.  For example, GPS
  collars often fail to locate the animal, so that the GPS data imported
  within R contain missing values.  It is recommended to *not remove*
  these missing values.
  
  However, sometimes, the data come without any information concerning
  the placement of these missing values.  If the trajectory is
  approximately regular (i.e. approximately constant time lag), it is
  possible to determine where these missing values should occur in the
  object of class \code{ltraj} (and in the optional attribute
  \code{infolocs}).  This is the role of the function \code{setNA}.

  The relocations in the object of class \code{ltraj} may not have been
  collected at exactly identical time lag (e.g. a relocation is
  collected at 17H57 instead of 18H00). The function \code{setNA}
  requires that the imprecision in the timing is at most equal to
  \code{tol}.  Because of this imprecision, it is necessary to pass a
  reference date as argument to the function \code{setNA}. This
  reference date is used to determine at which time the missing values
  should be placed.
  
  The reference date is chosen so that the rest of the division of
  (date.relocations - reference.date) by the time lag \code{dt} is equal
  to zero.  For example, if it is known that one of the relocations of the
  trajectory has been collected on January 16th 1996 at 18H00,
  and if the theoretical time lag between two relocations is of one
  hour, the date of reference could be (for example) the August 1st 2017
  at 05H00, because these two dates are separated by an exact number of
  hours (i.e. an exact number of \code{dt}). Therefore, any date
  fulfilling this condition could be passed as reference date.
  Alternatively, the August 1st 2007 at 05H30 is an uncorrect reference
  date, because the number of hours separating these two dates is not an
  integer.
}
\value{
  An object of class \code{ltraj}
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }
\seealso{\code{\link{ltraj}} for additional information about objects
of class \code{ltraj}. \code{\link{sett0}} (especially the examples
of this help page) and \code{\link{is.regular}} for additional
information about regular trajectories.}
\examples{

data(porpoise)
foc <- porpoise[1]

## the list foc does not contain any missing value:
foc
plotNAltraj(foc)

## we remove the second to tenth relocation
foc[[1]] <- foc[[1]][-c(2:10),]
foc <- rec(foc)

## The missing values are not visible:
foc
plotNAltraj(foc)

## The porpoise is located once a day.
## We use the first relocation as the reference date
foc2 <- setNA(foc, foc[[1]]$date[1], 24*3600)

## Missing values are now present
foc2
plotNAltraj(foc2)

}
\keyword{programming}
