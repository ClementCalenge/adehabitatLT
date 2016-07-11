\name{cutltraj}
\alias{cutltraj}
\alias{bindltraj}
\title{Split Trajectories into Several Bursts}
\description{

  The function \code{cutltraj} split the bursts in an object of class
  \code{ltraj} into several "sub-bursts", according to some specified
  criterion.
  
  The function \code{bindltraj} binds the bursts an object of class
  \code{ltraj} with the same attributes \code{"id"} into one unique
  burst.
}
\usage{
cutltraj(ltraj, criterion, value.NA = FALSE, nextr = TRUE, ...)
bindltraj(ltraj, ...)
}
\arguments{
  \item{ltraj}{an object of class \code{ltraj}}
  \item{criterion}{a character string giving any syntactically correct R
    logical expression implying the descriptive parameters in \code{x}}
  \item{value.NA}{logical.  The value that should be used to replace the
    missing values.}
  \item{nextr}{logical.  Whether the current "sub-burst" should stop
    after (\code{nextr = TRUE}) or before (\code{nextr = FALSE}) the
    first relocation matching \code{criterion}}
  \item{\dots}{additional arguments to be passed to other functions}
}
\details{
 
  Splitting a trajectory may be of interest in many situations. For
  example, if it is known that two kinds of activities of the monitored
  animals correspond to different properties of, say, the distance
  between successive relocations, it may be of interest to split the
  trajectory according to the values of these distances.

  The criterion used to cut the trajectory may imply any of the
  parameters describing a trajectory in the object \code{ltraj} (e.g.,
  \code{"dt"}, \code{"dist"}, \code{"dx"}, etc. see the help page of
  \code{as.ltraj}), as well as any variable stored in the attribute
  \code{"infolocs"} of the object.
  
  Two options are available in \code{cutltraj}, depending on the value
  of \code{nextr}. If \code{nextr = FALSE}, any sequence of
  successive relocations that *do not* match the criterion is considered
  as a new burst.  For example, if for a given burst, the criterion
  returns the vector \code{(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE,
    FALSE, FALSE, FALSE)}, then the function \code{cutltraj} creates
  two new bursts of relocations, the first one containing the first 3
  relocations and the second one the last 3 relocations.
  
  If \code{nextr = TRUE}, any sequence of successive relocations that
  *do not* match the criterion, *as well as the first relocation that
  does match it after this sequence* is considered as a new burst. 
  This option is available because many of the descriptive parameters
  associated to a given relocation in an object of class \code{ltraj}
  measure some specific feature concerning the position of the next
  relocation.  For example, one may want to consider as a burst any
  sequence of relocations for which the time lag is below one hour (the
  criterion is \code{"dt > 3600"}.  The first relocation for which this
  criterion is TRUE belong to the burst, and it is the next one which
  is excluded from the burst.  For example, if for a given burst, the
  criterion returns the vector \code{(FALSE, FALSE, FALSE, TRUE, TRUE,
    TRUE, FALSE, FALSE, FALSE)}, then the function \code{cutltraj}
  creates two new bursts of relocations, the first one containing the
  first 4 relocations and the second one the last 3 relocations.
  
}
\value{
  An object of class \code{ltraj}.
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }
\seealso{\code{\link{ltraj}} for additional information about objects
of class \code{ltraj} (and especially concerning the names of the
descriptive parameters that can be used in
\code{cutltraj}).  \code{\link{is.sd}} (especially the examples of this
help page) for other examples of use of this function}
\examples{
\dontrun{
#######################################################
##
## GPS monitoring of one bear

data(bear)

## We want to study the trajectory of the day at the scale
## of the day. We define one trajectory per day. The trajectory should begin
## at 22H00
## The following function returns TRUE if the date is comprised between
## 21H00 and 22H00 (i.e. correspond to the relocation taken at 21H30)

foo <- function(date) {
     da <- as.POSIXlt(date, "UTC")
     ho <- da$hour + da$min/60
     return(ho>21&ho<22)
}

## We cut the trajectory into bursts after the relocation taken at 21H30:

bea1 <- cutltraj(bear, "foo(date)", nextr = TRUE)
bea1

## Remove the first and last burst:
bea2 <- bea1[-c(1,length(bea1))]


#######################################################
##
## Bind the trajectories

bea3 <- bindltraj(bea2)
bea3
}
}
\keyword{programming}
