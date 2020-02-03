\name{residenceTime}
\alias{residenceTime}
\alias{print.resiti}
\alias{plot.resiti}
\title{Trajectory Analysis using the Residence Time Method}
\description{
  These functions can be used to apply the residence time method
  (Barraquand and Benhamou, 2008).
}
\usage{
residenceTime(lt, radius, maxt, addinfo = FALSE,
              units = c("seconds", "hours", "days"))

\method{print}{resiti}(x, \dots)

\method{plot}{resiti}(x, addpoints = FALSE, addlines = TRUE, \dots)
}
\arguments{
  \item{lt}{an object of class \code{ltraj}}
  \item{radius}{the radius of the patch (in units of the coordinates)}
  \item{maxt}{maximum time threshold that the animal is allowed to spend
    outside the patch before that we consider that the animal actually
    left the patch (see Details)} 
  \item{addinfo}{logical value.  If \code{TRUE}, then the residence time
    method is added as a variable in the \code{infolocs} component of
    the object \code{lt}.  If \code{FALSE} this function returns an
    object of class \code{resiti}}
  \item{units}{a character string indicating the time units of
    \code{maxt}}
  \item{x}{an object of class \code{"resiti"}}
  \item{addpoints}{logical.  Whether points should be added to the plot.}
  \item{addlines}{logical.  Whether lines should be added to the plot.}
  \item{\dots}{additionnal arguments to be passed to or from other methods}
}
\details{
  Barraquand and Benhamou (2008) proposed a new approach to identify the
  places where the animals spend the most of their time, relying on the
  calculation of their residence time in the various places where they
  have been relocated.  This approach is similar to the
  first passage time method: for a given value of \code{radius} and for
  a given relocation, the first passage time is defined as the time
  required by the animal to pass through a circle of given radius
  centred  on the relocation  (see the help page of the function
  \code{fpt} for additional details).  The residence time associated to
  a given relocation corresponds to the first passage time calculated at
  this place plus the passage times that occurred in this circle
  before or after the current relocation, *given* that 
  the animal did not spent a time greater than \code{maxt} before
  reentering the circle (see Barraquand and Benhamou, 2008, for
  details).  It is therefore computed by determining the various times
  at which the path intersects the perimeter of the circle
  centred on the current relocation, both forward and backward, and
  then by summing the durations associated with the
  various portions of the path occurring within the circle. The
  graphical examination of the changes with time allow to identify the
  dates and places where the animal spent most of its time.

  A partitionning method can be used to segment the series formed by the
  residence time into homogeneous segments.  Barraquand and Benhamou
  (2008) propose the method of Lavielle (1999, 2005). See the function
  \code{lavielle} for details about this method.
}
\value{
  If \code{addinfo = FALSE}, the function \code{residenceTime} returns a
  list of class \code{"resiti"} where each element corresponds to a
  burst of the object \code{lt}.  Each element is a \code{data.frame}
  with two columns: the date and the residence time associated with the
  date.
}
\references{
  Barraquand, F. and Benhamou, S. (2008) Animal movement in
  heterogeneous landscapes: identifying profitable places and
  homogeneous movement bouts. \emph{Ecology}, \bold{89}, 3336--3348. 
}
\author{Clement Calenge \email{clement.calenge@ofb.gouv.fr}}
\seealso{
  \code{\link{lavielle}} for the partitionning of the trajectory based
  on the residence time.
}
\examples{
\dontrun{
data(albatross)
ltr <- albatross[1]

## show the distances between successive relocations as a function
## of date
plotltr(ltr)

## focus on the first period
ltr <- gdltraj(ltr, as.POSIXct("2001-12-15", tz="UTC"),
               as.POSIXct("2003-01-10", tz="UTC"))

plot(ltr)
## We identify places that seem to be a patch and, with locator,
## we measure approximately their size.
## The approximate patch radius can be set equal to 100 km as a first try

plotltr(ltr, "dt")
## As a first try, we could set maxt equal to 15000 seconds, i.e.
## approximately 4 hours

## calculation of the residence time
res <- residenceTime(ltr, radius = 100000, maxt=4, units="hour")
plot(res)

## There seems to be about 10 segments. Let us try the method
## of Lavielle (1999, 2005) to segment this series:
## First calculate again the residence time as the infolocs attribute
## of the trajectory
res <- residenceTime(ltr, radius = 100000, maxt=4, addinfo = TRUE, units="hour")
res

## Note that the residence time is now an attribute of the infolocs
## component of res


## Now, use the Lavielle method, with Kmax set to 2-3 times the
## "optimal" number of segments, assessed visually according
## to the recommendations of Barraquand and Benhamou (2008)
## We set the minimum number of relocations in each segment to
## 10 observations (given that the relocations were theoretically
## taken every hour, this defines a patch as a place where the animal
## stays at least 10 hours: this also defines the scale of our study)

ii <- lavielle(res, which="RT.100000", Kmax=20, Lmin=10)

## Both the graphical method and the automated method to choose
## the optimal number of segments indicate 4 segments
## (see ?lavielle for a description of these methods):

chooseseg(ii)

## We identify the 4 segments: the method of Lavielle seems to do a good
## job:
(pa <- findpath(ii, 4))

## and we plot this partition:
plot(pa, perani=FALSE)


## Now, we could try a study at a smaller scale (patch = 50km):
res <- residenceTime(ltr, radius = 50000, maxt=4, addinfo = TRUE,
                     units="hour")
ii <- lavielle(res, which="RT.50000", Kmax=20, Lmin=10)

## 5 segments seem a good choice:
chooseseg(ii)

## There is more noise in the residence time, but
## the partition is still pretty clear:
(pa <- findpath(ii, 5))

## show the partition:
plot(pa, perani = FALSE)


## Now try at a larger scale (patch size=250 km)
res <- residenceTime(ltr, radius = 250000, maxt=4, addinfo = TRUE,
                     units="hour")
ii <- lavielle(res, which="RT.250000", Kmax=15, Lmin=10)

## 5 segments seem a good choice again:
chooseseg(ii)

## There is more noise in the residence time, but
## the partition is still pretty clear:
(pa <- findpath(ii, 5))

## show the partition:
plot(pa, perani = FALSE)

}
}
\keyword{chron}
\keyword{ts}
