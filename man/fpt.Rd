\name{fpt}
\alias{fpt}
\alias{plot.fipati}
\alias{varlogfpt}
\alias{meanfpt}
\title{Computation of the First Passage Time From Trajectories}
\description{
  These functions compute the first passage time using trajectories of
  class \code{"ltraj"} of type II (time recorded).
}
\usage{
fpt(lt, radii, units = c("seconds", "hours", "days"))
varlogfpt(f, graph = TRUE)
meanfpt(f, graph = TRUE)
\method{plot}{fipati}(x, scale, warn = TRUE, \dots)
}
\arguments{
  \item{lt}{an object of class \code{"ltraj"} of type II (time
    recorded)}
  \item{radii}{a numeric vector giving the radii of the circles}
  \item{units}{The time units of the results}
  \item{f,x}{an object of class \code{fipati} returned by the function
    \code{fpt}}
  \item{graph}{logical.  Whether the results should be plotted}
  \item{scale}{the value of the radius to be plotted}
  \item{warn}{logical. Whether the function should warn the user when
    the given scale does not correspond to possible radii available in
    the object of class \code{fipati}}
  \item{\dots}{additional arguments to be passed to the generic function
    \code{plot}}
}
\details{

  The first passage time (FPT) is a parameter often used to describe the
  scale at which patterns occur in a trajectory.  For a given scale r,
  it is defined as the time required by the animals to pass through a
  circle of radius r.  Johnson et al. (1992) indicated that the mean
  first passage time scales proportionately to the square of the radius
  of the circle for an uncorrelated random walk.  They used this
  property to differenciate facilitated diffusion and impeded diffusion,
  according to the value of the coefficient of the linear regression
  \code{log(FPT) = a * log(radius) + b}. Under the hypothesis of a
  random walk, \code{a} should be equal to 2 (higher for impeded
  diffusion, and lower for facilitated diffusion).  Note however, that
  the value of a converges to 2 only for large values of radius.
  
  Fauchald & Tveraa (2003) proposed another use of the FPT.  Instead of
  computing the mean of FPT, they propose the use of the variance of the
  log(FPT).  This variance should be high for scales at which patterns
  occur in the trajectory (e.g. area restricted search).  This method is
  often used to determine the scale at which an animal seaches for food.

}
\value{
  \code{fpt} computes the FPT for each relocation and each radius, and
  for each animals.  This function returns an object of class
  \code{"fipati"}, i.e. a list with one component per animal.  Each
  component is a data frame with each column corresponding to a value
  of \code{radii} and each row corresponding to a relocation.  An object
  of class \code{fipati} has an attribute named \code{"radii"}
  corresponding to the argument \code{radii} of the function
  \code{fpt}.

  \code{meanfpt} and \code{varlogfpt} return a data frame giving
  respectively the mean FPT and the variance of the log(FPT) for each
  animal (rows) and rach radius (column).  These objects also have an
  attribute \code{"radii"}.
}
\references{
  Johnson, A. R., Milne, B.T., & Wiens, J.A. (1992) Diffusion in
  fractal landscapes: simulations and experimental studies of
  tenebrionid beetle movements. \emph{Ecology} \bold{73}: 1968--1983.
  
  Fauchald, P. & Tveraa, T. (2003) Using first passage time in the
  analysis of area restricted search and habitat selection. \emph{Ecology}
  \bold{84}: 282--288. 
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}
\seealso{\code{\link{ltraj}} for additional information on objects of
  class \code{ltraj}}
\examples{

data(puechcirc)
i <- fpt(puechcirc, seq(300,1000, length=30))
plot(i, scale = 500, warn = FALSE)

toto <- meanfpt(i)
toto
attr(toto, "radii")


toto <- varlogfpt(i)
toto
attr(toto, "radii")


}
\keyword{spatial}
