\name{hist.ltraj}
\alias{hist.ltraj}
\title{Histogram of the Descriptive Parameters of a Trajectory}
\description{
  This function draws an histogram of any tranformation of the
  descriptive parameters of a trajectory in objects of class
  \code{ltraj}.
}
\usage{
\method{hist}{ltraj}(x, which = "dx/sqrt(dt)", ...)
}
\arguments{
  \item{x}{an object of class \code{ltraj}}
  \item{which}{a character string giving any syntactically correct R
    expression implying the descriptive elements in \code{x}, or the
    variables in the optional attribute \code{infolocs} .}
  \item{\dots}{parameters to be passed to the generic function
    \code{hist}.}
}
\value{
  a list of objects of class \code{"histogram"}  
}
\author{Clement Calenge \email{clement.calenge@ofb.gouv.fr}}
\seealso{\code{\link{hist}}, \code{\link{ltraj}} for additional
  information on the descriptive parameters of the trajectory,
  \code{\link{qqnorm.ltraj}} for examination of distribution. }
\examples{

## Simulation of a Brownian Motion
a <- simm.brown(c(1:300, seq(301,6000,by=20)))
plot(a, addpoints = FALSE)


## dx/sqrt(dt) and dy/sqrt(dt) are normally distributed (see
## ?qqchi) 
hist(a, "dx/sqrt(dt)", freq = FALSE)
lines(tutu <- seq(-5,5, length=50), dnorm(tutu), col="red")

hist(a, "dy/sqrt(dt)", freq = FALSE)
lines(tutu, dnorm(tutu), col="red")


## Look at the distribution of distances between
## successive relocations
hist(a, "dist/sqrt(dt)", freq = FALSE)
lines(tutu <- seq(0,5, length=50), dchi(tutu), col="red")
}
\keyword{hplot}

