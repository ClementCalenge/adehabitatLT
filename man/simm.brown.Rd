\name{simm.brown}
\alias{simm.brown}
\title{Simulate a Bivariate Brownian Motion}
\description{
  This function simulates a Bivariate Brownian Motion. 
}
\usage{
simm.brown(date = 1:100, x0 = c(0, 0), h = 1, id = "A1", burst = id)
}
\arguments{
  \item{date}{a vector indicating the date (in seconds) at which
    relocations should be simulated.  This vector can be of class
    \code{POSIXct}}
  \item{x0}{a vector of length 2 containing the coordinates of the
    startpoint of the trajectory}
  \item{h}{Scaling parameter for the brownian motion (larger values give
  smaller dispersion)}
  \item{id}{a character string indicating the identity of the simulated
    animal (see \code{help(ltraj)})}
  \item{burst}{a character string indicating the identity of the simulated
    burst (see \code{help(ltraj)})}
}
\details{
  A bivariate Brownian motion can be described by a vector
  \code{B2(t) = (Bx(t), By(t))}, where \code{Bx} and \code{By} are
  unidimensional Brownian motions.  Let \code{F(t)} the set of all
  possible realisations of the process \code{(B2(s), 0 < s < t)}.
  \code{F(t)} therefore corresponds to the known information at time
  \code{t}.  The properties of the bivariate Brownian motion are
  therefore the following: (i) \code{B2(0)= c(0,0)} (no uncertainty at
  time \code{t = 0}); (ii) \code{B2(t) - B2(s)} is independent of
  \code{F(s)} (the next increment does not depend on the present or past
  location); (iii) \code{B2(t) - B2(s)} follows a bivariate normal
  distribution with mean \code{c(0,0)} and with variance equal to
  \code{(t-s)}.

  Note that for a given parameter \code{h}, the process \code{1/h * B2(
    t * h^2 )} is a Brownian motion.  The function \code{simm.brown}
  simulates the process \code{B2(t * h^2)}.  Note that the function
  \code{hbrown} allows the estimation of this scaling factor from data. 
}
\value{
  An object of class \code{ltraj}
}
\references{ ~put references to the literature/web site here ~ }
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}\cr
  Stephane Dray \email{dray@biomserv.univ-lyon1.fr}\cr
  Manuela Royer \email{royer@biomserv.univ-lyon1.fr}\cr
  Daniel Chessel \email{chessel@biomserv.univ-lyon1.fr}
}
\seealso{\code{\link{ltraj}}, \code{\link{hbrown}}}
\examples{
plot(simm.brown(1:1000), addpoints = FALSE)

## Note the difference in dispersion:
plot(simm.brown(1:1000, h = 4), addpoints = FALSE)

}
\keyword{spatial}
