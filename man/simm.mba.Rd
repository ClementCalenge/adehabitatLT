\name{simm.mba}
\alias{simm.mba}
\title{Simulation of an Arithmetic Brownian Motion}
\description{
  This function simulates an Arithmetic Brownian Motion.
}
\usage{
simm.mba(date = 1:100, x0 = c(0, 0), mu = c(0, 0),
         sigma = diag(2), id = "A1", burst = id,
         proj4string=CRS())
}
\arguments{
  \item{date}{a vector indicating the date (in seconds) at which
    relocations should be simulated.  This vector can be of class
    \code{POSIXct}}
  \item{x0}{a vector of length 2 containing the coordinates of the
    startpoint of the trajectory}
  \item{mu}{a vector of length 2 describing the drift of the movement}
  \item{sigma}{a 2*2 positive definite matrix}
  \item{id}{a character string indicating the identity of the simulated
    animal (see \code{help(ltraj)})}
  \item{burst}{a character string indicating the identity of the simulated
    burst (see \code{help(ltraj)})}
  \item{proj4string}{a valid CRS object containing the projection
    information (see \code{?CRS} from the package \code{sp}).}
}
\details{

  The arithmetic Brownian motion (Brillinger et al. 2002) can be
  described by the stochastic differential equation:
  \deqn{d \mathbf{z}(t) = \mathbf{\mu} dt + \mathbf{\Sigma} d
    \mathbf{B}2(t)}{dz = mu * dt + Sigma dB2(t)}
  Coordinates of the animal at time t are contained in the vector
  \code{z(t)}.  \code{dz = c(dx, dy)} is the increment of the
  movement during dt.  \code{dB2(t)} is a bivariate brownian Motion (see
  \code{?simm.brown}).  The vector \code{mu} measures the drift of the
  motion.  The matrix \code{Sigma} controls for perturbations due to the
  random noise modeled by the Brownian motion. It can also be used to
  take into account a potential correlation between the components dx
  and dy of the animal moves during dt (see Examples).
  
}
\value{
  An object of class \code{ltraj}
}
\references{
Brillinger, D.R., Preisler, H.K., Ager, A.A. Kie, J.G. & Stewart,
B.S. (2002) Employing stochastic differential equations to
model wildlife motion. \emph{Bulletin of the Brazilian Mathematical
  Society} \bold{33}: 385--408.
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}\cr
  Stephane Dray \email{dray@biomserv.univ-lyon1.fr}\cr
  Manuela Royer \email{royer@biomserv.univ-lyon1.fr}\cr
  Daniel Chessel \email{chessel@biomserv.univ-lyon1.fr}
}
\seealso{ \code{\link{simm.brown}}, \code{\link{ltraj}},
  \code{\link{simm.crw}}, \code{\link{simm.mou}} }
\examples{

set.seed(253)
u <- simm.mba(1:1000, sigma = diag(c(4,4)), 
              burst = "Brownian motion")
v <- simm.mba(1:1000, sigma = matrix(c(2,-0.8,-0.8,2), ncol = 2),
              burst = "cov(x,y) > 0")
w <- simm.mba(1:1000, mu = c(0.1,0), burst = "drift > 0")
x <- simm.mba(1:1000, mu = c(0.1,0),
              sigma = matrix(c(2, -0.8, -0.8, 2), ncol=2),
              burst = "Drift and cov(x,y) > 0")
z <- c(u, v, w, x)
plot(z, addpoints = FALSE, perani = FALSE)

}
\keyword{spatial}
