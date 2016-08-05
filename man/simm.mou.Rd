\name{simm.mou}
\alias{simm.mou}
\title{Simulation of a Bivariate Ornstein-Uhlenbeck Process}
\description{
    This function simulates a bivariate Ornstein-Uhlenbeck process for
    animal movement.
}
\usage{
simm.mou(date = 1:100, b = c(0, 0),
         a = diag(0.5, 2), x0 = b,
         sigma = diag(2), id = "A1",
         burst = id, proj4string=CRS())
}
\arguments{
  \item{date}{a vector indicating the date (in seconds) at which
    relocations should be simulated.  This vector can be of class
    \code{POSIXct}}
  \item{b}{a vector of length 2 containing the coordinates of the
    attraction point}
  \item{a}{a 2*2 matrix}
  \item{x0}{a vector of length 2 containing the coordinates of the
    startpoint of the trajectory}
  \item{sigma}{a 2*2 positive definite matrix}
  \item{id}{a character string indicating the identity of the simulated
    animal (see \code{help(ltraj)})}
  \item{burst}{a character string indicating the identity of the simulated
    burst (see \code{help(ltraj)})}
  \item{proj4string}{a valid CRS object containing the projection
    information (see \code{?CRS} from the package \code{sp}).}
}
\details{
  The Ornstein-Uhlenbeck process can be used to take into account an
  "attraction point" into the animal movements (Dunn and Gipson 1977).
  This process can be simulated using the stochastic differential
  equation: 
  \deqn{d\mathbf{z} = \mathbf{a} ( \mathbf{b} - \mathbf{z}(t)) dt
    + \mathbf{\Sigma} d \mathbf{B2(t)}}{dz = a (b - z(t)) dt + Sigma
    dB2(t)}
  The vector \code{b} contains the coordinates of the attraction
  point. The matrix \code{a} (2 rows and 2 columns) contains
  coefficients controlling the force of the attraction.  The matrix
  \code{Sigma} controls the noise added to the movement (see
  \code{?simm.mba} for details on this matrix). 

}
\value{
  An object of class \code{ltraj}
}
\references{
  Dunn, J.E., & Gipson, P.S. (1977) Analysis of radio telemetry data in
  studies of home range. \emph{Biometrics} \bold{33}: 85--101. 
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}\cr
  Stephane Dray \email{dray@biomserv.univ-lyon1.fr}\cr
  Manuela Royer \email{royer@biomserv.univ-lyon1.fr}\cr
  Daniel Chessel \email{chessel@biomserv.univ-lyon1.fr}
}
\seealso{\code{\link{simm.brown}}, \code{\link{ltraj}},
  \code{\link{simm.crw}}, \code{\link{simm.mba}} }
\examples{
set.seed(253)
u <- simm.mou(1:50, burst="Start at the attraction point")
v <- simm.mou(1:50, x0=c(-3,3),
              burst="Start elsewhere")
w <- simm.mou(1:50, a=diag(c(0.5,0.1)), x0=c(-3,3),
              burst="Variable attraction")
x <- simm.mou(1:50, a=diag(c(0.1,0.5)), x0=c(-3,7),
              burst="Both")
z <- c(u,v,w,x)

plot(z, addpoints = FALSE, perani = FALSE)

}
\keyword{spatial}
