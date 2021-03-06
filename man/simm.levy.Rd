\name{simm.levy}
\alias{simm.levy}
\title{Simulates a Levy Walk}
\description{
  This function simulates a Levy walk
}
\usage{
simm.levy(date = 1:500, mu = 2, l0 = 1, x0 = c(0, 0),
          id = "A1", burst = id, typeII = TRUE,
          proj4string=CRS())
}
\arguments{
  \item{date}{a vector indicating the date (in seconds) at which
    relocations should be simulated.  This vector can be of class
    \code{POSIXct}. *Note that the time lag between two relocations
    should be constant* (regular trajectories required)}
  \item{mu}{The exponent of the Levy distribution}
  \item{l0}{The minimum length of a step}
  \item{x0}{a vector of length 2 containing the coordinates of the
    startpoint of the trajectory}
  \item{id}{a character string indicating the identity of the simulated
    animal (see \code{help(ltraj)})}
  \item{burst}{a character string indicating the identity of the simulated
    burst (see \code{help(ltraj)})}
  \item{typeII}{logical.  Whether the simulated trajectory should be of
    type II (\code{TRUE}, time recorded) or not (\code{FALSE}, time not
    recorded). See \code{help(ltraj)}.}
  \item{proj4string}{a valid CRS object containing the projection
    information (see \code{?CRS} from the package \code{sp}).}
}
\details{
  This function simulates a Levy flight with exponent \code{mu}.
  This is done by sampling a random relative angle from a uniform
  distribution (-pi, pi) for each step, and a step length generated by
  \code{dt * (l0 * (runif(1)^(1/(1 - mu))))}
}
\value{
  an object of class \code{ltraj}
}
\references{
  Bartumeus, F., da Luz, M.G.E., Viswanathan, G.M. Catalan, J. (2005)
  Animal search strategies: a quantitative random-walk
  analysis. \emph{Ecology}, \bold{86}: 3078--3087.
}
\author{Clement Calenge \email{clement.calenge@ofb.gouv.fr}}
\seealso{\code{\link{chi}}, \code{\link[CircStats]{rwrpnorm}},
  \code{\link{simm.brown}}, \code{\link{ltraj}},
  \code{\link{simm.crw}}, \code{\link{simm.mba}},
  \code{\link{simm.levy}}}
\examples{

suppressWarnings(RNGversion("3.5.0"))
set.seed(411)
w <- simm.levy(1:500, mu = 1.5, burst = "mu = 1.5")
u <- simm.levy(1:500, mu = 2, burst = "mu = 2")
v <- simm.levy(1:500, mu = 2.5, burst = "mu = 2.5")
x <- simm.levy(1:500, mu = 3, burst = "mu = 3")
par(mfrow=c(2,2))
lapply(list(w,u,v,x), plot, perani=FALSE)

}
\keyword{spatial}
