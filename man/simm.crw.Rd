\name{simm.crw}
\alias{simm.crw}
\title{Simulation of a Correlated Random Walk}
\description{
  This function simulates a correlated random walk
}
\usage{
simm.crw(date=1:100, h = 1, r = 0,
         x0=c(0,0), id="A1", burst=id,
         typeII=TRUE, proj4string=CRS())
}
\arguments{
  \item{date}{a vector indicating the date (in seconds) at which
    relocations should be simulated.  This vector can be of class
    \code{POSIXct}. *Note that the time lag between two relocations
    should be constant* (regular trajectories required)}
  \item{h}{the scaling parameter for the movement length}
  \item{r}{The concentration parameter for wrapped normal distribution
    of turning angles}
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
  Since the seminal paper of Kareiva and Shigesada (1983), most
  biologists describe the trajectories of an animal with the help of
  two distributions: the distribution of distances between successive
  relocations, and the distribution of turning angles between successive
  moves (relative angles in the class \code{ltraj}).  The CRW is
  built iteratively.  At each step of the simulation process,
  the orientation of the move is drawn from a wrapped normal
  distribution (with concentration parameter \code{r}).  The length of
  the move is drawn from a chi distribution, multiplied by \code{h *
    sqrt(dt)}. \code{h} is a scale parameter (the same as in the
  function \code{simm.brown()}, and the distribution is
  multiplied by sqrt(t) to make it similar to the discretized Brownian
  motion if \code{r == 0}.
}
\value{
  an object of class \code{ltraj}
}
\references{
  Kareiva, P. M. & Shigesada, N. (1983) Analysing insect movement as a
  correlated random walk. \emph{Oecologia}, \bold{56}: 234--238.
}
\author{Clement Calenge \email{clement.calenge@ofb.gouv.fr}\cr
  Stephane Dray \email{dray@biomserv.univ-lyon1.fr}\cr
  Manuela Royer \email{royer@biomserv.univ-lyon1.fr}\cr
  Daniel Chessel \email{chessel@biomserv.univ-lyon1.fr}
}
\note{
  This function requires the package \code{CircStats}.
}
\seealso{\code{\link{chi}}, \code{\link[CircStats]{rwrpnorm}},
  \code{\link{simm.brown}}, \code{\link{ltraj}},
  \code{\link{simm.crw}}, \code{\link{simm.mba}}}
\examples{

suppressWarnings(RNGversion("3.5.0"))
set.seed(876)
u <- simm.crw(1:500, r = 0.99, burst = "r = 0.99")
v <- simm.crw(1:500, r = 0.9, burst = "r = 0.9", h = 2)
w <- simm.crw(1:500, r = 0.6, burst = "r = 0.6", h = 5)
x <- simm.crw(1:500, r = 0, burst = "r = 0 (Uncorrelated random walk)",
              h = 0.1)
z <- c(u, v, w, x)
plot(z, addpoints = FALSE, perani = FALSE)

}
\keyword{spatial}

