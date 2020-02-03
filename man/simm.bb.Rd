\name{simm.bb}
\alias{simm.bb}
\title{Brownian bridge motion}
\description{
  This function simulates a brownian bridge motion
}
\usage{
simm.bb(date = 1:100, begin = c(0, 0), end = begin, id = "A1",
        burst = id, proj4string=CRS())
}
\arguments{
  \item{date}{a vector indicating the date (in seconds) at which
    relocations should be simulated.  This vector can be of class
    \code{POSIXct}}
  \item{begin}{a vector of length 2 giving the x and y coordinates of
    the location beginning of the trajectory}
  \item{end}{a vector of length 2 giving the x and y coordinates of
    the location ending the trajectory}
  \item{id}{ a character string indicating the identity of the simulated
    animal (see \code{help(ltraj)})}
  \item{burst}{ a character string indicating the identity of the
    simulated animal (see \code{help(ltraj)})}
  \item{proj4string}{a valid CRS object containing the projection
    information (see \code{?CRS} from the package \code{sp}).}
}
\value{
  An object of class \code{ltraj}.
}
\author{
  Clement Calenge \email{clement.calenge@ofb.gouv.fr}\cr
  Stephane Dray \email{dray@biomserv.univ-lyon1.fr}\cr
  Manuela Royer \email{royer@biomserv.univ-lyon1.fr}\cr
  Daniel Chessel \email{chessel@biomserv.univ-lyon1.fr}
}
\seealso{\code{\link{ltraj}}, \code{\link{hbrown}}}
\examples{

plot(simm.bb(1:1000, end=c(100,100)), addpoints = FALSE)

}
\keyword{spatial}

