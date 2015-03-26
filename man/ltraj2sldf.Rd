\name{ltraj2spdf}
\alias{ltraj2spdf}
\alias{ltraj2sldf}
\title{ Conversion of the class "ltraj" to the package "sp" }
\description{
  These functions convert the class "ltraj" available in adehabitatLT toward
  classes available in the package \code{sp}.

  \code{ltraj2spdf} converts an object of class \code{ltraj} into an
  object of class \code{SpatialPointsDataFrame}.

  \code{ltraj2sldf} converts an object of class \code{ltraj} into an
  object of class \code{SpatialLinesDataFrame}.
}
\usage{
ltraj2spdf(ltr)
ltraj2sldf(ltr, byid = FALSE)
}
\arguments{
  \item{ltr}{an object of class \code{ltraj}.}
  \item{byid}{logical.  If \code{TRUE}, one objects of class
    \code{Lines} correspond to one animal.  if \code{FALSE},
    one object of class \code{Lines} correspond to one burst.}
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}
\seealso{ \code{\link{ltraj}} for objects of class
  \code{ltraj}.}
\examples{
\dontrun{
if (require(sp)) {



#########################################
##
## Conversion ltraj -> SpatialPointsDataFrame
##

data(puechcirc)
plot(puechcirc)

toto <- ltraj2spdf(puechcirc)
plot(toto)



#########################################
##
## Conversion ltraj -> SpatialLinesDataFrame
##

toto <- ltraj2sldf(puechcirc)
plot(toto)


}
}
}
\keyword{hplot}
