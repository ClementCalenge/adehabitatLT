\name{buffalo}
\alias{buffalo}
\docType{data}
\title{
  GPS monitoring of a buffalo
}
\description{
  This data set contains the relocations of an African buffalo (Syncerus
  caffer) monitored in the W National Park (Niger) by D. Cornelis, as
  well as the habitat map of the study area.
}
\usage{data(buffalo)}
\format{
  This dataset is a list containing an object of class \code{ltraj} and
  a SpatialPixelsDataFrame.
}
\details{
  The "infolocs" component of the ltraj stores the proportion of the
  time duration between relocation i-1 and relocation i during which the
  animal was active.
}
\source{
  Cornelis D., Benhamou S., Janeau G., Morellet N., Ouedraogo M. & de Vissher
  M.-N. (submitted). The spatiotemporal segregation of limiting resources
  shapes space use patterns of West African savanna buffalo. Journal of
  Mammalogy.  
}
\examples{
data(buffalo)
plot(buffalo$traj, spixdf=buffalo$habitat)

}
\keyword{datasets}
