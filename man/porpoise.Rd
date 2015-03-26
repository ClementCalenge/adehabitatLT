\name{porpoise}
\alias{porpoise}
\docType{data}
\title{Argos monitoring of Porpoise Movements}
\description{
  This data set contains the relocations of 3 porpoises
}
\usage{data(porpoise)}
\format{
  This data set is a regular object of class \code{ltraj} (i.e. constant
  time lag of 24H).
}
\details{
    The coordinates are given in meters (UTM - zone 19).
}
\source{
  http://whale.wheelock.edu/
}
\examples{
data(porpoise)

plot(porpoise)

}
\keyword{datasets}
