\name{albatross}
\alias{albatross}
\docType{data}
\title{Argos Monitoring of Adult Albatross Movement}
\description{
  This data set contains the relocations of 6 adult albatross monitored
  in the Crozets Islands by the team of H. Weimerskirch from the
  CEBC-CNRS (Centre d'Etudes Biologiques de Chize, France).
}
\usage{data(albatross)}
\format{
  This data set is an object of class \code{ltraj}.
}
\details{
  The coordinates are given in meters (UTM - zone 42).
}
\source{
  http://suivi-animal.u-strasbg.fr/index.htm
}
\examples{
data(albatross)

plot(albatross)

}
\keyword{datasets}
