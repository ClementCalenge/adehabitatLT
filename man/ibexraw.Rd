\name{ibexraw}
\alias{ibexraw}
\docType{data}
\title{GPS Monitoring of Four Ibex in the Belledonne Mountain (irregular
data)}
\description{
  This dataset is an object of class "ltraj" (irregular trajectory,
  relocations roughly every 4 hours) containing the raw GPS relocations
  of four ibex during 15 days in the Belledonne mountain (French Alps).
}
\usage{data(ibexraw)}
\details{
  This dataset is nearly the same as the dataset \code{ibex}, except
  that the timing of relocations has not been rounded and the missing
  values have not been placed in the trajectory.
}
\source{
  Office national de la chasse et de la faune sauvage, CNERA Faune
  de Montagne, 95 rue Pierre Flourens, 34000 Montpellier, France.
}
\examples{
data(ibexraw)
plot(ibexraw)
}
\keyword{datasets}
