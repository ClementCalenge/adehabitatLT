\name{capreotf}
\alias{capreotf}
\docType{data}
\title{GPS Monitoring of one Roe Deer in Trois-Fontaines (France)}
\description{
  This dataset contains the relocations of a roe deer collected from
  May 1st to May 4th 2004 (every 5 minutes) using GPS collars in the
  wildlife reserve of Trois-Fontaines (Haute Marne, France) by the ONCFS
  (Office national de la chasse et de la faune sauvage).
}
\usage{data(capreotf)}
\format{
  This dataset is a regular object of class \code{ltraj} (i.e. constant
  time lag).
}
\source{
  Sonia Said, Office national de la chasse et de la faune sauvage,
  CNERA-CS, 1 place Exelmans, 55000 Bar-le-Duc (France).
}
\examples{
data(capreotf)

plot(capreotf)

}
\keyword{datasets}
