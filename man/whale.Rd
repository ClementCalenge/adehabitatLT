\name{whale}
\alias{whale}
\docType{data}
\title{Argos Monitoring of Whale Movement}
\description{
  This data set contains the relocations of one right whale.
}
\usage{data(whale)}
\format{
  This data set is a regular object of class \code{ltraj} (i.e. constant
  time lag of 24H00)
}
\details{
  The coordinates are given in decimal degrees (Longitude - latitude).
}
\source{
  http://whale.wheelock.edu/Welcome.html
}
\examples{
data(whale)

plot(whale)

}
\keyword{datasets}
