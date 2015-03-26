\name{hseal}
\alias{hseal}
\docType{data}
\title{Argos Monitoring of Hooded Seal}
\description{
  This data set contains the trajectory of one hooded seal.
}
\usage{data(hseal)}
\format{
  The dataset \code{hseal} is an object of class \code{ltraj}.  The
  coordinates are stored in meters (UTM - zone 21).
}
\source{
  Jonsen, I.D., Flemming, J.M. and Myers, R.A. (2005). Robust state-space
  modeling of animal movement data. \emph{Ecology}, \bold{86}, 2874--2880.  
}
\examples{
data(hseal)

plot(hseal)

}
\keyword{datasets}
