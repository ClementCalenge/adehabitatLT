\name{ld}
\alias{ld}
\alias{dl}
\title{
  Quick Conversion of Objects of Class ltraj from and to Dataframes
}
\description{
  The two functions \code{ld} and \code{dl} are useful to quickly
  convert objects of class ltraj from and to dataframes.
}
\usage{
ld(ltraj)
dl(x)
}
\arguments{
  \item{ltraj}{
    an object of class \code{ltraj}
  }
  \item{x}{
    an object of class \code{data.frame}, containing at least columns
    named \code{x,y,date}.
  }
}
\details{
  
  The function \code{ld} concatenates all bursts in an object of class
  \code{ltraj}, adds two columns named \code{id} and \code{burst}, and,
  when it is present, also adds the variables in the \code{infolocs}
  component.
  
  The function \code{dl} creates an object of class \code{ltraj} from a
  \code{data.frame}.  If no column named \code{id} exists, a random ID
  is generated.  If no column named \code{burst} exists, the ID is used
  as burst.  The columns named "dx", "dy", "dist", "dt", "R2n",
  "abs.angle" and "rel.angle" are recomputed by the function (see
  \code{?as.ltraj}).  Additional columns are used as the \code{infolocs}
  component.
}
\value{
  \code{ld} returns an object of class \code{data.frame}.\cr 
  \code{dl} returns an object of class \code{ltraj}. 
}
\author{
  Clement Calenge \email{clement.calenge@oncfs.gouv.fr}
}
\seealso{
  \code{\link{as.ltraj}} for additional information about objects of
  class \code{ltraj}
}
\examples{
 data(puechcirc)
 puechcirc ## class ltraj
 uu <- ld(puechcirc)
 head(uu)
 dl(uu)
}
\keyword{spatial}
