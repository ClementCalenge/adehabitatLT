\name{plotltr}
\alias{plotltr}
\title{Changes in Traject Parameters Over Time}
\description{
  This function allows a graphical examination of the changes in
  descriptive parameters in objects of class \code{ltraj}
}
\usage{
plotltr(x, which = "dist", pch = 16, cex = 0.7, addlines = TRUE,
        addpoints = TRUE,...)
}
\arguments{
  \item{x}{An object of class \code{ltraj}}
  \item{which}{a character string giving any syntactically correct R
    expression implying the descriptive elements in \code{x} or the
    variables in the optional attribute \code{infolocs}.}
  \item{pch}{the type of points on the plot (see \code{help(par)}).}
  \item{cex}{the size of points on the plot (see \code{help(par)}).}
  \item{addlines}{logical.  Indicates whether lines should be added to
    the plot.}
  \item{addpoints}{logical.  Indicates whether points should be added to
    the plot.}
  \item{\dots}{additional parameters to be passed to the generic
    function \code{plot}}
}
\author{Clement Calenge \email{clement.calenge@ofb.gouv.fr}}
\seealso{ \code{\link{ltraj}} for additional information about
  objects of class \code{ltraj}, and \code{\link{sliwinltr}} for a
  sliding window smoothing}
\examples{
data(puechcirc)

plotltr(puechcirc, "cos(rel.angle)")
plotltr(puechcirc, "dist")
plotltr(puechcirc, "dx")

}
\keyword{hplot}
