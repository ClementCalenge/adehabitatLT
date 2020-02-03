\name{trajdyn}
\alias{trajdyn}
\title{Interactive Display of Objects of Class 'ltraj'}
\description{
  This function provides an interactive version of \code{plot.ltraj}, 
  for the exploration of objects of class \code{ltraj}.  
}
\usage{
trajdyn(x, burst = attr(x[[1]], "burst"), hscale = 1, vscale = 1,
        recycle = TRUE, display = c("guess", "windows", "tk"), ...) 
}
\arguments{
  \item{x}{an object of class \code{ltraj}}
  \item{burst}{a character string indicating the burst identity to
    explore}
  \item{hscale}{passed to tkrplot}
  \item{vscale}{passed to tkrplot}
  \item{recycle}{logical.  Whether the trajectory should be recycled
  at the end of the display}
  \item{display}{type of display.  The default \code{guess} uses a
    windows graphics device if \code{getOption('device')=='windows'},
    otherwise it uses tk (requiring the \code{tkrplot} package).}
  \item{\dots}{additional arguments to be passed to the function
    \code{plot.ltraj}.}
}
\author{ Clement Calenge \email{clement.calenge@ofb.gouv.fr} }
\seealso{\code{\link{ltraj}} for further information on the class
    \code{ltraj}, and \code{\link{plot.ltraj}} for information on arguments
    that can be passed to this function.}
\examples{
\dontrun{

## Without map
data(puechcirc)
trajdyn(puechcirc)


## With map
data(puechabonsp)
trajdyn(puechcirc, spixdf = puechabonsp$map)

}
}
\keyword{hplot}
\keyword{spatial}
