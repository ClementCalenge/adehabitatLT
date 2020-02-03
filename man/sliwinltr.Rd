\name{sliwinltr}
\alias{sliwinltr}
\title{Apply a Function on an Object of Class "ltraj", Using a Sliding
  Window}
\description{
  This function applies a function on an object of class "ltraj", using
  a sliding window.
}
\usage{
sliwinltr(ltraj, fun, step, type = c("locs", "time"),
          units = c("sec", "min", "hour", "day"),
          plotit = TRUE, \dots)
}
\arguments{
  \item{ltraj}{an object of class \code{ltraj}}
  \item{fun}{the function to be applied, implying at least one of the
    descriptive  parameters in the object of class \code{ltraj} (see
    below)}
  \item{step}{the half-width of the sliding window.  If
    \code{type=="locs"}, it is a number of relocations.  If
    \code{type=="time"} it is a number described by \code{units}}
  \item{type}{character string.  If \code{type == "locs"}, \code{step}
    describes a number of relocations: if \code{type == "time"},
    \code{step} describes a time lag.}
  \item{units}{if \code{type == "time"}, the time units described by
    \code{step}.  Ignored otherwise}
  \item{plotit}{logical.  Whether the result should be plotted}
  \item{\dots}{additional arguments to be passed to the function
    \code{rec}}
}
\details{
  An object of class ltraj is a list with one component per burst of
  relocations.  The function \code{fun} is applied to each burst of
  relocations.  This burst of relocations should be refered as \code{x}
  in \code{fun}.  For example, to compute the mean of the distance
  between successive relocations, the function \code{fun} is equal to
  \code{function(x) mean(x$dist)}.

  Do not forget that some of the descriptive parameters in the object
  \code{ltraj} may contain missing values (see
  \code{help(ltraj)}).  The function should therefore specify how to
  manage these missing values.
}
\value{
  If \code{type=="locs"}, a list with one component per burst of
  relocation containing the smoothed values for each relocation.

  If \code{type=="locs"}, a list with one component per burst of
  relocation.  Each component is a data frame containing the time and
  the corresponding smoothed values for each date.
}
\author{Clement Calenge \email{clement.calenge@ofb.gouv.fr}}
\seealso{\code{\link{ltraj}} for additional information about objects
      of class \code{ltraj}}
\examples{
\dontrun{
data(capreotf)

## computes the average speed of the roe deer in a moving window of width
## equal to 60 minutes
toto <- sliwinltr(capreotf, function(x) mean(x$dist/x$dt, na.rm = TRUE),
                  step = 30, type = "time", units = "min")

## zoom before the peak
head(toto[[1]])
plot(toto[[1]][1:538,], ty="l")
}
}
\keyword{hplot}
