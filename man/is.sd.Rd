\name{is.sd}
\alias{is.sd}
\alias{sd2df}
\title{Handling of Trajectories of the Same Duration}
\description{
  \code{is.sd} tests whether the bursts of relocations in an object of
  class \code{ltraj} contain the same number of relocations, and cover
  the same duration ("sd" = "same duration").\cr
  \code{sd2df} gets one of the descriptive parameters of a regular "sd"
  trajectory (e.g. \code{"dt"}, \code{"dist"}, etc.) and returns a data
  frame with one relocation per row, and one burst per column. 
}
\usage{
is.sd(ltraj)
sd2df(ltraj, what)
}
\arguments{
  \item{ltraj}{an object of class \code{ltraj}}
  \item{what}{a character string indicating the descriptive parameter of
    the trajectory to be exported}
}
\value{
  \code{is.sd} returns a logical value.\cr
  \code{sd2df} returns a data frame with one column per burst of
  relocations, and one row per relocation.
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}
\seealso{\code{\link{set.limits}} for additional information about
  "sd" regular trajectories }
\examples{
\dontrun{
## Takes the example from the help page of cutltraj (bear):
data(bear)

## We want to study the trajectory of the animal at the scale
## of the day. We define one trajectory per day. The trajectory should begin
## at 22H00.
## The following function returns TRUE if the date is comprised between
## 21H00 and 22H00 and FALSE otherwise (i.e. correspond to the
## relocation taken at 21H30)

foo <- function(date) {
     da <- as.POSIXlt(date, "GMT")
     ho <- da$hour + da$min/60
     return(ho>21.1&ho<21.9)
}

## We cut the trajectory into bursts after the relocation taken at 21H30:

bea1 <- cutltraj(bear, "foo(date)", nextr = TRUE)
bea1

## Remove the first and last burst:
bea2 <- bea1[-c(1,length(bea1))]

## Is the resulting object "sd" ?
is.sd(bea2)

## Converts to data frame:
df <- sd2df(bea2, "dist")

## Plots the average distance per hour
meandi <- apply(df[-nrow(df),], 1, mean, na.rm = TRUE)
sedi <- apply(df[-nrow(df),], 1, sd, na.rm = TRUE) / sqrt(ncol(df))
plot(seq(0, 23.5, length = 47),
     meandi,
     ty = "b", pch = 16, xlab = "Hours (time 0 = 22H00)",
     ylab="Average distance covered by the bear in 30 mins",
     ylim=c(0, 500))
lines(seq(0, 23.5, length = 47),
      meandi+sedi, col="grey")
lines(seq(0, 23.5, length = 47),
      meandi-sedi, col="grey")


}
}
\keyword{programming}
