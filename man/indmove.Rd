\name{indmove}
\alias{indmove}
\alias{indmove.detail}
\alias{testang.ltraj}
\alias{testdist.ltraj}
\title{Testing Independence in Regular Trajectory Parameters}
\description{
  The function \code{indmove} tests for the independence between
  successive components \code{c(dx, dy)} for each burst in a regular
  object of class \code{ltraj}.
  
  The function \code{indmove.detail} tests for the independence between
  successive \code{dx} or \code{dy} for each burst in a regular object 
  of class \code{ltraj}.

  The function \code{testang.ltraj} tests for the independence between
  successive angles (relative or absolute) for each burst in a regular
  object of class \code{ltraj}.
  
  The function \code{testdist.ltraj} tests for the independence between
  successive distances between successive relocations for each burst in
  a regular object of class \code{ltraj}.
}
\usage{
indmove(ltr, nrep = 200, conflim = seq(0.95, 0.5, length=5),
        sep = ltr[[1]]$dt[1], units = c("seconds", "minutes",
                                        "hours", "days"),
        plotit = TRUE)


testang.ltraj(x, which = c("absolute", "relative"),
              nrep = 999, alter = c("two-sided","less","greater"))


testdist.ltraj(x, nrep = 999, alter = c("two-sided","less","greater"))


indmove.detail(x, detail=c("dx","dy"), nrep=999,
               alter = c("two-sided","less","greater"))
}
\arguments{
  \item{ltr,x}{an object of class \code{ltraj}}
  \item{conflim}{a vector giving the limits of the confidence intervals
    to be plotted}
  \item{nrep}{number of simulations}
  \item{units}{a character string indicating the time units for the result}
  \item{alter}{a character string specifying the alternative hypothesis,
  must be one of "greater", "less" or "two-sided" (default)}
  \item{which}{a character string indicating whether the absolute or
    relative angles are under focus}
  \item{detail}{a character string indicating whether \code{"dx"} or
    \code{"dy"} should be tested for independence}
  \item{plotit}{logical.  Whether the results should be plotted on a graph}
  \item{sep}{used in the case of variable time lag between
    relocations.  Indicates the theoretical time lag between two relocations}
}
\details{
  The function \code{indmove} randomises the order of the increments
  \code{c(dx, dy)} in a trajectory.  The criteria of the test is the
  Mean Squared Displacement (\code{R^2_n}) (Root & Kareiva 1984).
  
  The function \code{testang.ltraj} randomises the order of the angles in a
  trajectory.  The criteria of the test is \code{f^2 = sum_(i=1)^(n-1) 2*(1 -
    cos(angle[i+1] - angle[i]))}.  This measure corresponds to the
  mean squared length of the segment joining two successive angles on
  the trigonometric circle (see examples for an illustration).

  The function \code{testdist.ltraj} randomises the order of the
  distances between successive relocations in a trajectory.  The
  criteria of the test is \code{sum_(i=1)^(n-1) (dist[i+1] -
    dist[i])^2 } (Neuman 1941, Neuman et al. 1941).  The same criteria is
  used in \code{indmove.detail()}.

  Note that these functions require "regular" trajectories,
  i.e. trajectories for which the relocations are separated by a
  constant time lag.

  Finally, note that the functions \code{testang.ltraj} and
  \code{testdist.ltraj} are not affected by the presence of missing
  values in the bursts of relocations. The function \code{indmove} may
  be greatly affected by these missing values (they are removed prior to
  the test).
}
\value{
  \code{indmove()} returns a list with one component per burst.  Each
  component is a list of two data frames. The data frame \code{Time}
  contains the time points at which R2n is computed for the
  observation (first column) and the simulations (other ones).  The data
  frame \code{R2n} contains the values for the R2n (same dimensions).
  
  \code{testang.ltraj()}, \code{testdist.ltraj} and
  \code{indmove.detail} return lists of
  objects of class \code{randtest}.
}
\references{
  Root, R.B. & Kareiva, P.M. (1984) The search for resources by
  cabbage butterflies (Pieris Rapae): Ecological consequences and
  adaptive significance of markovian movements in a patchy
  environment. \emph{Ecology}, \bold{65}: 147--165.
  
  Neumann, J.V., Kent, R.H., Bellinson, H.R. & Hart, B.I. (1941) The mean
  square successive difference. \emph{Annals of Mathematical
    Statistics}, \bold{12}: 153--162.
  
  Neumann, J.V. (1941) Distribution of the ration of the mean square
  successive difference to the variance. \emph{The Annals of
    Mathematical Statistics}, \bold{12}: 367--395.
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}\cr
  Stephane Dray \email{dray@biomserv.univ-lyon1.fr}
}
\seealso{\code{\link{ltraj}}}
\examples{
\dontrun{
## theoretical independence between
br <- simm.brown(1:1000)
testang.ltraj(br)
testdist.ltraj(br)

indmove(br)
}

## Illustration of the statistic used for the test of the independence
## of the angles
opar <- par(mar = c(0,0,4,0))
plot(0,0, asp=1, xlim=c(-1, 1), ylim=c(-1, 1), ty="n", axes=FALSE,
main="Criteria f for the measure of independence between successive
angles at time i-1 and i") 
box()
symbols(0,0,circle=1, inches=FALSE, lwd=2, add=TRUE)
abline(h=0, v=0)
x <- c( cos(pi/3), cos(pi/2 + pi/4))
y <- c( sin(pi/3), sin(pi/2 + pi/4))
arrows(c(0,0), c(0,0), x, y)
lines(x,y, lwd=2, col="red")
text(0, 0.9, expression(f^2 == 2*sum((1 - cos(alpha[i]-alpha[i-1])),
i==1, n-1)), col="red")
foo <- function(t, alpha)
{
  xa <- sapply(seq(0, alpha, length=20), function(x) t*cos(x))
  ya <- sapply(seq(0, alpha, length=20), function(x) t*sin(x))
  lines(xa, ya)
}
foo(0.3, pi/3)
foo(0.1, pi/2 + pi/4)
foo(0.11, pi/2 + pi/4)
text(0.34,0.18,expression(alpha[i]), cex=1.5)
text(0.15,0.11,expression(alpha[i-1]), cex=1.5)
par(opar)
}
\keyword{htest}
