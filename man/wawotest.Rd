\name{wawotest}
\alias{wawotest}
\alias{wawotest.default}
\alias{wawotest.ltraj}
\title{Wald-Wolfowitz Test of Randomness}
\description{
  The function \code{wawotest.default} performs a Wald Wolfowitz test of
  the random distribution of the values in a vector.  The function
  \code{wawotest.ltraj} performs this tests for the descriptive
  parameters \code{dx}, \code{dy} and \code{dist} in an object of class
  \code{ltraj}.  The function \code{wawotest} is generic.
}
\usage{
wawotest(x, \dots)
\method{wawotest}{default}(x, alter = c("greater", "less"), \dots)
\method{wawotest}{ltraj}(x, \dots)
}
\arguments{
  \item{x}{for \code{wawotest.default}, a vector containing the
    successive observations building the series.
    For \code{wawotest.ltraj}, an object of class \code{ltraj}.}
  \item{alter}{a character string specifying the alternative hypothesis,
    must be one of "greater" (default), "less" or "two-sided"}
  \item{\dots}{additional arguments to be passed to other functions}
}
\details{
  The statistic of the test is equal to A = sum(y(i) y(i+1)),
  with y(N+1) = y(1).  Under the hypothesis of a random distribution of
  the values in the vector, this statistic is normally distributed, with
  theoretical means and variances given in Wald & Wolfowitz (1943).
}
\value{
  \code{wawotest.default} returns a vector containing the value of the
  statistic (\code{a}), its esperance (\code{ea}), its variance
  (\code{va}), the normed statistic (\code{za}) and the
  P-value. \code{wawotest.ltraj} returns a table giving these values for
  the descriptive parameters of the trajectory.
}
\references{
Wald, A. & Wolfowitz, J. (1943) An exact test for randomness in the
non-parametric case based on serial
correlation. \emph{Ann. Math. Statist.}, \bold{14}, 378--388.
}
\author{  Stephane Dray \email{dray@biomserv.univ-lyon1.fr} }
\seealso{ \code{\link{indmove}} and \code{\link{runsNAltraj}} for other
  tests of independence to be used with objects of class "ltraj" }
\examples{

data(puechcirc)
puechcirc
wawotest(puechcirc)

}
\keyword{htest}

