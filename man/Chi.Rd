\name{Chi}
\alias{chi}
\alias{dchi}
\alias{pchi}
\alias{qchi}
\alias{rchi}
\title{The Chi Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the chi distribution with \code{df} degrees of
  freedom.
}
\usage{
dchi(x, df = 2)
pchi(q, df = 2, lower.tail = TRUE, \dots)
qchi(p, df = 2, lower.tail = TRUE)
rchi(n, df = 2)
}

\arguments{
  \item{x,q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
  \item{df}{degrees of freedom (non-negative, but can be non-integer).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P[X <= x],
    otherwise, P[X > x].}
  \item{\dots}{additional arguments to be passed to the function
    \code{integrate}.}
}
\details{
  The chi distribution with \code{df} = n > 0 degrees of freedom
  has density
  
  \deqn{f_n (x) = 2^{1-n/2} x^{n-1} e^{\frac{-(x^2)}{2}} / \Gamma (n/2)}{
    f_n (x) = 2^(1-n/2) x^(n-1) e^(-(x^2)/2) / Gamma(n/2)}
  
  for x > 0.  This distribution is used to describe the square root of a
  variable distributed according to a chi-square distribution.
}
\value{
  \code{dchi} gives the density, \code{pchi} gives the distribution
  function, \code{qchi} gives the quantile function, and \code{rchi}
  generates random deviates.
}
\references{
  Evans, M., Hastings, N. and Peacock, B. (2000) Statistical
  Distributions, 3rd ed. Wiley, New York.
}
\author{Clement Calenge \email{clement.calenge@ofb.gouv.fr}}
\seealso{\code{\link{Chisquare}}}
\examples{

opar <- par(mfrow = c(2,2))

hist(rchi(100), ncla = 20, main="The Chi distribution")

plot(tutu <- seq(0, 5, length=20), dchi(tutu, df = 2), xlab = "x",
     ylab = "probability density", type = "l")

plot(tutu, pchi(tutu), xlab = "x", ylab = "Repartition function",
     type = "l")

par(opar)

}
\keyword{distribution}
