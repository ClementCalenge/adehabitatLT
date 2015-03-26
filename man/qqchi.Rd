\name{qqchi}
\alias{qqchi}
\alias{qqchi.default}
\alias{qqchi.ltraj}
\alias{qqnorm.ltraj}
\title{Quantile-Quantile Plots for Trajectories of Class 'ltraj'}
\description{
  The functions allow the examination of the distribution of trajectories
  descriptors (see Details).
}
\usage{

## Chi distribution of the increment length / sqrt(dt)
qqchi(y, ...)

\method{qqchi}{default}(y, df = 2, ylim, main = "Chi Q-Q Plot",
              xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
              plot.it = TRUE, datax = FALSE, \dots) 

\method{qqchi}{ltraj}(y, xlab = "Theoretical Quantiles", 
            ylab = "Sample Quantiles (Distances)", \dots)

## Normal Distribution of dx/sqrt(dt) or dy/sqrt(dt)
\method{qqnorm}{ltraj}(y, which=c("dx","dy"), ...)

}
\arguments{
  \item{y}{a vector containing the data sample for \code{qqchi.default}.
   an object of class \code{ltraj} for other functions.}
  \item{df}{the number of degrees of freedom of the Chi distribution
    (default to 2).}
  \item{xlab, ylab, main}{plot labels.}
  \item{plot.it}{logical. Should the result be plotted?}
  \item{datax}{logical. Should data values be on the x-axis?}
  \item{which}{a character string indicating the component (dx or dy) to
  be examined.}
  \item{ylim, \dots}{graphical parameters.}
}
\details{
  
  Among the numerous statistics that can be used to describe the
  movements of an animal, the length of the increment between two
  successive relocations is very common.  This increment can be
  described by a vector \code{i = c(dx, dy)}.  Under the hypothesis
  of a Brownian motion, dx and dy should be normally distributed with
  mean = 0 and variance = dt (where dt is the time interval between the
  two relocations).  Therefore, \code{dx/sqrt(dt)} and
  \code{dy/sqrt(dt)} should be normally distributed with mean = 0 and
  variance = 1.  The function \code{qqnorm.ltraj} performs a
  quantile-quantile plot of \code{dx/sqrt(dt)} or \code{dy/sqrt(dt)}
  vs. a normal distribution to verify wether the Brownian motion
  assumption is correct.
  
  Furthermore, the quantity \code{(dx^2 + dy^2)/dt} should be
  distributed according to a Chi-squared distribution with two degrees
  of freedom.  Thus, the quantity \code{distance / sqrt(dt)} should be
  distributed according to a Chi distribution with two degrees of
  freedom (where \code{distance} is the distance between the two
  relocations). The function \code{qqchi.ltraj} performs
  quantile-quantile plot of \code{distance/sqrt(dt)} vs. a Chi
  distribution to verify wether the Brownian motion
  assumption is correct.
}
\value{
  for functions dealing with objects of class \code{ltraj}, a list with
  components being themselves lists, with components:
  \item{x}{The x coordinates of the points that were/would be plotted}
  
  \item{y}{The original \code{y} vector, i.e., the corresponding y
    coordinates including 'NA's.}
  
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}
\seealso{\code{\link{chi}}, \code{\link{qqplot}}, \code{\link{ltraj}}.}
\examples{

## Example with an Arithmetic Brownian Process
toto <- simm.mba(1:500, sig = diag(c(5, 5)))
qqnorm(toto, "dx")
qqnorm(toto, "dy")
qqchi(toto)

## Example of wild boar
data(puechcirc)
qqnorm(puechcirc, "dx")
qqnorm(puechcirc, "dy")
qqchi(puechcirc)


}
\keyword{distribution}

