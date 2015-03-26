\name{lavielle}
\alias{lavielle}
\alias{lavielle.default}
\alias{lavielle.ltraj}
\alias{print.lavielle}
\alias{chooseseg}
\alias{findpath}
\title{Segmentation of a time series using the method of Lavielle (1999,
  2005)}
\description{
  These functions allow to perform a non-parametric segmentation of a
  time series using the penalized contrast method of Lavielle (1999,
  2005). The function \code{lavielle} computes the contrast matrix
  (i.e., the matrix used to segment the series) either from a series of
  observations or from an animal trajectory.  The function
  \code{chooseseg} can be used to estimate the number of segments
  building up the trajectory.  The function \code{findpath} can be used
  to find the limits of the segments (see Details). 
}
\usage{
lavielle(x, \dots)

\method{lavielle}{default}(x, Lmin, Kmax, ld = 1,
                           type = c("mean", "var", "meanvar"), \dots) 

\method{lavielle}{ltraj}(x, Lmin, Kmax, ld = 1, which = "dist",
                         type = c("mean", "var", "meanvar"), \dots) 

\method{print}{lavielle}(x, \dots)

chooseseg(lav, S = 0.75, output = c("full","opt"),
          draw = TRUE)

findpath(lav, K, plotit = TRUE)

}
\arguments{
  \item{x}{for \code{lavielle.default}, a vector containing the
    successive observations building up the series.  For
    \code{lavielle.ltraj}, an object of class \code{ltraj}.}
  \item{Lmin}{an integer value indicating the minimum number of
    observations in each segment.  Should be a multiple of \code{ld}.}
  \item{Kmax}{an integer value indicating the maximum number of segments
    expected in the series}
  \item{ld}{an integer value indicating the resolution for the
    calculation of the contrast function.  The contrast function will be
    evaluated for segments containing the observations \code{c(1:ld)},
    \code{c(1:(2*ld))}, \code{c(1:(3*ld))}, and all segments will
    necessarily contain a multiple of \code{ld} observations.  Note that
    \code{ld} should be set to values greater than 1 if memory problem
    occur}
  \item{type}{the type of contrast function to be used to segment the
    series (see Details)}
  \item{which}{a character string giving any syntactically correct R
    expression implying the descriptive elements in \code{x} or the
    variables in the optional attribute \code{infolocs}.}
  \item{lav}{an object of class \code{"lavielle"}}
  \item{S}{a value indicating the threshold in the second derivative of
    the contrast function}
  \item{output}{type of output expected (see the section value)}
  \item{draw}{a logical value indicating whether the decrease in the
    contrast function should be plotted}
  \item{K}{The number of segments}
  \item{plotit}{a logical value indicating whether the segmentation should
    be plotted}
  \item{\dots}{additional arguments to be passed from or to other
    functions}
}
\details{

  The method of Lavielle (1999, 2005) per se finds the best segmentation
  of a time series, given that it is built by \code{K} segments.  It
  searches the segmentation for which a contrast function (measuring the
  contrast between the actual series and the segmented series) is
  minimized.  Different contrast functions are available measuring
  different aspects of the variation of the series from one segment to
  the next: when \code{type = "mean"}, we suppose that only
  the mean of the segments varies between segments; when \code{type =
    "var"}, we suppose that only the variance of the segments varies
  between segments; when \code{type = "meanvar"}, we suppose that both
  the mean and the variance varies between segments.  It is required to
  specify a value for the minimum number of observations \code{Lmin} in
  a segment, as well as the maximum number of segments \code{Kmax} in
  the series.

  There are several approaches to estimate the best number of segments
  \code{K} to partition the time series.  One possible approach is
  the graphical examination of the decrease of the contrast function
  with the number of segments.  In theory, there should be a clear
  "break" in the decrease of this function after the optimal value of
  \code{K}.  Lavielle (2005) suggested an alternative way to estimate
  automatically the optimal number of segments, also relying on the
  presence of a "break" in the decrease of the contrast function.  He
  proposed to choose the last value of \code{K} for which the second
  derivative of a standardized constrast function is greater than a
  threshold \code{S} (see Lavielle, 2005 for details).  Based on
  numerical experiments, he proposed to choose the value \code{S =
  0.75}. Note, however, that for short time series (i.e. less than 500
  observations) some simulations indicated that this value may not be
  optimal and may depend on the value of \code{Kmax}, so that the
  graphical method is maybe more appropriate.
}
\note{
  The contrast matrix is a matrix of size \code{n*n} (with \code{n} the
  number of observations in the series).  If \code{n} is large, memory
  problems may occur.  In this case, setting \code{ld} to a value
  greater than one will allow to reduce the size of this matrix (i.e. it
  will be of size \code{k*k}, where \code{k = floor(n/ld)}).  However,
  this will also reduce the resolution of the segmentation, so that the
  segment limits will be less precisely estimated.
}
\value{
  The function \code{lavielle.default} returns a list of class
  \code{lavielle}, with an attribute \code{"typeseg"} set to
  \code{"default"}. This list contains the following elements:
  \item{contmat}{The contrast matrix}
  \item{sumcont}{The optimal contrast}
  \item{matpath}{The matrix of the paths from the first to the last
    observation}
  \item{Kmax}{The maximum number of segments}
  \item{Lmin}{The minimum number of observations in a segment}
  \item{ld}{the value of the resolution \code{ld}}
  \item{series}{The time series}

  The function \code{lavielle.ltraj} also returns a list of class
  \code{lavielle}, with an attribute \code{"typeseg"} set to
  \code{"ltraj"}. 
  
  The function \code{chooseseg} returns the optimal number of segments
  when \code{output = "opt"}, and a dataframe containing the value of
  the contrast function \code{Jk} and of the second derivative \code{D}
  of the standardized contrast function for each possible value of
  \code{K}, if \code{output = "full"}.

  The function \code{findpath} return a list containing vectors giving
  the index of the first and last observations in each segment, when the
  object of class \code{"lavielle"} passed as argument is characterized
  by an attribute \code{"typeseg"} set to \code{"default"}.  When the
  attribute \code{"typeseg"} is set to \code{"ltraj"}, this function
  returns an object of class ltraj where each burst correspond to a
  segment. 
}
\references{
  Lavielle, M. (1999) Detection of multiple changes in a sequence of
  dependent variables. \emph{Stochastic Processes and their
  Applications}, \bold{83}: 79--102.
  
  Lavielle, M. (2005) Using penalized contrasts for the change-point
  problem. Report number 5339, Institut national de recherche en
  informatique et en automatique.
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}.  The code
  is a C translation based on the Matlab code of M. Lavielle}
\examples{

#################################################
##
## A simulated series

set.seed(129)
seri <- c(rnorm(100), rnorm(100, mean=2),
          rnorm(100), rnorm(100, mean=-3),
          rnorm(100), rnorm(100, mean=2))
plot(seri, ty="l", xlab="time", ylab="Series")

## Segmentation:
(l <- lavielle(seri, Lmin=10, Kmax=20))

## choose the number of segments
chooseseg(l)

## There is a clear break in the
## decrease of the contrast function after K = 6
## Moreover, Jk(6) >> 0.75 and Jk(7) << 0.75
## We choose 6 segments:
fp <- findpath(l, 6)
fp

## This list gives the limits of the segments
## for example, to get the first segment:
seg <- 1
firstseg <- seri[fp[[seg]][1]:fp[[seg]][2]]

####################################################
##
## Now, changes of variance

## A simulated series
set.seed(129)
seri <- c(rnorm(100), rnorm(100, sd=2),
          rnorm(100), rnorm(100, sd=3),
          rnorm(100), rnorm(100, sd=2))
plot(seri, ty="l", xlab="time", ylab="Series")

## Segmentation:
(l <- lavielle(seri, Lmin=10, Kmax=20, type="var"))

## choose the number of segments
chooseseg(l)

## There is a clear break in the
## decrease of the contrast function after K = 6
## Moreover, Jk(6) >> 0.75 and Jk(7) << 0.75
## We choose 6 segments:
fp <- findpath(l, 6)
fp

## This list gives the limits of the segments
## for example, to get the first segment:
seg <- 1
firstseg <- seri[fp[[seg]][1]:fp[[seg]][1]]

#################################################
##
## Example of segmentation of a trajectory

## Show the trajectory
data(porpoise)
gus <- porpoise[1]
plot(gus)

## Show the changes in the distance between
## successive relocations with the time
plotltr(gus, "dist")

## Segmentation of the trajectory based on these distances
lav <- lavielle(gus, Lmin=2, Kmax=20)

## Choose the number of segments
chooseseg(lav)
## 4 segments seem a good choice

## Show the partition
kk <- findpath(lav, 4)
plot(kk)

}
\keyword{cluster}
\keyword{chron}
