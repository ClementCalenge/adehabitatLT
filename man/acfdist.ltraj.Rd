\name{acfdist.ltraj}
\alias{acfdist.ltraj}
\alias{acfang.ltraj}
\title{
Compute correlogram for angular and linear descriptors of a movement path
}
\description{
  The functions \code{acfdist.ltraj} and \code{acfang.ltraj} compute
  (and by default plot) a correlogram-like function .
}
\usage{
acfdist.ltraj(x, which = c("dist", "dx", "dy"), nrep = 999, lag = 1,
              plot = TRUE, xlab = "Lag", ylab = "autocorrelation")

acfang.ltraj(x, which = c("absolute", "relative"), nrep = 999, lag = 1,
             plot = TRUE, xlab = "Lag", ylab = "autocorrelation")
}

\arguments{
  \item{x}{an object of the class \code{ltraj}}
  \item{which}{ to select on which parameter the autocorrelation should
  be computed (see details).}
  \item{nrep}{ the number of repetitions used to test the significance
  of autocorrelation for each lag value.}
  \item{lag}{ maximum lag at which to calculate the autocorrelation.
  Default is 1.}
\item{plot}{ logical. If 'TRUE' (the default) the autocorrelation
  is plotted.}
  \item{xlab}{ a title for the x axis}
  \item{ylab}{ a title for the y axis}
}
\details{

  The function \code{acfdist.ltraj} is used to compute a correlogram for
  linear descriptors and \code{acfang.ltraj} for angular descriptors
  (see \code{as.ltraj} for a description of these descriptors).

  Statistics used are defined in Dray et al. (in press). They are based
  on squared differences between successive values. For angular
  descriptors, the statistic is based on the chord distance.

  In the case of missing data, the computation of the correlograms is
  restricted to the pairs of successive observed data and only observed
  data are permuted (i.e. the structure of the missing data is kept
  constant under permutation).

  The grey area represents a 95 \% interval obtained after permutation
  of the data. If the observed data is outside this region, it is
  considered as significant and represetend by a black symbol. Note that
  no multiple-comparison adjustement is performed.
  
}
\value{
  A list of matrices. Each matrix corresponds to a 'burst'. The matrix
  contains for each lag value (column), the values of autocorrelation
  (observed, and the 2.5 \%, 50 \% and 97.5 \% quantiles of for the set
  of \code{nrep} permutations of  values).
}
\references{
  Dray, S., Royer-Carenzi, M. and Calenge, C. The exploratory analysis
  of autocorrelation in animal movement studies. \emph{Ecological
    Research}, in press.\cr
  
  Calenge, C., Dray, S. and Royer-Carenzi, M. (2009) The concept of animals
  trajectories from a data analysis perspective.  \emph{Ecological
    Informatics}, \bold{4},34--41.
}
\author{
  Stephane Dray \email{dray@biomserv.univ-lyon1.fr}
}
\seealso{ \code{\link{as.ltraj}} for additional information on the class
  \code{ltraj},  \code{\link{wawotest}} for a simple test of the
  autocorrelation of the descriptive parameters on the trajectory. }
\examples{
\dontrun{
data(puechcirc)
puechcirc
acfang.ltraj(puechcirc, lag=5)
acfdist.ltraj(puechcirc, lag=5)
}
}
\keyword{ htest }

