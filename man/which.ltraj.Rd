\name{which.ltraj}
\alias{which.ltraj}
\title{Identify Relocations Fullfilling a Condition in an Object of
  Class "ltraj"}
\description{
  This function identifies the relocations fullfilling a condition in an
  object of class \code{ltraj}.
}
\usage{
which.ltraj(ltraj, expr)
}
\arguments{
  \item{ltraj}{an object of class \code{ltraj}}
  \item{expr}{a character string giving any syntactically correct R
    logical expression implying the descriptive elements in
    \code{ltraj} (or the name of a variable in the optional attribute
    \code{infolocs})}
}
\value{
  A data frame giving the ID, the Bursts and the relocations for which
  the condition described by \code{expr} is verified.
}
\author{Clement Calenge \email{clement.calenge@ofb.gouv.fr}}
\seealso{\code{\link{ltraj}} for additional information about objects
  of class \code{ltraj}}
\examples{
data(puechcirc)
puechcirc

## Identifies the relocations for which time lag is
## upper than one hour
which.ltraj(puechcirc, "dt>3600")
puechcirc[burst="CH930824"][[1]][27:28,]


## Identifies the speed between successive
## relocations upper than 0.8 meters/second
which.ltraj(puechcirc, "dist/dt > 0.8")

## This is the case for example for the
## relocations #28, 58, 59 and 60 of "CH930824"
puechcirc[burst="CH930824"][[1]][c(28,58,59,60),]

}
\keyword{programming}
