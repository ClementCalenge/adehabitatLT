\name{Extract.ltraj}
\alias{Extract.ltraj}
\alias{[.ltraj}
\alias{[<-.ltraj}
\title{Extract or Replace Parts of an Object of Class ltraj}
\description{
  Extract or replace subsets of objects of class \code{ltraj}.
}
\usage{
  \method{[}{ltraj}(x, i, id, burst)
  \method{[}{ltraj}(x, i, id, burst) <- value
}
\arguments{
\item{x}{an object of class \code{ltraj} }
\item{i}{numeric.  The elements to extract or replace}
\item{id}{a character vector indicating the identity of the animals to
  extract or replace}
\item{burst}{a character vector indicating the identity of the bursts of
  relocations to extract or replace}
\item{value}{an object of class \code{ltraj}}
}
\details{
  Objects of class \code{ltraj} contain several bursts of
  relocations.  This function subsets or replaces these bursts, based
  on their indices or on the attributes \code{id} *or* \code{burst}.\cr

  When replacement is done, it is required that \code{value} and
  \code{x} have the same variables in attribute \code{infolocs} (i.e.,
  both contain the same variables or both do not contain any variable,
  see the help page of \code{as.ltraj})
}
\value{
  An object of class \code{ltraj}.
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }
\seealso{ \code{\link{ltraj}}, \code{\link{gdltraj}} }
\examples{
  data(puechcirc)
  puechcirc

  ## Extract the second and third bursts
  (toto <- puechcirc[2:3])

  ## Extracts all bursts collected on the animal JE
  puechcirc[id = "JE93"]

  
  ## Replace one burst
  toto[2] <- puechcirc[1]
  toto
}
\keyword{spatial}
