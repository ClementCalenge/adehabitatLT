\name{burst}
\alias{burst}
\alias{id}
\alias{burst<-}
\alias{id<-}
\alias{removeinfo}
\alias{infolocs}
\alias{infolocs<-}
\title{ID, Bursts and infolocs of an Object of Class ltraj}
\description{
  Functions to get or set the attribute \code{"id"}, \code{"burst"}, and
  \code{"infolocs"} of the components of an object of class
  \code{ltraj}.  
}
\usage{
burst(ltraj)
burst(ltraj) <- value
id(ltraj)
id(ltraj) <- value
infolocs(ltraj, which)
infolocs(ltraj) <- value
removeinfo(ltraj)
}
\arguments{
  \item{ltraj}{an object of class \code{ltraj}}
  \item{value}{for the assignment functions \code{burst} and \code{id},
    a character vector of up to the same length as \code{ltraj}.  For
    \code{infolocs} a list of data frames of the same length of
    \code{ltraj} (with each component having the same number of rows as
    the corresponding element in \code{ltraj}).}
  \item{which}{an optional character vector containing the names of the
    variables in the \code{infolocs} attribute to be returned}
}
\details{
     The functions \code{id}, \code{burst} and \code{infolocs} are
     accessor functions, and \code{id<-} and \code{burst<-} are
     replacement function.  \code{removeinfo} removes the attribute
     \code{infolocs} from the object \code{ltraj} (see the help page of
     \code{as.ltraj}).
}
\value{
  For \code{id} and \code{burst}, a character vector of the same length
  as \code{ltraj}.  For \code{infolocs}, the data frame containing the
  information on the relocations.  \code{removeinfo} returns an object
  of class \code{ltraj}.\cr
  
  For \code{id<-} and \code{burst<-}, the updated object.  (Note that
  the value of \code{burst(x) <- value} is that of the assignment,
  \code{value}, not the return value from the left-hand side.)
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }
\seealso{\code{\link{ltraj}}, \code{\link{names}}}
\examples{
data(puechcirc)
puechcirc

## To see the ID and the burst
id(puechcirc)
burst(puechcirc)

## Change the burst
burst(puechcirc) <- c("glou", "toto", "titi")
puechcirc

burst(puechcirc)[2] <- "new name"
puechcirc

## Change the ID
id(puechcirc)[id(puechcirc)=="CH93"] <- "WILD BOAR"
puechcirc

## example of an object with an attribute "infolocs"
data(capreochiz)
head(capreochiz)

## Create an object of class "ltraj"
cap <- as.ltraj(xy = capreochiz[,c("x","y")], date = capreochiz$date,
                id = "Roe.Deer", typeII = TRUE,
                infolocs = capreochiz[,4:8])
cap
cap2 <- removeinfo(cap)
cap2

infolocs(cap)



}
\keyword{programming}
