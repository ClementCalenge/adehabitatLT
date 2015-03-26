\name{as.ltraj}
\alias{as.ltraj}
\alias{ltraj}
\alias{print.ltraj}
\alias{summary.ltraj}
\alias{rec}
\title{Working with Trajectories in 2D Space: the Class ltraj}
\description{

  The class \code{ltraj} is intended to store trajectories of
  animals.  Trajectories of type II correspond to trajectories for which the
  time is available for each relocation (mainly GPS and
  radio-tracking).  Trajectories of type I correspond to trajectories for which
  the time has not been recorded (e.g. sampling of tracks in the snow).
  
  \code{as.ltraj} creates an object of this class.
  
  \code{summary.ltraj} returns the number of relocations (and missing
  values) for each "burst" of relocations and each animal.
  
  \code{rec} recalculates the descriptive parameters of an object of
  class ltraj (e.g. after a modification of the contents of this object,
  see examples)
}
\usage{
as.ltraj(xy, date, id, burst = id, typeII = TRUE,
         slsp = c("remove", "missing"),
         infolocs = data.frame(pkey = paste(id, date, sep="."),
                               row.names=row.names(xy)))
\method{print}{ltraj}(x, \dots)
\method{summary}{ltraj}(object, \dots)
rec(x, slsp = c("remove", "missing"))
}
\arguments{
  \item{x, object}{an object of class \code{ltraj}}
  \item{xy}{a data.frame containing the x and y coordinates of the
    relocations} 
  \item{date}{for trajectories of type II, a vector of class \code{POSIXct}
    giving the date for each relocation. For trajectories of type I, this
    argument is not taken into account.}
  \item{id}{either a character string indicating the identity of the
    animal or a factor with length equal to \code{nrow(xy)} }
  \item{burst}{either a character string indicating the identity of the
    burst of relocations or a factor with length equal to
    \code{nrow(xy)}}
  \item{typeII}{logical. \code{TRUE} indicates a trajectory of type II
    (time recorded, e.g. radio-tracking), whereas \code{FALSE} indicates
    a trajectory of type I (time not recorded, e.g. sampling of tracks in
    the snow)}
  \item{slsp}{a character string used for the computation of the
    turning angles (see details)}
  \item{infolocs}{if not \code{NULL}, a data frame containing
    additionnal information on the relocations (e.g., precision).  By
    default, a primary key is generated.}
  \item{\dots}{For other functions, arguments to be passed to the generic
    functions \code{summary} and \code{print}}
}
\details{ 

  Objects of class \code{ltraj} allow the analysis of animal
  movements.  They contain the descriptive parameters of the moves
  generally used in such studies (coordinates of the relocations, date,
  time lag, relative and absolute angles, length of moves, increases
  in the x and y direction, and dispersion R2n, see below), as well as
  optionally metadata on the relocations (precision, etc.).
  
  The computation of turning angles may be problematic when successive
  relocations are located at the same place.  In such cases, at least
  one missing value is returned.  For example, let r1, r2, r3 and r4 be
  4 successive relocations of a given animal (with coordinates (x1,y1),
  (x2,y2), etc.).  The turning angle in r2 is computed between the moves
  r1-r2 and r2-r3.  If r2 = r3, then a missing value is returned for the
  turning angle at relocation r2.  The argument \code{slsp} controls the
  value returned for relocation r3 in such cases.  If \code{slsp ==
    "missing"}, a missing value is returned also for the relocation r3.
  If \code{slsp == "remove"}, the turning angle computed in r3 is the
  angle between the moves r1-r2 and r3-r4.

  For a given individual, trajectories are often sampled as "bursts"
  of relocations.  For example, when an animal is monitored using
  radio-tracking, the data may consist of several circuits of activity
  (two successive relocations on one circuit are often highly
  autocorrelated, but the data from two circuits may be sampled at long
  intervals in time).  These bursts are indicated by the attribute
  \code{burst}. Note that the bursts should be unique: do not use the
  same burst id for bursts collected on different animals.

  Two types of trajectories can be stored in objects of class \code{ltraj}:
  trajectories of type I correspond to trajectories where the time of
  relocations is not recorded.  It may be because it could not be noted
  at the time of sampling (e.g. sampling of animals' tracks in the snow)
  or because the analyst decided that he did not want to take it into
  account, i.e. to study only its geometrical properties. In this case,
  the variable \code{date} in each burst of the object contains a vector
  of integer giving the order of the relocations in the trajectory (i.e. 1,
  2, 3, ...). Trajectories of type II correspond to trajectories for which the
  time is available for each relocation. It is stored as a vector of
  class \code{POSIXct} in the column \code{date} of each burst of
  relocations. The type of trajectory should be defined when the object of
  class \code{ltraj} is defined, with the argument \code{typeII}.
  
  Concerning trajectories of type II, in theory, it is expected that the
  time lag between two relocations is constant in all the bursts and all
  the ids of one object of class 
  \code{ltraj} (don't mix animals located every 10 minutes and animals
  located every day in the same object).  Indeed, some of the
  descriptive parameters of the trajectory do not have any sense when
  the time lag varies.  For example, the  distribution of relative
  angles (angles between successive moves) depends on a given time
  scale; the angle between two during 10-min moves of a whitestork
  does not have the same biological meaning as the angle between two
  1-day move.  If the time lag varies, the underlying process varies
  too.  For this reason, most functions of adehabitatLT have been
  developed for "regular" trajectories, i.e. trajectories with a
  constant time lag (see \code{help(sett0)}).  Furthermore, several
  functions are intended to help the user to transform an object of
  class \code{ltraj} into a regular object (see for example
  \code{help(sett0)}, and particularly the examples to see how regular
  trajectories can be obtained from GPS data).
  
  Nevertheless, the class \code{ltraj} allows for variable time lag,
  which often occur with some modes of data collection (e.g. with Argos
  collars). But *we stress that their analysis is still an open
  question!!*
  
  Finally, the class \code{ltraj} deals with missing values in the
  trajectories.  Missing values are frequent in the trajectories of
  animals collected using telemetry: for example, GPS collar may not
  receive the signal of the satellite at the time of relocation. Most
  functions dealing with the class \code{ltraj} have a specified
  behavior in case of missing values.

  It is recommended to store the missing values in the data *before*
  the creation of the object of class \code{ltraj}.  For example, the
  GPS data imported within R contain missing values.  It is recommended
  to *not remove* these missing values before the creation of the
  object!!!  These missing values may present patterns (e.g. failure to
  locate the animal at certain time of the day or in certain habitat
  types), and *the analysis of these missing values should be part of the
  analysis of the trajectory* (e.g. see \code{help(runsNAltraj)} and
  \code{help(plotNAltraj)}.

  However, sometimes, the data come without any information concerning
  the location of these missing values.  If the trajectory is
  approximately regular (i.e. approximately constant time lag), it is
  possible to determine where these missing values should occur in the
  object of class \code{ltraj}.  This is the role of the function
  \code{setNA}.

  One word now about the attribute \code{infolocs} of the object of
  class \code{ltraj}.  This attribute is intended to store metadata
  concerning the relocations building the trajectory (the precision
  of the relocations, the value of environmental variables at this
  place, etc.).  There are constraints on the structure of this
  attribute. Although any variable can be stored in this attribute, it
  is required that: (i) all the relocations take a value (or a missing
  value) for all variables, (ii) all the variables are measured (or
  correspond to missing values) for all bursts and ids.  This means for
  example that the function \code{c.ltraj} cannot be used to combine an
  object where only variable "A" is stored as metadata and an object
  where only variable "B" is stored as metadata.  The function
  \code{removeinfo} can be used to remove this attribute.  Note also
  that the data.frames in the list infolocs should have the same
  row.names as the corresponding elements in the object of class
  \code{"ltraj"}. 
}
\value{
  \code{summary.ltraj} returns a data frame.\cr
  All other functions return objects of class \code{ltraj}.  An object
  of class \code{ltraj} is a list with one component per burst of
  relocations.  Each component is a data frame with two attributes and
  one optionnal attribute:
  the attribute \code{"id"} indicates the identity of the animal, and
  the attribute \code{"burst"} indicates the identity of the
  burst.  An optional attribute \code{"infolocs"} contains any
  additional information desired by the user (precision, etc.).
  Each main data frame stores the following columns:
\item{x}{the x coordinate for each relocation}
\item{y}{the y coordinate for each relocation}
\item{date}{the date for each relocation (type II) or a vector of
  integer giving the order of the relocations in the trajectory.}
\item{dx}{the increase of the move in the x direction.  At least two
  successive relocations are needed to compute \code{dx}.  Missing
  values are returned otherwise.}
\item{dy}{the increase of the move in the y direction.  At least two
  successive relocations are needed to compute \code{dy}.  Missing
  values are returned otherwise.}
\item{dist}{the length of each move.  At least two
  successive relocations are needed to compute \code{dist}.  Missing
  values are returned otherwise.}
\item{dt}{the time interval between successive relocations}
\item{R2n}{the squared net displacement between the current relocation
  and the first relocation of the trajectory}
\item{abs.angle}{the angle between each move and the x axis.  At least two
  successive relocations are needed to compute \code{abs.angle}.  Missing
  values are returned otherwise.}
\item{rel.angle}{the turning angles between successive moves.  At least
  three successive relocations are needed to compute \code{rel.angle}.
  Missing values are returned otherwise.}
}
\references{
  Calenge, C., Dray, S. and Royer, M. (2009) The concept of animals'
  trajectories from a data analysis perspective. \emph{Ecological
  Informatics}, \bold{4}: 34--41.
}
\author{
  Clement Calenge \email{clement.calenge@oncfs.gouv.fr}\cr
  Stephane Dray \email{dray@biomserv.univ-lyon1.fr}
}
\seealso{\code{\link{is.regular}} and \code{\link{sett0}} for additional
  information on "regular" trajectories.  \code{\link{setNA}} and
  \code{\link{runsNAltraj}} for additional information on missing values
  in trajectories.  \code{\link{c.ltraj}} to combine several objects of
  class \code{ltraj}, \code{\link{Extract.ltraj}} to extract or replace
  bursts of relocations, \code{\link{plot.ltraj}} and
  \code{\link{trajdyn}} for graphical displays, \code{\link{gdltraj}} to
  specify a time period. }
\examples{

data(puechabonsp)
locs <- puechabonsp$relocs
head(locs)
xy <- coordinates(locs)
df <- as.data.frame(locs)
id <- df[,1]

######################################################
##
## Example of a trajectory of type I (time not recorded)

(litrI <- as.ltraj(xy, id = id, typeII=FALSE))
plot(litrI)

## The components of the object of class "ltraj"
head(litrI[[1]])


######################################################
##
## Example of a trajectory of type II (time recorded)


### Conversion of the date to the format POSIX
da <- as.character(df$Date)
da <- as.POSIXct(strptime(as.character(df$Date),"\%y\%m\%d"))


### Creation of an object of class "ltraj", with for 
### example the first animal
(tr1 <- as.ltraj(xy[id=="Brock",],
                 date = da[id=="Brock"],
                 id="Brock"))

## The components of the object of class "ltraj"
head(tr1[[1]])

## With all animals
(litr <- as.ltraj(xy, da, id = id))

## Change something manually in the first burst:
head(litr[[1]])
litr[[1]][3,"x"] <- 700000

## Recompute the trajectory
litr <- rec(litr)
## Note that descriptive statistics have changed (e.g. dx)
head(litr[[1]])



######################################################
##
## Example of a trajectory of type II (time recorded)
## with an infolocs attribute:

data(capreochiz)
head(capreochiz)

## Create an object of class "ltraj"
cap <- as.ltraj(xy = capreochiz[,c("x","y")], date = capreochiz$date,
                id = "Roe.Deer", typeII = TRUE,
                infolocs = capreochiz[,4:8])
cap


}
\keyword{spatial}
