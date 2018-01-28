gdltraj <- function(x, min, max,
                    type=c("POSIXct","sec","min","hour",
                    "mday","mon","year","wday","yday"))
{
    ## Verifications
    if (!inherits(x, "ltraj"))
        stop("x should be of class \"ltraj\"")
    if (!attr(x, "typeII"))
        stop("x should be of type II (time recorded)")
    type <- match.arg(type)
    p4s <- .checkp4(x)
    ## gets the traj within the boundaries
    if (type=="POSIXct") {
        tz1 <- .checktz(x)
        tz2 <- .ctzda(min)
        tz3 <- .ctzda(max)
        tz <- c(tz2,tz3)
        if (any(tz!=tz1))
            stop("non consistent time zones")
        x <- lapply(x, function(y) {
            infol <- attr(y, "infolocs")
            if (!is.null(infol))
                infol <- infol[(y$date>min)&(y$date<max),,drop=FALSE]
            y <- y[(y$date>min)&(y$date<max),]
            if (!is.null(infol))
                attr(y, "infolocs") <- infol
            return(y)
        })
    } else {
        x <- lapply(x, function(y) {
            da <- unclass(as.POSIXlt(y$date))[[type]]
            infol <- attr(y, "infolocs")
            if (!is.null(infol))
                infol <- infol[(da>=min)&(da<max),,drop=FALSE]
            y <- y[(da>=min)&(da<max),]
            if (!is.null(infol))
                attr(y, "infolocs") <- infol
            return(y)
        })
    }
    if (all(sapply(x,nrow)==0))
        stop("No relocations within the specified interval")
    x[sapply(x, nrow)==0]<-NULL

    ## Output
    class(x) <- c("ltraj", "list")
    attr(x, "typeII") <-  TRUE
    attr(x, "regular") <-  is.regular(x)
    x <- rec(x)
    attr(x,"proj4string") <- p4s
    return(x)
  }
