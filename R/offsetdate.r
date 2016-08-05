offsetdate <- function(ltraj, offset, units=c("sec", "min", "hour", "day"))
{
    if (!inherits(ltraj,"ltraj"))
        stop("ltraj should be of class \"ltraj\"")
    if (!attr(ltraj,"typeII"))
        stop("ltraj should be of type II (time recorded)")
    p4s <- .checkp4(ltraj)
    units <- match.arg(units)
    offset <- .convtime(offset, units)

    res <- lapply(ltraj, function(x) {
        date <- as.numeric(x$date)
        date <- date - offset
        class(date) <- c("POSIXct","POSIXt")
        attr(date, "tzone") <- attr(x$date, "tzone")
        x$date <- date
        return(x)
    })
    class(res) <- c("ltraj", "list")
    attr(res,"typeII") <- TRUE
    attr(res,"regular") <- is.regular(res)
    attr(res,"proj4string") <- p4s
    return(res)
}
