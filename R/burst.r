burst <- function(ltraj)
{
    if (!inherits(ltraj, "ltraj"))
        stop("ltraj should be of class \"ltraj\"")
    return(unlist(lapply(ltraj, function(x) attr(x,"burst"))))
}##OK

"burst<-" <- function(ltraj, value)
{
    if (!inherits(ltraj, "ltraj"))
        stop("ltraj should be of class \"ltraj\"")
    p4s <- .checkp4(ltraj)
    res <- lapply(1:length(ltraj), function(i) {
        x <- ltraj[[i]]
        attr(x,"burst") <- value[i]
        return(x)
    })
    class(res) <- c("ltraj","list")
    attr(res, "typeII") <- attr(ltraj, "typeII")
    attr(res, "regular") <- is.regular(res)
    attr(res,"proj4string") <- p4s
    return(res)
}##OK


id <- function(ltraj)
{
    if (!inherits(ltraj, "ltraj"))
        stop("ltraj should be of class \"ltraj\"")
    return(unlist(lapply(ltraj, function(x) attr(x,"id"))))
}##OK

"id<-" <- function(ltraj, value)
{
    if (!inherits(ltraj, "ltraj"))
        stop("ltraj should be of class \"ltraj\"")
    p4s <- .checkp4(ltraj)
    res <- lapply(1:length(ltraj), function(i) {
        x <- ltraj[[i]]
        attr(x,"id") <- value[i]
        return(x)
    })
    class(res) <- c("ltraj","list")
    attr(res, "typeII") <- attr(ltraj, "typeII")
    attr(res, "regular") <- is.regular(res)
    attr(res,"proj4string") <- p4s
    return(res)
}##OK


"infolocs<-" <- function(ltraj, value)
{
    if (!inherits(ltraj, "ltraj"))
        stop("ltraj should be of class ltraj")
    if (length(value)!=length(ltraj))
        stop("the assignment should be a list of the same length as ltraj")
    if (length(value)>1) {
        toto <- apply(do.call(rbind,lapply(value, names)),2,function(x) all(x==x[1]))
        if (!toto) {
            stop("the names of the variables in infolocs differ between bursts")
        }
        toto <- apply(do.call(rbind,lapply(value, ncol)),2,function(x) all(x==x[1]))
        if (!toto) {
            stop("the number of variables in infolocs differ between bursts")
        }
        toto <- lapply(value, function(x) lapply(x, class))
        tutu <- all(sapply(1:length(toto[[1]]),
                           function(i) all(sapply(2:length(toto),
                                                  function(j) all(toto[[j]][[i]]==toto[[1]][[i]])))))
        if (!tutu) {
            stop("The class of the variables in the data.frames differ between bursts")
        }
    }

    for (i in (1:length(ltraj))) {
        df <- value[[i]]
        if (!inherits(df, "data.frame"))
            stop("value should be a list of data.frame")
        if (nrow(df)!=nrow(ltraj[[i]]))
            stop(paste("The burst", i,
                       "does not have the same number of elements in ltraj and value"))
        if (!all(row.names(df)==row.names(ltraj[[i]])))
            stop("The infolocs component should have the same row.names as the ltraj object")
        attr(ltraj[[i]], "infolocs") <- df
    }
    return(ltraj)
}
