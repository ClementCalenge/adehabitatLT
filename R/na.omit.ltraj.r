
na.omit.ltraj <- function(object, ...)
{
    if (!inherits(object, "ltraj"))
        stop("ltraj should be of class ltraj")
    p4s <- .checkp4(object)
    info <- infolocs(object)
    for (i in 1:length(object)) {
        x <- object[[i]]
        if (!is.null(info))
            info[[i]] <- info[[i]][(!is.na(x[,1]))&(!is.na(x[,2])),,drop=FALSE]
        object[[i]] <- x[(!is.na(x[,1]))&(!is.na(x[,2])),]
    }
    if (!is.null(info))
        infolocs(object) <- info
    res <- rec(object)
    attr(res, "proj4string") <- p4s
    return(res)
}
