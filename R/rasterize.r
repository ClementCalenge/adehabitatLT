rasterize.ltraj <- function(ltr, map)
{
    if (!inherits(ltr, "ltraj"))
        stop("ltr should be of class ltraj")
    if (!inherits(map, "SpatialPixels"))
        stop("map should inherit the class SpatialPixels")
#    ltr <- na.omit(ltr)
    pa <- gridparameters(map)
    pfs <- proj4string(map)

    xll <- pa[1,1]
    yll <- pa[2,1]
    cs <- pa[1,2]

    lapply(1:length(ltr), function(i) {
        x <- ltr[[i]][,c("x","y")]
        res <- .Call("RasterPas", x, xll, yll, cs, as.integer(0), PACKAGE="adehabitatLT")
        res <- as.data.frame(res)
        names(res) <- c("x","y","step")
        coordinates(res) <- c("x","y")
        if (!is.na(pfs))
            proj4string(res) <- CRS(pfs)
        return(res)
    })
}
