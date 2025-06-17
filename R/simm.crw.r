"simm.crw" <- function(date=1:100, h = 1, r = 0,
                       x0=c(0,0), id="A1", burst=id,
                       typeII=TRUE, proj4string=CRS())
{
    if (typeII) {
        if (!inherits(date, "POSIXct")) {
            class(date) <- c("POSIXct", "POSIXt")
            attr(date, "tzone") <- ""
        }
    }
    if (r < 0 | r > 1)
        stop("rho must be between 0 and 1")

    n <- length(date)
    dt <- c(diff(unclass(date)))
    if (all(dt-dt[1]>1e-7))
        stop("the time lag between relocations should be constant")

    if (r==0) {
        ang <- runif(n-2,0,2*pi)
    } else if (r==1) {
        ang <- rep(0, n-2)
    } else {
        sd <- sqrt(-2*log(r))
        ang<-rnorm(n-2,0,sd)%%(2*pi)
    }


    if (h>0) {
        v=sqrt(dt)*rchi(n-1) * h
    } else {
        v=-h
    }
    ang=cumsum(c(runif(1,0,2*pi),ang))
    si=c(x0[2], x0[2]+cumsum(v*sin(ang)))
    co=c(x0[1], x0[1]+cumsum(v*cos(ang)))
    res <- as.ltraj(data.frame(co,si),date, id, burst,
                    typeII=typeII, proj4string=proj4string)
    return(res)
}
