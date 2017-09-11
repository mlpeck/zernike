## Zernike coefficients of Wavefront error due to testing a conic at center of curvature
## this should be more accurate for numerical nulling than twice the height difference as in zconic

sconic <- function(D, rc, b = -1, lambda = 632.8, nmax = 6) {
    require(zernike)
    
    if ((nmax %%2) != 0) stop("nmax must be even")
    nz <- (nmax-2)/2
    zc <- numeric(nz)

    na <- D/2/rc
    p <- function(rho, na, b) {
        (rho*na)^2*(1 - 2/(1 + sqrt(1 - (1+b)*(rho*na)^2)) +
        (rho*na)^2/(1 + sqrt(1 - (1+b)*(rho*na)^2))^2)
    }
    q <- function(rho, na, b) 1 - sqrt(1 + p(rho, na, b))
    f <- function(rho, na, b, n) q(rho, na, b)*rzernike(rho, n, 0)*rho
    for (i in 1:nz) {
        n <- 2 + 2*i
        int <- integrate(f, lower=0, upper=1, na=na, b=b, n=n)
        zc[i] <- 4 * rc * sqrt(n+1) * int$value * 1.e6/lambda
        if (int$message != "OK") warning("coefficient may be inaccurate")
    }
    zc
}
