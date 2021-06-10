## Zernike coefficients for a conic

zconic <- function(D, rc, b = -1, lambda = 632.8, nmax = 6) {
    if ((nmax %%2) != 0) stop("nmax must be even")
    nterms <- nmax/2+1
    sa <- zc <- zs <- numeric(0)
    cz <- matrix(0, nterms, nterms)

    sd <- D/2
    na <- (sd/rc)^2
    sa[1] <- sa[2] <- 0
## 4th order conic term
    zs[3] <- (sd/rc)^3*sd*1e6/8/lambda
    zc[3] <- (1+b)*zs[3]
    sa[3] <- zc[3]-zs[3]
## binomial expansion of conic
    if (nterms >= 4) {
        for (k in 4:nterms) {
            mult <- (1-3/(2*(k-1)))*na
            zc[k] <- mult*(1+b)*zc[k-1]
            zs[k] <- mult*zs[k-1]
            sa[k] <- zc[k]-zs[k]
        }
    }	
## recurrence relation for coefficients of radial zernikes
    cz[1,] <- (-1)^(0:(nterms-1))
    cz[2,2] <- 2
    for (j in 3:nterms) {
        for (i in 2:j) {
            n <- 2*(j-1)
            cz[i, j] <- (4*(n-1)*cz[i-1, j-1]-2*(n-1)*cz[i,j-1]-(n-2)*cz[i,j-2])/n
        }
    }
    for (j in 1:nterms) {
        n <- 2*(j-1)
        cz[,j] <- cz[,j]*sqrt(n+1)
    }
    solve(cz,sa)[-(1:2)]
}

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



## returns cosine and sine components of Bath astigmatism

##
# D = diameter
# rc = radius of curvature
# s = beam separation
# lambda = source wavelength
# phi = angle from horizontal (degrees)
##

astig.bath <- function(D, rc, s, lambda=632.8, phi=0) {
    astig.tot <- D^2 * s^2/(32 * sqrt(6) * lambda* 1.e-6 * rc^3)
    astig.tot*c(cos(pi*phi/90), sin(pi*phi/90))
}


