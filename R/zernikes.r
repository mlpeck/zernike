## This is free software and is offered with no warranty whatsoever.

##utility function returns true if n-m is odd

.odd <- function(n,m) {
	if (((n-m)%%2) == 0) return(FALSE)
	else return(TRUE)
}


##derivative (wrt r) of zernike polynomials
## Iterative version of recurrence relation
## given in Eqn. 13 of Mathworld Zernike article.

drzernike <- function(rho, n, m) {
    if ((n<0) || (m<0) || (n<m) || (.odd(n,m))) stop("Bad argument to drzernike")
    if ((n==0) && (m==0))return(0*rho)
    if (n==m) return(n*rho^(n-1))
    if (m==0) dr <- 0*rho else dr <- m*rho^(m-1)
    for (nn in seq(m+2, n, by=2)) {
        dr <- dr + nn*(rzernike(rho,nn-1,abs(m-1))+rzernike(rho,nn-1,m+1))
    }
    dr
}


### Zernike polynomial

Zernike <- function(rho, theta, n, m, t) {
    sqrt(n+1) * rzernike(rho, n, m) * switch(t, n=1, c=sqrt(2)*cos(m*theta),
	 s=sqrt(2)*sin(m*theta))
}

### Radial derivative of zernike polynomial

DZernike <- function(rho, theta, n, m, t) {
    sqrt(n+1) * drzernike(rho, n, m) * switch(t, n=1, c=sqrt(2)*cos(m*theta),
	 s=sqrt(2)*sin(m*theta))
}

### Tangential derivative of zernike polynomial

DTZernike <- function(rho, theta, n, m, t) {
    if (m == 0) return (rho*0) else
    sqrt(n+1) * m * rzernike(rho, n, m) * 
                switch(t, n=0, c= -sqrt(2)*sin(m*theta), s = sqrt(2)*cos(m*theta))
}

## make a list of all orders up to maxorder

makezlist <- function(minorder=2, maxorder=14) {
    n <- numeric()
    m <- numeric()
    t <- character(length=0)

    for (order in seq(minorder, maxorder, by=2)) {
        mmax <- order/2
        mtemp <- numeric()
        for (j in mmax:1) mtemp <- c(mtemp, c(j, j))
        mtemp <- c(mtemp, 0)
        n <- c(n, order-mtemp)
        m <- c(m, mtemp)
        t <- c(t, rep(c("c", "s"), mmax), "n")
    }
    list(n=n, m=m, t=t)
}

## Augmented Fringe set

zlist.fr <- makezlist(2,12)


## Vector of factors from conversion between "normalized" and conventional Zernike definitions

zmult <- function(zlist = makezlist()) {
    mult <- sqrt(zlist$n+1)
    mult[zlist$m > 0] <- sqrt(2)*mult[zlist$m > 0]
    mult
}



## create a matrix of Zernike polynomial values


zpm.arb <- function(rho, theta, phi=0, zlist=makezlist()) {
	zm <- matrix(0, nrow=length(rho), ncol=length(zlist$n))
	for (i in (1:length(zlist$n))) {
		zm[,i] <- Zernike(rho,theta-pi*phi/180, zlist$n[i], zlist$m[i], zlist$t[i])
	}
	colnames(zm) <- colnames(zm, do.NULL=FALSE, prefix="Z")
	zm
}

## a faster Zernike matrix fill

zpm <- function(rho, theta, phi=0, maxorder = 14, nthreads=parallel::detectCores()/2) {
  if (phi != 0) theta <- theta - pi*phi/180
  if (nthreads == 1) {
    zm <- zpmC(rho, theta, maxorder)
  } else {
    RcppParallel::setThreadOptions(numThreads = nthreads)
    zm <- zpmCP(rho, theta, maxorder)
  }
  colnames(zm) <- paste("Z",0:(ncol(zm)-1), sep="")
  zm
}


## create a matrix of radial derivatives of Zernikes

filldzm <- function(rho, theta, phi=0, zlist=makezlist()) {
    dzm <- matrix(0, nrow=length(rho), ncol=length(zlist$n))
    for (i in (1:length(zlist$n))) {
            dzm[,i] <- DZernike(rho,theta-pi*phi/180, zlist$n[i], zlist$m[i], zlist$t[i])
    }
    colnames(dzm) <- colnames(dzm, do.NULL=FALSE, prefix="DZ")
    dzm
}

## create a matrix of gradient in polar coordinates of zernikes.
## radial derivative is on top, followed by tangential derivative.

fillgradientzm <- function(rho, theta, phi=0, zlist=makezlist()) {
    nr <- length(rho)
    drzm <- matrix(0, nrow=nr, ncol=length(zlist$n))
    dtzm <- matrix(0, nrow=nr, ncol=length(zlist$n))
    for (i in (1:length(zlist$n))) {
        drzm[,i] <- DZernike(rho, theta-pi*phi/180, zlist$n[i], zlist$m[i], zlist$t[i])
        dtzm[,i] <- DTZernike(rho, theta-pi*phi/180, zlist$n[i], zlist$m[i], zlist$t[i])/rho
        dtzm[rho==0,i] <- 0
    }
    gradzm <- rbind(drzm, dtzm)
    colnames(gradzm) <- colnames(gradzm, do.NULL=FALSE, prefix="Z")
    gradzm
}

## fit zernikes to data

fitzernikes <- function(wf, rho, theta, phi=0, maxorder = 14 , uselm=FALSE) {
    zm <- zpm(rho, theta, phi=phi, maxorder=maxorder)
    zm.names <- colnames(zm)
    if (uselm) {
        fmla <- as.formula(paste("wf ~ -1 + ", paste(zm.names, collapse="+")))
        dataf <- data.frame(cbind(wf, zm))
        fit <- lm(fmla, data=dataf)
    } else {
        fit <- qr.solve(crossprod(zm),crossprod(zm, wf))
    }
    fit
}

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

