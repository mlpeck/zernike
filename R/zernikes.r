## This is free software and is offered with no warranty whatsoever.

##utility function returns true if n-m is odd

.odd <- function(n,m) {
  if (((n-m)%%2) == 0) {
    FALSE
  } else {
    TRUE
  }
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

## make a list of radial and azimuthal indexes in extended Fringe sequence
## from order minorder to maxorder

makezlist <- function(minorder=0, maxorder=14) {
    n <- numeric()
    m <- numeric()
    t <- character(length=0)

    for (order in seq(minorder, maxorder, by=2)) {
        mmax <- order/2
        mtemp <- 0
        if (mmax > 0) {
          mtemp <- c(mtemp, rep(1:mmax, each=2))
        }
        mtemp <- rev(mtemp)
        n <- c(n, order-mtemp)
        m <- c(m, mtemp)
        t <- c(t, rep(c("c", "s"), mmax), "n")
    }
    list(n=n, m=m, t=t)
}

## Augmented Fringe set

zlist.fr <- makezlist(2,12)

## list of ZP indexes in ISO/ANSI sequence

#' Construct list of ZP indexes in ISO/ANSI sequence with sine terms first
#'
#' @param maxorder maximum radial and azimuthal order
#' @return a named list with n=radial indexes, m=azimuthal, and t indicating
#'  trig function to apply
#' @examples
#' zlist.iso <- makezlist.iso(maxorder=6)
#' zlist <- makezlist(0, 6)
makezlist.iso <- function(maxorder=12) {
  n <- numeric(0)
  m <- numeric(0)
  t <- "n"
  for (i in 0:maxorder) {
    n <- c(n, rep(i, i+1))
    m <- c(m, abs(seq(-i, i, by=2)))
    if (i >= 1) {
      t <- c(t, rep("s" , (i+1) %/% 2))
      if (i %% 2 == 0) {
        t <- c(t, "n")
      }
      t <- c(t, rep("c", (i+1) %/% 2))
    }
  }
  list(n=n, m=m, t=t)
}

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


## fit zernikes to data

fitzernikes <- function(wf, rho, theta, eps=0, phi=0, maxorder = 14, 
                        nthreads=parallel::detectCores()/2, uselm=FALSE, 
                        isoseq=FALSE, usecirc=FALSE) {
  if (isoseq) {
    theta <- theta - pi * phi/180
    if ((eps == 0.) | usecirc) {
      zm <- zpm_cart(x=rho*cos(theta), y=rho*sin(theta), maxorder=maxorder)
    } else {
      zm <- zapm_iso(rho, theta, eps, maxorder=maxorder)
    }
  } else {
    if ((eps == 0) | usecirc) {
      zm <- zpm(rho, theta, phi=phi, maxorder=maxorder, nthreads=nthreads)
    } else {
      zm <- zapm(rho, theta - pi * phi/180, eps, maxorder=maxorder)
    }
  }
  if (uselm) {
    zm.names <- paste("Z", 0:(ncol(zm)-1), sep="")
    fmla <- as.formula(paste("wf ~ -1 + ", paste(zm.names, collapse="+")))
    dataf <- data.frame(cbind(wf, zm))
    fit <- lm(fmla, data=dataf)
  } else {
    fit <- qr.solve(crossprod(zm),crossprod(zm, wf))
  }
  fit
}


