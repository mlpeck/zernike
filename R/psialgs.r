#############
## PSI routines
#############

## The following are low level routines to calculate the phase from various PSI algorithms.

## least squares fit of phase shifted interferograms with optional weights.

lspsi <- function(images, phases, wt=rep(1, length(phases))) {
  img.mat <- matrix(images, ncol=dim(images)[3])
  lspsiC(img.mat, phases, wt)
}

## Vargas et al.'s (2011) Principal Components method

pcapsi <- function(im.mat, BGSUB=TRUE, diagpos="v") {
    if (BGSUB) im.mat <- im.mat-rowMeans(im.mat)

	# svd of the crossproduct is faster!
	
    svd.cp <- svd(crossprod(im.mat))
    if (tolower(diagpos) == "v") {
        ph <- atan2(svd.cp$v[,2]*sqrt(svd.cp$d[2]),svd.cp$v[,1]*sqrt(svd.cp$d[1]))
        u <- im.mat %*% (svd.cp$u[,1:2] %*% diag(1/sqrt(svd.cp$d[1:2])))
    } else {
        ph <- atan2(svd.cp$v[,2],svd.cp$v[,1])
        u <- im.mat %*% svd.cp$u[,1:2]
    }
    ph <- wrap(ph-ph[1])
    phi <- atan2(-u[,2],u[,1])
    mod <- sqrt(u[,1]^2+u[,2]^2)
    r2 <- (svd.cp$d[1]+svd.cp$d[2])/sum(svd.cp$d)
    list(phi=phi, mod=mod/max(mod), phases=ph, snr=sqrt(r2/(1-r2)), eigen=svd.cp$d)
}

## Moore-Penrose generalized inverse (needed for genpca below)

mpinv <- function(X) {
    S <- svd(X)
    eps <- .Machine$double.eps * max(dim(X)) * S$d[1]
    dinv <- numeric(length(S$d))
    dinv[S$d >= eps] <- 1/S$d[S$d >= eps]
    tcrossprod(S$v %*% diag(dinv), S$u)
}

## my variation on a PCA based algorithm

gpcapsi <- function(im.mat, ptol=0.001, maxiter=20, trace=1) {
  gpcapsiC(im.mat, ptol, maxiter, trace)
}


## The "advanced iterative algorithm" of Wang & Han (2004)

aiapsi <- function(im.mat, phases,
		   ptol=0.001, maxiter=20, trace=1) {
  aiapsiC(im.mat, phases, ptol, maxiter, trace)
}


## A variation of Han & Kim's (1994) iterative algorithm

hkpsi <- function(im.mat, phases,
		maxiter=20, ptol=0.001, trace=1, plotprogress=TRUE) {
	.meps <- sqrt(.Machine$double.eps)
	M <- nrow(im.mat)
	nf <- ncol(im.mat)
	phases <- wrap(phases-phases[1])
	S <- rbind(rep(1,nf), cos(phases), sin(phases))
	sse <- numeric(maxiter+1)
	for (i in 1:maxiter) {
		phases.last <- phases
		
		# get the phase estimate from phase shifts
		
		SD <- svd(S)
		dp <- SD$d
		dp[dp>.meps] <- 1/dp[dp>.meps]
		dp[dp <= .meps] <- 0
		Phi <- tcrossprod(im.mat %*% (SD$v %*% diag(dp)), SD$u)
		sse[i] <- crossprod(as.vector(im.mat)-as.vector(Phi %*% S))
		if (i == 1) sse.1 <- sse[i]
		
		# get phase shifts from the phase. Note use crossproduct for speed.
		# Here is only difference from AIA - use full least squares estimate
		# at this step.
		
		PD <- svd(crossprod(Phi))
		dp <- PD$d
		dp[dp>.meps] <- 1/dp[dp>.meps]
		dp[dp <= .meps] <- 0
		S <- (PD$u %*% diag(dp) %*% t(PD$v)) %*% crossprod(Phi, im.mat)
		phases <- atan2(S[3,], S[2,])
		dphases <- sd(wrap(phases-phases.last))
		if (plotprogress) {
			if (i ==1)
				plot(1:maxiter, 1:maxiter, ylim=c(ptol, 1), type="n",
				  xlab="Iteration", ylab="", log="y")
			points(i, sse[i]/sse.1, pch=1)
			points(i, dphases, pch=2, col="green")
		}
		if ((trace > 0) && ((i-1)%%trace == 0)) {
		  cat(paste(i, ":", format(sse[i], digits=2), ":"), 
                      format(phases,digits=3), "\n")
		  flush.console()
		}
		if (dphases < ptol) break
		S <- rbind(rep(1,nf), cos(phases), sin(phases))
	}

	#final estimate of phase

	phases <- wrap(phases-phases[1])
	S <- rbind(rep(1,nf), cos(phases), sin(phases))
	SD <- svd(S)
	dp <- SD$d
	dp[dp>.meps] <- 1/dp[dp>.meps]
	dp[dp <= .meps] <- 0
	Phi <- tcrossprod(im.mat %*% (SD$v %*% diag(dp)), SD$u)
	sse[i+1] <- crossprod(as.vector(im.mat)-as.vector(Phi %*% S))
	phi <- atan2(-Phi[,3], Phi[,2])
	mod <- sqrt(Phi[,2]^2+Phi[,3]^2)
	list(phi=phi, mod=mod/max(mod), phases=phases, iter=i, sse=sse)
}

## PSI with variable tilt

tiltpsi <- function(im.mat, phases, x, y, tilts = NULL,
              nlpref=1, tlim=0.5, maxiter=20, ptol=0.001, trace=1, plotprogress=TRUE) {
  if (nlpref>1) require(Rsolnp)
  M <- nrow(im.mat)
  nf <- ncol(im.mat)
  phases <- wrap(phases-phases[1])
  gotfirst <- FALSE

  # if no guess at tilts set them to 0 and calculate the phase the fast way
  
  if (is.null(tilts)) {
    tilts <- matrix(0, nf, 2)
    X <- cbind(rep(1,nf), cos(phases), sin(phases))
    X <- X %*% solve(crossprod(X))
    B <- im.mat %*% X
    gotfirst <- TRUE
  }
  tilts[,1] <- tilts[,1]-tilts[1,1]
  tilts[,2] <- tilts[,2]-tilts[1,2]

  
  # this is the function minimized at the phase shift/tilt estimation step

  sseframe <- function(pt, im, phi, x, y, abar, bbar) {
    ph.xy <- pt[1] + 4*pi*(pt[2]*x+pt[3]*y)
    crossprod(im-abar-bbar*cos(phi+ph.xy))
  }

  # gradient of the above function. used by nlminb

  gradsse <- function(pt, im, phi, x, y, abar, bbar) {
    ph.xy <- phi+pt[1]+4*pi*(pt[2]*x+pt[3]*y)
    cpt <- cos(ph.xy)
    spt <- sin(ph.xy)
    grads <- numeric(3)
    grads[1] <- 2*bbar*(-abar*sum(spt)+crossprod(im,spt)-bbar*crossprod(cpt,spt))
    grads[2] <- 8*pi*bbar*(-abar*crossprod(x,spt)+crossprod((im*x),spt)-bbar*crossprod((x*cpt),spt))
    grads[3] <- 8*pi*bbar*(-abar*crossprod(y,spt)+crossprod((im*y),spt)-bbar*crossprod((y*cpt),spt))
    grads
  }


  pt.last <- cbind(phases, 2*pi*tilts)
  nltrace <- max(trace, 0)
  sse <- numeric(maxiter+1)
  sse.last <- 0

  for (i in 1:maxiter) {
    if ((i>1) || !gotfirst) B <- pxls(im.mat, phases, tilts, x, y)

    # first column of B contains the background estimate; next two the components of phase
    # note: either mean or median seem to work here

    abar <- mean(B[,1])
    bbar <- mean(sqrt(B[,2]^2+B[,3]^2))
    phi <- atan2(-B[,3], B[,2])
    sse[i] <- 0

    # use one of two nonlinear optimizers for the tilt and phase shift
    # estimation step
    
    for (n in 1:nf) {
      if (nlpref == 1) {
        smin <- nlminb(start=c(phases[n], tilts[n,]), objective=sseframe,
          gradient=gradsse,
          im=im.mat[,n], phi=phi, x=x, y=y, abar=abar, bbar=bbar,
          control = list(trace=nltrace, abs.tol=1e-20),
          lower=c(-2*pi, -tlim, -tlim), upper=c(+2*pi, tlim, tlim))
        sse[i] <- sse[i] + smin$objective
      } else {
        smin <- solnp(pars=c(phases[n], tilts[n,]), fun=sseframe,
          LB=c(-2*pi, -tlim, -tlim), UB=c(+2*pi, tlim, tlim),
          control = list(trace=nltrace),
          im=im.mat[,n], phi=phi, x=x, y=y, abar=abar, bbar=bbar)
        sse[i] <- sse[i] + smin$values[length(smin$values)]
      }
      phases[n] <- smin$par[1]
      tilts[n,] <- smin$par[2:3]
    }
    if(smin$convergence > 0) warning("Convergence reported failed")

    # the tilts and phase shifts are offset from frame 1
    
    tilts[,1] <- tilts[,1]-tilts[1,1]
    tilts[,2] <- tilts[,2]-tilts[1,2]
    phases <- wrap(phases-phases[1])
    pt <- cbind(phases, 2*pi*tilts)
    dp <- sqrt(sum(diag(crossprod(pt-pt.last))))/(3*nf)
    pt.last <- pt
    if (plotprogress) {
        if (i ==1) {
            sse.1 <- sse[1]
            plot(1:maxiter, 1:maxiter, ylim=c(ptol, 1), type="n",
              xlab="Iteration", ylab="", log="y")
        }
        points(i, sse[i]/sse.1, pch=1)
        points(i, dp, pch=2, col="green")
    }

    # this is slow, so it's best to print out some intermediate results so
    # the user knows something is happening.

    if (trace >= 0) print(paste("Iteration",i, "sse", format(sse[i], digits=2),
      "delta sse",format((sse[i]-sse.last)/sse.last, digits=2),
      "dp =", format(dp, digits=3)))
    if (dp < ptol) break
    sse.last <- sse[i]
  }
  B <- pxls(im.mat, phases, tilts, x, y)
  phi <- atan2(-B[,3], B[,2])
  mod <- sqrt(B[,2]^2+B[,3]^2)
  list(phi=phi,mod=mod/max(mod),phases=phases,tilts=tilts,iter=i,sse=sse)

}

