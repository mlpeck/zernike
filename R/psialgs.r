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

pcapsi <- function(im.mat, bgsub=TRUE, pcalg="pc1") {
    if (bgsub) im.mat <- im.mat-rowMeans(im.mat)

	# svd of the crossproduct is faster!
	
    svd.cp <- svd(crossprod(im.mat))
    if (tolower(pcalg) == "pc1") {
        ph <- atan2(svd.cp$v[,2]*sqrt(svd.cp$d[2]),svd.cp$v[,1]*sqrt(svd.cp$d[1]))
				wt <- sqrt(svd.cp$v[,2]^2*svd.cp$d[2] + svd.cp$v[,1]^2*svd.cp$d[1])
        u <- im.mat %*% (svd.cp$u[,1:2] %*% diag(1/sqrt(svd.cp$d[1:2])))
    } else {
        ph <- atan2(svd.cp$v[,2],svd.cp$v[,1])
				wt <- sqrt(svd.cp$v[,2]^2+svd.cp$v[,1]^2)
        u <- im.mat %*% svd.cp$u[,1:2]
    }
    if (tolower(pcalg) == "pc3") {
			rho <- sum(abs(u[,2]))/sum(abs(u[,1]))
			u[,1] <- rho*u[,1]
			ph <- atan2(svd.cp$v[,2], svd.cp$v[,1]/rho)
			wt <- sqrt(svd.cp$v[,2]^2+(svd.cp$v[,1]/rho)^2)
		}
    ph <- wrap(ph-ph[1])
		wt <- wt/mean(wt)
    phi <- atan2(-u[,2],u[,1])
    mod <- sqrt(u[,1]^2+u[,2]^2)
    r2 <- (svd.cp$d[1]+svd.cp$d[2])/sum(svd.cp$d)
    list(phi=phi, mod=mod/max(mod), phases=ph, wt=wt, snr=sqrt(r2/(1-r2)), eigen=svd.cp$d)
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

## tiltpsi (new version)

tiltpsi <- function(im.mat, phases, coords,
       		    ptol=0.01, maxiter=20, trace=1) {
  tiltpsiC(im.mat, phases, coords, ptol, maxiter, trace)
}

