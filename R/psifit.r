psfit_options <- function(refine=TRUE, puw_alg = "qual", fringescale=1,
                    wt=NULL, bgsub=TRUE,
                    maxiter=20, ptol=1.e-4, trace=1, nzcs = 2,
                    zc0=c(1:3, 6:7),
                    satarget=c(0,0), astig.bath=c(0,0),
                    maxorder=14, uselm=FALSE, sgs=1,
                    plots=TRUE, crop=FALSE, colors=topo.colors(256)) {
  list(refine=refine, puw_alg=puw_alg, fringescale=fringescale,
       wt=wt, bgsub=bgsub,
       maxiter=maxiter, ptol=ptol, trace=trace, nzcs=nzcs,
       zc0=zc0, satarget=satarget, astig.bath=astig.bath,
       maxorder=maxorder, uselm=uselm, sgs=sgs,
       plots=plots, crop=crop, colors=colors)
}

psifit <- function(images, phases, cp=NULL, satarget=NULL, psialg ="ls", options=psfit_options()) {
  dims <- dim(images)
  nr <- dims[1]
  nc <- dims[2]
  nf <- dims[3]
  im.mat <- matrix(images, ncol=nf)
  if (!is.null(satarget)) {
    options$satarget <- satarget
  }
  refine <- options$refine
  if (length(phases) < nf) {
    phases <- wrap((0:(nf-1)) * phases[1])
  }
  if (!is.null(cp)) {
    prt <- pupil.rhotheta(nr, nc, cp)
    mask <- as.vector(prt$rho)
    im.mat <- im.mat[!is.na(mask),]
    refine <- FALSE
  } else {
    mask <- numeric(nr*nc)
  }
  switch(psialg,
    ls = {
      if (is.null(options$wt)) {
        wt <- rep(1, nf)
      } else {
        wt <- options$wt
      };
      psfit <- lspsiC(im.mat, phases, wt);
      extras <- NULL
    },
    aia = {
      if (is.null(options$ptol)) {
        ptol <- 0.001
      } else {
        ptol <- options$ptol
      };
      if (is.null(options$maxiter)) {
        maxiter <- 20
      } else {
        maxiter <- options$maxiter
      };
      if (is.null(options$trace)) {
        trace <- 1
      } else {
        trace <- options$trace
      };
      psfit <- aiapsiC(im.mat, phases, ptol, maxiter, trace);
      phases <- psfit$phases;
    },
    pc1 = {
       if (is.null(options$bgsub)) {
        bgsub <- TRUE
      } else {
        bgsub <- options$bgsub
      };
      group_diag <- "v";
      psfit <- pcapsi(im.mat, bgsub, group_diag);
      phases <- psfit$phases
    },
    pc2 = {
       if (is.null(options$bgsub)) {
        bgsub <- TRUE
      } else {
        bgsub <- options$bgsub
      };
      group_diag <- "u";
      psfit <- pcapsi(im.mat, bgsub, group_diag);
      phases <- psfit$phases
    },
    gpc = ,
    gpcthentilt = {
      if (is.null(options$ptol)) {
        ptol <- 0.001
      } else {
        ptol <- options$ptol
      };
      if (is.null(options$maxiter)) {
        maxiter <- 20
      } else {
        maxiter <- options$maxiter
      };
      if (is.null(options$trace)) {
        trace <- 1
      } else {
        trace <- options$trace
      };
      psfit <- gpcapsiC(im.mat, ptol, maxiter, trace);
      phases <- psfit$phases
    },
    tilt = {
      if (is.null(cp)) stop("This algorithm must have measured interferogram outline");
      if (is.null(options$ptol)) {
        ptol <- 0.001
      } else {
        ptol <- options$ptol
      };
      if (is.null(options$maxiter)) {
        maxiter <- 20
      } else {
        maxiter <- options$maxiter
      };
      if (is.null(options$trace)) {
        trace <- 1
      } else {
        trace <- options$trace
      };
      if (is.null(options$nzcs)) {
        nzcs <- 2
      } else {
        nzcs <- min(options$nzcs, 8)
      };
      rho <- prt$rho;
      theta <- prt$theta;
      rho <- rho[!is.na(rho)];
      theta <- theta[!is.na(theta)];
      coords <- zpmC(rho, theta, maxorder=4);
      coords <- coords[, 2:(nzcs+1)]
      psfit <- tiltpsiC(im.mat, phases, coords,
                       maxiter=maxiter, ptol=ptol, trace=trace);
      phases <- psfit$phases;
    },
    { ## if no match to ls
      if (is.null(options$wt)) {
        wt <- rep(1, nf)
      } else {
        wt <- options$wt
      };
      psfit <- lspsiC(im.mat, phases, wt);
      extras <- NULL
    }
  )
  phi <- matrix(NA, nr, nc)
  mod <- matrix(0, nr, nc)
  phi[!is.na(mask)] <- psfit$phi
  mod[!is.na(mask)] <- psfit$mod
  if (is.null(cp)) {
    phi <- matrix(psfit$phi, ncol=nc)
    mod <- matrix(psfit$mod, ncol=nc)
    cp <- circle.pars(mod, plot=options$plots)
    prt <- pupil.rhotheta(nr, nc, cp)
  }
  if (refine || psialg=="gpcthentilt") {
    mask <- as.vector(prt$rho)
    im.mat <- matrix(images, ncol=nf)
    im.mat <- im.mat[!is.na(mask),]
    switch(psialg,
      ls = {
        mask <- numeric(nr*nc)
      },
      aia = {
        psfit <- aiapsiC(im.mat, phases, ptol, maxiter, trace);
        phases <- psfit$phases;
      },
      pc1 = {
        psfit <- pcapsi(im.mat, bgsub, group_diag);
        phases <- psfit$phases
      },
      pc2 = {
        psfit <- pcapsi(im.mat, bgsub, group_diag);
        phases <- psfit$phases
      },
      gpc = {
        psfit <- gpcapsiC(im.mat, ptol, maxiter, trace);
        phases <- psfit$phases
      },
      gpcthentilt = ,
      tilt = {
        if (is.null(options$ptol)) {
          ptol <- 0.001
        } else {
          ptol <- options$ptol
        };
        if (is.null(options$maxiter)) {
          maxiter <- 20
        } else {
          maxiter <- options$maxiter
        };
        if (is.null(options$trace)) {
          trace <- 1
        } else {
          trace <- options$trace
        };
        if (is.null(options$nzcs)) {
          nzcs <- 2
        } else {
          nzcs <- min(options$nzcs, 8)
        };
        rho <- prt$rho;
        theta <- prt$theta;
        rho <- rho[!is.na(rho)];
        theta <- theta[!is.na(theta)];
        coords <- zpmC(rho, theta, maxorder=4);
        coords <- coords[, 2:(nzcs+1)]
        psfit <- tiltpsiC(im.mat, phases, coords,
                          maxiter=maxiter, ptol=ptol, trace=trace);
        phases <- psfit$phases;
      }
    )
    phi <- matrix(NA, nr, nc)
    mod <- matrix(0, nr, nc)
    phi[!is.na(mask)] <- psfit$phi
    mod[!is.na(mask)] <- psfit$mod
  }
  phi[is.na(prt$rho)] <- NA
  class(phi) <- "pupil"
  cp.orig <- cp
  if (options$crop) {
    phi <- crop(phi, cp)
    mod <- crop(mod, cp)$im
    cp <- phi$cp
    phi <- phi$im
    nr <- nrow(phi)
    nc <- ncol(phi)
  }
  wf.raw <- switch(options$puw_alg,
                qual = qpuw(phi, mod),
                brcut = brcutpuw(phi),
                lp = lppuw::netflowpuw(phi, mod)
  )
  wf.raw <- options$fringescale * wf.raw
  class(wf.raw) <- "pupil"
  wfnets <- wf_net(wf.raw, cp, options)
  if(length(psfit) > 3) extras <- psfit[4:length(psfit)]
  list(phi=phi, mod=mod, phases=wrap(as.vector(phases)), cp=cp.orig,
       wf.net=wfnets$wf.net, wf.smooth=wfnets$wf.smooth,wf.residual=wfnets$wf.residual,
       fit=wfnets$fit, zcoef.net=wfnets$zcoef.net, extras=extras)
}

wf_net <- function(wf.raw, cp, options) {
  zlist <- makezlist(maxorder=options$maxorder)
  nr <- nrow(wf.raw)
  nc <- ncol(wf.raw)
  prt <- pupil.rhotheta(nr, nc, cp)
  rho <- prt$rho
  theta <- prt$theta
  if (options$sgs > 1) {
    subs <- matrix(FALSE, nr, nc)
    subs[seq(1, nr, by=options$sgs), seq(1, nc, by=options$sgs)] <- TRUE
    subs[is.na(rho)] <- FALSE
  } else {
    subs <- !is.na(rho)
  }
  wf.v <- wf.raw[subs]
  rho.v <- rho[subs]
  th.v <- theta[subs]
  rho.v <- rho.v[!is.na(wf.v)]
  th.v <- th.v[!is.na(wf.v)]
  wf.v <- wf.v[!is.na(wf.v)]
  fit <- fitzernikes(wf.v, rho.v, th.v, maxorder=options$maxorder, uselm=options$uselm)
  if (options$uselm) {
    cfit <- coef(fit)
  } else {
    cfit <- fit
  }
  if (sign(cfit[9])*sign(options$satarget[1]) < 0) {
    cfit <- -cfit
    wf.raw <- -wf.raw
    if (!options$uselm) fit <- -fit
  }
  zc.low <- rep(0,15)
  zc.low[options$zc0] <- cfit[options$zc0+1]
  zc.low[c(8,15)] <- zc.low[c(8,15)] + options$satarget
  zc.low[4:5] <- zc.low[4:5] + options$astig.bath
  wf.net <- wf.raw - pupil(zcoef=zc.low, zlist=makezlist(2,6), piston=cfit[1], 
                           nrow=nr, ncol=nc, cp=cp)
  if (options$plots) {
    if (tolower(.Platform$OS.type) == "windows") {
      windows(width=18, height=6)
    } else {
        x11(width=18, height=6)
    }
    split.screen(figs=c(1,3))
    screen(1)
    plot(wf.net, cp=cp, col=options$colors, addContours=FALSE)
    mtext(paste("RMS = ", format(pupilrms(wf.net),digits=3)))
  }
  zcoef.net <- cfit[-1]
  zcoef.net[1:15] <- zcoef.net[1:15] - zc.low
  wf.smooth <- pupil(zcoef=zcoef.net, zlist=zlist, cp=cp, nrow=nr, ncol=nc)
  if (options$plots) {
    screen(2)
    plot(wf.smooth, cp=cp, col=options$colors)
    mtext(paste("RMS = ", format(pupilrms(wf.smooth),digits=3)))
  }
  wf.residual <- wf.net - wf.smooth
  if (options$plots) {
    screen(3)
    plot(wf.residual, cp=cp, col=grey256, addContours=FALSE)
    mtext(paste("RMS = ", format(pupilrms(wf.residual),digits=3)))
    close.screen(all.screens=TRUE)
  }
  list(wf.net=wf.net, wf.smooth=wf.smooth, wf.residual=wf.residual, 
       fit=fit, zcoef.net=zcoef.net)
}
          
