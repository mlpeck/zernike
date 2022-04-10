## options for fringe analysis routines, wavefront fitting and wavefront display

psfit_options <- function(colors=topo.colors(256), refine=TRUE, puw_alg = "qual", fringescale=1,
                    wt=NULL, bgsub=TRUE,
                    maxiter=20, ptol=1.e-4, trace=1, nzcs = 2,
                    zc0=6:7,
                    satarget=c(0,0), astig.bath=c(0,0),
                    maxorder=14, uselm=FALSE, isoseq=FALSE, sgs=1,
                    nthreads=parallel::detectCores()/2,
                    plots=TRUE, crop=FALSE) {
  list(colors=colors, refine=refine, puw_alg=puw_alg, fringescale=fringescale,
       wt=wt, bgsub=bgsub,
       maxiter=maxiter, ptol=ptol, trace=trace, nzcs=nzcs,
       zc0=zc0, satarget=satarget, astig.bath=astig.bath,
       maxorder=maxorder, uselm=uselm, isoseq=isoseq, sgs=sgs,
       nthreads=nthreads,
       plots=plots, crop=crop)
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
                brcut = zernike::brcutpuw(phi),
                lpbrcut = lppuw::brcutpuw(phi),
                lp = lppuw::netflowpuw(phi, mod),
                qpuw(phi, mod)
  )
  wf.raw <- options$fringescale * wf.raw
  class(wf.raw) <- "pupil"
  wfnets <- wf_net(wf.raw, cp, options)
  if(length(psfit) > 3) extras <- psfit[4:length(psfit)]
  outs <- list(phi=phi, mod=mod, phases=wrap(as.vector(phases)), 
       cp=cp, cp.orig=cp.orig,
       wf.net=wfnets$wf.net, wf.smooth=wfnets$wf.smooth,wf.residual=wfnets$wf.residual,
       fit=wfnets$fit, zcoef.net=wfnets$zcoef.net, extras=extras)
  class(outs) <- append(class(outs), "wf_fitted")
  outs
}

