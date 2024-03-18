## options for fringe analysis routines, wavefront fitting and wavefront display

psfit_options <- function(colors=topo.colors(256), refine=TRUE, puw_alg = "qual", fringescale=1,
                    wt=NULL, bgsub=TRUE,
                    maxiter=20, ptol=1.e-4, trace=1, nzcs = 2,
                    zc0=6:7,
                    satarget=c(0,0), astig.bath=c(0,0),
                    maxorder=14, uselm=FALSE, sgs=1,
                    isoseq=FALSE, usecirc=FALSE, 
                    nthreads=parallel::detectCores()/2,
                    plots=TRUE, crop=FALSE) {
  list(colors=colors,             ## color palette for wavefront plots
       refine=refine,             ## do 2nd pass through psi algorithm?
       puw_alg=puw_alg,           ## phase unwrapping algorithm
       fringescale=fringescale,   ## waves per fringe in interferograms (usually 1 or 1/2)
       wt=wt,                     ## per frame weights in least squares PSI algorithm
       bgsub=bgsub,               ## subtract a background estimate in PC based PSI algorithms?
       maxiter=maxiter,           ## maximum no iterations for iterative PSI algorithms
       ptol=ptol,                 ## convergence tolerance for iterative PSI algorithms
       trace=trace,               ## some iterative algorithms can return info while working
       nzcs=nzcs,                 ## number of zernike coefficients to treat as variable in tiltpsiC
       zc0=zc0,                   ## coefficients to remove in net wavefront
       satarget=satarget,         ## target SA for numerical nulling
       astig.bath=astig.bath,     ## amount of astigmatism from Bath interferometer geometry
       maxorder=maxorder,         ## maximum order for Zernike polynomial fitting
       uselm=uselm,               ## use R function lm() to fit Zernikes to wavefront data?
       sgs=sgs,                   ## grid size for wavefront sampling
       isoseq=isoseq,             ## use ISO sequenced Zernikes?
       usecirc=usecirc,           ## use circular Zernikes even for obstructed apertures?
       nthreads=nthreads,         ## no threads to use with zpmCP
       plots=plots,               ## plot results in wavefront summaries?
       crop=crop)                 ## crop wavefront related matrixes?
}

## "high level" psi analysis function

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
      }
      psfit <- lspsiC(im.mat, phases, wt)
      extras <- NULL
    },
    aia = {
      if (is.null(options$ptol)) {
        ptol <- 0.001
      } else {
        ptol <- options$ptol
      }
      if (is.null(options$maxiter)) {
        maxiter <- 20
      } else {
        maxiter <- options$maxiter
      }
      if (is.null(options$trace)) {
        trace <- 1
      } else {
        trace <- options$trace
      }
      psfit <- aiapsiC(im.mat, phases, ptol, maxiter, trace)
      phases <- psfit$phases
    },
    pc1thenaia = ,
    pc1 = ,
    pc2 = ,
    pc3 = {
       if (is.null(options$bgsub)) {
        bgsub <- TRUE
      } else {
        bgsub <- options$bgsub
      }
      psfit <- pcapsi(im.mat, bgsub, pcalg = psialg)
      phases <- psfit$phases
    },
    gpc = ,
    gpcthentilt = {
      if (is.null(options$ptol)) {
        ptol <- 0.001
      } else {
        ptol <- options$ptol
      }
      if (is.null(options$maxiter)) {
        maxiter <- 20
      } else {
        maxiter <- options$maxiter
      }
      if (is.null(options$trace)) {
        trace <- 1
      } else {
        trace <- options$trace
      }
      psfit <- gpcapsiC(im.mat, ptol, maxiter, trace)
      phases <- psfit$phases
    },
    tilt = {
      if (is.null(cp)) stop("This algorithm must have measured interferogram outline")
      if (is.null(options$ptol)) {
        ptol <- 0.001
      } else {
        ptol <- options$ptol
      }
      if (is.null(options$maxiter)) {
        maxiter <- 20
      } else {
        maxiter <- options$maxiter
      }
      if (is.null(options$trace)) {
        trace <- 1
      } else {
        trace <- options$trace
      }
      if (is.null(options$nzcs)) {
        nzcs <- 2
      } else {
        nzcs <- min(options$nzcs, 8)
      }
      rho <- prt$rho
      theta <- prt$theta
      rho <- rho[!is.na(rho)]
      theta <- theta[!is.na(theta)]
      if (cp$obstruct == 0. || options$usecirc) {
        coords <- zpmC(rho, theta, maxorder=4)
      } else {
        coords <- zapm(rho, theta, eps=cp$obstruct, maxorder=4)
      }
      coords <- coords[, 2:(nzcs+1)]
      psfit <- tiltpsiC(im.mat, phases, coords,
                       maxiter=maxiter, ptol=ptol, trace=trace)
      phases <- psfit$phases
    },
    { ## if no match to ls
      if (is.null(options$wt)) {
        wt <- rep(1, nf)
      } else {
        wt <- options$wt
      }
      psfit <- lspsiC(im.mat, phases, wt)
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
  if (refine || psialg=="gpcthentilt" || psialg =="pc1thenaia") {
    mask <- as.vector(prt$rho)
    im.mat <- matrix(images, ncol=nf)
    im.mat <- im.mat[!is.na(mask),]
    switch(psialg,
      ls = {
        mask <- numeric(nr*nc)
      },
      pc1thenaia = ,
      aia = {
        psfit <- aiapsiC(im.mat, phases, ptol, maxiter, trace)
        phases <- psfit$phases
      },
      pc1 = ,
      pc2 = ,
      pc3 = {
        psfit <- pcapsi(im.mat, bgsub, pcalg = psialg)
        phases <- psfit$phases
      },
      gpc = {
        psfit <- gpcapsiC(im.mat, ptol, maxiter, trace)
        phases <- psfit$phases
      },
      gpcthentilt = ,
      tilt = {
        if (is.null(options$ptol)) {
          ptol <- 0.001
        } else {
          ptol <- options$ptol
        }
        if (is.null(options$maxiter)) {
          maxiter <- 20
        } else {
          maxiter <- options$maxiter
        }
        if (is.null(options$trace)) {
          trace <- 1
        } else {
          trace <- options$trace
        }
        if (is.null(options$nzcs)) {
          nzcs <- 2
        } else {
          nzcs <- min(options$nzcs, 8)
        }
        rho <- prt$rho
        theta <- prt$theta
        rho <- rho[!is.na(rho)]
        theta <- theta[!is.na(theta)]
        if (cp$obstruct == 0. || options$usecirc) {
          coords <- zpmC(rho, theta, maxorder=4)
        } else {
          coords <- zapm(rho, theta, eps=cp$obstruct, maxorder=4)
        }
        coords <- coords[, 2:(nzcs+1)]
        psfit <- tiltpsiC(im.mat, phases, coords,
                          maxiter=maxiter, ptol=ptol, trace=trace)
        phases <- psfit$phases
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
  if(length(psfit) > 3) {
    extras <- psfit[4:length(psfit)]
  }
  rundate <- date()
  algorithm <- paste(psialg, "with", nf, "frames")
  outs <- list(rundate=rundate, algorithm=algorithm,
       phi=phi, mod=mod, phases=wrap(as.vector(phases)), 
       cp=cp, cp.orig=cp.orig,
       wf.net=wfnets$wf.net, wf.smooth=wfnets$wf.smooth,wf.residual=wfnets$wf.residual,
       fit=wfnets$fit, zcoef.net=wfnets$zcoef.net, extras=extras)
  class(outs) <- c(class(outs), "wf_zfit")
  outs
}

