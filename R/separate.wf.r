## simple utility returns Euclidean length of a vector

hypot <- function(x) sqrt(crossprod(x))

## Zernike moments

zmoments <- function(zcoef, maxorder=14) {
  zcoef <- zcoef[-(1:3)]
  zlist <- makezlist(4, maxorder)
  nz <- length(zlist$n)
  sumstats <- NULL
  i <- 1
  repeat {
    if (i > nz) break
      if (zlist$m[i] == 0) {
        sumstats <- rbind(sumstats, c(zlist$n[i], zlist$m[i], zcoef[i]))
        i <- i+1
      }
      else {
        sumstats <- rbind(sumstats, c(zlist$n[i], zlist$m[i], hypot(zcoef[i:(i+1)])))
        i <- i+2
      }
  }
  sumstats <- data.frame(sumstats)
  names(sumstats) <- c("n", "m", "value")
  sumstats
}

## these two functions provide rudimentary project management. The first adds
## Zernike coefficients and rotation angles from one or more fits to a matrix.
## The second separates polished in from instrumental aberrations (if possible)

#' @describeIn separate.wf Add a wavefront fit to tracking list
#'
addfit <- function(..., th=0, zcm=NULL, theta=numeric(0)) {
  fits <- list(...)
  nt <- length(fits)
  if (length(th)==1 && nt>1) {
    th=rep(th, nt)
  }
  theta <- c(theta, th*pi/180)
  nz <- length(fits[[1]]$zcoef.net)
  for (i in 1:nt) {
    zcm <- cbind(zcm, fits[[i]]$zcoef.net[4:nz])
  }
  list(zcm=zcm, theta=theta)
}

 
#' Separate wavefronts
#'
#' Separate a set of wavefronts measured at different orientations into "polished in" and
#'   instrumental or "test stand" components.
#'
#' 
#' @param ... one or more wavefront fits as returned by [psifit], [vortexfit], [fftfit], or [wf_net].
#' @param th  the orientation angle(s) of object under test, in degrees.
#' @param zcm Zernike coefficient matrix from a previous call, or left empty to create a new one.
#' @param theta orientation angles in radians from a previous call, or leave empty to create a new one.
#'
#' @param zcm Zernike coefficient matrix from addfit.
#' @param theta vector of angles from addfit.
#' @param maxorder the maximum Zernike polynomial order of the fits.
#' @param nrow number of rows in the reconstructed wavefronts.
#' @param ncol number of columns in the reconstructed wavefronts.
#' @param cp list of values describing the location and size of the wavefront as returned by [circle.pars] or [pupil.pars].
#'
#' @return
#'
#' @param zcm Zernike coefficient matrix from the `zcoef.net` entry in the wavefront fits, minus the first 3 elements.
#` @param theta vector of orientation angles in radians, with the same number of elements as columns of `zcm`.
#'
#' @param zcb a data frame with zernike coefficients and standard errors of estimates for intrinsic and instrumental aberrations.
#' @param sumstats a data frame with summary statistics describing the fits to each set of coefficients.
#' @param wf.mirror the estimated derotated wavefront, stored in a matrix of size nrow x ncol with class `pupil`.
#'
#'
#' @details
#' The two functions `addfit()` and `separate.wf()` work together to provide a rudimentary project management capability when multiple
#'  optical tests have been run on an optical system rotated to one or more orientations.
#'  All `addfit()` does is extract the net Zernike coefficients in the list of values returned by
#'  [psifit], [vortexfit], [fftfit], or [wf_net]. It also accepts one or more orientation angles in degrees.
#'
#' The odd parameter list is to provide some flexibility in data entry. For example if several fits are available they could be entered as a group
#'  with the `theta`s entered as a vector the same length as the number of fits. Alternately fits could be entered one at a time.
#'  If the latter strategy is followed be sure to recycle the variable name that the return is assigned to. Also, the Zernike polynomial fits
#'  should be made to the same polynomial order. This isn't checked and will surely cause an error if different fit orders are used.
#'
#' `separate.wf()` makes use of the properties of Zernike polynomials under rotations to disentangle the contributions from the "mirror"
#'  and the instrument or test stand to the extent possible. Least squares fits are performed for each non-axysmmetric aberration and
#'  some possibly useful summary statistics from the fits are returned in the data frames `zcb` and `sumstats`.
#'
#' @describeIn separate.wf Separate wavefronts
separate.wf <- function(zcm, theta, maxorder=14,
                        nrow=nrow.default, ncol=ncol.default,
                        cp=cp.default) {
  
  nt <- length(theta)
  zlist <- makezlist(4, maxorder)
  nz <- length(zlist$n)
  zcb <- matrix(0, nz, 4)
  sumstats <- NULL
  cx <- c(rep(1,nt),rep(0,nt))
  cy <- c(rep(0,nt),rep(1,nt))
  i <- 1
  repeat {
    if (i > nz) break
      if (zlist$m[i] == 0) {
        zcb[i,1] <- mean(zcm[i,])
        zcb[i,3] <- sd(zcm[i,])/sqrt(nt)
        sumstats <- rbind(sumstats, c(zlist$n[i], zlist$m[i], zcb[i,1],
                                      NA, zcb[i, 3], rep(NA, 4)))
        i <- i+1
      } else {
        y <- c(zcm[i,],zcm[i+1,])
        rx <- c(cos(zlist$m[i]*theta),sin(zlist$m[i]*theta))
        ry <- c(-sin(zlist$m[i]*theta),cos(zlist$m[i]*theta))
        lsm <- lm(y ~ -1+rx+ry+cx+cy)
        cc <- coef(summary(lsm))[,1:2]
        zcb[i:(i+1),1] <- cc[1:2,1]
        zcb[i:(i+1),3] <- cc[1:2,2]
        if (nrow(cc)==4) {
          zcb[i:(i+1),2] <- cc[3:4,1]
          zcb[i:(i+1),4] <- cc[3:4,2]
        }
        lsm <- summary(lsm)
        sumstats <- rbind(sumstats,c(zlist$n[i],zlist$m[i],
                                     hypot(zcb[i:(i+1),1]), hypot(zcb[i:(i+1),2]),
                                     lsm$sigma, lsm$r.squared,lsm$fstatistic))
        i <- i+2
      }
  }
  zcb[is.na(zcb)] <- 0
  colnames(zcb) <- c("zc_mirror", "zc_inst", "se_zc_mirror", "se_zc_inst")
  colnames(sumstats)[1:7] <- c("n","m","rms_mirror", "rms_inst", "sigma","R2", "F")
  wf.mirror <- pupil.arb(zcoef=zcb[,1], zlist=zlist, nrow=nrow, ncol=ncol, cp=cp)
  wf.inst <- pupil.arb(zcoef=zcb[,2], zlist=zlist, nrow=nrow, ncol=ncol, cp=cp)
  list(zcb=data.frame(zcb),sumstats=data.frame(sumstats),
       wf.mirror=wf.mirror, wf.inst=wf.inst)
                        }
                        
