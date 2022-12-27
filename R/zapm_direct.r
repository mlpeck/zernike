############
#
# Radial Zernike annular polynomials from published formulas
#  
############ 

#' Radial Zernike annular polynomials from formulas
#'
#' The first 15 radial Zernike annular polynomials copied from
#' Table 1 of Mahajan (1994)
#'
#' @param rho a vector of radial coordinates.
#' @param eps the obstruction fraction 0 <= eps < 1.
#' @param j index of the value to fetch
#'
#' @return A vector of radial Zernike annular polynomial values at the coordinates `rho` for index `j`
#'
#' @details This function is included for testing and reference only.
#'
#' @references
#' Mahajan, V. N. (1994). Zernike Annular Polynomials and Optical Aberrations
#'   of Systems with Annular Pupils. *Supplement to Applied Optics* **33**, 8125-8127.
#'
#' @seealso This function is called by [zapm_direct()] and [zapm_iso_direct()].
#'   Complete sequences to arbitrary order can be calculated with [rzernike_ann()]
#'   and [rzernike_ann_128()].
#'
#' @md
rzernike_ann_direct <- function(rho, eps, j) {
  switch(j+1,
    rep(1, length(rho)),                                                                                ## n=0, m=0 (piston)
    rho/sqrt(1+eps^2),                                                                                  ## n=1, m=1 (tilt)
    rho^2/sqrt(1+eps^2+eps^4),                                                                          ## n=2, m=2 (astigmatism)
    2*rho^2/(1-eps^2) - (1+eps^2)/(1-eps^2),                                                            ## n=2, m=0 (defocus)
    rho^3/sqrt(1+eps^2+eps^4+eps^6),                                                                    ## n=3, m=3 (trefoil)
    (3*(1+eps^2)*rho^3-2*(1+eps^2+eps^4)*rho)/(1-eps^2)/sqrt((1+eps^2)*(1+4*eps^2+eps^4)),              ## n=3, m=1 (coma)
    rho^4/sqrt(1+eps^2+eps^4+eps^6+eps^8),                                                              ## n=4, m=4 (quadratic)
    (4*rho^4-3*((1-eps^8)/(1-eps^6))*rho^2)*sqrt(1-eps^2)/sqrt(16*(1-eps^10)-15*(1-eps^8)^2/(1-eps^6)), ## n=4, m=2 (secondary astigmatism)
    (6*rho^4-6*(1+eps^2)*rho^2+1+4*eps^2+eps^4)/(1-eps^2)^2,                                            ## n=4, m=0 (spherical aberration)
    rho^5/sqrt(1+eps^2+eps^4+eps^6+eps^8+eps^10),                                                       ## n=5, m=5 (5-fold)
    (5*rho^5-4*((1-eps^10)/(1-eps^8))*rho^3)/sqrt((25*(1-eps^12)-24*(1-eps^10)^2/(1-eps^8))/(1-eps^2)), ## n=5, m=3 (secondary trefoil)
    (10*(1+4*eps^2+eps^4)*rho^5-12*(1+4*eps^2+4*eps^4+eps^6)*rho^3+
        3*(1+4*eps^2+10*eps^4+4*eps^6+eps^8)*rho)/(1-eps^2)^2/
        sqrt((1+4*eps^2+eps^4)*(1+9*eps^2+9*eps^4+eps^6)),                                              ## n=5, m=1 (secodary coma)
    rho^6/sqrt(1+eps^2+eps^4+eps^6+eps^8+eps^10+eps^12),                                                ## n=6, m=6 (6-fold)
    (6*rho^6-5*((1-eps^12)/(1-eps^10))*rho^4)/
        sqrt((36*(1-eps^14)-35*(1-eps^12)^2/(1-eps^10))/(1-eps^2)),                                     ## n=6, m=4 (is there a name for this?)
    (15*(1+4*eps^2+10*eps^4+4*eps^6+eps^8)*rho^6-20*(1+4*eps^2+10*eps^4+10*eps^6+4*eps^8+eps^10)*rho^4+
        6*(1+4*eps^2+10*eps^4+20*eps^6+10*eps^8+4*eps^10+eps^12)*rho^2)/
        (1-eps^2)^2/
        sqrt((1+4*eps^2+10*eps^4+4*eps^6+eps^8)*(1+9*eps^2+45*eps^4+65*eps^6+45*eps^8+9*eps^10+eps^12)), ## n=6, m=2 (tertiary astigmatism)
    (20*rho^6-30*(1+eps^2)*rho^4+12*(1+3*eps^2+eps^4)*rho^2-(1+9*eps^2+9*eps^4+eps^6))/(1-eps^2)^3       ## n=6, m=0 (secondary spherical aberration)
  )
}

############
#
# matrix of zernike annular polynomials through order 6
#
############

#' Zernike Annular polynomials from formulas
#'
#' Create a matrix of Zernike Annular polynomial values
#'   for the counterparts of primary and secondary
#'   optical aberrations.
#'
#' @param rho a vector of radial coordinates with eps <= rho <= 1.
#' @param theta a vector of angular coordinates, in radians.
#' @param eps the obstruction fraction 0 <= eps < 1.
#'
#' @return a `length(rho) x 16` matrix of Zernike annular polynomial values.
#'
#' @details The values are from published formulas. This function is
#'    included for testing and reference only.
#'
#' @seealso Calls [rzernike_ann_direct()] for radial Zernike annular values. 
#'    Use [zapm()] or [zapm_128()] for complete sequences
#'    of values to arbitrary order.
#'
#' @md
zapm_direct <- function(rho, theta, eps) {
  if (length(rho) != length(theta)) {
    stop("coordinates have unequal length")
  }
  zm <- matrix(0, length(rho), 16)
  jvals <- c(0, 1, 3, 2, 5, 8, 4, 7, 11, 15)
  j <- 1
  k <- 1
  for (order in seq(0, 6, by=2)) {
    for (n in (order/2):order) {
      znorm <- sqrt(n+1)
      m <- order - n
      cat(paste("n =", n, "m =", m, "order =", order, "\n"))
      cat(paste("j =", j, "k =", k, "\n"))
      rz <- znorm * rzernike_ann_direct(rho, eps, jvals[j])
      if (m > 0) {
        zm[, k] <- sqrt(2) * rz * cos(m * theta)
        k <- k + 1
        zm[, k] <- sqrt(2) * rz * sin(m * theta)
        k <- k + 1
      } else {
        zm[, k] <- rz
        k <- k + 1
      }
      j <- j+1
    }
  }
  colnames(zm) <- paste("AZ", 0:15, sep="")
  zm
}

############
#
# matrix of zernike annular polynomials through order 6, ISO sequence
#
############


#' Zernike Annular polynomials from formulas
#'
#' Create a matrix of Zernike Annular polynomial values
#'   complete through radial and azimuthal orders 6
#'   in ISO/ANSI sequence.
#'
#' @param rho a vector of radial coordinates with eps <= rho <= 1.
#' @param theta a vector of angular coordinates, in radians.
#' @param eps the obstruction fraction 0 <= eps < 1.
#'
#' @return a `length(rho) x 28` matrix of Zernike annular polynomial values.
#'
#' @details The values are from published formulas. This function is
#'    included for testing and reference only.
#'
#' @seealso Calls [rzernike_ann_direct()] for radial Zernike annular values. 
#'    Use [zapm_iso()] or [zapm_iso_128()] for complete sequences
#'    of values to arbitrary order.
#'
#' @examples
#' sample_az_iso_direct <- function(eps=0.33, col=rev(zernike::rygcb(400)), addContours=TRUE, cscale=TRUE) {
#' 
#'  ## get coordinates for obstructed aperture
#'  cp <- cp.default
#'  cp$obstruct <- eps
#'  prt <- pupil.rhotheta(nrow.default,ncol.default,cp=cp)
#'  rho <- prt$rho[!is.na(prt$rho)]
#'  theta <- prt$theta[!is.na(prt$theta)]
#'  
#'  ## fill up matrixes of Annular Zernikes using direct and recursive formulas
#'  
#'  zam_dir <- zapm_iso_direct(rho, theta, eps=eps)
#'  zam_iso <- zapm_iso(rho, theta, eps=eps, maxorder=6)
#'  
#'  ## pick a column at random and look up its index pair
#'  
#'  zlist <- makezlist.iso(6)
#'  i <- sample(2:ncol(zam_iso), 1)
#'  n <- zlist$n[i]
#'  m <- zlist$m[i]
#'  
#'  ## fill up the wavefront representations and plot them
#'  
#'  wf_dir <- prt$rho
#'  wf_dir[!is.na(wf_dir)] <- zam_dir[,i]
#'  class(wf_dir) <- "pupil"
#'  
#'  wf_zapm <- wf_dir
#'  wf_zapm[!is.na(wf_zapm)] <- zam_iso[,i]
#'  class(wf_zapm) <- "pupil"
#'  
#'  plot(wf_dir, cp=cp, col=col, addContours=addContours, cscale=cscale)
#'  mtext(paste("Annular Zernike from formula, n =", n, " m =", m))
#'  
#'  x11()
#'  plot(wf_zapm, cp=cp, col=col, addContours=addContours, cscale=cscale)
#'  mtext(paste("Annular Zernike from `zapm_iso`, n =", n, " m =", m))
#'  
#'  cat(paste("Sample wavefronts: SD of diffs", format(sd(wf_dir-wf_zapm, na.rm=TRUE), digits=6), "; range of diffs", 
#'      format(min(wf_zapm-wf_dir, na.rm=TRUE), digits=6),":",
#'      format(max(wf_zapm-wf_dir, na.rm=TRUE), digits=6),"\n\n"))
#'  
#'  cat(paste("Zernike matrixes: SD of diffs", format(sd(zam_dir-zam_iso), digits=6), "; range of diffs", 
#'  format(min(zam_dir-zam_iso), digits=6),":",
#'  format(max(zam_dir-zam_iso), digits=6),"\n"))
#'  
#'  ## return Zernike matrices and wavefronts invisibly
#'  ## just in case user wants to do something with them
#'  
#'  invisible(list(zam_dir=zam_dir, zam_iso=zam_iso, wf_dir=wf_dir, wf_zapm=wf_zapm))
#' }
#'
#' sample_az_iso_direct()
#'
#' @md
zapm_iso_direct <- function(rho, theta, eps) {
  if (length(rho) != length(theta)) {
    stop("coordinates have unequal length")
  }
  ncol <- 28
  zm <- matrix(0, length(rho), ncol)
  i <- 1
  j0 <- 0
  for (n in 0:6) {
    j0 <- j0 + (n+2) %/% 2
    jbase <- j0 - (n+1) %% 2
    for (m in 0:n) {
      mp <- 2*m - n
      j = jbase - (abs(mp)+1) %/% 2
      znorm <- sqrt(n+1)
      zm[, i] <- znorm * rzernike_ann_direct(rho, eps, j)
      if (mp < 0) {
        zm[, i] <- -sqrt(2) * zm[, i] * sin(mp * theta)
      } else if (mp > 0) {
        zm[, i] <- sqrt(2) * zm[, i] * cos(mp * theta)
      }
      i <- i + 1
    }
  }
  colnames(zm) <- paste("AZ", 0:(ncol-1), sep="")
  zm
}    
