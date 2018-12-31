## Difference that maintains original array size by appending a row

.fdiff <- function(X) rbind(diff(X), rep(NA,ncol(X)))

## Ghiglia and Pritt's wrap operator. Thanks to Steve K.

wrap <- function(phase) (phase+pi)%%(2*pi) - pi


## Calculate residue map
## A "positive" charged residue has value +1
## A "negative" charge has value -1

rmap <- function(phase, dx=NULL, dy=NULL, plot=FALSE, ...) {
  nr <- nrow(phase)
  nc <- ncol(phase)

  if (is.null(dx)) dx <- wrap(.fdiff(phase))
  if (is.null(dy)) dy <- wrap(t(.fdiff(t(phase))))
  d2y <- rbind(dy[1:(nr-1),]-dy[2:nr,],rep(NA,nc))
  d2x <- cbind(dx[,2:nc]-dx[,1:(nc-1)],rep(NA,nr))
  residues <- round((d2y+d2x)/(2*pi))
  if (plot) {
  	## positive charges
	chp <- which(residues ==1, arr.ind=TRUE)
	## negative charges
	chm <- which(residues == -1, arr.ind=TRUE)
	image(1:nr, 1:nc, phase, col=grey256, asp=1, xlab="X", ylab="Y", useRaster=TRUE, ...)
	if(length(chp)>0) points(chp, col="green", pch=20)
	if(length(chm)>0) points(chm, col="red", pch=1)
	nr <- sum(abs(residues),na.rm=TRUE)
	return(nr)
  }
  else
  	return(residues)
}

## "Itoh's" method for phase unwrapping.

## integrate the wrapped phase differences
## This is called by brcutpuw, or it can be called
## directly by the user.

## Note: as of 15 May 2008 the unwrapper is written in C.
## It's fastest to bypass the computation of dx and dy in R
## and unwrap with a call to id_uw.
## netflowpuw (in lppuw) passes values of dx and dy with ucall==TRUE.
## brcutpuw passes dx and dy with ucall==FALSE.
## The call to id_dxy_uw is slower by about 80%.

idiffpuw <- function(phase, mask=phase, ucall=TRUE, dx=NULL, dy=NULL) {

    nr <- nrow(phase)
    nc <- ncol(phase)
	
    if (ucall) {
        if (is.null(dx)) {
            b <- .C(id_uw, as.integer(nr), as.integer(nc), as.double(phase/(2*pi)),
	           as.double(NA*phase), NAOK=TRUE)
	           puw <- matrix(b[[4]], nr, nc)
        } else {
            b <- .C(id_dxy_uw, as.integer(nr), as.integer(nc), as.double(phase/(2*pi)),
	         as.double(mask), as.double(dx/(2*pi)), as.double(dy/(2*pi)),
	         as.double(NA*phase), integer(nr*nc), 
	         NAOK=TRUE)
            puw <- matrix(b[[7]], nr, nc)
        }
        class(puw) <- "pupil"
        return(puw)
    } else {
        if (is.null(dx)) dx <- wrap(.fdiff(phase))
        if (is.null(dy)) dy <- wrap(t(.fdiff(t(phase))))
        b <- .C(id_dxy_uw, as.integer(nr), as.integer(nc), as.double(phase/(2*pi)), 
	      as.double(mask), as.double(dx/(2*pi)), as.double(dy/(2*pi)), 
	      as.double(NA*phase), integer(nr*nc), 
	      NAOK=TRUE)
        puw <- matrix(b[[7]], nr, nc)
        uw <- matrix(as.logical(b[[8]]), nr, nc)
        class(puw) <- "pupil"
        return(list(puw=2*pi*puw, uw=uw))
    }
}


## Quality guided unwrapper.
## This is mostly a wrapper for a call to the C function "q_uw".

qpuw <- function(phase, qual) {
    nr <- nrow(phase)
    nc <- ncol(phase)
	
    phase <- phase/(2*pi)
    qual[is.na(phase)] <- 0
	
    sol <- .C(q_uw, as.integer(nr), as.integer(nc), as.double(phase), as.double(qual),
		as.double(NA*phase), NAOK=TRUE)
    puw <- matrix(sol[[5]], nr, nc)
    class(puw) <- "pupil"
    puw
}

##
##  Branch cut algorithm.
## This now solves a variant of the assignment problem to minimize
## the total length of all branch cuts.
##
## note: parameter pen is a penalty added to distances from charges to edges.
##         I'm not sure this is really useful. Seems to work fine with the default
##         of 0.
##

    
brcutpuw <- function(phase, pen=0, details=FALSE) {

    require(zernike)
    require(clue)

    ## distance between points specified by their x,y coordinates.
    
    dmat <- function(p1, p2) {
        fn <- function(x,y) abs(y-x)
        xleg <- outer(p1[,1], p2[,1], FUN=fn)
        yleg <- outer(p1[,2], p2[,2], FUN=fn)
        xleg+yleg
    }
    
    ## Make a branch cut
    
    brcut <- function(startp, endp, mask) {
        xl <- endp[1]-startp[1]
        yl <- endp[2]-startp[2]
        if (abs(yl) <= abs(xl)) {
            slope <- yl/xl
            ix <- startp[1]:endp[1]
            iy <- startp[2] + round(slope*(ix-startp[1]))
        } else {
            slope <- xl/yl
            iy <- startp[2]:endp[2]
            ix <- startp[1] + round(slope*(iy-startp[2]))
        }
        mask[cbind(ix,iy)] <- NA
        mask
    }
    
    res <- rmap(phase)
    ## no residues, so a single call to idiffpuw is all we need
    if (sum(abs(res),na.rm=TRUE) == 0) {
        return(idiffpuw(phase,ucall=TRUE))
    }
    

    nr <- nrow(phase)
    nc <- ncol(phase)
    dx <- wrap(.fdiff(phase))
    dy <- wrap(t(.fdiff(t(phase))))
    
    mask <- matrix(1,nr,nc)
    mask[is.na(phase)] <- NA
    
    ## positive charges
    chp <- which(res==1, arr.ind=TRUE)
    ncp <- nrow(chp)
    ## negative charges
    chm <- which(res== -1, arr.ind=TRUE)
    ncm <- nrow(chm)
    ## get the boundary points - ones one pixel outside the area of defined phase.
    ## step in
    ptsb <- which(is.na(phase) & (!is.na(rbind(phase[-1,], rep(NA, nc)))), arr.ind=TRUE)
    ptsb <- rbind(ptsb, which(is.na(phase) & (!is.na(cbind(phase[,-1], rep(NA, nr)))), arr.ind=TRUE))
    ## step out
    ptsc <- which((!is.na(phase)) & is.na(rbind(phase[-1,], rep(NA, nc))), arr.ind=TRUE)
    ptsc[,1] <- ptsc[,1]+1
    ptsb <- rbind(ptsb, ptsc)
    ptsc <- which((!is.na(phase)) & is.na(cbind(phase[,-1], rep(NA, nr))), arr.ind=TRUE)
    ptsc[,2] <- ptsc[,2]+1
    ptsb <- rbind(ptsb, ptsc)
    
    
    if ((ncp>0) & (ncm>0)) {
    
      ## distances between residues
      dpm <- dmat(chp, chm)
      
      ## distance from each positive charge to nearest edge
      dce <- dmat(chp, ptsb)
      ipb <- apply(dce, 1, which.min)
      dpe <- apply(dce, 1, min)+pen
      pec <- matrix(ptsb[ipb,], ncp, 2)
      
      ## distance from each negative charge to nearest edge
      dce <- dmat(chm, ptsb)
      ipb <- apply(dce, 1, which.min)
      dme <- apply(dce, 1, min)+pen
      mec <- matrix(ptsb[ipb,], ncm, 2)
      
      ## cost matrix of assigning p to edge and m to edge for each (p,m). This has the same dimension as dpm
      cost.ec <- outer(dpe, dme, "+")
      isd <- (dpm <= cost.ec)
      
      ## elementwise minimum of dpm and cost.ec. 
      ## Idea is if dpm[p,m] < cost.ec[p,m] we should connect dipole, otherwise connect both to edge
      minc <- pmin(dpm, cost.ec)
      pecuts <- NULL
      mecuts <- NULL

      ## use Hungarian algorithm from package clue to get optimal assigments from p to m.
      ## Then figure out which are dipole cuts and which are edge cuts.
      ## If ncm > ncp there have to be some unassigned m's. Those connect to edge.
      
      if (ncm >= ncp) {
        yp <- solve_LSAP(minc)
        ass <- cbind(seq_along(yp), as.integer(as.vector(yp)))
        mecuts <- setdiff(1:ncm, ass[,2])
      } else {
        yp <- solve_LSAP(t(minc))
        ass <- cbind(as.integer(as.vector(yp)), seq_along(yp))
        pecuts <- setdiff(1:ncp, ass[,1])
      }
      dpcuts <- ass[which(isd[ass]),]
      edgecuts <- ass[which(!isd[ass]),]
      pecuts <- c(pecuts, edgecuts[,1])
      mecuts <- c(mecuts, edgecuts[,2])
      cost.dpcuts <- dpm[dpcuts]
      cost.pe <- dpe[pecuts]
      cost.me <- dme[mecuts]
      
      if (details) {
        cat(paste(nrow(dpcuts),"dipole connections, total cost =",sum(cost.dpcuts),"\n"))
        cat(paste(length(pecuts),"+ to edge, total cost =",sum(cost.pe),"\n"))
        cat(paste(length(mecuts),"- to edge, total cost =",sum(cost.me),"\n"))
      }
    
      cutlen <- sum(cost.dpcuts)+sum(cost.pe)+sum(cost.me)
      
      for (i in seq_along(dpcuts[,1])) {
          startp <- chp[dpcuts[i, 1],]
          endp <- chm[dpcuts[i, 2],]
          mask <- brcut(startp, endp, mask)
      }
      for (i in seq_along(pecuts)) {
          startp <- chp[pecuts[i],]
          endp <- pec[pecuts[i],]
          mask <- brcut(startp, endp, mask)
      }
      for (i in seq_along(mecuts)) {
          startp <- chm[mecuts[i],]
          endp <- mec[mecuts[i],]
          mask <- brcut(startp, endp, mask)
      }
    }
    ## if there are only plus or minus charges connect to nearest edge (should be rare case, but...)
    else if (ncp > 0) {
      ## distance from each positive charge to nearest edge
      dce <- dmat(chp, ptsb)
      ipb <- apply(dce, 1, which.min)
      dpe <- apply(dce, 1, min)+pen
      pe <- matrix(ptsb[ipb,],ncp,2)
      for (i in 1:ncp) {
          startp <- chp[i,]
          endp <- pe[i,]
          mask <- brcut(startp, endp, mask)
      }
    }
    else if (ncm > 0) {    
      ## distance from each negative charge to nearest edge
      dce <- dmat(chm, ptsb)
      ipb <- apply(dce, 1, which.min)
      dme <- apply(dce, 1, min)+pen
      me <- matrix(ptsb[ipb,],ncm,2)
      for (i in 1:ncm) {
          startp <- chm[i,]
          endp <- me[i,]
          mask <- brcut(startp, endp, mask)
      }
    }
    
    ## call idiffpuw to unwrap all unmasked pixels
    uwum <- idiffpuw(phase, mask, ucall=FALSE, dx=dx, dy=dy)
    puw <- uwum$puw
    uw <- uwum$uw
    uw[is.na(phase)] <- NA
    ## fill in whatever was left wrapped
    ## This inverts the logic of fill algorithm in idiffpuw
    ## todo initially contains a list of pixels that haven't been unwrapped
    ## that should be. We look for neighbors that have been unwrapped
    ## and unwrap from them, removing the newly unwrapped pixels from
    ## the todo list.
    ## This should work because branch cuts are just a pixel wide and should
    ## adjoin unmasked regions. It's possible that cuts could completely
    ## enclose an area, in which case it will still be filled
    ## with most likely erroneous values as we fill through the cuts.
    
    todo <- which(!uw)
    kn <- c(-1, 1, -nr, nr)
    while (length(todo) > 0) {
        for (i in 1:4) {
            b <- todo + kn[i]
            subs <- which(uw[b])
            puw[todo[subs]] <- puw[b[subs]] + switch(i, dx[b[subs]], -dx[todo[subs]],
             dy[b[subs]], -dy[todo[subs]])
            uw[todo[subs]] <- TRUE
            todo <- todo[-subs]
            if (length(todo)==0) break
        }
    }
    if (details) {
        bcuts <- matrix(NA, nr, nc)
        bcuts[is.na(mask) & (!is.na(phase))] <- 1
        list(puw=puw/(2*pi), cutlen=cutlen, bcuts=bcuts,
                    dpm=dpm, dpe=dpe, dme=dme,
                    dpcuts=dpcuts, pecuts=pecuts, mecuts=mecuts,
                    cost.dpcuts=cost.dpcuts, cost.pe=cost.pe, cost.me=cost.me
        )
    } else
        puw/(2*pi)
}
