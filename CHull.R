CHull <- function(loglik, nrpar, nsclust){
  # browser()
  Clus <- seq(nsclust[1], nsclust[2])
  fitMat <- matrix(data = c(Clus, loglik, nrpar), ncol = 3)
  k <- nrow(fitMat)
  
  screeratios <- matrix(NA, nsclust[2] - nsclust[1] + 1, 1)
  CHullcheck <- 0
  for(nclust in (nsclust[1] + 1):(nsclust[2] - 1)){
    LL_nclust <- fitMat[nclust-nsclust[1]+1, 2]
    npar_nclust <- fitMat[nclust-nsclust[1]+1, 3]
    LL_nclustmin1 <- fitMat[nclust-nsclust[1], 2]
    npar_nclustmin1 <- fitMat[nclust-nsclust[1], 3]
    LL_nclustplus1 <- fitMat[nclust-nsclust[1]+2, 2]
    npar_nclustplus1 <- fitMat[nclust-nsclust[1]+2, 3]
    # determine whether intermediate point is part of the convex hull
    slope <- (LL_nclustplus1-LL_nclustmin1)/(npar_nclustplus1-npar_nclustmin1)
    point_line <- LL_nclustmin1+slope*(npar_nclust-npar_nclustmin1)
    #diff_fit=(LL_nclust-LL_nclustmin1)/abs(LL_nclust)
    if(isTRUE(LL_nclust>=(point_line-.01))){ # && diff_fit>.0001)
      screeratios[nclust - nsclust[1] + 1] <- ((LL_nclust-LL_nclustmin1)/(npar_nclust-npar_nclustmin1))/((LL_nclustplus1-LL_nclust)/(npar_nclustplus1-npar_nclust))
    }
  }
  #screeratios[CHullcheck<0]=NA
  convexhull <- which(!is.na(screeratios))
  nrows <- nrow(fitMat)
  convexhull <- c(fitMat[1,1],convexhull,fitMat[nrows,1])
  nrhull <- length(convexhull)
  change <- 0
  if(nrhull < nrows){
    change <- 1
  }
  while(nrhull > 2 && change == 1){ # check again whether intermediate points are on the convex hull
    nsclusthull <- fitMat[convexhull, 1]
    change <- 0
    for(nclust in 2:(nrhull - 1)){
      if(!identical(convexhull[(nclust-1):(nclust+1)],c(convexhull[nclust]-1,convexhull[nclust],convexhull[nclust]+1))){
        LL_nclust <- fitMat[convexhull[nclust]-nsclust[1]+1,2]
        npar_nclust <- fitMat[convexhull[nclust]-nsclust[1]+1,3]
        LL_nclustmin1 <- fitMat[convexhull[nclust-1]-nsclust[1]+1,2]
        npar_nclustmin1 <- fitMat[convexhull[nclust-1]-nsclust[1]+1,3]
        LL_nclustplus1 <- fitMat[convexhull[nclust+1]-nsclust[1]+1,2]
        npar_nclustplus1 <- fitMat[convexhull[nclust+1]-nsclust[1]+1,3]
        # determine whether intermediate point is part of the convex hull
        slope <- (LL_nclustplus1-LL_nclustmin1)/(npar_nclustplus1-npar_nclustmin1)
        point_line <- LL_nclustmin1+slope*(npar_nclust-npar_nclustmin1)
        if(LL_nclust>=(point_line-.01)){
          # when the subset of three points spans across a point not on the hull, this is the corrected scree ratio (comparing with the previous and next point ON THE HULL)
          screeratios[convexhull[nclust]-nsclust[1]+1]=((LL_nclust-LL_nclustmin1)/(npar_nclust-npar_nclustmin1))/((LL_nclustplus1-LL_nclust)/(npar_nclustplus1-npar_nclust))
        } else {
          screeratios[convexhull[nclust]-nsclust[1]+1]=NA
          change=1
        }
      }
    }
    convexhull <- which(!is.na(screeratios))
    convexhull <- c(fitMat[1,1],convexhull,fitMat[nrows,1])
    nrhull <- length(convexhull)
  }
  fitMat <- cbind(fitMat, abs(screeratios))
  return(fitMat)
}
