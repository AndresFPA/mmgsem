EStep <- function(pi_ks, ngroup, nclus, loglik){

  max_g <-rep(0,ngroup)
  z_gks <- matrix(NA,nrow = ngroup,ncol = nclus)
  
  for(g in 1:ngroup){
    for(k in 1:nclus){
      z_gks[g,k] <- log(pi_ks[k])+loglik[g,k]
    }
    max_g[g] <- max(z_gks[g,]) # prevent arithmetic underflow 
    z_gks[g,] <- exp(z_gks[g,]-rep(max_g[g],nclus))
  }

  # divide by the rowwise sum of the above calculated part 
  z_gks <- diag(1/apply(z_gks,1,sum))%*%z_gks
  # z_gks <- round(z_gks, digits = 16)
  # z_gks <- diag(1/apply(z_gks,1,sum))%*%z_gks
  
  return(z_gks)
}
