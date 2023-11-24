#' Model selection wrapper for MMGSEM
#'
#' A wrapper to perform model selection of MMGSEM. It is used to answer the common clustering question: How many clusters should we use?
#' It calls the main function MMGSEM several times. One for each required model. Most of its arguments are the same as the MMGSEM main function.
#'
#' INPUT: Arguments required by the function
#' @param clusters Vector of length two, indicating the minimal and maximal number of clusters. Used for the model selection.
#' @param dat Observed data of interest for the MMGSEM model.
#' @param step1model Measurement model (MM). Used in step 1. Must be a string (like in lavaan). 
#'                   Can be a list of strings determing the number of measurement blocks (e.g., one string for the MM of 
#'                   factor 1, and a second string for the MM of factor 2)  
#' @param step2model Structural model. Used in step 2. Must be a string (like in lavaan).
#' @param group Name of the group variable. Must be a string.
#' @param seed Pre-defined seed for the random start in case a replication is needed in the future.
#' @param userStart Pre-defined start provided by the user. Must be a matrix of dimensions ngroups*nclus.
#'                  A 1 represents that the row group is in the column cluster. There must be only one 1
#'                  in each row. The function will return results only for this start.
#'                  Example for 6 groups and 2 clusters:
#'
#'                                   [,1] [,2]
#'                             [1,]    1    0
#'                             [2,]    1    0
#'                             [3,]    1    0
#'                             [4,]    0    1
#'                             [5,]    0    1
#'                             [6,]    0    1
#'
#' @param s1out Resulting lavaan object from a cfa() analysis. Can be used to directly input the results from the step 1
#'              if the user only wants to use MMGSEM() to estimate the structural model (Step 2). If not NULL, MMGSEM()
#'              will skip the estimation of Step 1 and use the s1out object as the input for Step 2.
#' @param max_it Maximum number of iterations for each start.
#' @param nstarts Number of random starts that will be performed (in an attempt to avoid local maxima).
#' @param partition Type of random partition used to start the EM algorithm. Can be "hard" and "soft".
#' @param constraints String containing the constrains of the model. Receives the same input as group.equal in lavaan.
#' @param NonInv String containing the non-invariances of the model. Similar to group.partial in lavaan.
#' @param Endo2Cov TRUE or FALSE argument to determine whether to allow or not covariance between endogenous 2 variables.
#'                 If TRUE (the default), the covariance between endogenous factors is allowed.
#' @param allG TRUE or FALSE. Determines whether the endogenous covariances are group (T) or cluster (F) specific.
#'             By default, is TRUE.
#' @param fit String argument. Either "factors" or "observed". Determines which loglikelihood will be used to fit the model.
#'
#' OUTPUT: The function will return a list with the following results:
#' @return z_gks: Posterior classification probabilities or cluster memberships (ngroups*nclus matrix).
#' @return final_fit: Lavaan fit of the best and final model (not useful, only for testing purposes).
#' @return psi_gks: Resulting psi matrices (i.e., residual factor covariances) for each group-cluster combination.
#' @return lambda: Resulting invariant lambda matrix (i.e., factor loadings) for all groups.
#' @return theta_gs: Resulting theta matrices (i.e., items' residual/unique variances) for each group.
#' @return beta_ks: Resulting beta matrices (i.e., regression coefficients) for each cluster.
#' @return loglikelihood: Total loglikelihood of the best model.
#' @return logliks_starts: Total loglikelihood of the models of each run (allows to check sensitivity to local maxima).
#' @return Obs_loglik: Reconstruction of the observed loglikelihood after the model has converged. Only valid when
#'                     fit is "factors".
#'
#' PLEASE NOTE: This function requires 'lavaan' package to work. It also requires the function 'EStep.R'.
#'

ModelSelection <- function(dat, step1model = NULL, step2model = NULL,
                           group, clusters, seed = NULL, userStart = NULL, s1out = NULL,
                           max_it = 10000L, nstarts = 1L, printing = FALSE,
                           partition = "hard", NonInv = NULL, constraints = "loadings",
                           Endo2Cov = TRUE, allG = TRUE, fit = "factors", est_method = "local") {
  
  # Prepare objects to compare in model selection
  nmodels   <- length(clusters[1]:clusters[2])
  nclus     <- clusters[1]:clusters[2]
  model_fit <- vector(mode = "list", length = nmodels)
  BIC_G     <- vector(mode = "list", length = nmodels)
  BIC_N     <- vector(mode = "list", length = nmodels)
  BIC_G_fac <- vector(mode = "list", length = nmodels)
  BIC_N_fac <- vector(mode = "list", length = nmodels)
  AIC3      <- vector(mode = "list", length = nmodels)
  AIC3_fac  <- vector(mode = "list", length = nmodels)
  AIC       <- vector(mode = "list", length = nmodels)
  AIC_fac   <- vector(mode = "list", length = nmodels)
  ICL       <- vector(mode = "list", length = nmodels)
  ICL_fac   <- vector(mode = "list", length = nmodels)
  LL        <- vector(mode = "list", length = nmodels)
  LL_fac    <- vector(mode = "list", length = nmodels)
  nrpar     <- vector(mode = "list", length = nmodels)
  nrpar_fac <- vector(mode = "list", length = nmodels)
  
  
  # Call MMGSEM using the arguments provided by the user k times (one per required model)
  for(k in clusters[1]:clusters[2]){
    # If the user provides output of the first step (s1out), use it for all models
    if(!is.null(s1out)){
      model_fit[[k]] <- MMGSEM(dat = dat, step1model = step1model, step2model = step2model, 
                               group = group, nclus = k, seed = seed, userStart = userStart, 
                               s1out = s1out, max_it = max_it, nstarts = nstarts, printing = printing, 
                               partition = partition, NonInv = NonInv, constraints = constraints, 
                               Endo2Cov = Endo2Cov, allG = allG, fit = fit, est_method = est_method)
    } else if (is.null(s1out)){
      # browser()
      if(k == clusters[1]){
        model_fit[[k]] <- MMGSEM(dat = dat, step1model = step1model, step2model = step2model, 
                                 group = group, nclus = k, seed = seed, userStart = userStart, 
                                 s1out = NULL, max_it = max_it, nstarts = nstarts, printing = printing, 
                                 partition = partition, NonInv = NonInv, constraints = constraints, 
                                 Endo2Cov = Endo2Cov, allG = allG, fit = fit, est_method = est_method)
        
        # Save the covariance matrix of step 1, so we do not run it for all models (it is the same for all of them).
        s1output <- model_fit[[k]]$MM
      } else {
        # browser()
        model_fit[[k]] <- MMGSEM(dat = dat, step1model = step1model, step2model = step2model, 
                                 group = group, nclus = k, seed = seed, userStart = userStart, 
                                 s1out = s1output, max_it = max_it, nstarts = nstarts, printing = printing, 
                                 partition = partition, NonInv = NonInv, constraints = constraints, 
                                 Endo2Cov = Endo2Cov, allG = allG, fit = fit, est_method = est_method)
      }
      print(paste("model", k, "finished"))
    }
    # Save the model selection criteria
    BIC_G[[k]]     <- model_fit[[k]]$model_sel$BIC$observed$BIC_G
    BIC_N[[k]]     <- model_fit[[k]]$model_sel$BIC$observed$BIC_N
    BIC_G_fac[[k]] <- model_fit[[k]]$model_sel$BIC$Factors$BIC_G
    BIC_N_fac[[k]] <- model_fit[[k]]$model_sel$BIC$Factors$BIC_N
    AIC[[k]]       <- model_fit[[k]]$model_sel$AIC$observed
    AIC_fac[[k]]   <- model_fit[[k]]$model_sel$AIC$Factors
    AIC3[[k]]      <- model_fit[[k]]$model_sel$AIC3$observed
    AIC3_fac[[k]]  <- model_fit[[k]]$model_sel$AIC3$Factors
    ICL[[k]]       <- model_fit[[k]]$model_sel$ICL$observed
    ICL_fac[[k]]   <- model_fit[[k]]$model_sel$ICL$Factors
    LL[[k]]        <- model_fit[[k]]$logLik$obs_loglik
    LL_fac[[k]]    <- model_fit[[k]]$logLik$loglik
    nrpar[[k]]     <- model_fit[[k]]$NrPar$Obs.nrpar
    nrpar_fac[[k]] <- model_fit[[k]]$NrPar$Fac.nrpar
    
  } # For loop ends here

  # Also do CHull
  Chull_res      <- CHull(loglik = unlist(LL), nrpar = unlist(nrpar), nsclust = clusters)
  Chull_res_fac  <- CHull(loglik = unlist(LL_fac), nrpar = unlist(nrpar_fac), nsclust = clusters)
  
  # Chull function already returns a matrix with LL and nrpar. Use it as a base for the rest of the results
  # browser()
  overview <- cbind(Chull_res, 
                    unlist(BIC_G), unlist(BIC_N),
                    unlist(AIC),
                    unlist(AIC3),
                    unlist(ICL),
                    Chull_res_fac[, c(2:4)],
                    unlist(BIC_G_fac), unlist(BIC_N_fac),
                    unlist(AIC_fac),
                    unlist(AIC3_fac),
                    unlist(ICL_fac)
                    )
  
  colnames(overview) <- c("Clusters", 
                          "LL", "nrpar", "Chull Scree", "BIC_G", "BIC_N", "AIC", "AIC3", "ICL",
                          "LL_fac", "nrpar_fac", "Chull Scree_fac", "BIC_G_fac", "BIC_N_fac", "AIC_fac", "AIC3_fac", "ICL_fac")
  
  overview <- as.data.frame(overview)
  
  return(list(Overview = overview, 
              Models   = model_fit)
         )
}


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

plot.MMGSEM <- function(ModelSel){
    overview <- ModelSel$Overview
    
    # Plots for visual inspection
    plot.BIC_G     <- ggplot(data = overview, mapping = aes(x = Clusters, y = BIC_G))     + geom_point() + geom_line() # plot(x = overview$Clusters, y = overview$BIC_G); lines(x = overview$Clusters, y = overview$BIC_G)
    plot.BIC_N     <- ggplot(data = overview, mapping = aes(x = Clusters, y = BIC_N))     + geom_point() + geom_line() # plot(x = overview$Clusters, y = overview$BIC_N); lines(x = overview$Clusters, y = overview$BIC_N)
    plot.BIC_G_fac <- ggplot(data = overview, mapping = aes(x = Clusters, y = BIC_G_fac)) + geom_point() + geom_line() # plot(x = overview$Clusters, y = overview$BIC_G_fac); lines(x = overview$Clusters, y = overview$BIC_G_fac)
    plot.BIC_N_fac <- ggplot(data = overview, mapping = aes(x = Clusters, y = BIC_N_fac)) + geom_point() + geom_line() # plot(x = overview$Clusters, y = overview$BIC_N_fac); lines(x = overview$Clusters, y = overview$BIC_N_fac)
    plot.AIC       <- ggplot(data = overview, mapping = aes(x = Clusters, y = AIC))       + geom_point() + geom_line() # plot(x = overview$Clusters, y = overview$AIC); lines(x = overview$Clusters, y = overview$AIC)
    plot.AIC_fac   <- ggplot(data = overview, mapping = aes(x = Clusters, y = AIC_fac))   + geom_point() + geom_line() # plot(x = overview$Clusters, y = overview$AIC_fac); lines(x = overview$Clusters, y = overview$AIC_fac)
    plot.AIC3      <- ggplot(data = overview, mapping = aes(x = Clusters, y = AIC3))      + geom_point() + geom_line() # plot(x = overview$Clusters, y = overview$AIC3); lines(x = overview$Clusters, y = overview$AIC3)
    plot.AIC3_fac  <- ggplot(data = overview, mapping = aes(x = Clusters, y = AIC3_fac))  + geom_point() + geom_line() # plot(x = overview$Clusters, y = overview$AIC3_fac); lines(x = overview$Clusters, y = overview$AIC3_fac)
    plot.ICL       <- ggplot(data = overview, mapping = aes(x = Clusters, y = ICL))       + geom_point() + geom_line() # plot(x = overview$Clusters, y = overview$ICL); lines(x = overview$Clusters, y = overview$ICL)
    plot.ICL_fac   <- ggplot(data = overview, mapping = aes(x = Clusters, y = ICL_fac))   + geom_point() + geom_line() # plot(x = overview$Clusters, y = overview$ICL_fac); lines(x = overview$Clusters, y = overview$ICL_fac)
    plot.Chull     <- ggplot(data = overview, mapping = aes(x = Clusters, y = LL))        + geom_point() + geom_line() # plot(x = overview$Clusters, y = overview$LL); lines(x = overview$Clusters, y = overview$LL)
    plot.Chull_fac <- ggplot(data = overview, mapping = aes(x = Clusters, y = LL_fac))    + geom_point() + geom_line() # plot(x = overview$Clusters, y = overview$LL_fac); lines(x = overview$Clusters, y = overview$LL_fac)
    
    plots <- list(BIC_G = plot.BIC_G, BIC_N = plot.BIC_N, BIC_G_fac = plot.BIC_G_fac,
                  BIC_N_fac = plot.BIC_N_fac, AIC = plot.AIC, AIC_fac = plot.AIC_fac, AIC3 = plot.AIC3, 
                  AIC3_fac = plot.AIC3_fac, ICL = plot.ICL, ICL_fac = plot.ICL_fac, 
                  Chull = plot.Chull, Chull_fac = plot.Chull_fac)
    
    return(plots)
}






