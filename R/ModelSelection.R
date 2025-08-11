#' Model selection wrapper for MMGSEM
#'
#' A wrapper to perform model selection of MMGSEM. It is used to answer the common clustering question: How many clusters should we use?
#' It calls the main function MMGSEM several times. One for each required model. Most of its arguments are the same as the MMGSEM main function.
#'
#' INPUT: Arguments required by the function
#' @param nclus Vector of length two, indicating the minimal and maximal number of clusters. Used for the model selection.
#' @param dat Observed data of interest for the MMGSEM model.
#' @param S1 Measurement model (MM). Used in step 1. Must be a string (like in lavaan).
#'                   Can be a list of strings determining the number of measurement blocks (e.g., one string for the MM of
#'                   factor 1, and a second string for the MM of factor 2)
#' @param S2 Structural model. Used in step 2. Must be a string (like in lavaan).
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
#' @param s1_fit Resulting lavaan object from a cfa() analysis. Can be used to directly input the results from the step 1
#'              if the user only wants to use MMGSEM() to estimate the structural model (Step 2). If not NULL, MMGSEM()
#'              will skip the estimation of Step 1 and use the s1out object as the input for Step 2. If NULL, MMGSEM() will
#'              estimate both step 1 and step 2.
#' @param max_it Maximum number of iterations for each start.
#' @param nstarts Number of random starts that will be performed (in an attempt to avoid local maxima).
#' @param partition Type of random partition used to start the EM algorithm. Can be "hard" and "soft".
#' @param constraints String containing the constrains of the model. Receives the same input as group.equal in lavaan.
#' @param NonInv String containing the non-invariances of the model. Similar to group.partial in lavaan.
#' @param endogenous_cov TRUE or FALSE argument to determine whether to allow or not covariance between purely endogeonus variables.
#'                       If TRUE (the default), the covariance between endogenous factors is allowed.
#' @param endo_group_specific TRUE or FALSE. Determines whether the endogenous covariances are group (T) or cluster (F) specific.
#'                            By default, it is TRUE.
#' @param sam_method either "local" or "global. Follows local and global approaches from the SAM method. GLOBAL NOT FUNCTIONAL YET.
#' @param rescaling Only used when data is ordered. By default, MMGSEM uses the marker variable scaling approach. But identification
#'                  constraints with ordinal data (by default) are handled by standardizing the factors' variance in the first step.
#'                  The rescaling argument (either T or F) rescales the factor variances and loadings to the marker variable scaling
#'                  before running step 2. It is set to T by default (rescaling happens). If set to F, the factor variances are kept fixed to 1.
#' @param ... MMGSEM relies on lavaan for the estimation of the first step (i.e., CFA). If needed, the users can pass any lavaan argument to MMGSEM
#'            and it will be considered when estimating the CFA. For instance, std.lv if users want standardized latent variables,
#'            group.equal for constraints, group.partial for non-invariances, etc.
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
#' PLEASE NOTE: This function requires 'lavaan' package to work.
#' @export
ModelSelection <- function(dat, S1 = NULL, S2 = NULL,
                           group, nclus, seed = NULL,
                           userStart = NULL, s1_fit = NULL,
                           max_it = 10000L, nstarts = 20L, printing = FALSE,
                           partition = "hard", endogenous_cov = TRUE,
                           endo_group_specific = TRUE,
                           sam_method = "local", meanstr = FALSE,
                           rescaling = F, ...) {

  # Prepare objects to compare in model selection
  nmodels   <- length(nclus[1]:nclus[2])
  clusters  <- nclus[1]:nclus[2]
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
  R2entropy <- vector(mode = "list", length = nmodels)


  # Call MMGSEM using the arguments provided by the user k times (one per required model)
  for(k in nclus[1]:nclus[2]){
    # If the user provides output of the first step (s1out), use it for all models
    if(!is.null(s1_fit)){
      model_fit[[k]] <- mmgsem::MMGSEM(dat = dat, S1 = S1, S2 = S2, group = group, nclus = k, seed = seed,
                                       userStart = userStart, s1_fit = s1_fit, max_it = max_it, nstarts = nstarts, printing = printing,
                                       partition = partition, endogenous_cov = endogenous_cov,
                                       endo_group_specific = endo_group_specific, sam_method = sam_method, ...)
    } else if (is.null(s1_fit)){
      # browser()
      if(k == nclus[1]){
        model_fit[[k]] <- mmgsem::MMGSEM(dat = dat, S1 = S1, S2 = S2, group = group, nclus = k, seed = seed,
                                         userStart = userStart, s1_fit = NULL, max_it = max_it, nstarts = nstarts, printing = printing,
                                         partition = partition, endogenous_cov = endogenous_cov,
                                         endo_group_specific = endo_group_specific, sam_method = sam_method, ...)

        # Save the covariance matrix of step 1, so we do not run it for all models (it is the same for all of them).
        s1_fit <- model_fit[[k]]$MM
      } else {
        # browser()
        model_fit[[k]] <- mmgsem::MMGSEM(dat = dat, S1 = S1, S2 = S2, group = group, nclus = k, seed = seed,
                                         userStart = userStart, s1_fit = s1_fit, max_it = max_it, nstarts = nstarts, printing = printing,
                                         partition = partition, endogenous_cov = endogenous_cov,
                                         endo_group_specific = endo_group_specific, sam_method = sam_method, ...)
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
    R2entropy[[k]] <- model_fit[[k]]$model_sel$R2_entropy
    LL[[k]]        <- model_fit[[k]]$logLik$obs_loglik
    LL_fac[[k]]    <- model_fit[[k]]$logLik$loglik
    nrpar[[k]]     <- model_fit[[k]]$NrPar$Obs.nrpar
    nrpar_fac[[k]] <- model_fit[[k]]$NrPar$Fac.nrpar

  } # For loop ends here
  # browser()
  # Also do CHull
  Chull_res      <- CHull(loglik = unlist(LL), nrpar = unlist(nrpar), nsclust = nclus)
  Chull_res_fac  <- CHull(loglik = unlist(LL_fac), nrpar = unlist(nrpar_fac), nsclust = nclus)

  # Chull function already returns a matrix with LL and nrpar. Use it as a base for the rest of the results
  # browser()
  overview <- cbind(unlist(R2entropy),
                    Chull_res,
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

  colnames(overview) <- c("R2entropy", "Clusters",
                          "LL", "nrpar", "Chull", "BIC_G", "BIC_N", "AIC", "AIC3", "ICL",
                          "LL_fac", "nrpar_fac", "Chull_fac", "BIC_G_fac", "BIC_N_fac", "AIC_fac", "AIC3_fac", "ICL_fac")

  overview         <- as.data.frame(overview)
  names(model_fit) <- paste0(nclus[1]:nclus[length(nclus)], "-cluster model")

  output <- (list(Overview = overview,
              Models   = model_fit)
             )

  class(output) <- "mmgsem"

  return(output)
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
        # nsclusthull <- fitMat[convexhull, 1]
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








