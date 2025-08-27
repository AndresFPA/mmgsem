#' Model selection wrapper for MMGSEM
#'
#' A wrapper to perform model selection of MMGSEM. It is used to answer the common clustering question: How many clusters should we use?
#' It calls the main function MMGSEM several times. One for each required model. Most of its arguments are the same as the MMGSEM main function.
#'
#' @param dat Data frame containing observed variables and a grouping variable..
#' @param S1 Step 1 (measurement model) specification using lavaan syntax. Can be a list of strings determing the number of measurement blocks (e.g., one string for the MM of factor 1, and a second string for the MM of factor 2).
#'           If s1_type = "mplus", S1 must be written in mplus syntax.
#' @param S2 Step 2 (structural model) specification using lavaan syntax.
#' @param group Name of the grouping variable (as a character string).
#' @param nclus Integer vector of length two. Determines the range of cluster solutions to estimate. It must contain the lower and upper boundaries of the models the user wants to estimate. For instance, nclus = c(1, 6) estimates models from 1 to 6 clusters.
#' @param seed Optional. Random seed for replicable results.
#' @param userStart Optional. A user-defined cluster membership matrix (dimensions: number of groups × number of clusters) with binary values (1 or 0; 1 indicates cluster membership). To be used when the user has prior knowledge about the data. There must be only one 1 for each row. Skips random starts.
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
#' @param s1_type String. Determines which model is used when estimating step 1. Can be "lavaan" (CFA), "blavaan" (BCFA), or "mplus" (ML-CFA).
#' @param s1_fit Optional. A fitted model object for Step 1. Can be from lavaan (CFA), blavaan (BCFA), or MplusAutomation (ML-CFA).
#' @param max_it Maximum number of iterations for Step 2 (default = 10000).
#' @param nstarts Number of random starts for Step 2 (default = 20)..
#' @param partition Initialisation partition for random starts: "hard" (default) or "soft".
#' @param endogenous_cov Logical. If TRUE (default), residual covariances among purely endogenous latent variables are estimated. If FALSE, they are fixed to 0 and only the residual variances are estimated.
#' @param endo_group_specific Logical. If TRUE (default), residual covariances are group- and cluster-specific (i.e., they are estimated for every group-cluster combination). If FALSE, they are cluster-specific (i.e., they are fixed to be equal across groups within a cluster). Note that, if FALSE, the residual covariances will also influence the clustering.
#' @param sam_method either "local" or "global. Follows local and global approaches from the SAM method. GLOBAL NOT FUNCTIONAL YET.
#' @param rescaling Only used when data is ordered. By default, MMGSEM uses the marker variable scaling approach. But identification
#'                  constraints with ordinal data (by default) are handled by standardizing the factors' variance in the first step.
#'                  The rescaling argument (either T or F) rescales the factor variances and loadings to the marker variable scaling
#'                  before running step 2. It is set to T by default (rescaling happens). If set to F, the factor variances are kept fixed to 1.
#' @param meanstr Logical. If TRUE, includes the mean structure in the model (e.g., for scalar invariance).
#' @param ordinal Logical. If TRUE, observed variables are treated as ordinal (default = FALSE)
#' @param ... MMGSEM relies on lavaan for the estimation of the first step (i.e., CFA). If needed, the users can pass any lavaan argument to MMGSEM
#'            and it will be considered when estimating the CFA. For instance, std.lv if users want standardized latent variables,
#'            group.equal for constraints, group.partial for non-invariances, etc.
#'
#' @returns Overview: Matrix. It contains the results of all the model selection measures for each model.
#' @returns Models: List containing the separate mmgsem() results for each of the fitted models.
#'
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
      model_fit[[k]] <- mmgsem(dat = dat, S1 = S1, S2 = S2, group = group, nclus = k, seed = seed,
                               userStart = userStart, s1_fit = s1_fit, max_it = max_it, nstarts = nstarts, printing = printing,
                               partition = partition, endogenous_cov = endogenous_cov,
                               endo_group_specific = endo_group_specific, sam_method = sam_method, ...)
    } else if (is.null(s1_fit)){
      # browser()
      if(k == nclus[1]){
        model_fit[[k]] <- mmgsem(dat = dat, S1 = S1, S2 = S2, group = group, nclus = k, seed = seed,
                                 userStart = userStart, s1_fit = NULL, max_it = max_it, nstarts = nstarts, printing = printing,
                                 partition = partition, endogenous_cov = endogenous_cov,
                                 endo_group_specific = endo_group_specific, sam_method = sam_method, ...)

        # Save the covariance matrix of step 1, so we do not run it for all models (it is the same for all of them).
        s1_fit <- model_fit[[k]]$MM
      } else {
        # browser()
        model_fit[[k]] <- mmgsem(dat = dat, S1 = S1, S2 = S2, group = group, nclus = k, seed = seed,
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








