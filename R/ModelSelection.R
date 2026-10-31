#' Model selection wrapper for MMGSEM
#'
#' A wrapper to select the number of clusters for MMG-SEM. By default it
#' estimates cross-sectional models with [mmgsem()]. If `time` is supplied, it
#' fits longitudinal models with [lmmgsem()] and reuses the measurement
#' invariance fit across candidate cluster counts.
#'
#' @param dat Data frame containing observed variables and a grouping variable..
#' @param S1 Step 1 (measurement model) specification. Can be omitted when
#'   `s1_fit` is supplied. In cross-sectional mode it may be lavaan syntax, or
#'   a list of strings for multiple measurement blocks.
#' @param S2 Step 2 (structural model) specification using lavaan syntax.
#' @param group Name of the grouping variable in `dat`.
#' @param time Optional name of a time variable. When supplied, fit longitudinal
#'   models with [lmmgsem()].
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
#' @param s1_fit Optional fitted Step 1 model. In longitudinal mode, this must
#'   be a lavaan multigroup fit with one group per group-time cell. In
#'   cross-sectional mode it may also be a blavaan or MplusAutomation fit.
#'   When omitted, `ModelSelection()` computes and reuses the fit.
#' @param max_it Maximum number of iterations for Step 2 (default = 10000).
#' @param nstarts Number of random starts for Step 2 (default = 20)..
#' @param printing Logical. If TRUE, prints progress messages during estimation.
#' @param partition Initialisation partition for random starts: "hard" (default) or "soft".
#' @param endogenous_cov Logical. If TRUE (default), residual covariances among purely endogenous latent variables are estimated. If FALSE, they are fixed to 0 and only the residual variances are estimated.
#' @param endo_group_specific Logical. If TRUE (default), residual covariances are group- and cluster-specific (i.e., they are estimated for every group-cluster combination). If FALSE, they are cluster-specific (i.e., they are fixed to be equal across groups within a cluster). Note that, if FALSE, the residual covariances will also influence the clustering.
#' @param sam_method either "local" or "global. Follows local and global approaches from the SAM method. GLOBAL NOT FUNCTIONAL YET.
#' @param rescaling Only used when data is ordered. By default, MMGSEM uses the marker variable scaling approach. But identification
#'                  constraints with ordinal data (by default) are handled by standardizing the factors' variance in the first step.
#'                  The rescaling argument (either T or F) rescales the factor variances and loadings to the marker variable scaling
#'                  before running step 2. It is set to T by default (rescaling happens). If set to F, the factor variances are kept fixed to 1.
#' @param meanstr Logical. If TRUE, includes the mean structure in the model (e.g., for scalar invariance).
#' @param only_slow Logical. If TRUE, uses the slower Step 2 estimation implementation.
#' @param s2_fast If `time` is supplied, use the fast longitudinal Step 2
#'   estimator when TRUE and the slow estimator when FALSE.
#' @param mi_method Score-test selection rule: `"group_violations"` or
#'   `"test_statistic"`; used when `time` is supplied.
#' @param mi_alpha Significance level for longitudinal score tests.
#' @param mi_d_cfi_threshold Maximum configural-to-metric CFI difference in the
#'   longitudinal invariance analysis.
#' @param mi_max_iter Maximum iterations per slice in the longitudinal
#'   invariance analysis.
#' @param mi_fit_indices Fit indices returned by the longitudinal invariance
#'   analysis.
#' @param mi_specific_noninv If `TRUE`, allow non-invariance for specific
#'   group-time pairs.
#'
#' @details
#' In longitudinal mode, every group must have data at all time points and at
#' least two time points must be present. The longitudinal information criteria
#' use the sequence-level factor likelihood and marginal state entropy. The
#' observed and factor likelihood columns are equal because a separate
#' indicator-level longitudinal likelihood is not calculated. Convex-hull
#' scores are unavailable when fewer than three cluster counts are fitted.
#' Parameter counts include free measurement parameters, group-time-specific
#' exogenous and endogenous covariance parameters, per-state regression
#' parameters, and initial/transition probabilities.
#' @param ... MMGSEM relies on lavaan for the estimation of the first step (i.e., CFA). If needed, the users can pass any lavaan argument to MMGSEM
#'            and it will be considered when estimating the CFA. For instance, std.lv if users want standardized latent variables,
#'            group.equal for constraints, group.partial for non-invariances, etc.
#'
#' @returns Overview: Matrix. It contains the results of all the model selection measures for each model.
#' @returns Models: List containing the separate `mmgsem()` results, or
#'   `lmmgsem()` results when `time` is supplied.
#'
#' @export
ModelSelection <- function(dat, S1 = NULL, S2 = NULL,
                           s1_type = "lavaan",
                           group, nclus, seed = NULL,
                           userStart = NULL, s1_fit = NULL,
                           max_it = 10000L, nstarts = 20L, printing = FALSE,
                           partition = "hard", endogenous_cov = TRUE,
                           endo_group_specific = TRUE,
                           sam_method = "local", meanstr = FALSE,
                           rescaling = F, only_slow = FALSE,
                           time = NULL, s2_fast = TRUE,
                           mi_method = c("group_violations", "test_statistic"),
                           mi_alpha = 0.05, mi_d_cfi_threshold = 0.01,
                           mi_max_iter = 20L,
                           mi_fit_indices = c("chisq", "df", "pvalue", "cfi",
                                              "rmsea", "srmr"),
                           mi_specific_noninv = FALSE,
                           ...) {

  if (length(nclus) != 2L || !is.numeric(nclus) ||
      any(!is.finite(nclus)) ||
      any(nclus != as.integer(nclus)) || nclus[1] < 1 ||
      nclus[2] < nclus[1]) {
    stop("'nclus' must contain increasing positive integer bounds.",
         call. = FALSE)
  }

  if (!is.null(time)) {
    mi_method <- match.arg(mi_method)
    if (!is.logical(s2_fast) || length(s2_fast) != 1L || is.na(s2_fast) ||
        !is.logical(only_slow) || length(only_slow) != 1L ||
        is.na(only_slow)) {
      stop("'s2_fast' and 'only_slow' must be TRUE or FALSE.",
           call. = FALSE)
    }
    if (!identical(s1_type, "lavaan") ||
        !is.null(userStart) || !identical(partition, "hard") ||
        !identical(endogenous_cov, TRUE) ||
        !identical(endo_group_specific, TRUE) ||
        !identical(sam_method, "local") || !identical(meanstr, FALSE) ||
        !identical(rescaling, FALSE)) {
      stop("Longitudinal model selection uses lavaan measurement models and ",
           "does not support the cross-sectional-only options.",
           call. = FALSE)
    }
    clusters <- seq.int(nclus[1], nclus[2])
    model_fit <- vector("list", length(clusters))
    for (i in seq_along(clusters)) {
      model_fit[[i]] <- lmmgsem(
        dat = dat,
        S1 = S1,
        S2 = S2,
        group = group,
        time = time,
        nclus = clusters[i],
        seed = seed,
        max_it = max_it,
        nstarts = nstarts,
        printing = printing,
        s2_fast = s2_fast && !only_slow,
        s1_fit = s1_fit,
        mi_method = mi_method,
        mi_alpha = mi_alpha,
        mi_d_cfi_threshold = mi_d_cfi_threshold,
        mi_max_iter = mi_max_iter,
        mi_fit_indices = mi_fit_indices,
        mi_specific_noninv = mi_specific_noninv,
        ...
      )
      if (is.null(s1_fit)) {
        s1_fit <- model_fit[[i]]$MM
      }
    }

    loglik <- vapply(model_fit, function(x) x$logLik$loglik, numeric(1))
    nrpar <- vapply(model_fit, function(x) x$NrPar$Obs.nrpar, numeric(1))
    nrpar_fac <- vapply(model_fit, function(x) x$NrPar$Fac.nrpar, numeric(1))
    chull_scores <- rep(NA_real_, length(clusters))
    chull_scores_fac <- rep(NA_real_, length(clusters))
    if (length(clusters) >= 3L) {
      chull_fit <- CHull(
        loglik = loglik,
        nrpar = nrpar,
        nsclust = c(1L, length(clusters))
      )
      chull_scores <- chull_fit[, 4]
      chull_fit_fac <- CHull(
        loglik = loglik,
        nrpar = nrpar_fac,
        nsclust = c(1L, length(clusters))
      )
      chull_scores_fac <- chull_fit_fac[, 4]
    }
    model_sel <- lapply(model_fit, `[[`, "model_sel")
    overview <- data.frame(
      R2entropy = vapply(model_sel, `[[`, numeric(1), "R2_entropy"),
      Clusters = clusters,
      LL = loglik,
      nrpar = nrpar,
      Chull = chull_scores,
      BIC_G = vapply(model_sel, function(x) x$BIC$observed$BIC_G, numeric(1)),
      BIC_N = vapply(model_sel, function(x) x$BIC$observed$BIC_N, numeric(1)),
      AIC = vapply(model_sel, function(x) x$AIC$observed, numeric(1)),
      AIC3 = vapply(model_sel, function(x) x$AIC3$observed, numeric(1)),
      ICL = vapply(model_sel, function(x) x$ICL$observed, numeric(1)),
      LL_fac = loglik,
      nrpar_fac = nrpar_fac,
      Chull_fac = chull_scores_fac,
      BIC_G_fac = vapply(model_sel, function(x) x$BIC$Factors$BIC_G, numeric(1)),
      BIC_N_fac = vapply(model_sel, function(x) x$BIC$Factors$BIC_N, numeric(1)),
      AIC_fac = vapply(model_sel, function(x) x$AIC$Factors, numeric(1)),
      AIC3_fac = vapply(model_sel, function(x) x$AIC3$Factors, numeric(1)),
      ICL_fac = vapply(model_sel, function(x) x$ICL$Factors, numeric(1))
    )
    names(model_fit) <- paste0(clusters, "-cluster model")
    output <- list(Overview = overview, Models = model_fit)
    class(output) <- c("lmmgsem", "mmgsem")
    return(output)
  }

  # Prepare objects to compare in model selection
  clusters  <- seq.int(nclus[1], nclus[2])
  nmodels   <- length(clusters)
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

  # Prepare object for the progress bar
  message("Model selection progress:\n")
  pb <- utils::txtProgressBar(min = 0, max = nmodels, style = 3)
  # Call MMGSEM once for each requested cluster count.
  for (i_model in seq_along(clusters)) {
    k <- clusters[i_model]
    # If the user provides output of the first step (s1out), use it for all models
    if(!is.null(s1_fit)){
      model_fit[[i_model]] <- mmgsem(dat = dat, S1 = S1, S2 = S2, s1_type = s1_type, group = group, nclus = k, seed = seed,
                               userStart = userStart, s1_fit = s1_fit, max_it = max_it, nstarts = nstarts, printing = printing,
                               partition = partition, endogenous_cov = endogenous_cov, only_slow = only_slow,
                               endo_group_specific = endo_group_specific, sam_method = sam_method, ...)
    } else if (is.null(s1_fit)){
      # browser()
      if(i_model == 1L){
        model_fit[[i_model]] <- mmgsem(dat = dat, S1 = S1, S2 = S2, s1_type = s1_type, group = group, nclus = k, seed = seed,
                                 userStart = userStart, s1_fit = NULL, max_it = max_it, nstarts = nstarts, printing = printing,
                                 partition = partition, endogenous_cov = endogenous_cov, only_slow = only_slow,
                                 endo_group_specific = endo_group_specific, sam_method = sam_method, ...)

        # Save the covariance matrix of step 1, so we do not run it for all models (it is the same for all of them).
        s1_fit <- model_fit[[i_model]]$MM
      } else {
        # browser()
        model_fit[[i_model]] <- mmgsem(dat = dat, S1 = S1, S2 = S2, s1_type = s1_type, group = group, nclus = k, seed = seed,
                                 userStart = userStart, s1_fit = s1_fit, max_it = max_it, nstarts = nstarts, printing = printing,
                                 partition = partition, endogenous_cov = endogenous_cov, only_slow = only_slow,
                                 endo_group_specific = endo_group_specific, sam_method = sam_method, ...)
      }
      # print(paste("model", k, "finished"))
    }
    
    # Update the progress bar
    utils::setTxtProgressBar(pb, i_model)
    
    # Save the model selection criteria
    BIC_G[[i_model]]     <- model_fit[[i_model]]$model_sel$BIC$observed$BIC_G
    BIC_N[[i_model]]     <- model_fit[[i_model]]$model_sel$BIC$observed$BIC_N
    BIC_G_fac[[i_model]] <- model_fit[[i_model]]$model_sel$BIC$Factors$BIC_G
    BIC_N_fac[[i_model]] <- model_fit[[i_model]]$model_sel$BIC$Factors$BIC_N
    AIC[[i_model]]       <- model_fit[[i_model]]$model_sel$AIC$observed
    AIC_fac[[i_model]]   <- model_fit[[i_model]]$model_sel$AIC$Factors
    AIC3[[i_model]]      <- model_fit[[i_model]]$model_sel$AIC3$observed
    AIC3_fac[[i_model]]  <- model_fit[[i_model]]$model_sel$AIC3$Factors
    ICL[[i_model]]       <- model_fit[[i_model]]$model_sel$ICL$observed
    ICL_fac[[i_model]]   <- model_fit[[i_model]]$model_sel$ICL$Factors
    R2entropy[[i_model]] <- model_fit[[i_model]]$model_sel$R2_entropy
    LL[[i_model]]        <- model_fit[[i_model]]$logLik$obs_loglik
    LL_fac[[i_model]]    <- model_fit[[i_model]]$logLik$loglik
    nrpar[[i_model]]     <- model_fit[[i_model]]$NrPar$Obs.nrpar
    nrpar_fac[[i_model]] <- model_fit[[i_model]]$NrPar$Fac.nrpar

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
  names(model_fit) <- paste0(clusters, "-cluster model")

  output <- (list(Overview = overview,
              Models   = model_fit)
             )

  class(output) <- "mmgsem"

  return(output)
}
