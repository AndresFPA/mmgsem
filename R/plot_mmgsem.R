#' Plot model selection results
#'
#' Plots the graphs often used to determine the optimal number of clusters (e.g., scree plots). Provides plots for all the model selection measures
#'
#' INPUT: Arguments required by the function
#' @param ModelSel Object returned from the ModelSel() function.
#'
#' @return plots List that contains the plots for all the model selection measures.
#'
#' @export
plot_MMGSEM <- function(ModelSel, criteria){
  if (!is.list(ModelSel)) {
    stop("Input must be a list, likely the output of ModelSelection().")
  }
  overview <- ModelSel$Overview

  # Plots for visual inspection
  plot.BIC_G     <- ggplot2::ggplot(data = overview, mapping = ggplot2::aes(x = Clusters, y = BIC_G))     + ggplot2::geom_point() + ggplot2::geom_line() # plot(x = overview$Clusters, y = overview$BIC_G); lines(x = overview$Clusters, y = overview$BIC_G)
  plot.BIC_N     <- ggplot2::ggplot(data = overview, mapping = ggplot2::aes(x = Clusters, y = BIC_N))     + ggplot2::geom_point() + ggplot2::geom_line() # plot(x = overview$Clusters, y = overview$BIC_N); lines(x = overview$Clusters, y = overview$BIC_N)
  plot.BIC_G_fac <- ggplot2::ggplot(data = overview, mapping = ggplot2::aes(x = Clusters, y = BIC_G_fac)) + ggplot2::geom_point() + ggplot2::geom_line() # plot(x = overview$Clusters, y = overview$BIC_G_fac); lines(x = overview$Clusters, y = overview$BIC_G_fac)
  plot.BIC_N_fac <- ggplot2::ggplot(data = overview, mapping = ggplot2::aes(x = Clusters, y = BIC_N_fac)) + ggplot2::geom_point() + ggplot2::geom_line() # plot(x = overview$Clusters, y = overview$BIC_N_fac); lines(x = overview$Clusters, y = overview$BIC_N_fac)
  plot.AIC       <- ggplot2::ggplot(data = overview, mapping = ggplot2::aes(x = Clusters, y = AIC))       + ggplot2::geom_point() + ggplot2::geom_line() # plot(x = overview$Clusters, y = overview$AIC); lines(x = overview$Clusters, y = overview$AIC)
  plot.AIC_fac   <- ggplot2::ggplot(data = overview, mapping = ggplot2::aes(x = Clusters, y = AIC_fac))   + ggplot2::geom_point() + ggplot2::geom_line() # plot(x = overview$Clusters, y = overview$AIC_fac); lines(x = overview$Clusters, y = overview$AIC_fac)
  plot.AIC3      <- ggplot2::ggplot(data = overview, mapping = ggplot2::aes(x = Clusters, y = AIC3))      + ggplot2::geom_point() + ggplot2::geom_line() # plot(x = overview$Clusters, y = overview$AIC3); lines(x = overview$Clusters, y = overview$AIC3)
  plot.AIC3_fac  <- ggplot2::ggplot(data = overview, mapping = ggplot2::aes(x = Clusters, y = AIC3_fac))  + ggplot2::geom_point() + ggplot2::geom_line() # plot(x = overview$Clusters, y = overview$AIC3_fac); lines(x = overview$Clusters, y = overview$AIC3_fac)
  plot.ICL       <- ggplot2::ggplot(data = overview, mapping = ggplot2::aes(x = Clusters, y = ICL))       + ggplot2::geom_point() + ggplot2::geom_line() # plot(x = overview$Clusters, y = overview$ICL); lines(x = overview$Clusters, y = overview$ICL)
  plot.ICL_fac   <- ggplot2::ggplot(data = overview, mapping = ggplot2::aes(x = Clusters, y = ICL_fac))   + ggplot2::geom_point() + ggplot2::geom_line() # plot(x = overview$Clusters, y = overview$ICL_fac); lines(x = overview$Clusters, y = overview$ICL_fac)
  plot.Chull     <- ggplot2::ggplot(data = overview, mapping = ggplot2::aes(x = Clusters, y = LL))        + ggplot2::geom_point() + ggplot2::geom_line() # plot(x = overview$Clusters, y = overview$LL); lines(x = overview$Clusters, y = overview$LL)
  plot.Chull_fac <- ggplot2::ggplot(data = overview, mapping = ggplot2::aes(x = Clusters, y = LL_fac))    + ggplot2::geom_point() + ggplot2::geom_line() # plot(x = overview$Clusters, y = overview$LL_fac); lines(x = overview$Clusters, y = overview$LL_fac)

  plots <- list(BIC_G = plot.BIC_G, BIC_N = plot.BIC_N, BIC_G_fac = plot.BIC_G_fac,
                BIC_N_fac = plot.BIC_N_fac, AIC = plot.AIC, AIC_fac = plot.AIC_fac, AIC3 = plot.AIC3,
                AIC3_fac = plot.AIC3_fac, ICL = plot.ICL, ICL_fac = plot.ICL_fac,
                Chull = plot.Chull, Chull_fac = plot.Chull_fac)

  if (is.null(criteria)){
    return(plots)
  } else {
    return_object <- plots[names(plots) == criteria]
    return(return_object)
  }

}
