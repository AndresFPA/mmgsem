#' Example data
#'
#' Example data used in the tutorial paper of mmgsem
#'
#' @format ##
#' A data frame with 100800 rows and 21 columns:
#' \describe{
#'   \item{V1-V20}{Observed variables of our example data}
#'   \item{group}{Grouping variable}
#' }
"data_example"


#' Standard errors (4K model)
#'
#' List object containing the standard errors of the 4-cluster model from the example data of the tutorial
#'
#' @format ##
#' List
#' \describe{
#'   \item{SE_vector}{List of vectors containing the standard errors per type of parameter.}
#'   \item{HESS}{Hessian matrix (second derivative) of the loglikelihood function with respect of all parameters.}
#'   \item{vcov_betas}{Variance-covariance matrix of the beta parameters. Used for hypothesis testing.}
#' }
"standard_errors"


#' Standard errors (5K model)
#'
#' List object containing the standard errors of the 5-cluster model from the example data of the tutorial
#'
#' @format ##
#' List
#' \describe{
#'   \item{SE_vector}{List of vectors containing the standard errors per type of parameter.}
#'   \item{HESS}{Hessian matrix (second derivative) of the loglikelihood function with respect of all parameters.}
#'   \item{vcov_betas}{Variance-covariance matrix of the beta parameters. Used for hypothesis testing.}
#' }
"five_K_se"
