% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MMG-SEM.R
\name{MMGSEM}
\alias{MMGSEM}
\title{Mixture Multi-Group Structural Equation Modelling (MMGSEM)}
\usage{
MMGSEM(
  dat,
  S1 = NULL,
  S2 = NULL,
  group,
  nclus,
  seed = NULL,
  userStart = NULL,
  s1_fit = NULL,
  max_it = 10000L,
  nstarts = 20L,
  printing = FALSE,
  partition = "hard",
  endogenous_cov = TRUE,
  endo_group_specific = TRUE,
  sam_method = "local",
  meanstr = FALSE,
  rescaling = F,
  only_slow = FALSE,
  ...
)
}
\arguments{
\item{dat}{Observed data of interest for the MMGSEM model.}

\item{S1}{Measurement model (MM). Used in step 1. Must be a string (like in lavaan).
Can be a list of strings determing the number of measurement blocks (e.g., one string for the MM of
factor 1, and a second string for the MM of factor 2)}

\item{S2}{Structural model. Used in step 2. Must be a string (like in lavaan).}

\item{group}{Name of the group variable. Must be a string.}

\item{nclus}{Pre-specified number of clusters.}

\item{seed}{Pre-defined seed for the random start in case a replication is needed in the future.}

\item{userStart}{Pre-defined start provided by the user. Must be a matrix of dimensions ngroups*nclus.
A 1 represents that the row group is in the column cluster. There must be only one 1
in each row. The function will return results only for this start.
Example for 6 groups and 2 clusters:

\if{html}{\out{<div class="sourceCode">}}\preformatted{                              [,1] [,2]
                        [1,]    1    0
                        [2,]    1    0
                        [3,]    1    0
                        [4,]    0    1
                        [5,]    0    1
                        [6,]    0    1
}\if{html}{\out{</div>}}}

\item{s1_fit}{Resulting lavaan object from a cfa() analysis. Can be used to directly input the results from the step 1
if the user only wants to use MMGSEM() to estimate the structural model (Step 2). If not NULL, MMGSEM()
will skip the estimation of Step 1 and use the s1out object as the input for Step 2. If NULL, MMGSEM() will
estimate both step 1 and step 2.}

\item{max_it}{Maximum number of iterations for each start.}

\item{nstarts}{Number of random starts that will be performed (in an attempt to avoid local maxima).}

\item{partition}{Type of random partition used to start the EM algorithm. Can be "hard" and "soft".}

\item{endogenous_cov}{TRUE or FALSE argument to determine whether to allow or not covariance between purely endogeonus variables.
If TRUE (the default), the covariance between endogenous factors is allowed.}

\item{endo_group_specific}{TRUE or FALSE. Determines whether the endogenous covariances are group (T) or cluster (F) specific.
By default, it is TRUE.}

\item{sam_method}{either "local" or "global. Follows local and global approaches from the SAM method. GLOBAL NOT FUNCTIONAL YET.}

\item{rescaling}{Only used when data is ordered. By default, MMGSEM uses the marker variable scaling approach. But identification
constraints with ordinal data (by default) are handled by standardizing the factors' variance in the first step.
The rescaling argument (either T or F) rescales the factor variances and loadings to the marker variable scaling
before running step 2. It is set to T by default (rescaling happens). If set to F, the factor variances are kept fixed to 1.}

\item{...}{MMGSEM relies on lavaan for the estimation of the first step (i.e., CFA). If needed, the users can pass any lavaan argument to MMGSEM
and it will be considered when estimating the CFA. For instance, std.lv if users want standardized latent variables,
group.equal for constraints, group.partial for non-invariances, etc.

OUTPUT: The function will return a list with the following results:}

\item{constraints}{String containing the constrains of the model. Receives the same input as group.equal in lavaan.}

\item{NonInv}{String containing the non-invariances of the model. Similar to group.partial in lavaan.}
}
\value{
z_gks: Posterior classification probabilities or cluster memberships (ngroups*nclus matrix).

final_fit: Lavaan fit of the best and final model (not useful, only for testing purposes).

psi_gks: Resulting psi matrices (i.e., residual factor covariances) for each group-cluster combination.

lambda: Resulting invariant lambda matrix (i.e., factor loadings) for all groups.

theta_gs: Resulting theta matrices (i.e., items' residual/unique variances) for each group.

beta_ks: Resulting beta matrices (i.e., regression coefficients) for each cluster.

loglikelihood: Total loglikelihood of the best model.

logliks_starts: Total loglikelihood of the models of each run (allows to check sensitivity to local maxima).

Obs_loglik: Reconstruction of the observed loglikelihood after the model has converged. Only valid when
fit is "factors".

PLEASE NOTE: This function requires 'lavaan' package to work.
}
\description{
Performs a mixture clustering based on the structural parameters (i.e., regressions) of a SEM model.
The estimation is done in a step-wise fashion and uses an expectation-maximization (EM) algorithm in the second step.
}
\details{
INPUT: Arguments required by the function
}
