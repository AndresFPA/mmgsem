[![R-CMD-check](https://github.com/AndresFPA/mmgsem/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/AndresFPA/mmgsem/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/AndresFPA/mmgsem/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/AndresFPA/mmgsem/actions/workflows/pkgdown.yaml)
[![License: MIT](https://img.shields.io/github/license/AndresFPA/mmgsem)](https://github.com/AndresFPA/mmgsem/blob/main/LICENSE)

# MMG-SEM: Mixture Multigroup SEM

R code for Mixture Multigroup SEM (MMG-SEM). MMG-SEM performs clustering of a SEM model based on the structural parameters while allowing for a combination of group-specific and invariant parameters in the measurement model. For more information (and citation) about the method, please see: [Mixture Multigroup Structural Equation Modeling: A Novel Method for Comparing Structural Relations Across Many Groups](https://doi.org/10.1037/met0000667)

## Installation

It is possible to download the (in development) mmgsem package from GitHub by:

```R
library("pak")
pak::pak("AndresFPA/mmgsem")
```

## Quick start

This example uses the package's example data to estimate a four-cluster model. The measurement model is estimated using `lavaan`; the structural model is then clustered by `mmgsem()`.

```R
data <- mmgsem::data_example

S1_model <- '
    F1 =~ V1 + V2 + V3 + V4 + V5
    F2 =~ V6 + V7 + V8 + V9 + V10
    F3 =~ V11 + V12 + V13 + V14 + V15
    F4 =~ V16 + V17 + V18 + V19 + V20
'

S2_model <- '
    F4 ~ F1 + F3
    F3 ~ F1 + F2
'

MM_fit <- lavaan::cfa(
  model         = S1_model, # Step 1 model syntax
  data          = data,     # Data
  group         = "group",  # Group variable
  group.equal   = "loadings", # Measurement model constraints (metric invariance)
  group.partial = c("F1 =~ V2", # Measurement non-invariances
                    "F2 =~ V7",
                    "F3 =~ V12",
                    "F4 =~ V17")
)

mmgsem_fit <- mmgsem::mmgsem(
  dat     = data,     # Data
  S1      = S1_model, # Step 1 (measurement) model syntax
  S2      = S2_model, # Step 2 (structural) model syntax
  group   = "group",  # Group variable
  nclus   = 4,        # Number of clusters
  seed    = 1,        # Seed for reproducible results
  s1_type = "lavaan", # Step 1 model type
  s1_fit  = MM_fit    # Optional fitted Step 1 lavaan model
)

summary(mmgsem_fit)
print(mmgsem_fit$posteriors)
```

To compare models with different numbers of clusters, use `ModelSelection()`. This fits each multiple models following the `nclus` argument (min, max). Its computation can take considerably longer than fitting a single model.

```R
modelSelection_fit <- mmgsem::ModelSelection(
  dat     = data,     # Data
  S1      = S1_model, # Step 1 (measurement) model syntax
  S2      = S2_model, # Step 2 (structural) model syntax
  group   = "group",  # Group variable
  nclus   = c(1, 6),  # Lower and upper cluster counts to fit
  seed    = 1,        # Seed for reproducible results
  s1_type = "lavaan", # Step 1 model type
  s1_fit  = MM_fit    # Optional fitted Step 1 lavaan model
)

summary(model = modelSelection_fit, model_selection = TRUE)
selected_fit <- mmgsem::extract(object = modelSelection_fit, nclus = 4)
```

Computing standard errors uses numerical differentiation and can take a long time, especially for larger models. Run this step when you need standard errors; the package also includes precomputed standard errors for this example data.

```R
# standard_errors <- mmgsem::standard_errors ## load the pre-computed standard errors for this example
standard_errors <- mmgsem::compute_se(object = selected_fit) # Fitted model for standard-error computation
summary(model = selected_fit, se = standard_errors)
```

## Tutorial

An in-depth tutorial on how to use this package has been published in _Structural Equation Modeling: A Multidisciplinary Journal_, please see: [Using Mixture Multigroup Structural Equation Modeling to Compare Structural Relations Across Many Groups: An Introduction to the mmgsem R Package](https://doi.org/10.1080/10705511.2026.2617170). An accompanying R script for the tutorial can also be found in this complementary [GitHub Repository](https://github.com/AndresFPA/Tutorial-mmgsem).

## (Main) Functions
+ `mmgsem()`: The main function performs the clustering based on the structural parameters.
+ `ModelSelection()`: Wrapper function. It runs MMG-SEM multiple times with different numbers of clusters and returns model selection measures.
+ `compute_se()`: Wrapper function. It computes the standard errors of the resulting MMG-SEM parameters.

### MMG-SEM
`mmgsem()` performs clustering at the group level (e.g., clusters of nations). It deals with regression parameter differences. For instance, one can ask: which countries show a similar relation between negative emotions and life satisfaction? For its estimation, we use a step-wise approach based on the Structural-After-Measurement approach (see [SAM](https://psycnet.apa.org/doi/10.1037/met0000503)).

### Model Selection
The `ModelSelection()` wrapper computes several model selection measures. Specifically, it computes the AIC, AIC3, the CHull, ICL, and the BIC. The BIC is computed in two different forms: (1) *BIC<sub>G</sub>* with the number of groups as sample size, and (2) *BIC<sub>N</sub>* with the number of observations as sample size. Performance of the different measures was studied in detail at [Selecting the number of clusters in Mixture Multigroup Structural Equation Modeling](https://doi.org/10.5964/meth.14931).

### Standard Errors
The `compute_se()` wrapper function computes standard errors from the Hessian matrix of the parameters. It uses numerical differentiation to compute the second derivative of the log-likelihood with respect to the parameters. Note that it also computes the *corrected* standard errors based on [Bakk et al](https://doi.org/10.1093/pan/mpu003). (2014).

# References
Bakk, Z., Oberski, D. L., & Vermunt, J. K. (2014). Relating Latent Class Assignments to External Variables: Standard Errors for Correct Inference. *Political Analysis, 22*(4), 520–540. [DOI](https://doi.org/10.1093/pan/mpu003)

Perez Alonso, A. F., Rosseel, Y., Vermunt, J. K., & De Roover, K. (2024). Mixture multigroup structural equation modeling: A novel method for comparing structural relations across many groups. *Psychological Methods*. [DOI](https://doi.org/10.1037/met0000667)

Perez Alonso, A. F., Vermunt, J. K., Rosseel, Y., & Roover, K. D. (2025). Selecting the number of clusters in Mixture Multigroup Structural Equation Modeling. *Methodology, 21*(1), 1–26. [DOI](https://doi.org/10.5964/meth.14931)

Perez Alonso, A. F., Vermunt, J. K., Rosseel, Y., Zhao, H., & De Roover, K. (2026). Using Mixture Multigroup Structural Equation Modeling to Compare Structural Relations Across Many Groups: An Introduction to the mmgsem R Package. *Structural Equation Modeling: A Multidisciplinary Journal*, 1–28. [DOI](https://doi.org/10.1080/10705511.2026.2617170)

Rosseel, Y., & Loh, W. W. (2022). A structural after measurement approach to structural equation modeling. *Psychological Methods*. [DOI](https://doi.org/10.1037/met0000503)

# Funding
The research behind this R package was funded by a Vidi grant [VI.Vidi.201.133](https://doi.org/10.61686/KLLGI61008) awarded to Kim De Roover by the Netherlands Organization for Scientific Research (NWO). 
