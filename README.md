# MMG-SEM: Mixture Multigroup SEM
R code for Mixture Multigroup SEM (MMG-SEM). MMG-SEM performs clustering of a SEM model based on the structural parameters while allowing for a combination of group-specific and invariant parameters in the measurement model. For more information (and citation) about the method, please see: [Mixture Multigroup Structural Equation Modeling: A Novel Method for Comparing Structural Relations Across Many Groups](https://doi.org/10.1037/met0000667)

## Installation
It is possible to download the (in development) mmgsem package from GitHub by:
```R
library("devtools")
devtools::install_github("AndresFPA/mmgsem")
```

# Functions
+ **MMG-SEM**: The main function is to perform the clustering based on the structural parameters.
+ **ModelSelection**: Wrapper function. It performs MMG-SEM as many times as the user requires with different numbers of clusters, returning model selection measures.
+ **SE**: Wrapper function. It computes the standard errors of the resulting parameters of an MMG-SEM analysis.

### MMG-SEM
MMG-SEM performs clustering at the group level (e.g., clusters of nations). It deals with regression parameter differences. For instance, one can answer the question: which countries present a similar relation between negative emotions and life satisfaction? For its estimation, we use a step-wise approach based on the Structural-After-Measurement approach (see [SAM](https://psycnet.apa.org/doi/10.1037/met0000503)).

### Model Selection
The model selection wrapper computes several model selection measures. Specifically, it computes the AIC, AIC3, the CHull, ICL, and the BIC. The BIC is computed in two different forms: (1) *BIC<sub>G</sub>* with the number of groups as sample size, and (2) *BIC<sub>N</sub>* with number of observations as sample size. 

### Standard Errors
The SE wrapper function computes the standard errors based on the Hessian matrix of the parameters. It relies on numerical derivation to obtain the second derivative of the loglikelihood function with respect to the parameters. Note that it also computes the *corrected* standard errors based on [Bakk et al](https://doi.org/10.1093/pan/mpu003). (2014).

## Funding
The research behind this R package was funded by a Vidi grant [VI.Vidi.201.133] awarded to Kim De Roover by the Netherlands Organization for Scientific Research (NWO). 
