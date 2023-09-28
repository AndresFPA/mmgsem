# MMG-SEM: Mixture Multigroup SEM
R code for Mixture Multigroup SEM (MMG-SEM). MMG-SEM performs clustering of a SEM model at the group level (i.e., clusters of countries). It performs the clustering based on the structural parameters while allowing for a combination of group-specific and invariant parameters in the measurement model. For more information about the method, please see: https://doi.org/10.31234/osf.io/mvd96

# Functions
MMG-SEM: Main function to perform the clustering based on the structural parameters.

E_Step: Necessary function for MMG-SEM. It performs the E-step of the EM algorithm.

ModelSelection: Wrapper function. It performs MMG-SEM as many times as the user requires with different numbers of clusters, returning model selection measures (i.e., BIC, AIC, Chull).

CHull: Necessary function for the model selection function. It performs the Convex Hull scree-ratio test for model selection.

SE: Wrapper function. It computes the standard errors of the resulting parameters of an MMG-SEM analysis.


