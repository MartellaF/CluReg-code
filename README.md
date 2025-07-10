# CluReg-code
MATLAB implementation of the ECM algorithm for Multivariate Regression Model Based on Latent Predictors
This repository contains the implementation of a novel multivariate regression model. The model builds latent predictors (LPs) as linear combinations of disjoint groups of explanatory variables, aiming to enhance interpretability and reduce dimensionality.
MRMoLP introduces a structured regression framework where:
- Explanatory variables are grouped into disjoint subsets.
- Each group defines a latent predictor (LP), which enters the regression model.
- Parameters are estimated using an Expectation Conditional Maximization (ECM) algorithm.
This approach is particularly useful when interpretability and grouped variable structure are important.
The method is described in detail in the paper under review **XXX (mettere archive)**

The **CluRegApplication** function applies the Multivariate Regression Model based on Latent Predictors to a real dataset, exploring multiple values for the number of latent groups (Q) and selecting the best solution for each Q based on the highest log-likelihood across several random initializations.
- Input:
    rangeQ: a vector of integers representing different values to try for the number of variable clusters 
    nrandstart: number of starting points to avoid local maxima
    X: an n x J matrix of covariates 
    Y: an n x M matrix of responses
- Output:
    The output is a cell array Results where each row (starting from row 2) corresponds to a value of Q and contains:
    Value of Q
    CPU time
    X and Y
    Estimated parameters: Copt, Wopt, Vopt, Sigmaeopt
    Log-likelihood, iteration count, AIC, BIC, number of parameters
    The number of starting points used
