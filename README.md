Readme
================

Introduction
------------

In the dimension reduction approach in the GJAM model, the species are clustered in their dependence behavior. In the extended model we can use the prior knowledge on the number of clusters in the GJAM model and obtain the cluster estimates.

Reproduce the results
---------------------

To reproduce the analysis results, there are two possibilites:

1.  Run the fit\_script.R, which will run and save the models in the folder "PCA analysis". This option could be time consuming.
2.  Download saved models from the following link [models](https://drive.google.com/open?id=1sRu7Q7rJ4aIYp-YCKOa_igA5Ofzc7tRH) and load the models in the folder "PCA analysis/r5\_models/chain\_1"(2) on your local repository

The analysis results could be reproduced by the "analysis.R" script.
**to run specifc functions for now use firstly script inlclude.R**

How to use the functions
------------------------

Firstly, we describe the initial GJAM model and then describe how to specify the parameters for the prior. Full description of the GJAM model is provided in the GJAM package documentation and [vignette](https://cran.r-project.org/web/packages/gjam/vignettes/gjamVignette.html)

Dimension reduction in GJAM model was proposed by \[Taylor et el 2018\]. Dimension reduction is used in two ways:

-   When dataset contains more species than can be fitted given sample size *n*.
-   When the number of species *S* is too large.

The first case is used with default parameters. For the second case we need to specify the parameters for dimension `reduction` in `reductList` and add this parameters in the `modelList`.

``` r
load("Bauges_plant.Rdata")
y<- Bauges_plant[,3:(ncol(Bauges_plant))]
Y_data<- gjamTrimY(y,20)$y  ## trim the species that occur less then in 20 cites
X_data <- Bauges_plant[,1:2] # environmental covariates
S<- ncol(Y_data)  # number of species
formula <- as.formula( ~   PC1  + PC2 + I(PC1^2) + I(PC2^2))
iterations=50# number of iterations
burn_period=10  # burn-in period
r_reduct = 5      # rank of the covariance matrix approximation
```

For the standard dimension reduction regime in GJAM model version, parameters for `reducltList` and `modelList` need to be specified. In the `reducltList`

-   N is maximum number of possible clusters. *N* ≤ *S*. Could be limited for *S* large to 150. (Technically: the truncation number in Dirichlet Multinomial process)
-   r is the number of latent factors. *r* ≤ *S*.

In the `'modelList'`

-   `typeNames` - specify the type of response data
-   `ng` is the number of iterations
-   `reductList` is the parameters for the dimension reduction

``` r
rl <- list(r =r_reduct, N = S)
ml   <- list(ng = iterations, burnin = burn_period, typeNames = 'PA', reductList = rl,PREDICTX = F)
fit<-gjam(formula, xdata = X_data, ydata = Y_data, modelList = ml)
```

To incorporate prior knowledge in the model, we firstly speciify, the Bayesian nonparametric prior in the variable `DRtype` (see Table 1)

Then `reducltList` in this case will be In the `reducltList`

-   `N` is maximum number of possible clusters. *N* ≤ *S*
-   `r` is the number of latent factors. *r* ≤ *S*.
-   `DRtype` is the type of the Baysian nonparametric prior.
-   `K` is prior number of clusters

For the Pitman--Yor process as there are two parameters PY(*α*, *σ*), we specify the desired value of *σ* ∈ (0, 1). (in the current implementation).

**Table 1. Dimension reduction types**

| `DRtype` | Process type | Parameters |           Priors on parameters           | Parameters to specify      |
|:--------:|:------------:|:----------:|:----------------------------------------:|----------------------------|
|   `'-'`  |   Dirichlet  |     *α*    |                   none                   | none                       |
|   `'1'`  |   Dirichlet  |     *α*    | *G**a*(*ν*<sub>1</sub>, *ν*<sub>2</sub>) | prior number of clusters K |
|   `'2'`  |  Pitman--Yor | (*α*, *σ*) |                   none                   | prior number of clusters K |

More precisely on the values of the parameter `DRtype`:

-   No value specified, then the model will use standard dimension reduction approach.
-   `1` is the dimension reduction with the Dirichlet process and prior distribution on the concentration parameter *α*
-   `2` is the dimension reduction with the Pitman--Yor process with fixed parameters *α* and *σ*.

For the example, in the **Bauges plant dataset ** , we take as a prior number of plant functional groups which is 16. If we want to take the Dirichlet process, then we use the following parameters for the `reductList`. Here, `DRtype =1`. The choice of *N* is equal to *S* is natural, as we consider any possible clustering and *K* is prior number of clusters.

``` r
rl1  <- list(r = r_reduct, N = S ,DRtype="1", K=16)  #prior is number of plant functional groups
```

For the Pitman--Yor process .Here, `DRtype =2` and we specify `K`,and value of *σ* = 0.5 and the choice of *N* is the same.

``` r
rl2  <- list(r = r_reduct, N = S ,DRtype="2", K=16)  #prior is number of plant functional groups
```

For instance, we fit the model with the first specification

``` r
ml1   <- list(ng = iterations, burnin = burn_period, typeNames = 'PA', reductList = rl1,PREDICTX = F)
fit1<-gjam(formula, xdata = X_data, ydata = Y_data, modelList = ml1)
```

We get the `fit1` gjam object. This object is the same as the one from the original model, except several additional parameters.

The function `gjjamCluster` estimates the optimal cluster,that summarize posterior cluster distribution. The function is build using the `GreedyEPL` package by
*R**a**s**t**e**l**l**i**e**t**a**l*.
. It takes as the input \* gjam object fitted model \* K possible number of clusters for initialization of greedy search algorithm
\* Clustering used for prior specification

-   last two options are used to set the initial clustering for the algorithm that estimate the optimal clustering.

This function return the values

-   VI\_list list of the estimated clusters, corresponding to the given starting points
-   EPL\_value list of the loss function values for the given starting points

The estimated cluster with the smallest EPL\_value best represent the posterior clustering distribution.
