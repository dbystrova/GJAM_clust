Readme
================

Introduction
------------

In the dimension reduction approach in the GJAM model, the species are clustered in their dependence behavior. In the extended model we can use the prior knowledgeon the number of groups in the dimension reduction in the GJAM model and obtain the cluster estimates.

Reproduce the analysis results
------------------------------

To reproduce the analysis results, there are two possibilites:

1.  Run the fit\_script.R, which will run and save the models in the folder "PCA analysis". This option could be time consuming
2.  Download saved models from the following link [models](https://drive.google.com/open?id=1sRu7Q7rJ4aIYp-YCKOa_igA5Ofzc7tRH) and load the models in the folder "PCA analysis/r5\_25" on your local repository

The analysis results could be reproduced by the "analysis.R" script.
**to run specifc functions for now use firstly script inlclude.R**

How to use the functions
------------------------

Firstly, we describe the initial GJAM model and then describe how to specify the parameters for the prior. Full description of the GJAM model is provided in the package documentation \[\] and [vignette](https://cran.r-project.org/web/packages/gjam/vignettes/gjamVignette.html)

Dimension reduction in GJAM model was proposed by \[Taylor et el 2018\]. Dimension reduction is used in two ways:

-   When dataset contains more species than can be fitted given sample size *n*.
-   When the number of species *S* is too large.

For the second case we need to specify the parameters for dimension `reduction` in `reductList` and add this parameters in the `modelList`.

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

Then we specify the `reducltList` and then `modelList` In the `reducltList`

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

We have described the standard version for the gjam Dimension reduction.

To incorporate prior knowledge in the model, we firstly speciify, the Bayesian nonparametric prior in the variable `DRtype` (see Table 1)

Then `reducltList` in this case will be In the `reducltList`

-   `N` is maximum number of possible clusters. *N* ≤ *S*, could be specified for only for `DRtype = 2`
-   `r` is the number of latent factors. *r* ≤ *S*.
-   `DRtype` is the type of the Baysian nonparametric prior.
-   `K` is prior number of clusters
-   `V` is the variance for the prior number of clusters, could be specified for only for `DRtype = 3`

For the Pitman--Yor process as there are two parameters PY(*α*, *σ*), for the case where only K is specified, we fix the value of *σ* = 0.25

**Table 1. Dimension reduction types**

<table style="width:100%;">
<colgroup>
<col width="8%" />
<col width="12%" />
<col width="16%" />
<col width="39%" />
<col width="22%" />
</colgroup>
<thead>
<tr class="header">
<th align="center"><code>DRtype</code></th>
<th align="center">Process type</th>
<th align="center">Parameters</th>
<th align="center">Priors on parameters</th>
<th>Parameters to specify</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="center"><code>'-'</code></td>
<td align="center">Dirichlet</td>
<td align="center"><span class="math inline"><em>α</em></span></td>
<td align="center">none</td>
<td>none</td>
</tr>
<tr class="even">
<td align="center"><code>'2'</code></td>
<td align="center">Dirichlet</td>
<td align="center"><span class="math inline"><em>α</em></span></td>
<td align="center"><span class="math inline"><em>G</em><em>a</em>(<em>ν</em><sub>1</sub>, <em>ν</em><sub>2</sub>)</span></td>
<td>prior number of clusters K</td>
</tr>
<tr class="odd">
<td align="center"><code>'3'</code></td>
<td align="center">Pitman--Yor</td>
<td align="center"><span class="math inline">(<em>α</em>, <em>σ</em>)</span></td>
<td align="center">none</td>
<td>prior number of clusters K / Variance around K</td>
</tr>
<tr class="even">
<td align="center"><code>'4'</code></td>
<td align="center">Pitman--Yor</td>
<td align="center"><span class="math inline">(<em>α</em>, <em>σ</em>)</span></td>
<td align="center">spike and slab</td>
<td>prior number of clusters K</td>
</tr>
</tbody>
</table>

More precisely on the values of the parameter `DRtype`:

-   No value specified, then the model will use standard dimension reduction approach.
-   `2` is the dimension reduction with the Dirichlet process and prior distribution on the concentration parameter *α*
-   `3` is the dimension reduction with the Pitman--Yor process with fixed parameters *α* and *σ*.

-   `4` is the dimension reduction with the Pitman--Yor process with spike and slab prior on *α* and *σ*.

For the example, in the **Bauges plant dataset ** , we take as a prior number of plant functional groups which is 16. If we want to take the Dirichlet process, then we use the following parameters for the `reductList`. Here, `DRtype =2`. The choice of *N* is equal to *S* is natural, as we consider any possible clustering and K is prior number of clusters.

``` r
rl2  <- list(r = r_reduct, N = S ,DRtype="2", K=16)  #prior is number of plant functional groups
```

For the Pitman--Yor process. If we want to take the Dirichlet process, then we use the following parameters for the `reductList`.Here, `DRtype =3` and we specify only `K`. The choice of *N* is the same.

``` r
rl3  <- list(r = r_reduct, N = S ,DRtype="3", K=16)  #prior is number of plant functional groups
```

For the Pitman--Yor process. If we want to take the Dirichlet process, then we use the following parameters for the `reductList`.Here, `DRtype =3` and we specify `K` and variance `V`. The choice of *N* is the same.

``` r
rl4  <- list(r = r_reduct, N = S ,DRtype="3", K=16, V=20)  #prior is number of plant functional groups
```

For instance, we fit the model with the first specification

``` r
ml2   <- list(ng = iterations, burnin = burn_period, typeNames = 'PA', reductList = rl2,PREDICTX = F)
fit2<-gjam(formula, xdata = X_data, ydata = Y_data, modelList = ml2)
```

We get the `fit2` gjam object. This object is the same as the one from the original model, except several additional parameters. Firstly, we can plot the prior distribution of the alpha parameter and get expected number of clusters for this prior distribution of *α* (We suggest to test it, due to stochastic nature of hyperparameter tunning). *α* ∼ Ga(*ν*<sub>1</sub>, *ν*<sub>2</sub>). Where *ν*<sub>1</sub> is shape and *ν*<sub>2</sub> rate parameters corresponingly in the modelList `otherpar` list in the `reductList` in `modelList` of gjam object `fit2`.

``` r
alpha <- concentation_parameter(fit2)
alpha$plot_prior
alpha$plot_all
```

The function `concentation_parameter` is made for visualization purposes, it takes as the input the object gjam and return the values

-   plot\_prior plot the prior distribution for *α* parameter
-   plot\_all plot for the prior and posterior for *α* parameter
-   alpha vector of samples from prior distribution of *α* parameter
-   E\_k expected number of clusters for this parameter of alpha

Then we need to estimate the clusters. The cluster estimation is done by the function gjam\_cluster, which would return optimal cluster, that summarize posterior cluster distribution. We use the method proposedby Wade et al., 2018 where we use the clusters

For building prior and posterior cluster distribution we can use the function...

We can also consider applcation of the dimension reduction with type 3 We can use the variance to have less peaky distribution
