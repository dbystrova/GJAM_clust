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

How to use the functions
------------------------

Firstly, we describe the initial GJAM model and then describe how to specify the parameters. Full description of the GJAM model is provided in the package documentation \[\] and [vignette](https://cran.r-project.org/web/packages/gjam/vignettes/gjamVignette.html)

Dimension reduction in GJAM model was proposed by \[Taylor et el 2018\]. Dimension reduction is used in two ways:

1.  When dataset contains more species than can be fitted given sample size n.\\

2.  When the number of species *S* is too large. \\

For the second case we need to specify the parameters for dimension reudction in reductList and add this parameters in the modelList.

**Table 2. Dimension reduction types**

| `DRtype` | Process type | Parameters |           Priors on parameters           | Parameters to specify                 |
|:--------:|:------------:|:----------:|:----------------------------------------:|---------------------------------------|
|   `'-'`  |   Dirichlet  |     *α*    |                   none                   | none                                  |
|   `'2'`  |   Dirichlet  |     *α*    | *G**a*(*ν*<sub>1</sub>, *ν*<sub>2</sub>) | prior number of clusters K            |
|   `'3'`  |  Pitman--Yor | (*α*, *σ*) |                   none                   | prior number of clusters K / Variance |
|   `'4'`  |  Pitman--Yor | (*α*, *σ*) |              spike and slab              | prior number of clusters K            |

**<sup>1</sup> For `'DA'` and `'CC'` data the second element of the partition is not zero, but rather depends on effort. There is thus zero-inflation. The default partition for each data type can be changed with the function `gjamCensorY` (see **specifying censored intervals**).
**<sup>2</sup> For `'CAT'` data species *s* has *K*<sub>*s*</sub> − 1 non-reference categories. The category with the largest *w*<sub>*i**s*, *k*</sub> is the '1', all others are zeros.
