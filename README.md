ReadMe
================
Daria Bystrova
30/04/2020

Results
-------

To reproduce the results, there are two possibilites:

1.  Run the fit\_script.R, which will save the models in the folder "PCA analysis"
2.  Download the models from the following link [linked models](https://drive.google.com/open?id=1sRu7Q7rJ4aIYp-YCKOa_igA5Ofzc7tRH)

Then, load the models in the folder "PCA analysis/r5\_25" on your local repository

The results are produced by the "analysis.R" script

How to use the functions
------------------------

In the dimension reduction approach in the GJAM model, the species are clustered in their dependence behavior. If we have prior knowledge of the number of groups of species that share the same behavior, we could use this prior knowledge in the dimension reduction in the GJAM model. For this, we need to specify the prior number of clusters and (possibly the variance around this number)
We can use one of two possible regimes, which are based on the clustering with the Dirichlet process and Pitman--Yor process. They differ in clustering flexibility.
