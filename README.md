# CancerFitness

This code implements methods for identifying interactions between mutations and CNAs within a cancer gene using log-linear regression model.

To determine the significance of the co-occurrence of a pair of somatic mutation and CNAs within a gene across cancer types, we generated a model using a log-linear regression with a Poisson function using the MASS package in R (Ripley, 2013). 

See Park et al for details.

Note: the R pacakge used to generate the results for the manuscript is here.
