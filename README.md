This code implements methods for identifying interactions between mutations and CNAs within a cancer gene using log-linear regression model.

To determine the significance of the co-occurrence of a pair of somatic mutation and CNAs within a gene across cancer types, we generated a model using a log-linear regression with a Poisson function using the MASS package in R (Ripley, 2013).

See Park et al, 2021 (submitted) for details.

Note 1: the R pacakge used to generate the results for the manuscript is here. Note 2 for the introduction of the input file:

Input file is required four types of sample types in each cancer type, including (i) Mut_CNAs, number of samples with both mutation and CNAs, (ii) NoMut_CNAs, number of samples with only CNAs, without mutation, (iii) Mut_WT, number of samples with only mutation without CNAs, and (iv) NoMut_WT, number of samples neither mutation nor CNAs. In the R code formulas, N indicates the number of samples across four conditions, including “mut” presents the number of samples with mutation event and “CNA gain (or loss)” presents the number of samples with CNA gain (or loss) event.
