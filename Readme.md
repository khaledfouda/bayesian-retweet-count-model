# Abstract
-----
Here we go over the process of building a hierarchical Bayesian model to model the future
number of retweets for a given tweet. We begin by giving a summary on the important
methods and theories used in our study. In section two we get introduced to the data and
choose appropriate likelihood distribution. In section three we build the graphical model
while defining the conjugate priors and computing the posterior's parameters. In the last
section we implement the model in R-language and sample from the conditional posterior
distributions. This study is accompanied with the code files for section 4 and the exploratory
data analysis done in section 2.

_____________________________________________________

- All explainations are model derivations as well as the results are explained in the report: [Probabilistic Graphical Models Case study: Hierarchical
Bayesian Model](https://github.com/khaledfouda/bayesian-retweet-count-model/blob/main/case_study__Hierarchical_Bayesian_Model_A4.pdf).  
- The project is based on [Tauhid Zaman, Emily B. Fox, and Eric T. Bradlow. A bayesian approach for predicting the
popularity of tweets](https://arxiv.org/abs/1304.6777v3).  
- Source code structure:
  - The notebook EDA_v2.ipynb reads, prepares and explores the data.
  - All the priors, posteriors and likelihoods are defined in models.R.
  - The MCMC sampling loop is defined in MCMC_trainer.R.
  - Gelman Rubin diagnostics are defined in gelmanRubinDiagnostic.R.

_________________________
