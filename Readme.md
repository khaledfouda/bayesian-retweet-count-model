# Probabilistic Graphical Models Case study: Hierarchical Bayesian Model  
-----
Here we go over the process of building a hierarchical Bayesian model to model the future
number of retweets for a given tweet. We begin by giving a summary of the important
methods and theories used in our study. In section two we get introduced to the data and
choose appropriate likelihood distribution. In section three we build the graphical model
while defining the conjugate priors and computing the posteriors' parameters. In the last
section, we implement the model in R-language and sample from the conditional posterior
distributions. This study is accompanied by the code files for section 4 and the exploratory
data analysis done in section 2.

_____________________________________________________

**The main report: [Probabilistic Graphical Models Case study: Hierarchical
Bayesian Model](https://github.com/khaledfouda/bayesian-retweet-count-model/blob/main/case_study__Hierarchical_Bayesian_Model_A4.pdf).**   

-------------

- The project is based on [Tauhid Zaman, Emily B. Fox, and Eric T. Bradlow. A bayesian approach for predicting the
popularity of tweets](https://arxiv.org/abs/1304.6777v3).  
- Source code structure:
  - The notebook [EDA_v2.ipynb](https://github.com/khaledfouda/bayesian-retweet-count-model/blob/main/src/EDA_v2.ipynb) reads, prepares and explores the data. A pdf version is avaiable  at [EDA_data_cleaning.pdf](https://github.com/khaledfouda/bayesian-retweet-count-model/blob/main/report/EDA_data_cleaning.pdf)
  - All the priors, posteriors and likelihoods are defined in [models.R](https://github.com/khaledfouda/bayesian-retweet-count-model/blob/main/src/models.R).
  - The MCMC sampling loop is defined in [MCMC_trainer.R](https://github.com/khaledfouda/bayesian-retweet-count-model/blob/main/src/MCMC_trainer.R).
  - Gelman Rubin diagnostics are defined in [gelmanRubinDiagnostic.R](https://github.com/khaledfouda/bayesian-retweet-count-model/blob/main/src/gelmanRubinDiagnostic.R).

_________________________
