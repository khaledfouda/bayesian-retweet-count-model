//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  int<lower=0> X;
  vector[X] J;
  vector[X] JCUM;
  vector[N] f;
  vector[N] d;
  vector[N] S;
  vector[N] M;
}
transformed data {
  real mu_beta = 0;
  real sigma_beta = 5;
  real a_sigma_b = .5;
  real b_sigma_b = .5;
  
  real mu_alpha = 0;
  real sigma_alpha= 100;
  real a_delta = .5;
  real b_delta = .5;
  real mu_a = 0;
  real sigma_a = 10 ;
  real k_b = 1;
  real theta_b = 500;
  
}
// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[3] beta;
  real<lower=0,upper=30> sigmaS_b;
  real<lower=0, upper=1> b_j_x;
  real M_j_x;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  beta ~ 
  y ~ normal(mu, sigma);
}

