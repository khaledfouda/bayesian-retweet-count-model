
//options(mc.cores = parallel::detectCores())
//rstan_options(auto_write = TRUE)
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

/*data {
  int<lower=0> J; // #schools
  real y[J]; // estim trt eff
  real<lower =0> sigma[J]; // stderr effec estim
}

parameters {
  real mu; // pop trt eff
  real<lower=0> tau; //  stderr trt eff
  vector[J] eta; //dev from mu by school
}

transformed parameters {
  vector[J] theta = mu + tau * eta; // school trt effects
}

model {
  target += normal_lpdf(eta|0,1); // prior log density
  target += normal_lpdf(y|theta, sigma); // log likelihood
}*/

data {
  int<lower=0> J;
  real y[J];
  real<lower=0> sigma[J];
}
parameters {
  real mu;
  real<lower=0> tau;
  real theta[J];
}
model {
  mu ~ normal(0,5);
  tau ~ cauchy(0, 5);
  theta ~ normal(mu, tau);
  y ~ normal(theta, sigma);
}
