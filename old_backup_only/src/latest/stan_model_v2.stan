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

functions{
 // real mujx(vector beta, real fj, real dj){
  //  return beta[1] + beta[2] * log(fj+1) + beta[3] * log(dj+1)
  //}
  //vector RAN(vector JCUM,int i){
 //   return(JCUM[i]+1):JCUM[i+1]);
  //}
  int NX(int n, int x, vector jc){
    for(i in 1:x){
      if(n <= jc[i]){
        return i;
      }
    }
    return 1;
  }
}

data {
  int<lower=0> N;
  int<lower=0> X;
  vector[X] J;
  vector[X+1] JCUM;
  int f[N];
  int M[N];
  int d[N];
  vector[N] S;
  int StoX[N];
}
transformed data {
  real mu_beta = 10;
  real sigma_beta = 100;
  real a_sigma_b = .5;
  real b_sigma_b = .5;
  
  real mu_alpha = 0;
  real sigma_alpha= 100;
  real a_delta = .5;
  real b_delta = .5;
  real mu_a = 0;
  real sigma_a =  10;
  real k_b = 1;
  real theta_b = 500;
}


parameters {
  vector[3] beta;
  real<lower=0,upper=30> sigmaS_b;
  vector<lower=0, upper=1>[X] logit_b;
  // Reaction time
  real<lower=0.01> a_t;
  real<lower=0> b_t;
  real<lower=0> sigmaS_delta;
  real alpha;
  vector[X] alphaX;
  vector<lower=0>[X] tauSX;
}
transformed parameters{
  vector<lower=0, upper=1>[X] b;
  for(i in 1:X){
    b[i] = inv_logit(logit_b[i]);
  }
}


model {
//  for(i in 1:3){
//    target += normal_lpdf(beta[i] | mu_beta, sigma_beta);
//  }
//  target += inv_gamma_lpdf(sigmaS_b | a_sigma_b, b_sigma_b);
//  for(i in 1:X){
//    target += normal_lpdf(logit_b[i]|beta[1] + beta[2] * log(f[i]+1) + beta[3] * log(f[i]+1)^2, sqrt(sigmaS_b));
    //target += log_inv_logit(b[i]) + log1m_inv_logit(b[i]);
//  }
//  for(i in 1:X){
//    target += binomial_lpmf(M[i]| f[i],b[i]);
//  }
  // Reaction time modeling
  a_t ~ lognormal(mu_a, sigma_a);
  b_t ~ gamma(k_b, theta_b);
  alpha ~ normal(mu_alpha, sigma_alpha);
  sigmaS_delta ~ inv_gamma(a_delta,b_delta);
  
  for(i in 1:X){
    alphaX[i] ~ normal(alpha,sqrt(sigmaS_delta));
    tauSX[i] ~ inv_gamma(a_t,b_t);
  }
  for(i in 1:N){
    S[i] ~ normal(alphaX[StoX[i]],sqrt(tauSX[StoX[i]]));
  }
}
