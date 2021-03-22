require(logitnorm)
require(lognorm)
require(statmod)
require(mvtnorm)

#------------------------------------
# INPUT DATA from observations
model_input = readRDS('../data/model_input.rds')
X = model_input$X # the total number of root users
J = model_input$J # where J_x is the number of retweets user x got including
# them.
JCUM = model_input$JCUM
N = model_input$N # the total number of observed reaction time for all training and
# predictions tweets.
f = model_input$f # the number of followers of user
d = model_input$d
S = model_input$S
M = model_input$M

#------------------------------------------
# Retweet graph constant parameters. ******
#------------------------------------------
# Each of the betas has a normal distribution with the following
#     means and standard deviations
mu.beta =0
sigma.beta = 100
# The variance of the retweet probability(sigmaS.b) has IG distribution with
# the following parameters
a.sigma.b = .5
b.sigma.b = .5
#------------------------------------------
# Reaction time constant parameters. ******
#------------------------------------------
# Alpha has a normal prior with the following param
mu.alpha = 0
sigma.alpha = 100
# Sigma squared delta has IG prior with the folowwing
a.delta = .5
b.delta = .5
# a_t has logNorm prior with the following param
mu.a = 0
sigma.a = 10 
# b_t has gamma prior the following parameters
k.b = 1
theta.b = 500
#------------------------------------------------
# Retweet graph prior distributions ***********
#------------------------------------------------
prior_beta = function(){
  c(rmvnorm(1, rep(mu.beta,3), diag(rep(sigma.beta,3))))
}
prior_sigmaS.b = function() rinvgauss(1,a.sigma.b,b.sigma.b)
mu.j.x = function(beta, f.j.x, d.j.x) {
  return(beta[1] + beta[2]*log(f.j.x+1)+beta[3]*log(d.j.x+1))
}
prior_b.j.x = function(mu.j.x, sigmaS.b) rlogitnorm(1, mu.j.x,sqrt(sigmaS.b)) 
prior_M.j.x = function(f.j.x, b.j.x) rbinom(1,f.j.x,b.j.x)
#------------------------------------------
# Reaction time prior distribution. ******
#------------------------------------------
prior_alpha = function() rnorm(1,mu.alpha, sigma.alpha)
prior_sigmaS.delta = function() rinvgauss(1,a.delta,b.delta)
prior_a.t = function() rlnorm(1, mu.a,sigma.a)
prior_b.t = function() rgamma(1, k.b, theta.b)
prior_alpha.x = function(alpha, sigmaS.delta) rnorm(1, alpha, sqrt(sigmaS.delta))
prior_tauS.x = function(a.t, b.t) rinvgauss(1, a.t, b.t)

prior_S.j.x = function(alpha.x, tauS.x){
  rlnorm(1, alpha.x, sqrt(tauS.x))
}
#------------------------------------------------
# Retweet graph conditional posteriors ***********
#------------------------------------------------
post_beta = function(sigmaS.b, b){
  # The conditional joint posterior of the three betas.
  # The follow a multivariate normal distribution with parameters mu and C 
  # as defined below.
  N1 = N + sigmaS.b * sigma.beta^(-2)
  E = sum(log(f+1)*log(d+1))
  D = sum(log(d+1))
  D2 = sum(log(d+1)^2) + sigmaS.b * sigma.beta^(-2)
  W = sum(log(f+1))
  W2 = sum(log(f+1)^2) + sigmaS.b * sigma.beta^(-2)
  Y0 = sum(log(b+1))
  YF = sum(log(b+1)*log(f+1))
  Yd = sum((log(b+1)^2)*log(d+1)) + sigmaS.b * sigma.beta^(-2)
  C = matrix(c(
    N1, W, D,
    W, W2, E,
    D, E, D2),nrow=3,byrow = TRUE)
  C = sigmaS.b * solve(C)
  mu = C %*% t(c(Y0, YF, Yd))
  return(rmvnorm(1, mu, C))
}
post_sigmaS.b = function(b, beta){
  ad = a.sigma.b + N/2
  mu = c()
  for(i in 1:N){
    mu = c(mu, mu.j.x(beta, f[i],d[i]))
  }
  bd = b.sigma.b + .5 * sum( (logit(b)-mu)^2)
  return(rinvgauss(1, ad, bd))
}
# Note: b.j.x is sample by MH
# the posterior of it returns a density instead of a sample.
# the sample is taken from the proposal
post_b.j.x = function(x, f.j.x, b.j.x, d.j.x, beta, sigmaS.b){
  mu = mu.j.x(beta, f.j.x, d.j.x)
  den =  dbinom(x,f.j.x,b.j.x) *  dlogitnorm(x, mu,sqrt(sigmaS.b))
  return(den)
}
propos_b.j.x = function(beta, f.j.x, d.j.x, sigmaS.b){
  mu = mu.j.x(beta, f.j.x, d.j.x)
  return(rlogitnorm(1, mu, sqrt(sigmaS.b)))
}
#---------

#------------------------------------------
# Reaction time conditional posteriors. ******
#------------------------------------------
post_alpha = function(sigmaS.delta, alpha.X){
  mu = (X + sigmaS.delta/sigma.alpha)^(-2) *sum(alpha.X)
  sig = (X + sigmaS.delta/sigma.alpha)^(-2) * sigmaS.delta
  return(rnorm(1, mu, sqrt(sig)))
}
post_sigmaS.delta = function(alpha, alpha.X){
  ad = a.delta + X/2
  bd = b.delta + .5 * sum((alpha.X-alpha)^2)
  return(rinvgauss(1, ad, bd))
}
# Note a.t is sampled using MH
# The posterior returns a density instead of a sample.
# the sample is taken from the proposal
post_a.t = function(x, a.t, b.t, tauS.X){
  den = dlnorm(x, mu.a,sigma.a) * prod(dinvgauss(sqrt(tauS.X), a.t, b.t))
  return(den)
}
propos_a.t = function(prev){
  return(rlnorm(1,prev, .2))
}
post_b.t = function(a.t, tauS.X){
  kd = k.b + X * a.t
  thetad = (theta.b^(-1) + sum(tau.X^(-2)) )^(-1)
  return(rgamma(1, kd, thetad))
}
post_alpha.x = function(M.x, tauS.x, sigmaS.delta, S.x){
  mud = (M.x + tauS.x / sigmaS.delta)^(-1) * sum(log(S.x))
  sigd = (M.x + tauS.x / sigmaS.delta)^(-1) * tauS.x
  return(rnorm(1, mud, sqrt(sigd)))
}
post_tauS.x = function(a.t, b.t, M.x, S.x, alpha.x){
  ad = a.t + M.x/2
  bd = b.t + .5 * sum((log(S.x)-alpha.x)^2)
  return(rinvgauss(1, ad, bd))
}
#----------------------------------------------------------
#----------------------------------------------------------