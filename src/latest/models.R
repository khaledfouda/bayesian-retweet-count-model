require(logitnorm)
require(lognorm)
require(statmod)
require(mvtnorm)
require(invgamma)
#------------------------------------
# INPUT DATA from observations
model_input = readRDS('../../data/model_input.rds')
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
W = matrix(c(rep(1,N),log(f+1),log(d+1)),N,3,byrow = FALSE)

#S = log(S)
#f = log(f+1)
#------------------------------------------
# Retweet graph constant parameters. ******
#------------------------------------------
# Each of the betas has a normal distribution with the following
#     means and standard deviations
mu.beta = c(0,0,0)
S.beta = diag(c(5^2,5^2,5^2),3,3)
WTWS.inv = solve(t(W)%*%W + solve(S.beta))
# The variance of the retweet probability(sigmaS.b) has IG distribution with
# the following parameters
a.sigma.b = .5
b.sigma.b = .5
#------------------------------------------
# Reaction time constant parameters. ******
#------------------------------------------
# Alpha has a normal prior with the following param
mu.alpha = 0
sigma.alpha = 5
# Sigma squared delta has IG prior with the folowwing
a.delta = .5
b.delta = .5
# a_t has logNorm prior with the following param
mu.a = 0
sigma.a = 5 
# b_t has gamma prior the following parameters
k.b = 1
theta.b = 1/5
eps = 1e-323
#--------------------------
# Beta
prior_beta = function(){
  c(rmvnorm(1, mu.beta, S.beta))
}
post_beta = function(sigmaS.b, b){
  #print("Nah")
  covar = sigmaS.b * WTWS.inv
  mu = c((covar %*% t(W)) %*% logit(b)) / sigmaS.b
  return(c(rmvnorm(1, mu, covar)))
}
#-------------------------------------------------
# Sigma Squared b
prior_sigmaS.b = function() rinvgamma(1,a.sigma.b,rate=b.sigma.b)
post_sigmaS.b = function(b, beta){
  #print("DONT")
  ad = a.sigma.b + N/2
  mu = c(beta %*% t(W))
  bd = b.sigma.b + .5 * sum( (logit(b)-mu)^2)
  return(rinvgamma(1, ad,rate= bd))
}
#--------------------------
# b_j^x
prior_b = function(beta, sigmaS.b) {
  mu = c(beta %*% t(W))
  bj = rep(2,N)
  for(i in 1:N){
    while(bj[i]<eps|bj[i]>=1)
    bj[i] = rlogitnorm(1,mu[i],sqrt(sigmaS.b))
  }
  return(bj)
}

trans_b = function(beta, sigmaS.b){
  #print("HUH")
  mu = c(beta %*% t(W))
  bJ = sapply(mu, function(d){ r=0;while(r<eps) r=rlogitnorm(1,d,sigmaS.b);return(r)} )
  return(bJ)
}
post_b = function(beta, sigmaS.b,b){
  #print("POST?")
  mu = c(beta %*% t(W))
  den = (b^M) * (1-b)^(f-M) * exp( (-1/(2*sigmaS.b)) * (logit(b)-mu)^2)
  den = sapply(den, function(d)max(d,eps))
  return(den)
}
trans_b.j.x = function(mui, sigmaS.b){
  bj = rlogitnorm(1,mu[i],sqrt(sigmaS.b))
  ifelse(bj>eps, return(bj), trans_b.j.x(mui, sigmaS.b))
}
post_b.j.x = function(bi, mui, mi,fi, sigmaS.b){
  den = (bi^mi) * (1-bi)^(fi-mi) * exp( (-1/(2*sigmaS.b)) * (logit(bi)-mui)^2)
  #den = replace(den, den<eps, eps) 
  return(max(den, eps))
}
#---------
# a_t
prior_a.t = function() {
  a = rlnorm(1,mu.a,sigma.a)
  ifelse(a>.01&a<10, return(a), prior_a.t())
}
post_a.t = function(x, b.t, tauS.X){
  den = exp(-(log(x)^2)/(2*sigma.a^2)) * prod(sqrt(tauS.X)^(-x) * ((b.t^x)/gamma(x)))
  return(den)
}
trans_a.t = function(prev){
  a = rnorm(1,prev, .1)
  ifelse(a>.01&a<10, return(a), trans_a.t(prev))
}
#------------------------
# b_t
prior_b.t = function() 
  {
  bt = rgamma(1, k.b, rate=theta.b)
  ifelse(bt>0&bt<10, return(bt), prior_b.t())
}
post_b.t = function(a.t, tauS.X){
  kd = k.b + X * a.t
  thetad = theta.b + sum(1/sqrt(tauS.X))
  bt = rgamma(1, kd, rate=thetad)
  while(bt < 0|bt>15){bt = rgamma(1, kd, rate=thetad)}
  return(bt)
}
#-----------------------------
# alpha
prior_alpha = function() rnorm(1,mu.alpha, sigma.alpha)
post_alpha = function(sigmaS.delta, alpha.X){
  mu = (X + sigmaS.delta* (sigma.alpha)^(-2))^(-1) *sum(alpha.X)
  sig =(X + sigmaS.delta* (sigma.alpha)^(-2))^(-1) * sigmaS.delta
  return(rnorm(1, mu, sqrt(sig)))
}
#-----------------------------------------
# sigma sq delta
prior_sigmaS.delta = function() {
  sd = rinvgamma(1,a.delta,rate=b.delta)
  while(sd < 0){sd=rinvgamma(1,a.delta,rate=b.delta)}
  return(sd)
}
post_sigmaS.delta = function(alpha, alpha.X){
  ad = a.delta + X/2
  bd = b.delta + .5 * sum((alpha.X-alpha)^2)
  sd = rinvgamma(1, ad, rate=bd)
  while(sd <0){sd=rinvgamma(1, ad, rate=bd)}
  return(sd)
}
#---------------------------------
# alpha x
prior_alpha.x = function(alpha, sigmaS.delta) rnorm(X, alpha, sqrt(sigmaS.delta))

post_alpha.x = function(tauS.X, sigmaS.delta){
  alx = rep(NA,X)
  
  for(i in 1:X){
  mu = (J[i]-1 + tauS.X[i]*sigmaS.delta^-1)^(-1) * sum(S[(JCUM[i]+1):JCUM[i+1]])
  sigs = (J[i]-1 + tauS.X[i]*sigmaS.delta^-1)^(-1) *  tauS.X[i]
  alx[i] = rnorm(1, mu, sqrt(sigs))
  }
  return(alx)
}
#-----------------------------------------------
# TauSq
prior_tauS.x = function(a.t, b.t) {
  rinvgamma(X, a.t, b.t)
}

post_tauS.x = function(a.t, b.t, alpha.x){
  ad = a.t + (J-1)/2
  ta = rep(NA,X)
  for(i in 1:X){
    bd = b.t + .5 * sum((S[(JCUM[i]+1):JCUM[i+1]] - alpha.x[i])^2)
    ta[i] = rinvgamma(1, ad[i],bd)
    }
  return(ta)
}

#------------------------------------------------
# Retweet graph prior distributions ***********
#------------------------------------------------


prior_M.j.x = function(f.j.x, b.j.x) rbinom(1,f.j.x,b.j.x)
#------------------------------------------
# Reaction time prior distribution. ******
#------------------------------------------

prior_S.j.x = function(alpha.x, tauS.x){
  rnorm(1, alpha.x, sqrt(tauS.x))
}
#------------------------------------------------
# Retweet graph conditional posteriors ***********
#------------------------------------------------

# Note: b.j.x is sample by MH
# the posterior of it returns a density instead of a sample.
# the sample is taken from the proposal


#------------------------------------------
# Reaction time conditional posteriors. ******
#------------------------------------------

# Note a.t is sampled using MH
# The posterior returns a density instead of a sample.
# the sample is taken from the proposal


#----------------------------------------------------------
#----------------------------------------------------------