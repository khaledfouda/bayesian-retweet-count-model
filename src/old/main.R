require(logitnorm)
#------------------------------------------
# Retweet graph constant parameters. ******
#------------------------------------------
# Each of the betas has a normal distribution with the following
#     means and variances
mu.beta0 =0
mu.betaf =0
mu.betad =0
sigmaS.beta0 = 100
sigmaS.betaf = 100
sigmaS.betad = 100
# The variance of the retweet probability has IG distribution with
# the following parameters
a.sigma.b = .5
b.sigma.b = .5
#------------------------------------------------
# Retweet graph prior distributions ***********
#------------------------------------------------
prior_M.j.x = function(f.j.x, b.j.x) rbinom(1,f.j.x,b.j.x)
prior_b.j.x = function(mu.j.x, sigmaSb) rlogitnorm(1, mu.j.x,sigmaSb) 
prior_mu.j.x = function(b0, bf, bd, f.j.x, d.j.x) {
  return(b0 + bf*log(f.j.x+1)+bd*log(d.j.x+1))
}
prior_beta