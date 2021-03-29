source('./models.R')

# params
#------------
betas = c(NA,NA,NA)
sigmaS.b = NA
b.J.X = vector(length = sum(J))
#---------------
alpha = NA
sigmaS.delta = NA
a.t = NA
b.t = NA
alpha.X = vector(length = X)
tauS.X = vector(length = X)
#-------------
n = 500 # number of samples to be taken from each parameter
#---------
sam_betas = array(NA, dim=c(n,3))
sam_sigmaS.b = vector(length=n)
sam_b = array(NA, dim=c(n,sum(J)))
#----------------------
sam_alpha= vector(length = n)
sam_sigmaS.delta = vector(length = n)
sam_a.t = vector(length = n)
sam_b.t = vector(length = n)
sam_alpha = array(NA, dim=c(n,X))
sam_tauS = array(NA, dim=c(n,X))
#-------------------------------------------------
#------------------------------------------------
# Initialization ::
# Intialize each according to its prior distribution. 
betas = prior_beta()
sigmaS.b = prior_sigmaS.b()
for(i in 1:N){
  muj = mu.j.x(betas, f[i], d[i])
  b.J.X[i] = prior_b.j.x(muj, sigmaS.b)
}
alpha = prior_alpha()
sigmaS.delta = prior_sigmaS.delta()
a.t = prior_a.t()
b.t = prior_b.t()
for(i in 1:X){
  alpha.X[i] = prior_alpha.x(alpha, sigmaS.delta)
}
for(i in 1:X){
  tauS.X[i] = prior_tauS.x(a.t, b.t)
}
# add each of the sampled values to the parameter's list.
sam_betas[1,] = betas
sam_sigmaS.b[1] = sigmaS.b
sam_b[1,] = b.J.X
sam_alpha[1] = alpha
sam_sigmaS.delta[1] = sigmaS.delta
sam_a.t[1] = a.t
sam_b.t[1] = b.t
sam_alpha[1,] = alpha.X
sam_tauS[1,] = tauS.X
#--------------------------------------------------
for(t in 1:n){
  # Here for each time step, we sample each parameter once,
  # Note that a.t and b.j.x are sampled using MH
  betas = post_beta(sigmaS.b, b.J.X)
  sigmaS.b = post_sigmaS.b(b.J.X, betas)
  for(i in 1:N){
  	prob = 0
  	while(runif(1) > prob){
  	  # newp ~ Q(p|p')
  		newp = propos_b.j.x(betas, f[i], d[i], sigmaS.b)
  		# probability of acceptance = p(p)/p(p')
  		prob = post_b.j.x(newp,M[i], f[i], d[i], betas, sigmaS.b) /
  		post_b.j.x(b.J.X[i],M[i], f[i], d[i], betas, sigmaS.b)
  	}
  	b.J.X[i] = newp
  }
  alpha = post_alpha(sigmaS.delta, alpha.X)
  sigmaS.delta= post_sigmaS.delta(alpha, alpha.X)
  prob = 2
  while(runif(1) > prob){
    newa = propos_a.t(a.t)
    prob = post_a.t(newa,a.t, b.t, tauS.X) / post_a.t(a.t,a.t,b.t,tauS.X)
  }
  a.t = newa 
  b.t = post_b.t(a.t,tauS.X)
  alpha.X[1] = post_alpha.x(M[1], tauS.X[1], sigmaS.delta, S.x[1:JCUM[1]])
  for(i in 2:X){
    alpha.X[i] = post_alpha.x(M[i], tauS.X[i], sigmaS.delta,
                            S.x[(JCUM[i-1]+1):JCUM[i]])
  }
  tauS.X[1] = post_tauS.x(a.t, b.t, M[1],S.x[1:JCUM[1]],
                          alpha.X[1])
  for(i in 2:X){
    tauS.X[i] = post_tauS.x(a.t, b.t, M[i],S.x[(JCUM[i-1]+1):JCUM[i]],
                            alpha.X[i])
  }
  #----------------
  # we append the new samples to the arrays.
  sam_betas[t,] = betas
  sam_sigmaS.b[t] = sigmaS.b
  sam_b[t,] = b.J.X
  sam_alpha[t] = alpha
  sam_sigmaS.delta[t] = sigmaS.delta
  sam_a.t[t] = a.t
  sam_b.t[t] = b.t
  sam_alpha[t,] = alpha.X
  sam_tauS[t,] = tauS.X
}





