source('./models.R')
require(dplyr)
# params
#------------
# betas = c(NA,NA,NA)
# sigmaS.b = NA
# b.J.X = vector(length = sum(J))
# #---------------
# alpha = NA
# sigmaS.delta = NA
# a.t = NA
# b.t = NA
# alpha.X = vector(length = X)
# tauS.X = vector(length = X)
#-------------
n = 8000 # number of samples to be taken from each parameter
#---------
sam_beta = array(NA, dim=c(n/2,3))
sam_sigmaS.b = vector(length=n/2)
sam_b = array(NA, dim=c(n/2,sum(J)))
#----------------------
sam_alpha= vector(length = n/2)
sam_sigmaS.delta = vector(length = n/2)
sam_a.t = vector(length = n/2)
sam_b.t = vector(length = n/2)
sam_alpha.X = array(NA, dim=c(n/2,X))
sam_tauS.X = array(NA, dim=c(n/2,X))
#-------------------------------------------------
#------------------------------------------------
# Initialization ::
# Intialize each according to its prior distribution. 
C = 5
summ_beta = array(NA, c(3,C,3))
summ_sigmaS.b = array(NA,c(C,3))
summ_alpha= array(NA,c(C,3))
summ_sigmaS.delta = array(NA,c(C,3))
summ_a.t = array(NA,c(C,3))
summ_b.t = array(NA,c(C,3))
summ_alpha.X = array(NA,c(X,C,3))
summ_tauS.X = array(NA,c(X, C,3))

#RNGkind(sample.kind = "Rounding")

#-----------------------------------
for(c in 1:C){
  #set.seed(12)
  #-----1
  #beta = prior_beta()
  beta = runif(3, -5,5)
  #cat(beta,'\n')
  #sigmaS.b = prior_sigmaS.b()
  sigmaS.b = runif(1,.01,2)
  #cat(sigmaS.b,'\n')
  #b = prior_b(beta, sigmaS.b)
  b = runif(N,.4,.6)
  #cat(max(b),min(b),'\n')
  
  #------2
  alpha = prior_alpha()
  sigmaS.delta = prior_sigmaS.delta()
  a.t = prior_a.t()
  b.t = prior_b.t()
  alpha.X = prior_alpha.x(alpha, sigmaS.delta)
  tauS.X = prior_tauS.x(a.t, b.t)
  #------1
  # add each of the sampled values to the parameter's list.
  # sam_betas[1,] = beta
  # sam_sigmaS.b[1] = sigmaS.b
  # sam_b[1,] = b
  # #---------2
  # sam_alpha[1] = alpha
  # sam_sigmaS.delta[1] = sigmaS.delta
  # sam_a.t[1] = a.t
  # sam_b.t[1] = b.t
  # sam_alpha[1,] = alpha.X
  # sam_tauS[1,] = tauS.X
  #--------------------------------------------------
  for(t in 1:n){
    # Here for each time step, we sample each parameter once,
    # Note that a.t and b.j.x are sampled using MH
    #--------1
    beta = post_beta(sigmaS.b, b)
    #cat(beta,'\n')
    sigmaS.b = post_sigmaS.b(b, beta)
    #cat(sigmaS.b,'\n')
    mu = c(beta %*% t(W))
    #cat(mu[1:4],'\n')
    ### b
    newp = trans_b(beta, sigmaS.b)
    prob = post_b(beta, sigmaS.b, newp)/
      post_b(beta, sigmaS.b, b)
    Accepted = runif(N) < prob
    b[Accepted] = newp[Accepted]
    # 
    # for(i in 1:N){
    # 	prob = 0
    # 	counter = 0
    # 	while(runif(1) > prob & counter < 20){
    # 	  # newp ~ Q(p|p')
    # 	  counter = counter + 1
    # 	  newp = trans_b.j.x(mu[i], sigmaS.b)
    # 		#newp = propos_b.j.x(betas, f[i], d[i], sigmaS.b)
    # 		# probability of acceptance = p(p)/p(p')
    # 		prob = post_b.j.x(newp,mu[i],M[i],f[i], sigmaS.b) / 
    # 		  post_b.j.x(b[i],mu[i],M[i],f[i], sigmaS.b)
    # 	}
    # 	print(counter)
    # 	if(counter!=20) {b[i] = newp}
    # }
  
    #----------------------2
    alpha = post_alpha(sigmaS.delta, alpha.X)
    sigmaS.delta= post_sigmaS.delta(alpha, alpha.X)
    prob = 0
    while(runif(1) > prob){
      newa = trans_a.t(a.t)
      prob = post_a.t(newa, b.t, tauS.X) / post_a.t(a.t,b.t,tauS.X)
    }
    a.t = newa
    b.t = post_b.t(a.t,tauS.X)
    alpha.X = post_alpha.x(tauS.X, sigmaS.delta)

    tauS.X = post_tauS.x(a.t, b.t, alpha.X)
    if(t>n/2){
      j = t-n/2
      #----------------1
      #we append the new samples to the arrays.
      sam_beta[j,] = beta
      sam_sigmaS.b[j] = sigmaS.b
      sam_b[j,] = b
      #-------------------------2
      sam_alpha[j] = alpha
      sam_sigmaS.delta[j] = sigmaS.delta
      sam_a.t[j] = a.t
      sam_b.t[j] = b.t
      sam_alpha.X[j,] = alpha.X
      sam_tauS.X[j,] = tauS.X
      }
}
  summ_beta[,c,] = c(apply(sam_beta,2,mean), apply(sam_beta,2,var),
                     apply(sam_beta,2,sum))
  summ_sigmaS.b[c,] = c(mean(sam_sigmaS.b),var(sam_sigmaS.b),sum(sam_sigmaS.b))
  summ_alpha[c,] = c(mean(sam_alpha),var(sam_alpha),sum(sam_alpha))
  summ_sigmaS.delta[c,] = c(mean(sam_sigmaS.delta),var(sam_sigmaS.delta),
                            sum(sam_sigmaS.delta))
  summ_a.t[c,] = c(mean(sam_a.t),var(sam_a.t),sum(sam_a.t))
  summ_b.t[c,] = c(mean(sam_b.t),var(sam_b.t),sum(sam_b.t))
  summ_alpha.X[,c,] = c(apply(sam_alpha.X,2,mean), apply(sam_alpha.X,2,var),
                        apply(sam_alpha.X,2,sum))
  summ_tauS.X[,c,] = c(apply(sam_tauS.X,2,mean), apply(sam_tauS.X,2,var),
                       apply(sam_tauS.X,2,sum))
  cat("Chain ",c, " is complete.", '\n')
  
}
#sam_a.t

save.image(file = "../../output/workspace.RData")