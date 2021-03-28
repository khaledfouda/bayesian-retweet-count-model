source('./models.R')

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
C = 5
summ_a.t = array(NA,c(C,3))
summ_b.t = array(NA,c(C,3))
RNGkind(sample.kind = "Rounding")

#-----------------------------------
for(c in 1:C){
  set.seed(12)
  #-----1
  #beta = prior_beta()
  beta = c(1,1,1)
  cat(beta,'\n')
  #sigmaS.b = prior_sigmaS.b()
  sigmaS.b = .1
  cat(sigmaS.b,'\n')
  #b = prior_b(beta, sigmaS.b)
  b = rep(.5,N)
  cat(max(b),min(b),'\n')
  
  #------2
  # alpha = prior_alpha()
  # sigmaS.delta = prior_sigmaS.delta()
  # a.t = prior_a.t()
  # b.t = prior_b.t()
  # alpha.X = prior_alpha.x(alpha, sigmaS.delta)
  # tauS.X = prior_tauS.x(a.t, b.t)
  #------1
  # add each of the sampled values to the parameter's list.
  sam_betas[1,] = beta
  sam_sigmaS.b[1] = sigmaS.b
  sam_b[1,] = b
  #---------2
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
    b[Accepted] = newp[Accepted]}
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
  }
    #----------------------2
    # alpha = post_alpha(sigmaS.delta, alpha.X)
    # sigmaS.delta= post_sigmaS.delta(alpha, alpha.X)
    # prob = 0
    # while(runif(1) > prob){
    #   newa = trans_a.t(a.t)
    #   prob = post_a.t(newa, b.t, tauS.X) / post_a.t(a.t,b.t,tauS.X)
    # }
    # a.t = newa 
    # b.t = post_b.t(a.t,tauS.X)
    # alpha.X = post_alpha.x(tauS.X, sigmaS.delta)
    # 
    # tauS.X = post_tauS.x(a.t, b.t, alpha.X)
    #----------------1
    #we append the new samples to the arrays.
    #sam_betas[t,] = beta
    #sam_sigmaS.b[t] = sigmaS.b
    #sam_b[t,] = b
    #-------------------------2
    # sam_alpha[t] = alpha
    # sam_sigmaS.delta[t] = sigmaS.delta
    # sam_a.t[t] = a.t
    # sam_b.t[t] = b.t
    # sam_alpha[t,] = alpha.X
    # sam_tauS[t,] = tauS.X
  }
  summ_a.t[c,] = c(mean(sam_a.t[(n/2):n]), var(sam_a.t[(n/2):n]),
                   sum(sam_a.t[(n/2):n]))
  summ_b.t[c,] = c(mean(sam_b.t[(n/2):n]), var(sam_b.t[(n/2):n]),
                   sum(sam_b.t[(n/2):n]))
}
#sam_a.t

#length(sam_a.t)
#summ_a.t = array(NA,c(5,3))
#summ_a.t[1,] = c(mean(sam_a.t[250:500]), var(sam_a.t[250:500]),
#                 sum(sam_a.t[250:500]))

#summ_b.t = array(NA,c(5,3))
#summ_b.t[1,] = c(mean(sam_b.t[250:500]), var(sam_b.t[250:500]),
#                 sum(sam_b.t[250:500]))


#require(coda)


#summ_a.t

#ovrall_a.t = sum(summ_a.t[,3])/(5*250)
# ovrall_b.t = sum(summ_b.t[,3])/(5*250)
# 
# B_a.t = (250/4) * sum( (summ_a.t[,1]-ovrall_a.t)^2 ) 
# W_a.t = (1/5) * sum(summ_a.t[,2])
# 
# varEs_a.t = (249/250) * W_a.t + (1/250) * B_a.t
# 
# R_a.t = sqrt(varEs_a.t/W_a.t)
# 

