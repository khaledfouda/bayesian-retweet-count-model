source('./models.R')
require(dplyr)

#------------
# n is the number of samples kept at each chain and C is the number of chains.
# m is the number of samples to be discarded.
n = 80000
m = 5000
C = 20
#-------------
# The following arrays will hold the samples drawn for each parameter per chain.
# We only start storing samples after reaching half the samples
sam_beta = array(NA, dim=c(n,3))
sam_sigmaS.b = vector(length=n)
sam_b = array(NA, dim=c(n,sum(J)))
#----------------------
sam_alpha= vector(length = n)
sam_sigmaS.delta = vector(length = n)
sam_a.t = vector(length = n)
sam_b.t = vector(length = n)
sam_alpha.X = array(NA, dim=c(n,X))
sam_tauS.X = array(NA, dim=c(n,X))
#-------------------------------------------------
#------------------------------------------------
# The following will hold the summaries across chains for each parameter
# i'll be storing the mean and variance of each set
MC = 2*C # The number of resulting MC after dividing the sample in two
summ_beta = array(NA, c(3,MC,2))
summ_sigmaS.b = array(NA,c(MC,2))
summ_alpha= array(NA,c(MC,2))
summ_sigmaS.delta = array(NA,c(MC,2))
summ_a.t = array(NA,c(MC,2))
summ_b.t = array(NA,c(MC,2))
summ_alpha.X = array(NA,c(X,MC,2))
summ_tauS.X = array(NA,c(X, MC,2))
summ_b = array(NA, c(N, MC, 2))
#-----------------------------------
mc = 1 # Index for the markov chains
for(c in 1:C){ # for each chain
  cat("Chain ", c, "is starting", "\n")
  # 1. we begin by giving initial values for our parameters.
  #-----##
  #beta = prior_beta()
  beta = runif(3, -5,5)
  #sigmaS.b = prior_sigmaS.b()
  sigmaS.b = runif(1,.01,2)
  #b = prior_b(beta, sigmaS.b)
  b = runif(N,.4,.6)
  #------##
  alpha = prior_alpha()
  sigmaS.delta = prior_sigmaS.delta()
  a.t = prior_a.t()
  b.t = prior_b.t()
  alpha.X = prior_alpha.x(alpha, sigmaS.delta)
  tauS.X = prior_tauS.x(a.t, b.t)
 #--------------------------------------------------
  for(t in 1:(n+m)){ # draw one sample from each parameter
    
    # Note that a.t and b.j.x are sampled using MH
    #-------#
    beta = post_beta(sigmaS.b, b)
    sigmaS.b = post_sigmaS.b(b, beta)
    mu = c(beta %*% t(W))
    for(i in 1:N)
    {
      prob = 0
      counter = 0
      while(runif(1) > prob & counter < 30){
        newp = trans_bj(mu[i], sigmaS.b)
        prob = post_bj(newp, mu[i], M[i], f[i], sigmaS.b) / 
          post_bj(b[i], mu[i], M[i], f[i], sigmaS.b)
        counter = counter + 1
      }
      if(counter != 30) b[i] = newp
    }
    
    #newp = trans_b(beta, sigmaS.b)
    #prob = post_b(beta, sigmaS.b, newp)/ post_b(beta, sigmaS.b, b)
    #Accepted = runif(N) < prob
    #b[Accepted] = newp[Accepted]
    #----------##
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
    #---------------------
    if(t%%500==0) cat( t," ... ") 
    # for the second half of the samples, we store them
    if(t>m){
      if(t == m+1) cat("\nBurn-up period for the current chain is over.\n")
      j = t-m
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
  #--------------------------------------------------
  # Updating Markov chains summaries for Diagonsis.
  #-----------------------------------------------------
  summ_beta[,mc,] = c(apply(sam_beta[1:(n/2),],2,mean),
                      apply(sam_beta[1:(n/2),],2,var))
  summ_beta[,(mc+1),] = c(apply(sam_beta[(n/2+1):(n),],2,mean),
                        apply(sam_beta[(n/2+1):(n),],2,var))
  
  summ_sigmaS.b[mc,] = c(mean(sam_sigmaS.b[1:(n/2)]),var(sam_sigmaS.b[1:(n/2)]))
  summ_sigmaS.b[(mc+1),] = c(mean(sam_sigmaS.b[(n/2+1):(n)]),
                           var(sam_sigmaS.b[(n/2+1):(n)]))

  summ_alpha[mc,] = c(mean(sam_alpha[1:(n/2)]),var(sam_alpha[1:(n/2)]))
  summ_alpha[(mc+1),] = c(mean(sam_alpha[(n/2+1):(n)]),
                        var(sam_alpha[(n/2+1):(n)]))

  summ_sigmaS.delta[mc,] = c(mean(sam_sigmaS.delta[1:(n/2)]),
                             var(sam_sigmaS.delta[1:(n/2)]))
  summ_sigmaS.delta[(mc+1),] = c(mean(sam_sigmaS.delta[(n/2+1):(n)]),
                               var(sam_sigmaS.delta[(n/2+1):(n)]))

  summ_a.t[mc,] = c(mean(sam_a.t[1:(n/2)]),var(sam_a.t[1:(n/2)]))
  summ_a.t[(mc+1),] = c(mean(sam_a.t[(n/2+1):(n)]),
                      var(sam_a.t[(n/2+1):(n)]))

  summ_b.t[mc,] = c(mean(sam_b.t[1:(n/2)]),var(sam_b.t[1:(n/2)]))
  summ_b.t[(mc+1),] = c(mean(sam_b.t[(n/2+1):(n)]),
                      var(sam_b.t[(n/2+1):(n)]))

  summ_alpha.X[,mc,] = c(apply(sam_alpha.X[1:(n/2),],2,mean),
                         apply(sam_alpha.X[1:(n/2),],2,var))
  summ_alpha.X[,(mc+1),] = c(apply(sam_alpha.X[(n/2+1):(n),],2,mean),
                           apply(sam_alpha.X[(n/2+1):(n),],2,var))

  summ_tauS.X[,mc,] = c(apply(sam_tauS.X[1:(n/2),],2,mean),
                        apply(sam_tauS.X[1:(n/2),],2,var))
  summ_tauS.X[,(mc+1),] = c(apply(sam_tauS.X[(n/2+1):(n),],2,mean),
                        apply(sam_tauS.X[(n/2+1):(n),],2,var))
  
  summ_b[,mc,] = c(apply(sam_b[1:(n/2),],2,mean),
                        apply(sam_b[1:(n/2),],2,var))
  summ_b[,(mc+1),] = c(apply(sam_b[(n/2+1):(n),],2,mean),
                            apply(sam_b[(n/2+1):(n),],2,var))
  
  cat("\nChain ",c, "out of ", C, " is complete. Saving to disk...", '\n')
  save(summ_beta, summ_sigmaS.b, summ_sigmaS.delta , summ_alpha, summ_a.t,
       summ_b.t, summ_alpha.X, summ_tauS.X, summ_b,
       file=paste0('../output/tues',c,'.RData'))
  mc = mc+2
  cat("Chain's summary is saved.\n")
}
#-----------------------------------------------------
save.image(file = "../output/workspace_tuesFin.RData")
#-------------------------------------------------------