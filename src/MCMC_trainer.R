source('./models.R')
require(dplyr)

#------------
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
C = 4 # the number of Markov chains to simulate
MC = 2*C # The number of resulting MCs
summ_beta = array(NA, c(3,MC,2))
summ_sigmaS.b = array(NA,c(MC,2))
summ_alpha= array(NA,c(MC,2))
summ_sigmaS.delta = array(NA,c(MC,2))
summ_a.t = array(NA,c(MC,2))
summ_b.t = array(NA,c(MC,2))
summ_alpha.X = array(NA,c(X,MC,2))
summ_tauS.X = array(NA,c(X, MC,2))

#RNGkind(sample.kind = "Rounding")

#-----------------------------------
mc = 1 # Index for the markov chains
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
  #--------------------------------------------------
  # Updating Markov chains summaries for Diagonsis.
  #-----------------------------------------------------
  summ_beta[,mc,] = c(apply(sam_beta[1:(n/4),],2,mean),
                      apply(sam_beta[1:(n/4),],2,var))
  summ_beta[,(mc+1),] = c(apply(sam_beta[(n/4+1):(n/2),],2,mean),
                        apply(sam_beta[(n/4+1):(n/2),],2,var))
  
  summ_sigmaS.b[mc,] = c(mean(sam_sigmaS.b[1:(n/4)]),var(sam_sigmaS.b[1:(n/4)]))
  summ_sigmaS.b[(mc+1),] = c(mean(sam_sigmaS.b[(n/4+1):(n/2)]),
                           var(sam_sigmaS.b[(n/4+1):(n/2)]))
  
  summ_alpha[mc,] = c(mean(sam_alpha[1:(n/4)]),var(sam_alpha[1:(n/4)]))
  summ_alpha[(mc+1),] = c(mean(sam_alpha[(n/4+1):(n/2)]),
                        var(sam_alpha[(n/4+1):(n/2)]))
  
  summ_sigmaS.delta[mc,] = c(mean(sam_sigmaS.delta[1:(n/4)]),
                             var(sam_sigmaS.delta[1:(n/4)]))
  summ_sigmaS.delta[(mc+1),] = c(mean(sam_sigmaS.delta[(n/4+1):(n/2)]),
                               var(sam_sigmaS.delta[(n/4+1):(n/2)]))
  
  summ_a.t[mc,] = c(mean(sam_a.t[1:(n/4)]),var(sam_a.t[1:(n/4)]))
  summ_a.t[(mc+1),] = c(mean(sam_a.t[(n/4+1):(n/2)]),
                      var(sam_a.t[(n/4+1):(n/2)]))
  
  summ_b.t[mc,] = c(mean(sam_b.t[1:(n/4)]),var(sam_b.t[1:(n/4)]))
  summ_b.t[(mc+1),] = c(mean(sam_b.t[(n/4+1):(n/2)]),
                      var(sam_b.t[(n/4+1):(n/2)]))
  
  summ_alpha.X[,mc,] = c(apply(sam_alpha.X[1:(n/4),],2,mean),
                         apply(sam_alpha.X[1:(n/4),],2,var))
  summ_alpha.X[,(mc+1),] = c(apply(sam_alpha.X[(n/4+1):(n/2),],2,mean),
                           apply(sam_alpha.X[(n/4+1):(n/2),],2,var))
  
  summ_tauS.X[,mc,] = c(apply(sam_tauS.X[1:(n/4),],2,mean),
                        apply(sam_tauS.X[1:(n/4),],2,var))
  summ_tauS.X[,(mc+1),] = c(apply(sam_tauS.X[(n/4+1):(n/2),],2,mean),
                          apply(sam_tauS.X[(n/4+1):(n/2),],2,var))
  
  cat("Chain ",c, "out of ", C, " is complete.", '\n')
  
  mc = mc+2
  
}
#sam_a.t

save.image(file = "../output/workspace.RData")