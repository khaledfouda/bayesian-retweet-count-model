require(knitr)
require(dplyr)
#--------------------------
gelmanR = function(summ, n){
  m = nrow(summ)
  overallmean = sum(summ[,1])/(m)
  B = (n/(m-1)) * sum((summ[,1]-overallmean)^2)
  We = (1/m) * sum(summ[,2])
  varE = ((n-1)/n) * We + (1/n)* B
  Rh = sqrt(varE/We)
  meansd = mean(sqrt(summ[,2]))
  Lbound = overallmean - 1.96 * meansd / sqrt(n)
  Ubound = overallmean + 1.96 * meansd / sqrt(n)
  neff = m*n*varE/B
  MCSE = sqrt(varE/neff)
  cv = paste0((round(mean( sqrt(summ[,2])/summ[,1]),3)*100),'%')
  return(c(overallmean,varE,Lbound, Ubound,Rh,neff,MCSE,cv))
}
#-------------------------------------
n.samples = n/2

#-------------------------------------
output = matrix( c(
  'Beta 0', gelmanR(summ_beta[1,,],n.samples),
  'beta f', gelmanR(summ_beta[2,,],n.samples),
  'beta d', gelmanR(summ_beta[3,,],n.samples),
  'sigmaS.b',gelmanR(summ_sigmaS.b,n.samples),
  'alpha', gelmanR(summ_alpha,n.samples),
  'sigmaS.delta', gelmanR(summ_sigmaS.delta,n.samples),
  'a.t', gelmanR(summ_a.t,n.samples),
  'b.t', gelmanR(summ_b.t,n.samples),
  'alpha.X.1', gelmanR(summ_alpha.X[1,,],n.samples),
  'alpha.X.2', gelmanR(summ_alpha.X[2,,],n.samples),
  'alpha.X.3', gelmanR(summ_alpha.X[3,,],n.samples),
  'tauS.X.1', gelmanR(summ_tauS.X[1,,],n.samples),
  'tauS.X.2', gelmanR(summ_tauS.X[2,,],n.samples),
  'tauS.X.3', gelmanR(summ_tauS.X[3,,],n.samples),
  'b.1', gelmanR(summ_b[1,,],n.samples),
  'b.2', gelmanR(summ_b[23,,],n.samples),
  'b.3', gelmanR(summ_b[47,,],n.samples)
  ),ncol=9,nrow=17,byrow = TRUE)
#---------------------------------------------
output = as.data.frame(output)
colnames(output) = c("Parameter", 'mean', 'Estim var', "low.bound",
                     'upp.bound', 'R_hat','n_eff', 'MC SE', 'CV')
#-------------------------------------------
output %>% select(-c(Parameter,CV)) %>% mutate_all(as.numeric) %>%
  mutate_if(is.numeric, round, digits=5) %>% 
   mutate(R_hat=round(R_hat,3), n_eff=round(n_eff))  -> output[,-c(1,9)]
#------------------------------------------------------
View(output)
kable(output, "latex")




#--------------------------------------------------
#load('../output/last2.RData')
#--------------------------------------------------
# require(abind)
# sam_b = NA
# sam_alpha.X = NA
# sam_alpha = NA
# sam_sigmaS.b = NA
# sam_beta = NA
# sam_tauS.X = NA
# sam_sigmaS.delta= NA
# sam_a.t = NA
# sam_b.t = NA
# #------------------------------------------------
# load('../output/tues1.RData')
# #---------------------------------------------
# summ_a.t.a = summ_a.t[1:2,]
# summ_alpha.a = summ_alpha[1:2,]
# summ_b.t.a = summ_b.t[1:2,]
# summ_sigmaS.b.a = summ_sigmaS.b[1:2,]
# summ_sigmaS.delta.a = summ_sigmaS.delta[1:2,]
# summ_alpha.X.a = summ_alpha.X[,1:2,]
# summ_b.a = summ_b[,1:2,]
# summ_beta.a = summ_beta[,1:2,]
# summ_tauS.X.a = summ_tauS.X[,1:2,]
# #-------------------------------------
# load('../output/tues2_1.RData')
# #--------------------------------------
# summ_a.t.a = rbind(summ_a.t.a, summ_a.t[1:4,])
# summ_alpha.a = rbind(summ_alpha.a, summ_alpha[1:4,])
# summ_b.t.a = rbind(summ_b.t.a, summ_b.t[1:4,])
# summ_sigmaS.b.a = rbind(summ_sigmaS.b.a, summ_sigmaS.b[1:4,])
# summ_sigmaS.delta.a = rbind(summ_sigmaS.delta.a, summ_sigmaS.delta[1:4,])
# summ_alpha.X.a = abind(summ_alpha.X.a, summ_alpha.X[,1:4,],along=2)
# summ_b.a = abind(summ_b.a, summ_b[,1:4,],along=2)
# summ_beta.a = abind(summ_beta.a, summ_beta[,1:4,], along=2)
# summ_tauS.X.a = abind(summ_tauS.X.a, summ_tauS.X[,1:4,],along=2)
# #-----------------------------------------------------------
# summ_a.t = summ_a.t.a
# summ_alpha = summ_alpha.a
# summ_b.t = summ_b.t.a
# summ_sigmaS.b = summ_sigmaS.b.a
# summ_sigmaS.delta = summ_sigmaS.delta.a
# summ_alpha.X = summ_alpha.X.a
# summ_b = summ_b.a
# summ_beta = summ_beta.a
# summ_tauS.X = summ_tauS.X.a
# #------------------------------------------------------