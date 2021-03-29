gelmanR = function(summ, n){
  m = nrow(summ)
  #eachsum = summ[,1] * n
  overallmean = sum(summ[,3]/(n*m))
  B = (n/(m-1)) * sum((summ[,1]-overallmean)^2)
  We = (1/m) * sum(summ[,2])
  varE = ((n-1)/n) * We + (1/n)* B
  Rh = sqrt(varE/We)
  meansd = mean(sqrt(summ[,2]))
  Lbound = overallmean - 1.96 * meansd / sqrt(n)
  Ubound = overallmean + 1.96 * meansd / sqrt(n)
  return(c(overallmean,varE,Lbound, Ubound,Rh) )
}


gelmanR(summ_b.t, n/2)
gelmanR(summ_a.t, n/2)
gelmanR(summ_alpha, n/2)

output = matrix( c(
  'Beta 0', gelmanR(summ_beta[1,,],n/2),
  'beta f', gelmanR(summ_beta[2,,],n/2),
  'beta d', gelmanR(summ_beta[3,,],n/2),
  'sigmaS.b',gelmanR(summ_sigmaS.b,n/2),
  'alpha', gelmanR(summ_alpha,n/2),
  'sigmaS.delta', gelmanR(summ_sigmaS.delta,n/2),
  'a.t', gelmanR(summ_a.t,n/2),
  'b.t', gelmanR(summ_b.t,n/2),
  'alpha.X.1', gelmanR(summ_alpha.X[1,,],n/2),
  'alpha.X.2', gelmanR(summ_alpha.X[2,,],n/2),
  'alpha.X.3', gelmanR(summ_alpha.X[3,,],n/2),
  'tauS.X.1', gelmanR(summ_tauS.X[1,,],n/2),
  'tauS.X.2', gelmanR(summ_tauS.X[2,,],n/2),
  'tauS.X.3', gelmanR(summ_tauS.X[3,,],n/2)
  ),ncol=6,nrow=14,byrow = TRUE)

output = as.data.frame(output)
colnames(output) = c("Parameter", 'mean', 'Estim var', "low.bound",
                     'upp.bound', 'R-hat')

output %>% mutate_if(is.numeric, round, digits=2)
output %>% select(-Parameter) %>% mutate_all(as.numeric) %>%
  mutate_if(is.numeric, round, digits=6) -> output[,-1]

View(output)

kable(output, "latex")