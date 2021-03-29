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
  return(c(overallmean,varE,Lbound, Ubound,Rh,neff,MCSE) )
}

n.samples = n/4
gelmanR(summ_b.t, n.samples)
gelmanR(summ_a.t, n.samples)
gelmanR(summ_alpha, n.samples)

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
  'tauS.X.3', gelmanR(summ_tauS.X[3,,],n.samples)
  ),ncol=8,nrow=14,byrow = TRUE)

output = as.data.frame(output)
colnames(output) = c("Parameter", 'mean', 'Estim var', "low.bound",
                     'upp.bound', 'R-hat','n_eff', 'MC SE')

output %>% select(-Parameter) %>% mutate_all(as.numeric) %>%
  mutate_if(is.numeric, round, digits=6) -> output[,-1]

View(output)

kable(output, "latex")

#load('../output/workspace.RData')

summ_a.t