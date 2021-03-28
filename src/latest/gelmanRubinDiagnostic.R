gelmanR = function(summ, n){
  m = nrow(summ)
  overall = sum(summ[,3]/(n*m))
  B = (n/(m-1)) * sum((summ[,1]-overall)^2)
  We = (1/m) * sum(summ[,2])
  varE = ((n-1)/n) * We + (1/n)* B
  Rh = sqrt(varE/We)
  return(Rh)
}


gelmanR(summ_b.t, n/2)
gelmanR(summ_a.t, n/2)

