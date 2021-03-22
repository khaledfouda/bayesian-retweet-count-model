setwd("D:/concordia/winter21/graphical models/project")


pts = matrix(NA, nrow=10000,ncol=2)
x = rnorm(2,0,.2)
pts[1,] = x
plot(0,0)
rejected = list()
for(i in 2:10000){
  y = rnorm(2,0,.2)#Q
  prob = min(1, dnorm(y)/dnorm(x))
  if(runif(1) <= prob){
    x = y
    points(x[1],x[2], col= rgb(red = 1, green = 0, blue = 0, alpha = 0.2))
  }else{points(x[1],x[2], col='black')}
  pts[i,] = x
}
#pts
prob