require(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


data = readRDS('../data/model_input.rds')
model_input = list(
X = data$X, # the total number of root users
J = data$J, # where J_x is the number of retweets user x got including
# them.
JCUM = data$JCUM,
N = data$N, # the total number of observed reaction time for all training and
# predictions tweets.
f = data$root_f, # the number of followers of user
d = data$d,
S = data$S,
M = data$M)

fit = stan('./stan_model.stan', data=model_input,
               iter = 5000, warmup=1000, chains=3)

la <- extract(fit, permuted = TRUE)
lab <- la
print(fit)
plot(fit)


rnorm(1,0,100) -> a1;a1
rinvgamma(1,.5,.5) -> a2; a2

which(data$M==67)
data$f[11]
data$M[11]

ind =which(data$root_f<data$M)
data$root_f[ind]
data$M[ind]