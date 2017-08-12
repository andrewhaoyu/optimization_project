rm(list=ls());

library(rstan)
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())
expit <- function(x){
  exp(x)/(1+exp(x))
}
a=0
b=0
c=1

##simulation the data with link f = a+b*w+c*w^2
simu.data <- function(ntrain,ntest,n, a=0, b=0, c=1) {
  w <- runif(n);
  y <- a+b*w+c*w^2+rnorm(n,0,0.1)
  list(Ntrain=ntrain,Ntest=ntest, N=n, y=y, w=w)
}

stan.model <-'
functions {
real kappa_f(real x, real xprime, real lambda, int covfn) {
real d;
real tmp;
real rst;

d = (x-xprime)^2;
d = lambda * sqrt(d);

if (0 == covfn) {
rst = exp(-d^2/2);
} else if (1 == covfn) {
tmp = sqrt(3)*d;
rst = (1 + tmp) * exp(-tmp);
} else if (2 == covfn) {
tmp = sqrt(5)*d;
rst = (1 + tmp + tmp^2/3) * exp(-tmp);
}
return(rst);
}

matrix kxx(vector x, int ntot, int covfn, real lambda, real tau2, real sigma2) {
matrix[ntot, ntot] rst;
for (i in 1:ntot) {
# nugget 0.00001 for stability
rst[i,i] = tau2 + sigma2 + 0.00001;
for (j in (i+1):ntot) {
rst[i,j] = tau2 * kappa_f(x[i], x[j], lambda, covfn);
rst[j,i] = rst[i,j];
}
}

return(rst);
}
}

data {
//define data
int<lower=0> Ntrain;
int<lower=0> Ntest;
int<lower=0> N;
int          y[Ntrain];
vector[N]    w;
int          COVFN;
}

parameters{
//define prior
real           beta;
real<lower=0>  tau2;
real<lower=0>  sigma2;
real<lower=0>  lambda;
vector[N]      m;
}

transformed parameters{
//declare mu,Sigma
vector[N]   mu;
matrix[N,N] Sigma;
vector[Ntrain]   f;

//define mu
for(i in 1:N)
mu[i] = beta*w[i];

Sigma = kxx(w, N, COVFN, lambda, tau2, sigma2);

for(k in 1:Ntrain)
f[k] = normal_cdf(m[k],0,1);
}

model{
//prior
beta   ~ normal(0,1000);
lambda ~ cauchy(0, 2.5);
sigma2 ~ lognormal(0, 2);
tau2   ~ lognormal(0, 2);

//relation between f and m, since z|m~normal(m,1)
target += multi_normal_lpdf(m | mu, Sigma);
target += binomial_lpmf(y | 1, f);
}'

smodel <- stan_model(model_code = stan.model);
n_train = 150
n_test = 10
data   <- simu.data(n_train,n_test,n_train+n_test);
data$y = data$y[1:n_train]

fit1   <- sampling(smodel,
                   data    = c(data, list(COVFN=1)),
                   warmup  = 4000,
                   iter    = 8000,
                   control = list(adapt_delta = 0.95),
                   chains  = 4)


 m= extract(fit1,"m")[[1]]
traceplot(fit1,'beta')
traceplot(fit1,'sigma2')
traceplot(fit1,'tau2')
traceplot(fit1,'lambda')

fit_summary = summary(fit1)

p = apply(m,2,pnorm)
festimate = colMeans(p)
w= data$w
plot(expit(a+b*w+c*w^2)[(n_train+1):(n_train+n_test)],festimate[(n_train+1):(n_train+n_test)],xlab="true_value",ylab="estimate_value")




save.image(file="/data/zhangh20/Dan/simuation_prediction/result.Rdata")
