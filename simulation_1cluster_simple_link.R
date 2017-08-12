library(rstan)
expit = function(x){
  exp(x)/(1+exp(x))
}
n = 100
w = rnorm(n)
a =0
b = 0
c = 1
#simulation the data with link f = a+b*w+c*w^2
y = rbinom(n,1,expit(a+b*w+c*w^2))



data = list(N=n,
            y=y,
            w=w)

stan.model = '
data{
//define data
int<lower=0> N;
int y[N];
vector[N] w;
}
parameters{
//define prior
real beta;
real tau;
real eta;
}
transformed parameters{
//declare m,mu,Sigma
vector[N] m;
vector[N] mu;
matrix[N,N] Sigma;



//declare f;
vector[N] f;


//define mu
for(i in 1:N)
mu[i] = beta*w[i];

for(i in 1:N)
for(j in 1:N)
Sigma[i,j] = tau^2*exp(-(w[i]-w[j])^2/eta^2);

for(i in 1:N)
Sigma[i,i] = tau^2;

for(k in 1:N)
f[k] =  normal_cdf(m[k],0,1);


}
model{
//prior
beta ~ normal(0,1000);
tau ~ cauchy(0,2.5);
eta ~ cauchy(0,2.5);
//relation between m and w;
m ~ multi_normal(mu,Sigma);


//relation between f and m, since z|m~normal(m,1).



for(k in 1:N)
y[k] ~ binomial(1,f[k]);
}

'

fit1 = stan(model_code=stan.model,data=data,
            warmup = 100,
            iter = 5000,
            chains = 4)





