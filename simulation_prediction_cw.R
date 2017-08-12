
rm(list=ls());
library(rstan);
rstan_options(auto_write=TRUE);
options(mc.cores=parallel::detectCores(), error=recover);


##--------------------------------------------------------------
##               TOOLS FUNCTION
##--------------------------------------------------------------
expit <- function(x) {
    rst <- exp(x)/(1+exp(x));
    rst[which(is.nan(rst))] <- 1;
    rst
}

get.fw <- function(w, a=-2, b=4, c=5, d=-6) {
    expit(a+b*w+c*w^2+d*w^3);
}

##simulation the data with link f = a+b*w+c*w^2+d*w^3
simu.data <- function(ntrain, ntest, ...) {
    n <- ntrain + ntest;
    w <- runif(n);
    y <- rbinom(ntrain, 1, get.fw(w[1:ntrain], ...));
    list(Ntrain=ntrain, Ntest=ntest, N=n, y=y, w=w)
}


##get nearst neighbor for a single point
get.nearst <- function(x12, exi.x12, n.nearst=10) {
    dists <- apply(exi.x12, 1, function(x) {
        rst <- sum((x-x12)^2);
        if (0 == rst)
            rst <- Inf;
        rst
    });

    nids <- order(dists)[1:min(n.nearst, length(dists))];
    nids
}


##spatial matrix
get.kappa.f <- function(x, xprime, lambda, covfn) {
    rst.all <- apply(xprime, 1, function(y) {
        d  <- sum((x-y)^2);
        d  <- lambda * sqrt(d);
        if (0 == covfn) {
            rst <- exp(-d^2/2);
        } else if (1 == covfn) {
            tmp <- sqrt(3)*d;
            rst <- (1 + tmp) * exp(-tmp);
        } else if (2 == covfn) {
            tmp <- sqrt(5)*d;
            rst <- (1 + tmp + tmp^2/3) * exp(-tmp);
        }
        rst
    });
    rbind(rst.all);
}

get.kxx <- function(x12, lambda, tau2, sigma2, covfn) {
    ntot <- nrow(x12);
    rst  <- diag(tau2 + sigma2 + 0.00001, ntot);
    for (i in 1:(ntot-1)) {
        rst[i, (i+1):ntot] <- tau2 * get.kappa.f(x12[i,],
                                                 x12[(i+1):ntot,,drop=FALSE],
                                                 lambda,
                                                 covfn);
        rst[(i+1):ntot, i] <- rst[i, (i+1):ntot];
    }

    rst;
}

## Exponential/Matern Model
get.post <- function(rst.bayes, newdata.loc, existing.loc, covfn=0,
                     n.nearst=50, npost=2000, para=TRUE) {
    ##check overlap
    neighbor <- apply(newdata.loc,
                      1,
                      function(x) {
        get.nearst(x, existing.loc, n.nearst);
    })

    ##
    pars  <- rstan::extract(rst.bayes,
                            pars = c("tau2", "sigma2", "lambda", "m"));

    n.iter <- length(pars$tau2);
    subsmp <- sample(1:n.iter, min(npost, n.iter));
    lmu    <- parallel::mclapply(1:length(subsmp),
                                 function(ti) {
                            print(ti);
                            i      <- subsmp[ti];
                            new.uk <- sapply(1:nrow(newdata.loc),
                                               function(j) {
                                cur.uk  <- pars$m[i, neighbor[,j]];
                                cur.exi <- existing.loc[neighbor[,j],,drop=FALSE];
                                cur.kxx <- get.kxx(cur.exi,
                                                   pars$lambda[i],
                                                   pars$tau2[i],
                                                   pars$sigma2[i],
                                                   covfn);

                                cur.sxx <- solve(cur.kxx);
                                x       <- newdata.loc[j,,drop=FALSE];
                                v12     <- pars$tau2[i] * get.kappa.f(x, cur.exi,  pars$lambda[i], covfn);
                                B       <- cur.kxx[1,1];
                                C       <- v12;
                                CDinv   <- C %*% cur.sxx;
                                cMu     <- 0 + CDinv %*% cur.uk;
                                cVar    <- B - CDinv %*% t(C);
                                rnorm(1, cMu, sqrt(cVar));
                            })

                            cur.lmu <- pnorm(new.uk);
                        }, mc.cores=ifelse(para,
                                           parallel::detectCores()-1,
                                           1));

    lmu <- simplify2array(lmu);
}


##--------------------------------------------------------------
##               STAN MODEL
##--------------------------------------------------------------
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

transformed data {
    vector[N] mu0;
    for (k in 1:N) {
        mu0[k] = 0;
    }
}

parameters{
real<lower=0>  tau2;
real<lower=0>  sigma2;
real<lower=0>  lambda;
vector[N]      m;
}

transformed parameters{
matrix[N,N]    Sigma;
vector[Ntrain] f;

Sigma = kxx(w, N, COVFN, lambda, tau2, sigma2);
for(k in 1:Ntrain) {
  f[k] = normal_cdf(m[k],0,1);
}
}

model{
  //prior
  lambda ~ cauchy(0, 2.5);
  sigma2 ~ lognormal(0, 2);
  tau2   ~ lognormal(0, 2);

  //relation between f and m, since z|m~normal(m,1)
  //target += multi_normal_lpdf(m|mu, Sigma);
  //target += binomial_lpmf(y|1, f);
  m  ~  multi_normal_cholesky(mu0, cholesky_decompose(Sigma));
  y  ~ binomial(1, f);
}'
smodel <- stan_model(model_code = stan.model);

##--------------------------------------------------------------
##               STAN
##--------------------------------------------------------------
n_train <- 50;
n_test  <- 2;
COVFN   <- 1;

data    <- simu.data(n_train,n_test);
fit1    <- sampling(smodel,
                    data    = c(data, list(COVFN=COVFN)),
                    warmup  = 5000,
                    iter    = 10000,
                    control = list(adapt_delta = 0.95),
                    chains  = 4)

traceplot(fit1, pars=c('sigma2', "tau2", "lambda"));

#m= extract(fit1,"m")[[1]]
##
## traceplot(fit1,'sigma2')
## traceplot(fit1,'tau2')
## traceplot(fit1,'lambda')
## fit_summary <- summary(fit1)

##--------------------------------------------------------------
##               POST
##--------------------------------------------------------------
new.w  <- sort(runif(100));
pred.y <- get.post(fit1, cbind(new.w), cbind(data$w), covfn = COVFN);

pdf(file="tmp.pdf");
plot(new.w, get.fw(new.w), type="l", ylim=c(0,1), col="red");
lines(new.w, apply(pred.y, 1, mean));
lines(new.w, apply(pred.y, 1, quantile, 0.025), col="gray");
lines(new.w, apply(pred.y, 1, quantile, 0.975), col="gray");
legend("bottomright", legend=c("True", "GP Mean", "CI"), col=c("red", "black", "gray"), lty=1);
dev.off();

## compare prediction in and out STAN
pred.inx <- n_train+(1:n_test);
new.w2   <- data$w[pred.inx];
pred.y2  <- get.post(fit1, cbind(new.w2), cbind(data$w), covfn=COVFN, para=TRUE);
fit.m    <- pnorm(extract(fit1, "m")$m[, pred.inx]);

pdf(file="tmp.pdf");
i <- 2;
plot(density(pred.y2[i,]));
lines(density(fit.m[,i]), col="red");
dev.off();




#p = apply(m,2,pnorm)
#festimate = colMeans(p)
#w= data$w
#plot(expit(a+b*w+c*w^2),festimate)
## save.image(file="/data/zhangh20/Dan/simuation_prediction/result.Rdata")

library(GPfit)
n = 5; d = 1;
computer_simulator <- function(x){
  x = 2*x+0.5;
  y = sin(10*pi*x)/(2*x) + (x-1)^4;
  return(y)
}
set.seed(3);
library(lhs);
x = maximinLHS(n,d);
y = computer_simulator(x);
GPmodel = GP_fit(x,y);
print(GPmodel)


n = 7; d = 1;
computer_simulator <- function(x) {
  y <- log(x+0.1)+sin(5*pi*x);
  return(y)
}
set.seed(1);
library(lhs);
x = maximinLHS(n,d);
y = computer_simulator(x);
GPmodel = GP_fit(x,y);
print(GPmodel, digits = 4)
plot.GP(GPmodel)


n = 7; d = 1;
computer_simulator <- function(x) {
  y = log(x+0.1)+sin(5*pi*x);
  return(y)
}
set.seed(1);
library(lhs);
x = maximinLHS(n,d);
y = computer_simulator(x);
GPmodel = GP_fit(x,y);
## Plotting with changes from the default line type and characters
plot.GP(GPmodel, resolution = 100, line_type = c(6,2), pch = 5)

n = 200; d = 1;
computer_simulator <- function(x){
  x = 2*x+0.5;
  y = sin(10*pi*x)/(2*x) + (x-1)^4;
  return(y)
}
set.seed(3);
library(lhs);
x = runif(n)
y = computer_simulator(x);
xvec <- seq(from=0,to=1,length.out=200);
GPmodel = GP_fit(x,y);
GPprediction = predict.GP(GPmodel,xvec);
yhat = GPprediction$Y_hat;
mse = GPprediction$MSE;
completedata = GPprediction$complete_data;
completedata;
plot.GP(GPmodel)
plot(xvec,computer_simulator(xvec))
