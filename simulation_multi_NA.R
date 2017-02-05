
rm(list=ls())
#options(error = recover);
commandarg <- commandArgs(trailingOnly=T)
i1 <- as.numeric(commandarg[1])
library(rstan);
library(NlcOptim)
library(MASS)
library(mvtnorm)
CI95_empircal <- function(post.bounds){
  x <- post.bounds[,1]
  y <- post.bounds[,2]
  max_x <- quantile(x,probs = 0.05)
  x.seq <- seq(min(x),max_x,0.001)
  y.seq <- seq(min(y),max(y),0.001)
  n <- nrow(post.bounds)
  x.y <- expand.grid(x.seq,y.seq)
  pdf <- rep(0,nrow((x.y)))
  colnames(x.y) <- c("x","y")
  for( i in 1:nrow(x.y)){
    idx <- which(post.bounds[,1]>=x.y[i,1]&post.bounds[,2]<=x.y[i,2])
    pdf[i] <- length(idx)/n
  }
  result <- data.frame(x.y,pdf=pdf)
  
  idx <- which(abs(result$pdf-0.95)<1e-03)
  differ <- result$y[idx]-result$x[idx]
  jdx <- which.min(differ)
  x.y.95 <- cbind(result$x[idx],result$y[idx])
  if(length(jdx)>1){
    jdx <- jdx[1]
  }
  final <- x.y.95[jdx,]
  return(final)
}

computefunction <- function(a,b){
  
  mu1 <- mu[1]
  sigma1 <- sqrt(Sigma[1,1])
  result <- pmvnorm(lower=c(a,-Inf),upper = c(Inf,b),mean=mu,sigma=Sigma)
  return((result[1]-0.95))
}

fn <- function(x){return(x[2]-x[1])}
gr <- function(x){
  g <- rep(NA,2)
  g[1] <- -1
  g[2] <- 1
  return(g)
}

CI_95 <- function(post.bounds){
  
  
  p0 <- c(0.5,0.6)
  ans <- auglag(par=p0,fn=fn,gr=gr,heq=heq,heq.jac = heq.jac,hin = hin,hin.jac = hin.jac)
  return(ans$par)
}






stan.model <- "
data {

int<lower=0> N[6];
}

transformed data {

vector[10] ones;

for(i in 1:10){
ones[i] = 1;
}
}

parameters {

simplex[10] p;

}

transformed parameters {

simplex[6] op;



op[1] = p[1];
op[2] = p[2];
op[3] = p[3];
op[4] = p[4];
op[5] = p[5];
op[6]= p[6]+p[7]+p[8]+p[9]+p[10];

}

model {
p  ~ dirichlet(ones);
N  ~ multinomial(op);
}

generated quantities {

real bounds[2];
bounds[1] = p[1]+2*p[2]+3*p[3]+4*p[4]+5*p[5]+(p[6]+p[7]+
p[8]+p[9]+p[10]);
bounds[2] = p[1]+2*p[2]+3*p[3]+4*p[4]+5*p[5]+5*(p[6]+p[7]+
p[8]+p[9]+p[10]);
}

"


gernating.function <- function(r0,r1){
  result <- rep(0,length(r0))
  idx <- which(r0==0&r1==0)
  result[idx] <- t(rmultinom(length(idx),1,prob=c(0.6,0.1,0.1,0.1,0.1)))%*%c(1,2,3,4,5)
  idx <- which(r0==0&r1==1)
  result[idx] <- t(rmultinom(length(idx),1,prob=c(0.1,0.1,0.6,0.1,0.1)))%*%c(1,2,3,4,5)
  idx <- which(r0==1&r1==1)
  result[idx] <- t(rmultinom(length(idx),1,prob=c(0.1,0.1,0.1,0.1,0.6)))%*%c(1,2,3,4,5)
  return(result)
  
}



simulationtimes <- 1
post.95 <- matrix(0,simulationtimes,2)
post.b <- matrix(0,simulationtimes,2)
true.y.result <- rep(0,simulationtimes)
set.seed(i1)
for(simulation in 1:simulationtimes){
  n <-3000# total number of people
  K <- 5 # number of levels for Likert outcome
  #result <- matrix(0,nrow=simulationtimes,ncol=(6+2*3*K*n))
  final_result <- NULL
  #simulation <- simulationtimes
  
  #print(paste0("we are in",simulation,"th run"))
  # n <-25# total number of people
  # K <- 2 # number of levels for Likert outcome
  # 
  # Start simulation
  
  #y <- sample(1:K,n,replace = T)
  r1samp <- rbinom(n,1,0.8)
  r2samp <- rbinom(n,1,0.8)
  r0 <- pmin(r1samp,r2samp) # whether the individual will respond under low incentive
  r1 <- pmax(r1samp,r2samp) # whether the individual will respond under high incentive
  z <- rbinom(n,1,0.2) # the individual receive high or low incentive
  r <- ifelse(z==1,r1,r0) # observed response indicator
  
  y <- gernating.function(r0,r1)
  completedata <- cbind(z,r,y)
  
  true.y <- mean(completedata[,3])
  true.y.result[simulation] <- true.y
  
  
  
  idx <- which(r==0)
  y[idx] <- NA
  data <- cbind(z,r,y)
  
  data_count <- table(r,y,useNA = "ifany")
  n11 <- data_count[2,1]
  n12 <- data_count[2,2]
  n13 <- data_count[2,3]
  n14 <- data_count[2,4]
  n15 <- data_count[2,5]
  m <- data_count[1,6]
  
  n.obs <- c(n11,n12,n13,n14,n15,m)
  # p_int<- c(0.024,rep(0.004,4),0.032,0.032,0.192,0.032,0.032,rep(0.064,4),0.384)
  # q_int <- 0.2
  # int <- list(list(p=p_int,q=q_int+runif(1,0,0.01)),
  #             list(p=p_int,q=q_int+runif(1,0,0.01)),
  #             list(p=p_int,q=q_int+runif(1,0,0.01)),
  #             list(p=p_int,q=q_int+runif(1,0,0.01)),
  #             list(p=p_int,q=q_int+runif(1,0,0.01))
  # )
  #n.obs <- c(100, 100, 100, 100, 100, 100);
  #names(n.obs) <- c("n010", "n011", "m00",
  #                 "n110", "n111", "m10");
  
  rst <- stan(model_code = stan.model,
              data = list(N=n.obs),
              iter =35000, warmup = 6000, chains = 5,
              thin = 5, verbose = TRUE);
  print(rst);
  
  post.p <- extract(rst, "p")$p;
  post.bounds <- extract(rst, "bounds")$bounds;
  

  a <- CI95_empircal(post.bounds)
  
  post.95[simulation,1] <- a[1]
  post.95[simulation,2] <- a[2]
  post.b[simulation,] <- colMeans(post.bounds)
  
  
}

result <- list(post.tp = post.tp, post.bounds=post.bounds,post.95 = post.95)

save(result,file=paste0("/users/hzhang1/R/Dan/simulation2_bayesian_multi_na/",i1,".RData"))
