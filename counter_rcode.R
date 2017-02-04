
rm(list=ls())
options(error = recover);

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
hin <- function(x){
  h <- rep(NA,3)
  h[1] <- x[2]-x[1]
  h[2] <- x[1]
  h[3] <- 1-x[2]
  return(h)
}
hin.jac <- function(x){
  j <- matrix(NA,3,2)
  j[1,] <- c(-1,1)
  j[2,] <- c(1,0)
  j[3,] <- c(0,-1)
  return(j)
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
heq <- function(x){
  h <- (computefunction(x[1],x[2]))
  return(h)
}
heq.jac <- function(x){
  j <- matrix(NA,1,length(x))
  j[1,1] <- computefunction_grad_1(x[1],x[2])
  j[1,2] <- computefunction_grad_2(x[1],x[2])
  return(j)
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
  vector[6] ones;
for (i in 1:6) {
ones[i] = 1;
}
}

  parameters {
    simplex[6] p;
    real<lower=0, upper=1> q;
  }

  transformed parameters {
    simplex[6] op;

    op[1] = (1-q)*p[5];
    op[2] = (1-q)*p[6];
    op[3] = (1-q)*(p[1]+p[2]+p[3]+p[4]);
    op[4] = q*(p[3] + p[5]);
    op[5] = q*(p[4] + p[6]);
    op[6] = q*(p[1] + p[2]);
   }

  model {
    p  ~ dirichlet(ones);
    q  ~ uniform(0,1);
    N  ~ multinomial(op);
  }

  generated quantities {
    real<lower=0, upper=1> tp[7];
    real bounds[2];

    tp[1] = p[6]/(p[5] + p[6]);
    tp[2] = p[5] + p[6];
    tp[3] = p[1] + p[2];
    tp[4] = p[3] + p[4];
    tp[5] = (p[4] + p[6])/(p[3]+p[4]+p[5]+p[6]);
    tp[6] = p[4]/(p[3]+p[4]);
    tp[7] = p[2] + p[4] + p[6];

    bounds[1] = tp[5] * (tp[2] + tp[4]);
    bounds[2] = tp[5] * (tp[2] + tp[4]) + tp[3];
  }

  "


gernating.function <- function(r0,r1){
  result <- rep(0,length(r0))
  idx <- which(r0==0&r1==0)
  result[idx] <- (rbinom(length(idx),1,0.3)+1)
  idx <- which(r0==0&r1==1)
  result[idx] <- (rbinom(length(idx),1,0.5)+1)
  idx <- which(r0==1&r1==1)
  result[idx] <- (rbinom(length(idx),1,0.7)+1)
  return(result)
  
}



simulationtimes <- 500
post.95 <- matrix(0,simulationtimes,2)
post.b <- matrix(0,simulationtimes,2)
true.y.result <- rep(0,simulationtimes)
set.seed(1234)
for(simulation in 1:simulationtimes){
  n <-1500# total number of people
  K <- 2 # number of levels for Likert outcome
  result <- matrix(0,nrow=simulationtimes,ncol=(6+2*3*K*n))
  final_result <- NULL
  #simulation <- simulationtimes
  
  #print(paste0("we are in",simulation,"th run"))
  # n <-25# total number of people
  # K <- 2 # number of levels for Likert outcome
  # 
  # Start simulation
  
  #y <- sample(1:K,n,replace = T)
  r1samp <- rbinom(n,1,0.8)
  r2samp <- rbinom(n,1,0.5)
  r0 <- pmin(r1samp,r2samp) # whether the individual will respond under low incentive
  r1 <- pmax(r1samp,r2samp) # whether the individual will respond under high incentive
  z <- rbinom(n,1,0.5) # the individual receive high or low incentive
  r <- ifelse(z==1,r1,r0) # observed response indicator
  
  y <- gernating.function(r0,r1)
  c <- cbind(sapply(1:K,function(sk) as.integer(r0==1 & y==sk)),sapply(1:K,function(sk) as.integer(r0==0 & r1==1 & y==sk)),sapply(1:K,function(sk) as.integer(r0==0 & r1==0 & y==sk)))
  tmpid <- expand.grid(1:3,1:K)
  colnames(c) <- apply(tmpid[order(tmpid[,1],tmpid[,2]),],1,paste0,collapse=".")
  y <- y-1

  completedata <- cbind(z,r,y)
  
  true.y <- mean(completedata[,3])
  true.y.result[simulation] <- true.y
  
  
  
  idx <- which(r==0)
  y[idx] <- NA
  data <- cbind(z,r,y)
  
  data_count <- table(z,r,y,useNA = "ifany")
  n010 <- data_count[[3]]
  n110 <- data_count[[4]]
  n011 <- data_count[[7]]
  n111 <- data_count[[8]]
  m00 <- data_count[[9]]
  m10 <- data_count[[10]]
  n.obs <- c(n010,n011,m00,n110,n111,m10)
 
  p_int <- c(0.07,0.03,0.25,0.25,0.12,0.28)
  q_int <- 0.5
  # p6_int <- n.obs[2]/((1-q_int)*n)
  # p5_int <- n.obs[1]/((1-q_int)*n)
  # p3_int <- n.obs[4]/(q_int*n)-p5_int
  # p4_int <- n.obs[5]/(q_int*n)-p6_int
  # p2_int <- m10/(q_int*n/2)
  # p1_int <- m10/(q_int*n/2)
  # p_int_la <- c(p2_int,p3_int,p4_int,p5_int,p6_int)
  # p_int_la1 <- p_int_la+runif(5,-0.01,0.01)
  # p_int_la2 <- p_int_la+runif(5,-0.01,0.01)
  # p_int_la3 <- p_int_la+runif(5,-0.01,0.01)
  # p_int_la4 <- p_int_la+runif(5,-0.01,0.01)
  # p_int_la5 <- p_int_la+runif(5,-0.01,0.01)
  # p_int1 <- c(1-sum(p_int_la1),p_int_la1)
  # p_int2 <- c(1-sum(p_int_la2),p_int_la2)
  # p_int3 <- c(1-sum(p_int_la3),p_int_la3)
  # p_int4 <- c(1-sum(p_int_la4),p_int_la4)
  # p_int5 <- c(1-sum(p_int_la5),p_int_la5)
  # 
  # 
  # p_int1 <- c(1-sum(p_int_la),p_int_la)
  int <- list(list(p=p_int,q=q_int+runif(1,0,0.01)),
              list(p=p_int,q=q_int+runif(1,0,0.01)),
              list(p=p_int,q=q_int+runif(1,0,0.01)),
              list(p=p_int,q=q_int+runif(1,0,0.01)),
              list(p=p_int,q=q_int+runif(1,0,0.01))
              )
  
  #n.obs <- c(100, 100, 100, 100, 100, 100);
  #names(n.obs) <- c("n010", "n011", "m00",
  #                 "n110", "n111", "m10");
  
  rst <- stan(model_code = stan.model,
              data = list(N=n.obs),
              iter =20000, warmup = 6000, chains = 5,
              thin = 5, verbose = TRUE,init = int);
  print(rst);
  
  post.p <- extract(rst, "p")$p;
  post.q <- extract(rst, "q")$q;
  post.tp <- extract(rst, "tp")$tp;
  post.bounds <- extract(rst, "bounds")$bounds;
 
 mu <- colMeans(post.bounds)
  Sigma <- var(post.bounds)
  mu1 <- mu[1]
  mu2 <- mu[2]
  sigma1 <- sqrt(2)
  sigma2 <- sqrt(4)
  p <- Sigma[1,2]/(sigma1*sigma2)
  
  
  ##check if P(Y=1) is always in the bounds
  mean(post.tp[,7] > post.bounds[,1] & post.tp[,7] < post.bounds[,2])
  a <- CI95_empircal(post.bounds)
  
  post.95[simulation,1] <- a[1]
  post.95[simulation,2] <- a[2]
  post.b[simulation,] <- colMeans(post.bounds)
  
  
}


idx <- which(post.95[,1]<=0.53&post.95[,2]>=0.63)
jdx <- which(post.95[,1]<=true.y.result&post.95[,2]>=true.y.result)
post.95 <- as.data.frame(post.95)
colnames(post.95) <- c("lower","upper")
library(ggplot2)
plotdata <- data.frame(post.95,post.b)
colnames(plotdata) <- c("lower","upper","l","u")
ggplot(plotdata,aes(lower,upper))+geom_point()+geom_vline(xintercept=0.53)+geom_hline(yintercept = 0.63)+geom_point(aes(l,u),col="blue")

plotdata2 <- data.frame(post.bounds)
colnames(plotdata2) <- c("l","u")
ggplot(plotdata2,aes(l))+geom_density()
