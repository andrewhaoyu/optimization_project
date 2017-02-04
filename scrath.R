library(MASS)
library(mvtnorm)
library(NlcOptim)
library(Rsolnp)
library(alabama)
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

CI_95 <- function(post.bounds){
  
 
  p0 <- c(0.5,0.6)
  ans <- auglag(par=p0,fn=fn,gr=gr,heq=heq,heq.jac = heq.jac,hin = hin,hin.jac = hin.jac)
  return(ans$par)
}




mu <- c(0.5272844,0.6289225)
mu1 <- mu[1]
mu2 <- mu[2]
Sigma <- matrix(c(0.0003443427,0.0002887762,0.0002887762, 0.0003390882),2,2)
sigma1 <- sqrt(2)
sigma2 <- sqrt(4)

p <- Sigma[1,2]/(sigma1*sigma2)

derivative_function <- function(par){
  a <- par[1]
  b <- par[2]
  mu1 <- mu[1]
  mu2 <- mu[2]
  sigma1 <- sqrt(Sigma[1,1])
  sigma2 <- sqrt(Sigma[2,2])
  p <- Sigma[1,2]/(sigma1*sigma2)
  muy_x <- mu2+p*sigma2*(a-mu1)/sigma1
  sigmay_x <- sigma2^2*(1-p^2)
  mux_y <- mu1+p*sigma1*(b-mu2)/sigma2
  sigmax_y <- sigma1^2*(1-p^2)
  
  f_a <- dnorm(a,mean=mu1,sd=sigma1)-dnorm(a,mean=mu1,sd=sigma1)*pnorm(b,mean=muy_x,sd=sigmay_x)
  f_b <- -dnorm(b,mean=mu2,sd=sigma2)*pnorm(a,mean=mux_y,sd=sigmay_x)
  
  tanget <- (-f_b/f_a-1)
  return(tanget^2)
}

parstart <- c(2,3)
optim(parstart,derivative_function)
nleqslv(parstart,derivative_function)
n <- 4800
data <- mvrnorm(n,mu,Sigma)

computefunction <- function(a,b){
 
  mu1 <- mu[1]
  sigma1 <- sqrt(Sigma[1,1])
  result <- pmvnorm(lower=c(a,-Inf),upper = c(Inf,b),mean=mu,sigma=Sigma)
  return((result[1]-0.95))
}
computefunction_grad_1<- function(a,b){
 
  mu1 <- mu[1]
  mu2 <- mu[2]
  sigma1 <- sqrt(Sigma[1,1])
  sigma2 <- sqrt(Sigma[2,2])
  p <- Sigma[1,2]/(sigma1*sigma2)
  muy_x <- mu2+p*sigma2*(a-mu1)/sigma1
  sigmay_x <- sqrt(sigma2^2*(1-p^2))
  mux_y <- mu1+p*sigma1*(b-mu2)/sigma2
  sigmax_y <- sqrt(sigma1^2*(1-p^2))
  
  f_a <- -dnorm(a,mean=mu1,sd=sigma1)*pnorm(b,mean=muy_x,sd=sigmay_x)
  #f_b <- dnorm(b,mean=mu2,sd=sigma2)-dnorm(b,mean=mu2,sd=sigma2)*pnorm(a,mean=mux_y,sd=sigmay_x)
  return(f_a)
}

computefunction_grad_2<- function(a,b){
  
  mu1 <- mu[1]
  mu2 <- mu[2]
  sigma1 <- sqrt(Sigma[1,1])
  sigma2 <- sqrt(Sigma[2,2])
  p <- Sigma[1,2]/(sigma1*sigma2)
  muy_x <- mu2+p*sigma2*(a-mu1)/sigma1
  sigmay_x <- sqrt(sigma2^2*(1-p^2))
  mux_y <- mu1+p*sigma1*(b-mu2)/sigma2
  sigmax_y <- sqrt(sigma1^2*(1-p^2))
  
  f_b <- dnorm(b,mean=mu2,sd=sigma2)-dnorm(b,mean=mu2,sd=sigma2)*pnorm(a,mean=mux_y,sd=sigmay_x)
  return(f_b)
}
a <- seq(min(data[,1]),max(data[,1]),0.001)
value <- rep(0,length(a))
b_result <- rep(0,length(a))
for(i in 1:length(a)){
  try <- optim(0,computefunction,a=a[i],method = "Brent",lower = -10,upper = 10)
  value[i] <- try$value
  b_result[i] <- try$par-a[i]
}
optim(0,computefunction,a=0.58,method = "Brent",lower = -10,upper = 10)
idx <- which(value<=1e-4)
value <- value[idx]
a <- a[idx]
b_result <- b_result[idx]


# eval_f <- function(x){
#   return(list("objective"=x[2]-x[1],
#               "gradient"=c(-1,1)))
# }
# eval_g_ineq <- function(x){
#   constr <- c(x[1]-x[2])
#   #grad <- c(1,-1)
#   return(list("constrains"=constr))
# }
# eval_g_eq <- function(x){
#   constr <- computefunction(x[1],x[2])
#  # grad <- c(computefunction_grad_1(x[1],x[2]),
#           #  computefunction_grad_2(x[1],x[2]))
#   #return(list("constrains"=constr,
#               #"jacobian"=grad))
#   return(list("constrains"=constr))
# }
# x0 <- c(0.4,0.6)
# local_opts <-  list( "algorithm" = "NLOPT_LD_MMA",
#                       "xtol_rel" = 1.0e-7 )
# opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
#               "xtol_rel" = 1.0e-7,
#               "maxeval" = 1000,
#               "local_opts" = local_opts )
# res <- nloptr( x0=x0,
#                eval_f=eval_f,
#                 eval_g_eq=eval_g_eq,
#                opts=opts)



#  objfun <- function(x){
#    return(x[2]-x[1])
#  }
#  confun <- function(x){
#    f <- NULL
#    f <- rbind(f,computefunction(x[1],x[2]))
#    return(list(ceq=f,c=NULL))
#  }
#  x0 <- c(0.5,0.7)
#  NlcOptim(x0,objfun=objfun,confun=confun)
# 
# fn1 <-  function(x){
#      return(x[2]-x[1])
# }
# eqn1 <- function(x){
#   return(computefunction(x[1],x[2]))
# }
# eqB <- 0
# x0 <- c(0.4,0.7)
# result <- solnp(x0,fun=fn1,eqfun=eqn1,eqB=eqB)

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
 }

 p0 <- c(0.5,0.7)
 ans <- auglag(par=p0,fn=fn,gr=gr,heq=heq)
 
 
 x1 = runif (10)
 x2 = runif (10)
 x = cbind (x1, x2)
 m = mecdf (x)

library(data.table)
library(mltools)
 
 
 x.seq <- seq(min(x),max(x),0.001)
 y.seq <- seq(min(y),max(y),0.001)
 n <- nrow(post.bounds)
 pdf <- rep(0,nrow((x.y)))
 x.y <- expand.grid(x.seq,y.seq)
 for( i in 1:nrow(x.y)){
   if(i%%1000==0){
     print(i)
   }
   idx <- which(post.bounds[,1]>=x.y[i,1]&post.bounds[,2]<=x.y[i,2])
   pdf[i] <- length(idx)/n
 }
 
 
 result <- data.frame(x.y,pdf=pdf)
 
 
 CI_empircal_2 <- function(post.bounds){
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
 x <- post.bounds[,1]
 y <- post.bounds[,2]
dt <- data.table(x=x,y=y)
colnames(dt) <- c("x","y")
result <- empirical_cdf(dt,ubounds=list(x=seq(min(x),max(x),0.001),y=seq(min(y),max(y),0.001)))


