library(readstata13)
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
  
  idx <- which(abs(result$pdf-0.95)<1e-06)
  differ <- result$y[idx]-result$x[idx]
  jdx <- which.min(differ)
  x.y.95 <- cbind(result$x[idx],result$y[idx])
  if(length(jdx)>1){
    jdx <- jdx[1]
  }
  final <- x.y.95[jdx,]
  return(final)
}


load("sm.bd.rdata")

computebounds <- function(variable){
  z <- data$randomsample
  r <- rep(1,nrow(data))
  idx <- which(is.na(variable))
  r[idx] <- 0
  y <- variable
  data <- cbind(z,r,y)
  data_count <- table(z,r,y,useNA = "ifany")
  n011 <- data_count[[3]]
  n111 <- data_count[[4]]
  n012 <- data_count[[7]]
  n112 <- data_count[[8]]
  n013 <- data_count[[11]]
  n113 <- data_count[[12]]
  n014 <- data_count[[15]]
  n114 <- data_count[[16]]
  n015 <- data_count[[19]]
  n115 <- data_count[[20]]
  m00 <- data_count[[21]]
  m10 <- data_count[[22]]
  
  n.obs <- c(n011,n012,n013,n014,n015,m00,n111,n112,n113,n114,n115,m10)
  
  rst <- sampling(sm.bd,
              data = list(N=n.obs),
              iter =35000, warmup = 6000, chains = 5,
              thin = 5, verbose = TRUE);
  print(rst);
  
  post.p <- extract(rst, "p")$p;
  post.q <- extract(rst, "q")$q;
  post.tp <- extract(rst, "tp")$tp;
  post.bounds <- extract(rst, "bounds")$bounds;
  a <- CI95_empircal(post.bounds)
  return(a)
}
worst.bounds <- function(variable){
  y <- variable
  y[is.na(variable)] <- min(variable,na.rm = T)
  a <- rep(0,2)
  a[1] <- mean(y)
  y[is.na(variable)] <- max(variable,na.rm = T)
  a[2] <- mean(y)
  return(a)
}

setwd("../realdata")
data <- read.dta13("crisis_study_dataset.dta")

y.interest <- data[,24:33]
z <- data$randomsample
createtable <- function(y){
  table(z,y,useNA ="ifany")
}

summary.result <- NULL
for(i in 1:ncol(y.interest)){
  summary.result <- rbind(summary.result,createtable(y.interest[,i]))
}
summary.result <- cbind(rep(c(0,1),ncol(y.interest)),summary.result)


Bd.result <- NULL

for(i in 1:ncol(y.interest)){
  Bd.result <- rbind(Bd.result,computebounds(y.interest[,i]))
}





######
realdata.Result <- cbind(Na.result,Bd.result)
places <- 2
