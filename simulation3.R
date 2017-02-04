# Parameters to be tuned
# commandarg <- commandArgs(trailingOnly = T)
# i1 <- as.integer(commandarg[[1]])
# print(i1)
# set.seed(i1)

library(Matrix)

set.seed(12345)
# simulationtimes <- 1
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






quadmax <- function(z,r,y){
  idx_1 <- which(z==1)
  idx_0 <- which(z==0)
  idx_1_NA <- which(z==1&r==0)
  idx_0_NA <- which(z==0&r==0)
  n1_NA <- length(idx_1_NA)
  n0_NA <- length(idx_0_NA)
  n1 <- sum(data[,1])
  Q <- matrix(-1/((n-1)*n1*(n-n1)),n,n)
  diag(Q) <- 1/(n1*(n-n1))
  
  y_fit <- y
  
  
  
  
  for(k in 0:(n1_NA+n0_NA)){
    
    temp_y_1 <- rep(2,length(idx_1_NA))
    temp_y_0 <- rep(2,length(idx_0_NA))
    if(k==0){
      y_fit[idx_1_NA] <- temp_y_1
      y_fit[idx_0_NA] <- temp_y_0
      test_result <- test(y_fit,idx_1,idx_0,Q)
      if(test_result==1){
        #print(k)
        return(list(quadmax=mean(y_fit),y_fit=y_fit))
      }
      
    }else{
      pcom <- generate(k,n1_NA,n0_NA)
      if(class(pcom)=="numeric"){
        pcom <- matrix(pcom,1,2)
      }
      for(i in 1:nrow(pcom)){
        y_fit <- y
        n1_change <- pcom[i,1]
        n0_change <- pcom[i,2]
        temp_y_1 <- rep(2,length(idx_1_NA))
        temp_y_0 <- rep(2,length(idx_0_NA))
        if(n1_change!=0){
          temp_y_1[1:n1_change] <- 1 
        }
        if(n0_change!=0){
          temp_y_0[1:n0_change] <- 1
        }
        
        y_fit[idx_1_NA] <- temp_y_1
        y_fit[idx_0_NA] <- temp_y_0
        test_result <- test(y_fit,idx_1,idx_0,Q)
        if(test_result==1){
          #print(c(k,i))
          return(list(quadmax=mean(y_fit),y_fit=y_fit))
        }
      }
    }
    
  }
}

quadmin <- function(z,r,y){
  idx_1 <- which(z==1)
  idx_0 <- which(z==0)
  idx_1_NA <- which(z==1&r==0)
  idx_0_NA <- which(z==0&r==0)
  n1_NA <- length(idx_1_NA)
  n0_NA <- length(idx_0_NA)
  n1 <- sum(data[,1])
  Q <- matrix(-1/((n-1)*n1*(n-n1)),n,n)
  diag(Q) <- 1/(n1*(n-n1))
  
  y_fit <- y
  
  
  
  
  for(k in 0:(n1_NA+n0_NA)){
    
    temp_y_1 <- rep(1,length(idx_1_NA))
    temp_y_0 <- rep(1,length(idx_0_NA))
    if(k==0){
      y_fit[idx_1_NA] <- temp_y_1
      y_fit[idx_0_NA] <- temp_y_0
      test_result <- test(y_fit,idx_1,idx_0,Q)
      if(test_result==1){
        #print(k)
        return(list(quadmin=mean(y_fit),y_fit=y_fit))
      }
      
    }else{
      pcom <- generate(k,n1_NA,n0_NA)
      if(class(pcom)=="numeric"){
        pcom <- matrix(pcom,1,2)
      }
      for(i in 1:nrow(pcom)){
        y_fit <- y
        n1_change <- pcom[i,1]
        n0_change <- pcom[i,2]
        temp_y_1 <- rep(1,length(idx_1_NA))
        temp_y_0 <- rep(1,length(idx_0_NA))
        if(n1_change!=0){
          temp_y_1[1:n1_change] <- 2
        }
        if(n0_change!=0){
          temp_y_0[1:n0_change] <- 2
        }
        
        y_fit[idx_1_NA] <- temp_y_1
        y_fit[idx_0_NA] <- temp_y_0
        test_result <- test(y_fit,idx_1,idx_0,Q)
        if(test_result==1){
          #print(c(k,i))
          return(list(quadmin=mean(y_fit),y_fit=y_fit))
        }
      }
    }
    
  }
}


test <- function(y_fit,idx_1,idx_0,Q){
  differ <- abs(mean(y_fit[idx_1])-mean(y_fit[idx_0]))
  sd_differ <- sqrt(t(y_fit)%*%Q%*%y_fit)
  test <- differ/sd_differ
  result <- ifelse(test<=qnorm(0.975),1,0)
  return(result)
}


generate <- function(k,n1_NA,n0_NA){
  if(k==0){
    return(matrix(0,1,2))
  }else{
    result <- matrix(0,k+1,2)
    result[,1] <- c(0:k)
    result[,2] <- c(k:0)
    idx <- which(result[,1]>n1_NA|result[,2]>n0_NA)
    if(length(idx)==0){
      return(result)
    }else{
      result <- result[-idx,]}
    return(result)
  }
}


generate.c <- function(y.fit,z,r){
  result <- rep(0,3*K*length(z))
  for(i in 1:length(z)){
    result[(6*i-5):(6*i)] <- generate.c.indi(y.fit[i],z[i],r[i])
  }
  return(result)
}

generate.c.indi <- function(y.fit.ind,z,r){
  result <- NULL
  if(z==1&r==0){
    if(y.fit.ind==1){
      result <- c(0,0,0,0,1,0)
    }else{
      result <- c(0,0,0,0,0,1)
    }
  }else if(z==0&r==1){
    if(y.fit.ind==1){
      result <- c(1,0,0,0,0,0)
    }else{
      result <- c(0,1,0,0,0,0)
    }
  }else{
    if(y.fit.ind==1){
      result <- c(0,0,1,0,0,0)
    }else{
      result <- c(0,0,0,1,0,0)
    }
  }
  return(result)
}

explain.c <- function(x){
  y_ex <- rep(0,length(x)/(3*K))
  for(i in 1:(length(x)/(3*K))){
    temp <- x[(6*i-5):(6*i)]
    if(sum(temp[c(1,3,5)])==1){
      y_ex[i] <- 1
    }else{
      y_ex[i] <- 2
    }
  }
  return(y_ex)
}



direct.est <- function(z,r,y){
  n1 <- sum(z==1)
  n <- length(z)
  n0 <- n-n1
  idx.z1.r1 <- which(z==1&r==1)
  e.z1.r1 <- mean(y[idx.z1.r1])
  p.r1_z1 <- length(idx.z1.r1)/n1
  b <- e.z1.r1*p.r1_z1
  p.r0_z1 <- 1-p.r1_z1
  a <- p.r0_z1
  
  idx.z0.r1 <- which(z==0&r==1)
  p.r1_z0 <- length(idx.z0.r1)/n0
  e.r1_z0 <- mean(y[idx.z0.r1])
  d <- p.r1_z0*e.r1_z0
  p.r0_z0 <- 1-p.r1_z0
  c <- p.r0_z0
  
  if(((a+b-d)/c)>=1){
    lx <- 1
  }else{
    lx <- (c+d-b)/a
  }
  
  if(((a+b-d)/c)<=2){
    ux <- 2
  }else{
    ux <- (2*c+d-b)/a
  }
  
  obl <- lx*a+b
  obu <- ux*a+b
  return(c(obl,obu))
  
}






###main funciton
simulationtimes <- 300
result_test <- matrix(0,simulationtimes,4)
for(simulation in 1:simulationtimes){
  n <-500# total number of people
  K <- 2 # number of levels for Likert outcome
  result <- matrix(0,nrow=simulationtimes,ncol=(6+2*3*K*n))
  final_result <- NULL
  #simulation <- simulationtimes
  
  print(paste0("we are in",simulation,"th run"))
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
  
  # Check conditions
  # table(rowSums(c)) # check whether each row of c sum up to 1
  # # condition 1 should be 1, 0
  # unique(rowSums(c[z==1 & r==0,grep("3.",colnames(c))]))
  # unique(rowSums(c[z==1 & r==0,!grepl("3.",colnames(c))]))
  # # condition 2 should be 1, 0
  # table(rowSums(c[z==0 & r==0,grepl("2.",colnames(c)) | grepl("3.",colnames(c))]))
  # table(rowSums(c[z==0 & r==0,grepl("1.",colnames(c))]))
  # # condition 3 should be 1, 0
  # sapply(1:K, function(sk) {
  #   c(unique(rowSums(c[z==1 & r==1 & y==sk,c(paste0(c(1,2),".",sk))])),
  #     unique(rowSums(c[z==1 & r==1 & y==sk,!colnames(c) %in% c(paste0(c(1,2),".",sk))])))
  # })
  # # condition 4 should be 1, 0
  # sapply(1:K, function(sk) {
  #   c(unique(c[z==0 & r==1 & y==sk,paste0("1.",sk)]),
  #     unique(rowSums(c[z==0 & r==1 & y==sk,!colnames(c) %in% paste0("1.",sk)])))
  # })
  completedata <- cbind(z,r,y)
  
  true.y <- mean(completedata[,3])
  
  
  
  idx <- which(r==0)
  y[idx] <- NA
  data <- cbind(z,r,y)
  idx_0<- which(z==0)
  idx_1 <- which(z==1)
  idx_0_NA <- which(z==0&r==0)
  idx_1_NA <- which(z==1&r==0)
  n1 <- sum(data[,1])
  Q <- matrix(-1/((n-1)*n1*(n-n1)),n,n)
  diag(Q) <- 1/(n1*(n-n1))
  
  temp_y_0 <- rep(1,length(idx_0_NA))
  temp_y_1 <- rep(1,length(idx_1_NA))
  
  result_mean <- rep(1,length(idx_1_NA))
  
  ### make both groups as all 1
  temp_y_0 <- rep(2,length(idx_0_NA))
  
  temp_y_1 <- rep(2,length(idx_1_NA))
  
  #y_fit <- y
  #y_fit[idx_0_NA] <- temp_y_0
  # y_fit[idx_1_NA] <- temp_y_1
  # differ <- abs(mean(y_fit[idx_0])-mean(y_fit[idx_1]))
  
  #sd_dif <- sqrt(t(y_fit)%*%Q%*%y_fit)
  result_test[simulation,1] <- quadmin(z,r,y)$quadmin
  result_test[simulation,2] <- quadmax(z,r,y)$quadmax
  result_test[simulation,c(3:4)] <- direct.est(z,r,y)
  
}

lower <- result_test[,1]
upper <- result_test[,2]
ln <- result_test[,3]
un <- result_test[,4]
library(ggplot2)
result_test <- data.frame(result_test)
colnames(result_test) <- c("lower","upper","ln","un")
f <- ggplot(result_test,aes(ln,un))+geom_point()+geom_vline(xintercept = 1.53,col="red")+geom_hline(yintercept = 1.63,col="red")+geom_point(aes(lower,upper),col="blue")

