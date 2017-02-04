# Parameters to be tuned
commandarg <- commandArgs(trailingOnly = T)
i1 <- as.integer(commandarg[[1]])
set.seed(i1)

library(gurobi)
simulationtimes <- 1
result <- matrix(0,nrow=simulationtimes,ncol=2)
simulation <- simulationtimes
  n <- 100 # total number of people
  K <- 2 # number of levels for Likert outcome
  
  # Start simulation
  
  #y <- sample(1:K,n,replace = T)
  r1samp <- rbinom(n,1,0.5)
  r2samp <- rbinom(n,1,0.5)
  r0 <- pmin(r1samp,r2samp) # whether the individual will respond under low incentive
  r1 <- pmax(r1samp,r2samp) # whether the individual will respond under high incentive
  z <- rbinom(n,1,0.5) # the individual receive high or low incentive
  r <- ifelse(z==1,r1,r0) # observed response indicator
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
  y <- gernating.function(r0,r1)
   c <- cbind(sapply(1:(K),function(sk) as.integer(r0==1 & y==sk)),sapply(1:(K),function(sk) as.integer(r0==0 & r1==1 & y==sk)),sapply(1:(K),function(sk) as.integer(r0==0 & r1==0 & y==sk)))
   tmpid <- expand.grid(1:3,1:(K))
   colnames(c) <- apply(tmpid[order(tmpid[,1],tmpid[,2]),],1,paste0,collapse=".")

  # Check conditions
  table(rowSums(c)) # check whether each row of c sum up to 1
  # condition 1 should be 1, 0
  unique(rowSums(c[z==1 & r==0,grep("3.",colnames(c))]))
  unique(rowSums(c[z==1 & r==0,!grepl("3.",colnames(c))]))
  # condition 2 should be 1, 0
  table(rowSums(c[z==0 & r==0,grepl("2.",colnames(c)) | grepl("3.",colnames(c))]))
  table(rowSums(c[z==0 & r==0,grepl("1.",colnames(c))]))
  # condition 3 should be 1, 0
  sapply(1:(K), function(sk) {
    c(unique(rowSums(c[z==1 & r==1 & y==sk,c(paste0(c(1,2),".",sk))])),
      unique(rowSums(c[z==1 & r==1 & y==sk,!colnames(c) %in% c(paste0(c(1,2),".",sk))])))
  })
  # condition 4 should be 1, 0
  sapply( 1:(K), function(sk) {
    c(unique(c[z==0 & r==1 & y==sk,paste0("1.",sk)]),
      unique(rowSums(c[z==0 & r==1 & y==sk,!colnames(c) %in% paste0("1.",sk)])))
  })
  completedata <- cbind(z,r,y)
  true.y <- mean(completedata[,3])
  c.true <- as.vector(t(c))
  
  idx <- which(r==0)
  y[idx] <- NA
  data <- cbind(z,r,y)
  b <- rep(c(1:K),3*n)
  
  total <- length(b)
  cvec <- 1/n*(b)
  ub <- 1
  lb <- 0
  Amat <- matrix(0,nrow=n,ncol=total)
  temp <- rep(1,3*K)
  K3 <- 3*K
  for(i in 1:n){
    Amat[i,(1+(i-1)*length(temp)):(i*length(temp))] <- temp
  }
  bvec <- rep(1,n)
  sense <- rep('E',n)
  idx <- which(data[,1]==1&data[,2]==0)
  Amat_temp <- matrix(0,nrow=length(idx),ncol=total)
  for(i in 1:length(idx)){
    Amat_temp[i,((idx[i]-1)*K3+2*K+1):((idx[i]-1)*K3+K3)] <- rep(1,K)
  }
  sense_temp <- rep('E',length(idx))
  bvec_temp <- rep(1,length(idx))
  Amat <- rbind(Amat,Amat_temp)
  sense <- c(sense,sense_temp)
  bvec <- c(bvec,bvec_temp)
  idx <- which(data[,1]==0&data[,2]==0)
  Amat_temp <- matrix(0,nrow=length(idx),ncol=total)
  for(i in 1:length(idx)){
    Amat_temp[i,((idx[i]-1)*K*3+K+1):((idx[i]-1)*K3+K3)] <- rep(1,2*K)
  }
  sense_temp <- rep('E',length(idx))
  bvec_temp <- rep(1,length(idx))
  Amat <- rbind(Amat,Amat_temp)
  sense <- c(sense,sense_temp)
  bvec <- c(bvec,bvec_temp)
  
  jdx <- which(data[,1]==1&data[,2]==1)
  Amat_temp <- matrix(0,nrow=length(jdx),ncol=total)
  for(j in 1:length(jdx)){
    y.jdx <- data[jdx[j],3]
    Amat_temp[j,(jdx[j]-1)*K3+y.jdx] <- 1
    Amat_temp[j,(jdx[j]-1)*K3+K+y.jdx] <- 1
  }
  sense_temp <- rep('E',length(jdx))
  bvec_temp <- rep(1,length(jdx))
  Amat <- rbind(Amat,Amat_temp)
  sense <- c(sense,sense_temp)
  bvec <- c(bvec,bvec_temp)
  
  jdx <- which(data[,1]==0&data[,2]==1)
  Amat_temp <- matrix(0,nrow=length(jdx),ncol=total)
  for(j in 1:length(jdx)){
    y.jdx <- data[jdx[j],3]
    Amat_temp[j,(jdx[j]-1)*K3+y.jdx] <- 1
  }
  sense_temp <- rep('E',length(jdx))
  bvec_temp <- rep(1,length(jdx))
  Amat <- rbind(Amat,Amat_temp)
  sense <- c(sense,sense_temp)
  bvec <- c(bvec,bvec_temp)
  n1 <- sum(data[,1])
  a <- data[,1]/n1-(1-data[,1])/(n-n1)
  Q <- matrix(-1/((n-1)*n1*(n-n1)),n,n)
  diag(Q) <- 1/(n1*(n-n1))
  B <- kronecker(diag(n),rep(1:K,3))
  Qmat <- B%*%(a%*%t(a)-3.84*Q)%*%t(B)
  c.true%*%Qmat%*%c.true
  eigen_decom <- eigen(Qmat)
  u <- try$u
  eigen_adjust <- eigen_decom$values
  eigen_adjust[abs(eigen_adjust)<1e-15] <- 0
  D <- diag(eigen_adjust)
  delta <- eigen_decom$vectors
  Q_eigen <- delta%*%D%*%solve(delta)
  c.true%*%Qmat%*%c.true
  c.true%*%Q_eigen%*%c.true
  ub <- rep(1,total)
  lb <- rep(0,total)
  vtypes <- rep("B",total)
  sense <- rep("=",length(bvec))
  model <- list()
  model$A <- Amat
  model$modelsense <- "max"
  model$obj <- cvec
  model$rhs <- bvec
  model$sense <- sense
  model$vtypes <- vtypes
  # model$ub <- ub
  # model$lb <- lb
  
  qc1 <- list()
  qc1$Qc <- Q_eigen
  qc1$rhs <- 0.0
  model$quadcon <- list(qc1)
  res <- gurobi(model)
  c_gurobi <- res$x
  result[simulation,1] <- res$objval
  model$modelsense <- "min"
  res <- gurobi(model)
  result[simulation,2] <- res$objval

save(result,file=paste0("/users/hzhang1/R/Dan/simulation1/result",i1,".Rdata"))








