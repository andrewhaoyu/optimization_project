# Parameters to be tuned
commandarg <- commandArgs(trailingOnly = T)
i1 <- as.integer(commandarg[[1]])
print(i1)
set.seed(i1)

library(Matrix)

library(gurobi)
simulationtimes <- 1
n <-50# total number of people
K <- 2 # number of levels for Likert outcome
result <- matrix(0,nrow=simulationtimes,ncol=(6+2*3*K*n))
final_result <- NULL
simulation <- simulationtimes
  
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
  b <- rep(c(1:K),3*n)
  
  total <- length(b)
  cvec <- c(1/n*(b))
  
  ub <- 1
  lb <- 0
  Amat <- matrix(0,nrow=n,ncol=total)
  temp <- rep(1,3*K)
  K3 <- 3*K
  for(i in 1:n){
    Amat[i,(1+(i-1)*length(temp)):(i*length(temp))] <- temp
  }
  bvec <- rep(1,n)
  sense <- rep('=',n)
  idx <- which(data[,1]==1&data[,2]==0)
  if(length(idx)!=0){
    Amat_temp <- matrix(0,nrow=length(idx),ncol=total)
    for(i in 1:length(idx)){
      Amat_temp[i,((idx[i]-1)*K3+2*K+1):((idx[i]-1)*K3+K3)] <- rep(1,K)
    }
    sense_temp <- rep('=',length(idx))
    bvec_temp <- rep(1,length(idx))
    Amat <- rbind(Amat,Amat_temp)
    sense <- c(sense,sense_temp)
    bvec <- c(bvec,bvec_temp)
  }
  
  idx <- which(data[,1]==0&data[,2]==0)
  if(length(idx)!=0){
    Amat_temp <- matrix(0,nrow=length(idx),ncol=total)
    for(i in 1:length(idx)){
      Amat_temp[i,((idx[i]-1)*K*3+K+1):((idx[i]-1)*K3+K3)] <- rep(1,2*K)
    }
    sense_temp <- rep('=',length(idx))
    bvec_temp <- rep(1,length(idx))
    Amat <- rbind(Amat,Amat_temp)
    sense <- c(sense,sense_temp)
    bvec <- c(bvec,bvec_temp)
  }
  
  
  jdx <- which(data[,1]==1&data[,2]==1)
  if(length(jdx)!=0){
    Amat_temp <- matrix(0,nrow=length(jdx),ncol=total)
    for(j in 1:length(jdx)){
      y.jdx <- data[jdx[j],3]
      Amat_temp[j,(jdx[j]-1)*K3+y.jdx] <- 1
      Amat_temp[j,(jdx[j]-1)*K3+K+y.jdx] <- 1
    }
    sense_temp <- rep('=',length(jdx))
    bvec_temp <- rep(1,length(jdx))
    Amat <- rbind(Amat,Amat_temp)
    sense <- c(sense,sense_temp)
    bvec <- c(bvec,bvec_temp)
  }
  
  
  jdx <- which(data[,1]==0&data[,2]==1)
  if(length(jdx)!=0){
    Amat_temp <- matrix(0,nrow=length(jdx),ncol=total)
    for(j in 1:length(jdx)){
      y.jdx <- data[jdx[j],3]
      Amat_temp[j,(jdx[j]-1)*K3+y.jdx] <- 1
    }
    sense_temp <- rep('=',length(jdx))
    bvec_temp <- rep(1,length(jdx))
    Amat <- rbind(Amat,Amat_temp)
    sense <- c(sense,sense_temp)
    bvec <- c(bvec,bvec_temp)
  }
  
  
  
  
  n1 <- sum(data[,1])
  a <- data[,1]/n1-(1-data[,1])/(n-n1)
  Q <- matrix(-1/((n-1)*n1*(n-n1)),n,n)
  diag(Q) <- 1/(n1*(n-n1))
  B <- kronecker(diag(n),rep(1:K,3))
  Qmat <- B%*%(a%*%t(a)-qnor*Q)%*%t(B)
  eigen_decom <- eigen(Qmat)
  eigen_adjust <- eigen_decom$values
  eigen_adjust[abs(eigen_adjust)<1e-15] <- 0
  D <- diag(eigen_adjust)
  delta <- eigen_decom$vectors
  Q_eigen <- delta%*%D%*%solve(delta)
  
 
  
  
  
  #ub <- rep(1,total)
  #lb <- rep(0,total)
  vtypes <- rep("B",ncol(Amat))
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
  print("first step is okay")
  params <- list(MIPGap=0.0059,TimeLimit = 3600)
  #res <- gurobi(model,params)
  res1 <- gurobi(model,params)
  print("gurobi can be run out for max!!!!!!")
  if (is.null(res1$objval)){
    result[simulation,1] <- NA
    result[simulation,7:(6+3*K*n)] <- NA
  }else{
    result[simulation,1] <- res1$objval
    result[simulation,7:(6+3*K*n)] <- res1$x
  }
  result[simulation,4] <- res1$runtime
  
  model$modelsense <- "min"
  res2 <- gurobi(model,params)
  print("gurobi can run out for min!!!!!")
  if(is.null(res2$objval)){
    result[simulation,2] <- NA
    result[simulation,(7+3*K*n):(6+3*K*2*n)] <- NA
  }else{
    result[simulation,2] <- res2$objval
    result[simulation,(7+3*K*n):(6+3*K*2*n)] <- res2$x
  }
  # result[simulation,2] <- res$objval
  result[simulation,3] <- true.y
  result[simulation,5] <- res2$runtime
  result[simulation,6] <- i1
  
  
  
  final_result_temp <- list(result=result,res1=res1,res2=res2,completedata=completedata,
                            data=data)
  
  final_result[[simulation]] <- final_result_temp
  
  

save(result,file=paste0("/users/hzhang1/R/Dan/simulation1/result_final_adj",i1,".Rdata"))











