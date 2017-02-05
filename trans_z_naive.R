trans.z.naive <- function(p,pmat,pvec){
  K <- ncol(pmat)
  pn <- p*(pvec[2]+pvec[3])+(1-p)*pvec[3]
  pnmat <- matrix(0,2,K)
  for(k in 1:K){
    pnmat[1,k] <- pmat[1,k]*pvec[1]/(1-pn)+pmat[2,k]*(1-p)*pvec[2]/(1-pn)
    pnmat[2,k] <- pmat[3,k]*pvec[3]/pn+pmat[2,k]*(p)*pvec[2]/pn
  }
  
  ####mean for bayesian method
  Mean.B <- Mean.by.Cate(pmat)
  L.B <- pvec[3]*Mean.B[3]+pvec[2]*Mean.B[2]+pvec[1]
  U.B <- pvec[3]*Mean.B[3]+pvec[2]*Mean.B[2]+pvec[1]*K
  ####mean for Naive method
  Mean.N <- Mean.by.Cate(pnmat)
  L.N <- pn*Mean.N[2]+(1-pn)
  U.N <- pn*Mean.N[2]+K*(1-pn)
  bd.B <- c(L.B,U.B)
  bd.N <- c(L.N,U.N)
  return(list(pnmat=pnmat,pn=pn,bd.B = bd.B,bd.N = bd.N))
}

Mean.by.Cate <- function(pmat){
  K <- ncol(pmat)
  R <- nrow(pmat)
  MeanResult <- rep(0,R)
  for(i in 1:K){
    MeanResult <- i*pmat[,i]+MeanResult
  }
  return(MeanResult)
}

p <- 0.2
pvec <- c(0.04,0.32,0.64)
pmat <- matrix(c(c(0.6,.1,.1,.1,.1),c(.1,.1,.6,.1,.1),c(.1,.1,.1,.1,.6)),3,5,
               byrow=T)

result <- trans.z.naive(p,pmat,pvec)
pnmat <- result[[1]]
pn <- result[[2]]
