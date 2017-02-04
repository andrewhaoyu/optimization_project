setwd("/users/hzhang1/R/Dan/simulation2_bayesian_multi")
filesDir <- getwd()
files <- dir(filesDir,pattern=".RData")
post.95 <- NULL
for(file in files){
  load(file)
  post.95 <- rbind(post.95,result$post.95)
}
idx <- which(post.95[,1]<=3.56&post.95[,2]>=3.72)