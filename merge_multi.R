simulationtimes <- 1000
setwd("/users/hzhang1/R/Dan/simulation2_bayesian_multi")
filesDir <- getwd()
files <- dir(filesDir,pattern=".RData")
post.95 <- NULL
for(i in 1:16){
	file <- paste0(i,".RData")
        load(file)
        post.95 <- rbind(post.95,result$post.95)
}
for(i in 18:742){
	file <- paste0(i,".RData")
        load(file)
        post.95 <- rbind(post.95,result$post.95)
}
for(i in 744:simulationtimes){
	file <- paste0(i,".RData")
        load(file)
        post.95 <- rbind(post.95,result$post.95)
}
idx <- which(post.95[,1]<3.56&post.95[,2]>3.72)
leng <- post.95[,2]-post.95[,1]
setwd("/users/hzhang1/R/Dan/simulation2_bayesian_multi_na")
filesDir <- getwd()
files <- dir(filesDir,pattern=".RData")
post.95.na <- NULL
for(i in 1:16){
	file <- paste0(i,".RData")
        load(file)
        post.95.na <- rbind(post.95.na,result$post.95)
}
for(i in 18:742){
	file <- paste0(i,".RData")
        load(file)
        post.95.na <- rbind(post.95.na,result$post.95)
}
for(i in 744:simulationtimes){
	file <- paste0(i,".RData")
        load(file)
        post.95.na <- rbind(post.95.na,result$post.95)
}
idx.na <- which(post.95.na[,1]<3.56&post.95.na[,2]>3.72)
idx.na.self <- which(post.95.na[,1]<3.048&post.95.na[,2]>4.232)
leng.na <-post.95.na[,2]-post.95.na[,1]

ratio <- leng/leng.na