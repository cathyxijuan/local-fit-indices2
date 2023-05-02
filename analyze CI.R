source("functions.R")
load("pop.indices.RData")
load("pop.fit.mod.orig.RData")
load("rmsea.ci.mod.orig.n1000.RData")


pop.rmsea <- pop.indices["pop.rmsea",]
pop.rmsea


rmsea.ci.check <- list()

for(i in 1:5000){
  ci <- rmsea.ci.mod.orig.n1000[[i]]
  rmsea.ci.default <- ci[1,] <= pop.rmsea &pop.rmsea <=ci[2,]
  rmsea.ci.adj.str.exp <- ci[3,] <= pop.rmsea &pop.rmsea <=ci[4,]
  rmsea.ci.adj.str.exp.tri <- ci[5,] <= pop.rmsea &pop.rmsea <=ci[6,]
  rmsea.ci.check[[i]] <- rbind(rmsea.ci.default,
                                rmsea.ci.adj.str.exp,
                               rmsea.ci.adj.str.exp.tri)
}


list.mean(rmsea.ci.check)



str(rmsea.ci.check)
str(rmsea.ci.mod.orig.n200)
str(simplify2array(rmsea.ci.mod.orig.n200))
list.mean(rmsea.ci.mod.orig.n200)
apply(simplify2array(rmsea.ci.mod.orig.n200), 1:2, mean, na.rm = T)

apply(simplify2array(rmsea.ci.check), 1:2, mean, na.rm = T)

str(simplify2array(rmsea.ci.check))
####Usage: find a matrix of means based on a list of matrices
#Argument: lis: a list of matrices
list.mean <- function(lis ){
  if(sum(sapply(lis, is.null))==0){
    apply(simplify2array(lis), 1:2, mean, na.rm = T)
  } else{
    lis <- lis[!(sapply(lis, is.null))]
    apply(simplify2array(lis), 1:2, mean, na.rm = T)
  }
}