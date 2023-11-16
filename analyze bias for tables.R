source("functions.R")
load("pop.indices.RData")

load("fit.mod.orig.n150.RData")
load("fit.mod.orig.n200.RData")
load("fit.mod.orig.n300.RData")
load("fit.mod.orig.n500.RData")
load("fit.mod.orig.n800.RData")
load("fit.mod.orig.n1000.RData")






means.orig.mod <- list("n150"=list.mean(fit.mod.orig.n150),
                       "n200"=list.mean(fit.mod.orig.n200),
                       "n300"=list.mean(fit.mod.orig.n300),
                       "n500"=list.mean(fit.mod.orig.n500),
                       "n800"=list.mean(fit.mod.orig.n800),
                       "n1000"=list.mean(fit.mod.orig.n1000))

means.orig.mod

bias.orig.mod <- vector(mode="list")

for(i in 1:length(means.orig.mod)){
  bias.orig.mod[[names(means.orig.mod)[i]]] <- 
    round(means.orig.mod[[i]] - pop.indices, 9)
}

bias.orig.mod 


rmsea200 <- bias.orig.mod$n200[1:7,]

apply(rmsea200,2,  function(x) which.max(abs(x)) )
apply(rmsea200,2,  function(x) which.min(abs(x))-1 )

rmsea500 <- bias.orig.mod$n500[1:7,]
apply(rmsea500,2,  function(x) which.min(abs(x))-1)

rmsea800 <- bias.orig.mod$n800[1:7,]
apply(rmsea800,2,  function(x) which.min(abs(x))-1 )
apply(rmsea800,2,  function(x) which.max(abs(x)) )





srmr200 <- bias.orig.mod$n200[15:23,]

apply(srmr200,2,  function(x) which.min(abs(x))-3 )
apply(srmr200,2,  function(x) which.max(abs(x))-3 )

srmr500 <- bias.orig.mod$n500[15:23,]

apply(srmr500,2,  function(x) which.min(abs(x))-3 )
apply(srmr500,2,  function(x) which.max(abs(x))-3 )

srmr800 <- bias.orig.mod$n800[15:23,]

apply(srmr800,2,  function(x) which.min(abs(x))-3 )
apply(srmr800,2,  function(x) which.max(abs(x))-3 )





cfi200 <- bias.orig.mod$n200[8:14,]
apply(cfi200,2,  function(x) which.max(abs(x)) )

apply(cfi200,2,  function(x) which.min(abs(x))-1 )

cfi500 <- bias.orig.mod$n500[8:14,]
apply(cfi500,2,  function(x) which.max(abs(x)) )


apply(cfi500,2,  function(x) which.min(abs(x))-1 )

cfi800 <- bias.orig.mod$n800[8:14,]

apply(cfi800,2,  function(x) which.min(abs(x))-1 )
apply(cfi800,2,  function(x) which.max(abs(x)))


###high reliability
load("fit.mod.high.n150.RData")
load("fit.mod.high.n200.RData")
load("fit.mod.high.n300.RData")
load("fit.mod.high.n500.RData")
load("fit.mod.high.n800.RData")
load("fit.mod.high.n1000.RData")






means.high.mod <- list("n150"=list.mean(fit.mod.high.n150),
                       "n200"=list.mean(fit.mod.high.n200),
                       "n300"=list.mean(fit.mod.high.n300),
                       "n500"=list.mean(fit.mod.high.n500),
                       "n800"=list.mean(fit.mod.high.n800),
                       "n1000"=list.mean(fit.mod.high.n1000))


means.high.mod

bias.high.mod <- vector(mode="list")

for(i in 1:length(means.high.mod)){
  bias.high.mod[[names(means.high.mod)[i]]] <- 
    round(means.high.mod[[i]] - pop.indices, 9)
}

rmsea200 <- bias.high.mod$n200[1:7,]

apply(rmsea200,2,  function(x) which.max(abs(x)) )
apply(rmsea200,2,  function(x) which.min(abs(x))-1 )

rmsea500 <- bias.high.mod$n500[1:7,]

apply(rmsea500,2,  function(x) which.max(abs(x)) )

apply(rmsea500,2,  function(x) which.min(abs(x))-1 )

rmsea800 <- bias.high.mod$n800[1:7,]
apply(rmsea800,2,  function(x) which.max(abs(x)) )

apply(rmsea800,2,  function(x) which.min(abs(x))-1 )


cfi200 <- bias.high.mod$n200[8:14,]
apply(cfi200,2,  function(x) which.max(abs(x)) )

apply(cfi200,2,  function(x) which.min(abs(x))-1 )

cfi500 <- bias.high.mod$n500[8:14,]

apply(cfi500,2,  function(x) which.max(abs(x)))
apply(cfi500,2,  function(x) which.min(abs(x))-1 )

cfi800 <- bias.high.mod$n800[8:14,]
apply(cfi800,2,  function(x) which.max(abs(x)) )

apply(cfi800,2,  function(x) which.min(abs(x))-1 )






srmr200 <- bias.high.mod$n200[15:23,]

apply(srmr200,2,  function(x) which.min(abs(x))-3 )
apply(srmr200,2,  function(x) which.max(abs(x))-3 )


srmr500 <- bias.high.mod$n500[15:23,]

apply(srmr500,2,  function(x) which.min(abs(x))-3 )
apply(srmr500,2,  function(x) which.max(abs(x))-3 )

srmr800 <- bias.high.mod$n800[15:23,]

apply(srmr800,2,  function(x) which.min(abs(x))-3 )
apply(srmr800,2,  function(x) which.max(abs(x))-3 )



