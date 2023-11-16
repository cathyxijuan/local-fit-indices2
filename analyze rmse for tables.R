source("functions.R")
load("pop.indices.RData")

load("fit.mod.orig.n150.RData")
load("fit.mod.orig.n200.RData")
load("fit.mod.orig.n300.RData")
load("fit.mod.orig.n500.RData")
load("fit.mod.orig.n800.RData")
load("fit.mod.orig.n1000.RData")




rmse.orig.mod <- list("n150"=rmse(fit.mod.orig.n150,pop.indices),
                      "n200"=rmse(fit.mod.orig.n200,pop.indices),
                      "n300"=rmse(fit.mod.orig.n300,pop.indices),
                      "n500"=rmse(fit.mod.orig.n500,pop.indices),
                      "n800"=rmse(fit.mod.orig.n800,pop.indices),
                      "n1000"=rmse(fit.mod.orig.n1000,pop.indices))

rmsea <- rmse.orig.mod$n200[1:7, ]
apply(rmsea,2,  function(x) which.min(abs(x))-1 )
apply(rmsea,2,  function(x) which.max(abs(x))-1 )

cfi <- rmse.orig.mod$n200[8:14, ]
apply(cfi,2,  function(x) which.min(abs(x))-1 )
apply(cfi,2,  function(x) which.max(abs(x))-1 )


srmr <- rmse.orig.mod$n200[15:nrow(rmse.orig.mod$n200), ]
apply(srmr,2,  function(x) which.min(abs(x))-3 )
apply(srmr,2,  function(x) which.max(abs(x)) -3)



load("fit.mod.high.n150.RData")
load("fit.mod.high.n200.RData")
load("fit.mod.high.n300.RData")
load("fit.mod.high.n500.RData")
load("fit.mod.high.n800.RData")
load("fit.mod.high.n1000.RData")






rmse.high.mod <- list("n150"=rmse(fit.mod.high.n150,pop.indices),
                      "n200"=rmse(fit.mod.high.n200,pop.indices),
                      "n300"=rmse(fit.mod.high.n300,pop.indices),
                      "n500"=rmse(fit.mod.high.n500,pop.indices),
                      "n800"=rmse(fit.mod.high.n800,pop.indices),
                      "n1000"=rmse(fit.mod.high.n1000,pop.indices))
rmse.high.mod 
rmsea <- rmse.high.mod$n200[1:7, ]
apply(rmsea,2,  function(x) which.min(abs(x))-1 )
apply(rmsea,2,  function(x) which.max(abs(x))-1 )

cfi <- rmse.high.mod$n200[8:14, ]
apply(cfi,2,  function(x) which.min(abs(x))-1 )
apply(cfi,2,  function(x) which.max(abs(x))-1 )


srmr <- rmse.high.mod$n200[15:nrow(rmse.high.mod$n200), ]
apply(srmr,2,  function(x) which.min(abs(x))-3 )
apply(srmr,2,  function(x) which.max(abs(x))-3 )


