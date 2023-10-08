source("functions.R")
load("pop.indices.RData")

load("fit.mod.orig.n150.RData")
load("fit.mod.orig.n200.RData")
load("fit.mod.orig.n300.RData")
load("fit.mod.orig.n500.RData")
load("fit.mod.orig.n800.RData")
load("fit.mod.orig.n1000.RData")






sd.orig.mod <- list("n150"=list.sd(fit.mod.orig.n150),
                    "n200"=list.sd(fit.mod.orig.n200),
                    "n300"=list.sd(fit.mod.orig.n300),
                    "n500"=list.sd(fit.mod.orig.n500),
                    "n800"=list.sd(fit.mod.orig.n800),
                    "n1000"=list.sd(fit.mod.orig.n1000))


rmsea <-  sd.orig.mod$n500[2:7,]

apply(rmsea,2,  function(x) which.min(abs(x)) )

apply(rmsea,2,  function(x) which.max(abs(x)) )



cfi <-  sd.orig.mod$n500[9:14,]

apply(cfi,2,  function(x) which.min(abs(x)) )

apply(cfi,2,  function(x) which.max(abs(x)) )

srmr <-  sd.orig.mod$n500[18:23,]

apply(srmr ,2,  function(x) which.min(abs(x)))

apply(srmr ,2,  function(x) which.max(abs(x)))
###high reliability
load("fit.mod.high.n150.RData")
load("fit.mod.high.n200.RData")
load("fit.mod.high.n300.RData")
load("fit.mod.high.n500.RData")
load("fit.mod.high.n800.RData")
load("fit.mod.high.n1000.RData")






sd.high.mod <- list("n150"=list.sd(fit.mod.high.n150),
                    "n200"=list.sd(fit.mod.high.n200),
                    "n300"=list.sd(fit.mod.high.n300),
                    "n500"=list.sd(fit.mod.high.n500),
                    "n800"=list.sd(fit.mod.high.n800),
                    "n1000"=list.sd(fit.mod.high.n1000))

rmsea <-  sd.high.mod$n500[2:7,]

apply(rmsea,2,  function(x) which.min(abs(x)) )

apply(rmsea,2,  function(x) which.max(abs(x)) )


cfi <-  sd.high.mod$n500[9:14,]

apply(cfi,2,  function(x) which.min(abs(x)))

apply(cfi,2,  function(x) which.max(abs(x)))

srmr <-  sd.high.mod$n500[18:23,]

apply(srmr ,2,  function(x) which.min(abs(x)))

apply(srmr ,2,  function(x) which.max(abs(x)))

