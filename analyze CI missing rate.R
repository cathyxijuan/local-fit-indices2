source("functions.R")
load("pop.indices.RData")
load("pop.fit.mod.orig.RData")
load("rmsea.ci.mod.orig.n150.RData")
load("rmsea.ci.mod.orig.n200.RData")
load("rmsea.ci.mod.orig.n300.RData")
load("rmsea.ci.mod.orig.n500.RData")
load("rmsea.ci.mod.orig.n800.RData")
load("rmsea.ci.mod.orig.n1000.RData")
load("rmsea.ci.mod.high.n150.RData")
load("rmsea.ci.mod.high.n200.RData")
load("rmsea.ci.mod.high.n300.RData")
load("rmsea.ci.mod.high.n500.RData")
load("rmsea.ci.mod.high.n800.RData")
load("rmsea.ci.mod.high.n1000.RData")
load("srmr.ci.mod.orig.n150.RData")
load("srmr.ci.mod.orig.n200.RData")
load("srmr.ci.mod.orig.n300.RData")
load("srmr.ci.mod.orig.n500.RData")
load("srmr.ci.mod.orig.n800.RData")
load("srmr.ci.mod.orig.n1000.RData")
load("srmr.ci.mod.high.n150.RData")
load("srmr.ci.mod.high.n200.RData")
load("srmr.ci.mod.high.n300.RData")
load("srmr.ci.mod.high.n500.RData")
load("srmr.ci.mod.high.n800.RData")
load("srmr.ci.mod.high.n1000.RData")
load("cfi.ci.mod.orig.n150.RData")
load("cfi.ci.mod.orig.n200.RData")
load("cfi.ci.mod.orig.n300.RData")
load("cfi.ci.mod.orig.n500.RData")
load("cfi.ci.mod.orig.n800.RData")
load("cfi.ci.mod.orig.n1000.RData")
load("cfi.ci.mod.high.n150.RData")
load("cfi.ci.mod.high.n200.RData")
load("cfi.ci.mod.high.n300.RData")
load("cfi.ci.mod.high.n500.RData")
load("cfi.ci.mod.high.n800.RData")
load("cfi.ci.mod.high.n1000.RData")

pop.rmsea <- pop.indices["pop.rmsea",]
rmsea.ci.orig.mod.missing.below <- 
  list("n150"=missing.rate.below.rmsea(pop.rmsea, rmsea.ci.mod.orig.n150),
       "n200"=missing.rate.below.rmsea(pop.rmsea, rmsea.ci.mod.orig.n200),
       "n300"=missing.rate.below.rmsea(pop.rmsea, rmsea.ci.mod.orig.n300),
       "n500"=missing.rate.below.rmsea(pop.rmsea, rmsea.ci.mod.orig.n500),
       "n800"=missing.rate.below.rmsea(pop.rmsea, rmsea.ci.mod.orig.n800),
       "n1000"=missing.rate.below.rmsea(pop.rmsea, rmsea.ci.mod.orig.n1000))

rmsea.ci.orig.mod.missing.below
save(rmsea.ci.orig.mod.missing.below, file="rmsea.ci.orig.mod.missing.below.RData")


rmsea.ci.orig.mod.missing.above <- 
  list("n150"=missing.rate.above.rmsea(pop.rmsea, rmsea.ci.mod.orig.n150),
       "n200"=missing.rate.above.rmsea(pop.rmsea, rmsea.ci.mod.orig.n200),
       "n300"=missing.rate.above.rmsea(pop.rmsea, rmsea.ci.mod.orig.n300),
       "n500"=missing.rate.above.rmsea(pop.rmsea, rmsea.ci.mod.orig.n500),
       "n800"=missing.rate.above.rmsea(pop.rmsea, rmsea.ci.mod.orig.n800),
       "n1000"=missing.rate.above.rmsea(pop.rmsea, rmsea.ci.mod.orig.n1000))

rmsea.ci.orig.mod.missing.above
save(rmsea.ci.orig.mod.missing.above, file="rmsea.ci.orig.mod.missing.above.RData")

rmsea.ci.high.mod.missing.below <- 
  list("n150"=missing.rate.below.rmsea(pop.rmsea, rmsea.ci.mod.high.n150),
       "n200"=missing.rate.below.rmsea(pop.rmsea, rmsea.ci.mod.high.n200),
       "n300"=missing.rate.below.rmsea(pop.rmsea, rmsea.ci.mod.high.n300),
       "n500"=missing.rate.below.rmsea(pop.rmsea, rmsea.ci.mod.high.n500),
       "n800"=missing.rate.below.rmsea(pop.rmsea, rmsea.ci.mod.high.n800),
       "n1000"=missing.rate.below.rmsea(pop.rmsea, rmsea.ci.mod.high.n1000))

rmsea.ci.high.mod.missing.below
save(rmsea.ci.high.mod.missing.below, file="rmsea.ci.high.mod.missing.below.RData")



rmsea.ci.high.mod.missing.above <- 
  list("n150"=missing.rate.above.rmsea(pop.rmsea, rmsea.ci.mod.high.n150),
       "n200"=missing.rate.above.rmsea(pop.rmsea, rmsea.ci.mod.high.n200),
       "n300"=missing.rate.above.rmsea(pop.rmsea, rmsea.ci.mod.high.n300),
       "n500"=missing.rate.above.rmsea(pop.rmsea, rmsea.ci.mod.high.n500),
       "n800"=missing.rate.above.rmsea(pop.rmsea, rmsea.ci.mod.high.n800),
       "n1000"=missing.rate.above.rmsea(pop.rmsea, rmsea.ci.mod.high.n1000))

rmsea.ci.high.mod.missing.above
save(rmsea.ci.high.mod.missing.above, file="rmsea.ci.high.mod.missing.above.RData")






pop.srmr <- pop.indices["pop.srmr",]
srmr.ci.orig.mod.missing.below <- 
  list("n150"=missing.rate.below.srmr(pop.srmr, srmr.ci.mod.orig.n150),
       "n200"=missing.rate.below.srmr(pop.srmr, srmr.ci.mod.orig.n200),
       "n300"=missing.rate.below.srmr(pop.srmr, srmr.ci.mod.orig.n300),
       "n500"=missing.rate.below.srmr(pop.srmr, srmr.ci.mod.orig.n500),
       "n800"=missing.rate.below.srmr(pop.srmr, srmr.ci.mod.orig.n800),
       "n1000"=missing.rate.below.srmr(pop.srmr, srmr.ci.mod.orig.n1000))

srmr.ci.orig.mod.missing.below
save(srmr.ci.orig.mod.missing.below, file="srmr.ci.orig.mod.missing.below.RData")


srmr.ci.orig.mod.missing.above <- 
  list("n150"=missing.rate.above.srmr(pop.srmr, srmr.ci.mod.orig.n150),
       "n200"=missing.rate.above.srmr(pop.srmr, srmr.ci.mod.orig.n200),
       "n300"=missing.rate.above.srmr(pop.srmr, srmr.ci.mod.orig.n300),
       "n500"=missing.rate.above.srmr(pop.srmr, srmr.ci.mod.orig.n500),
       "n800"=missing.rate.above.srmr(pop.srmr, srmr.ci.mod.orig.n800),
       "n1000"=missing.rate.above.srmr(pop.srmr, srmr.ci.mod.orig.n1000))

srmr.ci.orig.mod.missing.above
save(srmr.ci.orig.mod.missing.above, file="srmr.ci.orig.mod.missing.above.RData")

srmr.ci.high.mod.missing.below <- 
  list("n150"=missing.rate.below.srmr(pop.srmr, srmr.ci.mod.high.n150),
       "n200"=missing.rate.below.srmr(pop.srmr, srmr.ci.mod.high.n200),
       "n300"=missing.rate.below.srmr(pop.srmr, srmr.ci.mod.high.n300),
       "n500"=missing.rate.below.srmr(pop.srmr, srmr.ci.mod.high.n500),
       "n800"=missing.rate.below.srmr(pop.srmr, srmr.ci.mod.high.n800),
       "n1000"=missing.rate.below.srmr(pop.srmr, srmr.ci.mod.high.n1000))

srmr.ci.high.mod.missing.below
save(srmr.ci.high.mod.missing.below, file="srmr.ci.high.mod.missing.below.RData")



srmr.ci.high.mod.missing.above <- 
  list("n150"=missing.rate.above.srmr(pop.srmr, srmr.ci.mod.high.n150),
       "n200"=missing.rate.above.srmr(pop.srmr, srmr.ci.mod.high.n200),
       "n300"=missing.rate.above.srmr(pop.srmr, srmr.ci.mod.high.n300),
       "n500"=missing.rate.above.srmr(pop.srmr, srmr.ci.mod.high.n500),
       "n800"=missing.rate.above.srmr(pop.srmr, srmr.ci.mod.high.n800),
       "n1000"=missing.rate.above.srmr(pop.srmr, srmr.ci.mod.high.n1000))

srmr.ci.high.mod.missing.above
save(srmr.ci.high.mod.missing.above, file="srmr.ci.high.mod.missing.above.RData")




pop.cfi <- pop.indices["pop.cfi",]
cfi.ci.orig.mod.missing.below <- 
  list("n150"=missing.rate.below.cfi(pop.cfi, cfi.ci.mod.orig.n150),
       "n200"=missing.rate.below.cfi(pop.cfi, cfi.ci.mod.orig.n200),
       "n300"=missing.rate.below.cfi(pop.cfi, cfi.ci.mod.orig.n300),
       "n500"=missing.rate.below.cfi(pop.cfi, cfi.ci.mod.orig.n500),
       "n800"=missing.rate.below.cfi(pop.cfi, cfi.ci.mod.orig.n800),
       "n1000"=missing.rate.below.cfi(pop.cfi, cfi.ci.mod.orig.n1000))

cfi.ci.orig.mod.missing.below
save(cfi.ci.orig.mod.missing.below, file="cfi.ci.orig.mod.missing.below.RData")


cfi.ci.orig.mod.missing.above <- 
  list("n150"=missing.rate.above.cfi(pop.cfi, cfi.ci.mod.orig.n150),
       "n200"=missing.rate.above.cfi(pop.cfi, cfi.ci.mod.orig.n200),
       "n300"=missing.rate.above.cfi(pop.cfi, cfi.ci.mod.orig.n300),
       "n500"=missing.rate.above.cfi(pop.cfi, cfi.ci.mod.orig.n500),
       "n800"=missing.rate.above.cfi(pop.cfi, cfi.ci.mod.orig.n800),
       "n1000"=missing.rate.above.cfi(pop.cfi, cfi.ci.mod.orig.n1000))

cfi.ci.orig.mod.missing.above
save(cfi.ci.orig.mod.missing.above, file="cfi.ci.orig.mod.missing.above.RData")

cfi.ci.high.mod.missing.below <- 
  list("n150"=missing.rate.below.cfi(pop.cfi, cfi.ci.mod.high.n150),
       "n200"=missing.rate.below.cfi(pop.cfi, cfi.ci.mod.high.n200),
       "n300"=missing.rate.below.cfi(pop.cfi, cfi.ci.mod.high.n300),
       "n500"=missing.rate.below.cfi(pop.cfi, cfi.ci.mod.high.n500),
       "n800"=missing.rate.below.cfi(pop.cfi, cfi.ci.mod.high.n800),
       "n1000"=missing.rate.below.cfi(pop.cfi, cfi.ci.mod.high.n1000))

cfi.ci.high.mod.missing.below
save(cfi.ci.high.mod.missing.below, file="cfi.ci.high.mod.missing.below.RData")



cfi.ci.high.mod.missing.above <- 
  list("n150"=missing.rate.above.cfi(pop.cfi, cfi.ci.mod.high.n150),
       "n200"=missing.rate.above.cfi(pop.cfi, cfi.ci.mod.high.n200),
       "n300"=missing.rate.above.cfi(pop.cfi, cfi.ci.mod.high.n300),
       "n500"=missing.rate.above.cfi(pop.cfi, cfi.ci.mod.high.n500),
       "n800"=missing.rate.above.cfi(pop.cfi, cfi.ci.mod.high.n800),
       "n1000"=missing.rate.above.cfi(pop.cfi, cfi.ci.mod.high.n1000))

cfi.ci.high.mod.missing.above
save(cfi.ci.high.mod.missing.above, file="cfi.ci.high.mod.missing.above.RData")




