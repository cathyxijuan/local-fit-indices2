source("functions.R")
load("pop.indices.RData")
load("pop.fit.mod.orig.RData")
load("rmsea.ci.mod.orig.n200.RData")
load("rmsea.ci.mod.orig.n300.RData")
load("rmsea.ci.mod.orig.n500.RData")
load("rmsea.ci.mod.orig.n800.RData")
load("rmsea.ci.mod.orig.n1000.RData")

pop.rmsea <- pop.indices["pop.rmsea",]
rmsea.ci.orig.mod <- list(
                          "n200"=ci.coverage(pop.rmsea, rmsea.ci.mod.orig.n200),
                          "n300"=ci.coverage(pop.rmsea, rmsea.ci.mod.orig.n300),
                          "n500"=ci.coverage(pop.rmsea, rmsea.ci.mod.orig.n500),
                          "n800"=ci.coverage(pop.rmsea, rmsea.ci.mod.orig.n800),
                          "n1000"=ci.coverage(pop.rmsea, rmsea.ci.mod.orig.n1000))

rmsea.ci.orig.mod
save(rmsea.ci.orig.mod, file="rmsea.ci.orig.mod.RData")




ci.coverage(pop.rmsea, rmsea.ci.mod.orig.n200)
