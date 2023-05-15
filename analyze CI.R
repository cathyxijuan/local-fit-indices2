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
rmsea.ci.orig.mod <- list("n150"=ci.coverage(pop.rmsea, rmsea.ci.mod.orig.n200),
                          "n200"=ci.coverage(pop.rmsea, rmsea.ci.mod.orig.n200),
                          "n300"=ci.coverage(pop.rmsea, rmsea.ci.mod.orig.n300),
                          "n500"=ci.coverage(pop.rmsea, rmsea.ci.mod.orig.n500),
                          "n800"=ci.coverage(pop.rmsea, rmsea.ci.mod.orig.n800),
                          "n1000"=ci.coverage(pop.rmsea, rmsea.ci.mod.orig.n1000))

rmsea.ci.orig.mod
save(rmsea.ci.orig.mod, file="rmsea.ci.orig.mod.RData")




pop.rmsea <- pop.indices["pop.rmsea",]
rmsea.ci.high.mod <- list("n150"=ci.coverage(pop.rmsea, rmsea.ci.mod.high.n200),
                          "n200"=ci.coverage(pop.rmsea, rmsea.ci.mod.high.n200),
                          "n300"=ci.coverage(pop.rmsea, rmsea.ci.mod.high.n300),
                          "n500"=ci.coverage(pop.rmsea, rmsea.ci.mod.high.n500),
                          "n800"=ci.coverage(pop.rmsea, rmsea.ci.mod.high.n800),
                          "n1000"=ci.coverage(pop.rmsea, rmsea.ci.mod.high.n1000))

rmsea.ci.high.mod
save(rmsea.ci.high.mod, file="rmsea.ci.high.mod.RData")



pop.rmsea <- pop.indices["pop.rmsea",]
rmsea.ci.orig.mod <- list("n150"=ci.coverage(pop.rmsea, rmsea.ci.mod.orig.n200),
                          "n200"=ci.coverage(pop.rmsea, rmsea.ci.mod.orig.n200),
                          "n300"=ci.coverage(pop.rmsea, rmsea.ci.mod.orig.n300),
                          "n500"=ci.coverage(pop.rmsea, rmsea.ci.mod.orig.n500),
                          "n800"=ci.coverage(pop.rmsea, rmsea.ci.mod.orig.n800),
                          "n1000"=ci.coverage(pop.rmsea, rmsea.ci.mod.orig.n1000))

rmsea.ci.orig.mod
save(rmsea.ci.orig.mod, file="rmsea.ci.orig.mod.RData")




pop.rmsea <- pop.indices["pop.rmsea",]
rmsea.ci.high.mod <- list("n150"=ci.coverage(pop.rmsea, rmsea.ci.mod.high.n200),
                          "n200"=ci.coverage(pop.rmsea, rmsea.ci.mod.high.n200),
                          "n300"=ci.coverage(pop.rmsea, rmsea.ci.mod.high.n300),
                          "n500"=ci.coverage(pop.rmsea, rmsea.ci.mod.high.n500),
                          "n800"=ci.coverage(pop.rmsea, rmsea.ci.mod.high.n800),
                          "n1000"=ci.coverage(pop.rmsea, rmsea.ci.mod.high.n1000))

rmsea.ci.high.mod
save(rmsea.ci.high.mod, file="rmsea.ci.high.mod.RData")





pop.srmr <- pop.indices["pop.srmr",]
srmr.ci.orig.mod <- list("n150"=ci.coverage(pop.srmr, srmr.ci.mod.orig.n200),
                         "n200"=ci.coverage(pop.srmr, srmr.ci.mod.orig.n200),
                         "n300"=ci.coverage(pop.srmr, srmr.ci.mod.orig.n300),
                         "n500"=ci.coverage(pop.srmr, srmr.ci.mod.orig.n500),
                         "n800"=ci.coverage(pop.srmr, srmr.ci.mod.orig.n800),
                         "n1000"=ci.coverage(pop.srmr, srmr.ci.mod.orig.n1000))

srmr.ci.orig.mod
save(srmr.ci.orig.mod, file="srmr.ci.orig.mod.RData")




pop.srmr <- pop.indices["pop.srmr",]
srmr.ci.high.mod <- list("n150"=ci.coverage(pop.srmr, srmr.ci.mod.high.n200),
                         "n200"=ci.coverage(pop.srmr, srmr.ci.mod.high.n200),
                         "n300"=ci.coverage(pop.srmr, srmr.ci.mod.high.n300),
                         "n500"=ci.coverage(pop.srmr, srmr.ci.mod.high.n500),
                         "n800"=ci.coverage(pop.srmr, srmr.ci.mod.high.n800),
                         "n1000"=ci.coverage(pop.srmr, srmr.ci.mod.high.n1000))

srmr.ci.high.mod
save(srmr.ci.high.mod, file="srmr.ci.high.mod.RData")

pop.cfi <- pop.indices["pop.cfi",]
cfi.ci.orig.mod <- list("n150"=ci.coverage.cfi(pop.cfi, cfi.ci.mod.orig.n200),
                        "n200"=ci.coverage.cfi(pop.cfi, cfi.ci.mod.orig.n200),
                        "n300"=ci.coverage.cfi(pop.cfi, cfi.ci.mod.orig.n300),
                        "n500"=ci.coverage.cfi(pop.cfi, cfi.ci.mod.orig.n500),
                        "n800"=ci.coverage.cfi(pop.cfi, cfi.ci.mod.orig.n800),
                        "n1000"=ci.coverage.cfi(pop.cfi, cfi.ci.mod.orig.n1000))

cfi.ci.orig.mod
save(cfi.ci.orig.mod, file="cfi.ci.orig.mod.RData")




pop.cfi <- pop.indices["pop.cfi",]
cfi.ci.high.mod <- list("n150"=ci.coverage.cfi(pop.cfi, cfi.ci.mod.high.n200),
                        "n200"=ci.coverage.cfi(pop.cfi, cfi.ci.mod.high.n200),
                        "n300"=ci.coverage.cfi(pop.cfi, cfi.ci.mod.high.n300),
                        "n500"=ci.coverage.cfi(pop.cfi, cfi.ci.mod.high.n500),
                        "n800"=ci.coverage.cfi(pop.cfi, cfi.ci.mod.high.n800),
                        "n1000"=ci.coverage.cfi(pop.cfi, cfi.ci.mod.high.n1000))

cfi.ci.high.mod
save(cfi.ci.high.mod, file="cfi.ci.high.mod.RData")
