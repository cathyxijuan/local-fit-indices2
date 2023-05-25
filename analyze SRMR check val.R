source("functions.R")
load("pop.indices.RData")
load("srmr.ci.check.val.mod.orig.n150.RData")
load("srmr.ci.check.val.mod.orig.n200.RData")
load("srmr.ci.check.val.mod.orig.n300.RData")
load("srmr.ci.check.val.mod.orig.n500.RData")
load("srmr.ci.check.val.mod.orig.n800.RData")
load("srmr.ci.check.val.mod.orig.n1000.RData")


srmr.check.val.orig.mod <- list("n150"=list.mean(srmr.ci.check.val.mod.orig.n200),
                                "n200"=list.mean(srmr.ci.check.val.mod.orig.n200),
                                "n300"=list.mean(srmr.ci.check.val.mod.orig.n300),
                                "n500"=list.mean(srmr.ci.check.val.mod.orig.n500),
                                "n800"=list.mean(srmr.ci.check.val.mod.orig.n800),
                                "n1000"=list.mean(srmr.ci.check.val.mod.orig.n1000))

save(srmr.check.val.orig.mod, file="srmr.check.val.orig.mod.RData")