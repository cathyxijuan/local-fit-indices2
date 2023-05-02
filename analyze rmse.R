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

rmse.orig.mod
save(rmse.orig.mod, file="rmse.orig.mod.RData")







###high reliability
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
save(rmse.high.mod, file="rmse.high.mod.RData")




##low reliability

load("fit.mod.low.n150.RData")
load("fit.mod.low.n200.RData")
load("fit.mod.low.n300.RData")
load("fit.mod.low.n500.RData")
load("fit.mod.low.n800.RData")
load("fit.mod.low.n1000.RData")






rmse.low.mod <- list("n150"=rmse(fit.mod.low.n150,pop.indices),
                   "n200"=rmse(fit.mod.low.n200,pop.indices),
                   "n300"=rmse(fit.mod.low.n300,pop.indices),
                   "n500"=rmse(fit.mod.low.n500,pop.indices),
                   "n800"=rmse(fit.mod.low.n800,pop.indices),
                   "n1000"=rmse(fit.mod.low.n1000,pop.indices))

rmse.low.mod
save(rmse.low.mod, file="rmse.low.mod.RData")
