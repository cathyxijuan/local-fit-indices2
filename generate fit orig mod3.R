source("models original.R")
source("functions3.R")





set.seed(777)
fit.mod.orig.v3.n150 <- simu.fit(pop.mod, sat.struct, path.mod.list, 
                              sample.size=150, rep.num=5000)
save(fit.mod.orig.v3.n150, file="fit.mod.orig.v3.n150.RData")





set.seed(777)
fit.mod.orig.v3.n200 <- simu.fit(pop.mod, sat.struct, path.mod.list, 
                              sample.size=200, rep.num=5000)
save(fit.mod.orig.v3.n200, file="fit.mod.orig.v3.n200.RData")





set.seed(777)
fit.mod.orig.v3.n300 <- simu.fit(pop.mod, sat.struct, path.mod.list, 
                              sample.size=300, rep.num=5000)
save(fit.mod.orig.v3.n300, file="fit.mod.orig.v3.n300.RData")





set.seed(777)
fit.mod.orig.v3.n500 <- simu.fit(pop.mod, sat.struct, path.mod.list, 
                              sample.size=500, rep.num=5000)
save(fit.mod.orig.v3.n500, file="fit.mod.orig.v3.n500.RData")



set.seed(777)
fit.mod.orig.v3.n800 <- simu.fit(pop.mod, sat.struct, path.mod.list, 
                              sample.size=800, rep.num=5000)
save(fit.mod.orig.v3.n800, file="fit.mod.orig.v3.n800.RData")


set.seed(777)
fit.mod.orig.v3.n1000 <- simu.fit(pop.mod, sat.struct, path.mod.list, 
                               sample.size=1000, rep.num=5000)
save(fit.mod.orig.v3.n1000, file="fit.mod.orig.v3.n1000.RData")



