source("models original.R")
source("functions.R")





set.seed(777)
fit.mod.orig.n150 <- simu.fit(pop.mod, sat.struct, path.mod.list, 
                              sample.size=150, rep.num=5000)
save(fit.mod.orig.n150, file="fit.mod.orig.n150.RData")

sum(sapply(fit.mod.orig.n150, is.null))



set.seed(777)
fit.mod.orig.n200 <- simu.fit(pop.mod, sat.struct, path.mod.list, 
                              sample.size=200, rep.num=5000)
save(fit.mod.orig.n200, file="fit.mod.orig.n200.RData")

sum(sapply(fit.mod.orig.n200, is.null))




set.seed(777)
fit.mod.orig.n300 <- simu.fit(pop.mod, sat.struct, path.mod.list, 
                              sample.size=300, rep.num=5000)
save(fit.mod.orig.n300, file="fit.mod.orig.n300.RData")

sum(sapply(fit.mod.orig.n300, is.null))




set.seed(777)
fit.mod.orig.n500 <- simu.fit(pop.mod, sat.struct, path.mod.list, 
                              sample.size=500, rep.num=5000)
save(fit.mod.orig.n500, file="fit.mod.orig.n500.RData")
sum(sapply(fit.mod.orig.n500, is.null))



set.seed(777)
fit.mod.orig.n800 <- simu.fit(pop.mod, sat.struct, path.mod.list, 
                              sample.size=800, rep.num=5000)
save(fit.mod.orig.n800, file="fit.mod.orig.n800.RData")
sum(sapply(fit.mod.orig.n800, is.null))


set.seed(777)
fit.mod.orig.n1000 <- simu.fit(pop.mod, sat.struct, path.mod.list, 
                               sample.size=1000, rep.num=5000)
save(fit.mod.orig.n1000, file="fit.mod.orig.n1000.RData")
sum(sapply(fit.mod.orig.n1000, is.null))



