source("models low reliable.R")
source("functions.R")


set.seed(7)
fit.mod.low.n150 <- simu.fit(pop.mod, sat.struct, path.mod.list, 
                             sample.size=150, rep.num=5000)
save(fit.mod.low.n150, file="fit.mod.low.n150.RData")
sum(sapply(fit.mod.low.n150, is.null))


set.seed(777)
fit.mod.low.n200 <- simu.fit(pop.mod, sat.struct, path.mod.list, 
                              sample.size=200, rep.num=5000)
save(fit.mod.low.n200, file="fit.mod.low.n200.RData")
sum(sapply(fit.mod.low.n200, is.null))




set.seed(777)
fit.mod.low.n300 <- simu.fit(pop.mod, sat.struct, path.mod.list, 
                              sample.size=300, rep.num=5000)
save(fit.mod.low.n300, file="fit.mod.low.n300.RData")

sum(sapply(fit.mod.low.n300, is.null))




set.seed(777)
fit.mod.low.n500 <- simu.fit(pop.mod, sat.struct, path.mod.list, 
                              sample.size=500, rep.num=5000)
save(fit.mod.low.n500, file="fit.mod.low.n500.RData")
sum(sapply(fit.mod.low.n500, is.null))



set.seed(777)
fit.mod.low.n800 <- simu.fit(pop.mod, sat.struct, path.mod.list, 
                              sample.size=800, rep.num=5000)
save(fit.mod.low.n800, file="fit.mod.low.n800.RData")
sum(sapply(fit.mod.low.n800, is.null))


set.seed(777)
fit.mod.low.n1000 <- simu.fit(pop.mod, sat.struct, path.mod.list, 
                               sample.size=1000, rep.num=5000)
save(fit.mod.low.n1000, file="fit.mod.low.n1000.RData")
sum(sapply(fit.mod.low.n1000, is.null))
