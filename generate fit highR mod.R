source("models high reliable.R")
source("functions.R")






set.seed(777)
fit.mod.high.n150 <- simu.fit(pop.mod, sat.struct, path.mod.list, 
                              sample.size=150, rep.num=5000)
save(fit.mod.high.n150, file="fit.mod.high.n150.RData")
sum(sapply(fit.mod.high.n150, is.null))


set.seed(777)
fit.mod.high.n200 <- simu.fit(pop.mod, sat.struct, path.mod.list, 
                              sample.size=200, rep.num=5000)
save(fit.mod.high.n200, file="fit.mod.high.n200.RData")
sum(sapply(fit.mod.high.n200, is.null))




set.seed(777)
fit.mod.high.n300 <- simu.fit(pop.mod, sat.struct, path.mod.list, 
                              sample.size=300, rep.num=5000)
save(fit.mod.high.n300, file="fit.mod.high.n300.RData")

sum(sapply(fit.mod.high.n300, is.null))




set.seed(777)
fit.mod.high.n500 <- simu.fit(pop.mod, sat.struct, path.mod.list, 
                              sample.size=500, rep.num=5000)
save(fit.mod.high.n500, file="fit.mod.high.n500.RData")
sum(sapply(fit.mod.high.n500, is.null))



set.seed(777)
fit.mod.high.n800 <- simu.fit(pop.mod, sat.struct, path.mod.list, 
                              sample.size=800, rep.num=5000)
save(fit.mod.high.n800, file="fit.mod.high.n800.RData")
sum(sapply(fit.mod.high.n800, is.null))


set.seed(777)
fit.mod.high.n1000 <- simu.fit(pop.mod, sat.struct, path.mod.list, 
                               sample.size=1000, rep.num=5000)
save(fit.mod.high.n1000, file="fit.mod.high.n1000.RData")
sum(sapply(fit.mod.high.n1000, is.null))
