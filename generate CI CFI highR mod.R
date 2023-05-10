source("models high reliable.R")
source("functions.R")




set.seed(777)
cfi.ci.mod.high.n150 <- simu.cfi.ci(pop.mod, sat.struct, path.mod.list, 
                                        sample.size=150, rep.num=5000)

save(cfi.ci.mod.high.n150, file="cfi.ci.mod.high.n150.RData")





set.seed(777)
cfi.ci.mod.high.n200 <- simu.cfi.ci(pop.mod, sat.struct, path.mod.list, 
                                        sample.size=200, rep.num=5000)

save(cfi.ci.mod.high.n200, file="cfi.ci.mod.high.n200.RData")




set.seed(777)
cfi.ci.mod.high.n300 <- simu.cfi.ci(pop.mod, sat.struct, path.mod.list, 
                                        sample.size=300, rep.num=5000)

save(cfi.ci.mod.high.n300, file="cfi.ci.mod.high.n300.RData")




set.seed(777)
cfi.ci.mod.high.n500 <- simu.cfi.ci(pop.mod, sat.struct, path.mod.list, 
                                        sample.size=500, rep.num=5000)

save(cfi.ci.mod.high.n500, file="cfi.ci.mod.high.n500.RData")



set.seed(777)
cfi.ci.mod.high.n800 <- simu.cfi.ci(pop.mod, sat.struct, path.mod.list, 
                                        sample.size=800, rep.num=5000)

save(cfi.ci.mod.high.n800, file="cfi.ci.mod.high.n800.RData")


set.seed(777)
cfi.ci.mod.high.n1000 <- simu.cfi.ci(pop.mod, sat.struct, path.mod.list, 
                                         sample.size=1000, rep.num=5000)

save(cfi.ci.mod.high.n1000, file="cfi.ci.mod.high.n1000.RData")