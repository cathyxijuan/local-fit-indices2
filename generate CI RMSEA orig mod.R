source("models original.R")
source("functions.R")




set.seed(777)
rmsea.ci.mod.orig.n150 <- simu.rmsea.ci(pop.mod, sat.struct, path.mod.list, 
                                        sample.size=150, rep.num=5000)

save(rmsea.ci.mod.orig.n150, file="rmsea.ci.mod.orig.n150.RData")





set.seed(777)
rmsea.ci.mod.orig.n200 <- simu.rmsea.ci(pop.mod, sat.struct, path.mod.list, 
                                        sample.size=200, rep.num=5000)

save(rmsea.ci.mod.orig.n200, file="rmsea.ci.mod.orig.n200.RData")




set.seed(777)
rmsea.ci.mod.orig.n300 <- simu.rmsea.ci(pop.mod, sat.struct, path.mod.list, 
                                        sample.size=300, rep.num=5000)

save(rmsea.ci.mod.orig.n300, file="rmsea.ci.mod.orig.n300.RData")




set.seed(777)
rmsea.ci.mod.orig.n500 <- simu.rmsea.ci(pop.mod, sat.struct, path.mod.list, 
                                        sample.size=500, rep.num=5000)

save(rmsea.ci.mod.orig.n500, file="rmsea.ci.mod.orig.n500.RData")



set.seed(777)
rmsea.ci.mod.orig.n800 <- simu.rmsea.ci(pop.mod, sat.struct, path.mod.list, 
                                        sample.size=800, rep.num=5000)

save(rmsea.ci.mod.orig.n800, file="rmsea.ci.mod.orig.n800.RData")


set.seed(777)
rmsea.ci.mod.orig.n1000 <- simu.rmsea.ci(pop.mod, sat.struct, path.mod.list, 
                                         sample.size=1000, rep.num=5000)

save(rmsea.ci.mod.orig.n1000, file="rmsea.ci.mod.orig.n1000.RData")