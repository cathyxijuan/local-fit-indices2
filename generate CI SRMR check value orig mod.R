source("models high reliable.R")
source("functions.R")

set.seed(777)
srmr.ci.check.val.mod.high.n150 <- simu.srmr.ci.check.val(pop.mod, sat.struct, path.mod.list, 
                                                          sample.size=150, rep.num=5000)

save(srmr.ci.check.val.mod.high.n150, file="srmr.ci.check.val.mod.high.n150.RData")





set.seed(777)
srmr.ci.check.val.mod.high.n200 <- simu.srmr.ci.check.val(pop.mod, sat.struct, path.mod.list, 
                                                          sample.size=200, rep.num=5000)

save(srmr.ci.check.val.mod.high.n200, file="srmr.ci.check.val.mod.high.n200.RData")




set.seed(777)
srmr.ci.check.val.mod.high.n300 <- simu.srmr.ci.check.val(pop.mod, sat.struct, path.mod.list, 
                                                          sample.size=300, rep.num=5000)

save(srmr.ci.check.val.mod.high.n300, file="srmr.ci.check.val.mod.high.n300.RData")




set.seed(777)
srmr.ci.check.val.mod.high.n500 <- simu.srmr.ci.check.val(pop.mod, sat.struct, path.mod.list, 
                                                          sample.size=500, rep.num=5000)

save(srmr.ci.check.val.mod.high.n500, file="srmr.ci.check.val.mod.high.n500.RData")



set.seed(777)
srmr.ci.check.val.mod.high.n800 <- simu.srmr.ci.check.val(pop.mod, sat.struct, path.mod.list, 
                                                          sample.size=800, rep.num=5000)

save(srmr.ci.check.val.mod.high.n800, file="srmr.ci.check.val.mod.high.n800.RData")


set.seed(777)
srmr.ci.check.val.mod.high.n1000 <- simu.srmr.ci.check.val(pop.mod, sat.struct, path.mod.list, 
                                                           sample.size=1000, rep.num=5000)

save(srmr.ci.check.val.mod.high.n1000, file="srmr.ci.check.val.mod.high.n1000.RData")

