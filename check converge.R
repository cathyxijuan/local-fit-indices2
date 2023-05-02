source("functions.R")
source("models original.R")



load("fit.mod.orig.n150.RData")
load("fit.mod.orig.n200.RData")
load("fit.mod.orig.n300.RData")
load("fit.mod.orig.n500.RData")
load("fit.mod.orig.n800.RData")
load("fit.mod.orig.n1000.RData")

count.null(fit.mod.orig.n150) #207
count.null(fit.mod.orig.n200) #56
count.null(fit.mod.orig.n300) #6
count.null(fit.mod.orig.n500) #0
count.null(fit.mod.orig.n800) #0
count.null(fit.mod.orig.n1000) #0

set.seed(777)
nonconverge.mod.orig.n150 <- nonconverge.count(pop.mod, sat.struct, 
                                               sample.size=150, rep.num=5000)
nonconverge.mod.orig.n150 #0


   
set.seed(777)
nonconverge.mod.orig.n200 <- nonconverge.count(pop.mod, sat.struct, 
                                               sample.size=200, rep.num=5000)
nonconverge.mod.orig.n200 #0



#####low reliability
source("functions.R")
source("models low reliable.R")
load("fit.mod.low.n150.RData")
load("fit.mod.low.n200.RData")
load("fit.mod.low.n300.RData")
load("fit.mod.low.n500.RData")
load("fit.mod.low.n800.RData")
load("fit.mod.low.n1000.RData")



count.null(fit.mod.low.n150) #4077
count.null(fit.mod.low.n200) #3424
count.null(fit.mod.low.n300) #2205
count.null(fit.mod.low.n500) #894
count.null(fit.mod.low.n800) #272
count.null(fit.mod.low.n1000) #118


set.seed(777)
nonconverge.mod.low.n150 <- nonconverge.count(pop.mod, sat.struct, 
                                              sample.size=150, rep.num=5000)
nonconverge.mod.low.n150 #345


set.seed(777)
nonconverge.mod.low.n200 <- nonconverge.count(pop.mod, sat.struct, 
                                              sample.size=200, rep.num=5000)
nonconverge.mod.low.n200 #112


set.seed(777)
nonconverge.mod.low.n300 <- nonconverge.count(pop.mod, sat.struct, 
                                              sample.size=300, rep.num=5000)
nonconverge.mod.low.n300 #14


set.seed(777)
nonconverge.mod.low.n500 <- nonconverge.count(pop.mod, sat.struct, 
                                              sample.size=500, rep.num=5000)
nonconverge.mod.low.n500 #0


set.seed(777)
nonconverge.mod.low.n800 <- nonconverge.count(pop.mod, sat.struct, 
                                              sample.size=800, rep.num=5000)
nonconverge.mod.low.n800 #0


set.seed(777)
nonconverge.mod.low.n1000 <- nonconverge.count(pop.mod, sat.struct, 
                                               sample.size=1000, rep.num=5000)
nonconverge.mod.low.n1000 #0




##### high reliability
source("functions.R")
source("models high reliable.R")
load("fit.mod.high.n150.RData")
load("fit.mod.high.n200.RData")
load("fit.mod.high.n300.RData")
load("fit.mod.high.n500.RData")
load("fit.mod.high.n800.RData")
load("fit.mod.high.n1000.RData")



count.null(fit.mod.high.n150) #0
count.null(fit.mod.high.n200) #0
count.null(fit.mod.high.n300) #0
count.null(fit.mod.high.n500) #0
count.null(fit.mod.high.n800) #0
count.null(fit.mod.high.n1000) #0
