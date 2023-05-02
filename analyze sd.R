source("functions.R")
load("pop.indices.RData")

load("fit.mod.orig.n150.RData")
load("fit.mod.orig.n200.RData")
load("fit.mod.orig.n300.RData")
load("fit.mod.orig.n500.RData")
load("fit.mod.orig.n800.RData")
load("fit.mod.orig.n1000.RData")






sd.orig.mod <- list("n150"=list.sd(fit.mod.orig.n150),
                       "n200"=list.sd(fit.mod.orig.n200),
                       "n300"=list.sd(fit.mod.orig.n300),
                       "n500"=list.sd(fit.mod.orig.n500),
                       "n800"=list.sd(fit.mod.orig.n800),
                       "n1000"=list.sd(fit.mod.orig.n1000))

sd.orig.mod
save(sd.orig.mod, file="sd.orig.mod.RData")







###high reliability
load("fit.mod.high.n150.RData")
load("fit.mod.high.n200.RData")
load("fit.mod.high.n300.RData")
load("fit.mod.high.n500.RData")
load("fit.mod.high.n800.RData")
load("fit.mod.high.n1000.RData")






sd.high.mod <- list("n150"=list.sd(fit.mod.high.n150),
                       "n200"=list.sd(fit.mod.high.n200),
                       "n300"=list.sd(fit.mod.high.n300),
                       "n500"=list.sd(fit.mod.high.n500),
                       "n800"=list.sd(fit.mod.high.n800),
                       "n1000"=list.sd(fit.mod.high.n1000))
sd.high.mod 
save(sd.high.mod, file="sd.high.mod.RData")




##low reliability

# load("fit.mod.low.n150.RData")
# load("fit.mod.low.n200.RData")
# load("fit.mod.low.n300.RData")
# load("fit.mod.low.n500.RData")
# load("fit.mod.low.n800.RData")
# load("fit.mod.low.n1000.RData")
# 
# 
# 
# 
# 
# 
# sd.low.mod <- list("n150"=list.sd(fit.mod.low.n150),
#                       "n200"=list.sd(fit.mod.low.n200),
#                       "n300"=list.sd(fit.mod.low.n300),
#                       "n500"=list.sd(fit.mod.low.n500),
#                       "n800"=list.sd(fit.mod.low.n800),
#                       "n1000"=list.sd(fit.mod.low.n1000))
# 
# sd.low.mod
# save(sd.low.mod, file="sd.low.mod.RData")
