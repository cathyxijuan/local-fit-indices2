source("functions.R")
load("pop.indices.RData")

load("fit.mod.orig.n150.RData")
load("fit.mod.orig.n200.RData")
load("fit.mod.orig.n300.RData")
load("fit.mod.orig.n500.RData")
load("fit.mod.orig.n800.RData")
load("fit.mod.orig.n1000.RData")






means.orig.mod <- list("n150"=list.mean(fit.mod.orig.n150),
                       "n200"=list.mean(fit.mod.orig.n200),
                       "n300"=list.mean(fit.mod.orig.n300),
                       "n500"=list.mean(fit.mod.orig.n500),
                       "n800"=list.mean(fit.mod.orig.n800),
                       "n1000"=list.mean(fit.mod.orig.n1000))

means.orig.mod

bias.orig.mod <- vector(mode="list")

for(i in 1:length(means.orig.mod)){
  bias.orig.mod[[names(means.orig.mod)[i]]] <- 
    round(means.orig.mod[[i]] - pop.indices, 4)
}

bias.orig.mod


save(bias.orig.mod, file="bias.orig.mod.RData")




###high reliability
load("fit.mod.high.n150.RData")
load("fit.mod.high.n200.RData")
load("fit.mod.high.n300.RData")
load("fit.mod.high.n500.RData")
load("fit.mod.high.n800.RData")
load("fit.mod.high.n1000.RData")






means.high.mod <- list("n150"=list.mean(fit.mod.high.n150),
                       "n200"=list.mean(fit.mod.high.n200),
                       "n300"=list.mean(fit.mod.high.n300),
                       "n500"=list.mean(fit.mod.high.n500),
                       "n800"=list.mean(fit.mod.high.n800),
                       "n1000"=list.mean(fit.mod.high.n1000))


means.high.mod

bias.high.mod <- vector(mode="list")

for(i in 1:length(means.high.mod)){
  bias.high.mod[[names(means.high.mod)[i]]] <- 
    round(means.high.mod[[i]] - pop.indices, 4)
}

bias.high.mod

save(bias.high.mod, file="bias.high.mod.RData")



# ##low reliability
# 
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
# means.low.mod <- list("n150"=list.mean(fit.mod.low.n150),
#                       "n200"=list.mean(fit.mod.low.n200),
#                       "n300"=list.mean(fit.mod.low.n300),
#                       "n500"=list.mean(fit.mod.low.n500),
#                       "n800"=list.mean(fit.mod.low.n800),
#                       "n1000"=list.mean(fit.mod.low.n1000))
# 
# bias.low.mod <- vector(mode="list")
# 
# for(i in 1:length(means.low.mod)){
#   bias.low.mod[[names(means.low.mod)[i]]] <- 
#     round(means.low.mod[[i]] - pop.indices, 4)
# }
# 
# bias.low.mod
# save(bias.low.mod, file="bias.low.mode.RData")
