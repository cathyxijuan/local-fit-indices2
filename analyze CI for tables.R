load("rmsea.ci.orig.mod.RData")
load("rmsea.ci.high.mod.RData")
load("cfi.ci.orig.mod.RData")
load("cfi.ci.high.mod.RData")
load("srmr.ci.orig.mod.RData")
load("srmr.ci.high.mod.RData")

rmsea200 <- rmsea.ci.orig.mod$n200

apply(rmsea200, 2, function(x) which.max(abs(x-0.9)))
apply(rmsea200, 2, function(x) which.min(abs(x-0.9)))

rmsea200 <- rmsea.ci.high.mod$n200

apply(rmsea200, 2, function(x) which.max(abs(x-0.9)))
apply(rmsea200, 2, function(x) which.min(abs(x-0.9)))


cfi200 <- cfi.ci.orig.mod$n200

apply(cfi200, 2, function(x) which.max(abs(x-0.9)))
apply(cfi200, 2, function(x) which.min(abs(x-0.9)))

cfi200 <- cfi.ci.high.mod$n200

apply(cfi200, 2, function(x) which.max(abs(x-0.9)))
apply(cfi200, 2, function(x)  which.min(abs(x-0.9)))


which.min(rowMeans(abs(cfi.ci.high.mod$n200-0.9)))
which.min(rowMeans(cfi.ci.high.mod$n200)-0.9)
which.max(rowMeans(cfi.ci.high.mod$n200)-0.9)
which.max(rowMeans(cfi.ci.high.mod$n200-0.9))
0.9-rowMeans(abs(cfi.ci.high.mod$n200-0.9))


srmr200 <- srmr.ci.orig.mod$n200

apply(srmr200, 2, function(x) which.max(abs(x-0.9)))
apply(srmr200, 2, function(x) which.min(abs(x-0.9)))

srmr200 <- srmr.ci.high.mod$n200

apply(srmr200, 2, function(x) which.max(abs(x-0.9)))
apply(srmr200, 2, function(x) which.min(abs(x-0.9)))





rmsea500 <- rmsea.ci.orig.mod$n500

apply(rmsea500, 2, function(x) which.max(abs(x-0.9)))
apply(rmsea500, 2, function(x) which.min(abs(x-0.9)))

rmsea500 <- rmsea.ci.high.mod$n500

apply(rmsea500, 2, function(x) which.max(abs(x-0.9)))
apply(rmsea500, 2, function(x) which.min(abs(x-0.9)))


cfi500 <- cfi.ci.orig.mod$n500

apply(cfi500, 2, function(x) which.max(abs(x-0.9)))
apply(cfi500, 2, function(x) which.min(abs(x-0.9)))

cfi500 <- cfi.ci.high.mod$n500

apply(cfi500, 2, function(x) which.max(abs(x-0.9)))
apply(cfi500, 2, function(x) which.min(abs(x-0.9)))


srmr500 <- srmr.ci.orig.mod$n500

apply(srmr500, 2, function(x) which.max(abs(x-0.9)))
apply(srmr500, 2, function(x) which.min(abs(x-0.9)))

srmr500 <- srmr.ci.high.mod$n500

apply(srmr500, 2, function(x) which.max(abs(x-0.9)))
apply(srmr500, 2, function(x) which.min(abs(x-0.9)))


